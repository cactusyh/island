"""Query one actual assignment or summarize the versioned polymer benchmark."""

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from hashlib import sha256
from pathlib import Path

from island.exceptions import IslandError

from .assignment import assign_polymer_parameters, assign_radonpy
from .emc import FIELDS, _field_sources
from .radonpy import FAMILIES


def benchmark(evidence):
    """Summarize archived observations; preserve each library and stereo policy."""
    root = Path(evidence)
    panel = root / "reference_assignment.csv.gz"
    fields = defaultdict(Counter)
    published = Counter()
    gaps = defaultdict(Counter)
    with gzip.open(panel, "rt") as handle:
        for r in csv.DictReader(handle):
            field = r["field"]
            success = r["assigned"] == "True"
            fields[field].update(
                total=1,
                assigned=int(success),
                strict=int(success and r["bond_stereo_omitted"] != "True"),
            )
            if r["published_positive"] == "True":
                published.update(total=1, assigned=int(success))
            if not success:
                pairs = {
                    ",".join(sorted(pair.split(", ")))
                    for pair in re.findall(
                        r"increment pair \{([^}]+)\}", r["diagnostic"]
                    )
                }
                gaps[field].update(pairs)
    fixed = json.loads((root / "phase_4j94/summary.json").read_text())
    cure = json.loads((root / "phase_4j95_j97/cure-study.json").read_text())
    return {
        "panel_sha256": sha256(panel.read_bytes()).hexdigest(),
        "linear": dict(fields),
        "published_positive": dict(published),
        "missing_increments": dict(gaps),
        "fixed_graph_md": fixed,
        "md_runs_by_architecture": dict(
            Counter(
                "branched_block"
                if r["run"].startswith("branched_styrene_block")
                else "branch"
                if r["run"].startswith("ethyl_branch")
                else "block"
                for r in fixed
            )
        ),
        "periodic_epoxy": [
            {k: c[k] for k in ("dgeba", "pacm", "seed", "statistics")}
            for c in cure["cases"]
        ],
    }


def probe(
    structure,
    field,
    *,
    dp=None,
    stereo_policy="strict",
    field_dir=None,
    emc_mode="graph",
    cache_dir=None,
    head_end_group=None,
    tail_end_group=None,
):
    """Accept a repeat string or an already built graph through its native backend."""
    report = {
        "field": field,
        "input_kind": "repeat" if isinstance(structure, str) else "graph",
        "psmiles": structure if isinstance(structure, str) else None,
        "dp": dp,
        "stereo_policy": stereo_policy,
        "emc_mode": emc_mode,
        "hashes": {},
        "assignment": "unattempted",
        "single_point": "unmeasured",
        "dynamics": "unmeasured",
        "properties": "unmeasured",
    }
    if field not in {*FIELDS, *FAMILIES}:
        raise ValueError(f"Unknown field: {field}")
    if not isinstance(structure, str) and any(
        v is not None for v in (dp, head_end_group, tail_end_group)
    ):
        raise ValueError("A supplied graph already has fixed size and end groups")
    try:
        if field in FIELDS:
            family = FIELDS[field]
            report["actual_field"] = family
            report["hashes"] = _field_sources(family, None, field_dir)[3]
            from island.forcefields.assignment_core import PolymerSpec, normalize
            from island.forcefields.assignment_diagnostics import preflight_assignment

            graph = (
                normalize(
                    PolymerSpec(
                        (("A", structure),),
                        dp,
                        end_groups=(head_end_group, tail_end_group),
                    )
                )
                if isinstance(structure, str)
                else structure
            )
            report["preflight"] = preflight_assignment(
                graph, field, field_dir=field_dir
            )
            if emc_mode == "graph":
                if stereo_policy != "strict":
                    raise ValueError(
                        "Complete graphs preserve stereo; select emc_mode='repeat' for polyjarvis relaxation"
                    )
                result = assign_polymer_parameters(
                    graph, field, field_dir=field_dir, cache_dir=cache_dir
                )
            else:
                if cache_dir is not None:
                    raise ValueError(
                        "Assignment caching applies only to complete graphs"
                    )
                result = assign_polymer_parameters(
                    structure,
                    field,
                    dp=dp,
                    field_dir=field_dir,
                    emc_mode=emc_mode,
                    stereo_policy=stereo_policy,
                    head_end_group=head_end_group,
                    tail_end_group=tail_end_group,
                )
            report.update(
                engine=result.emc_version,
                actual_field=result.field,
                bond_stereo_omitted=getattr(result, "bond_stereo_omitted", False),
            )
        else:
            if (
                field_dir is not None
                or stereo_policy != "strict"
                or emc_mode != "graph"
                or cache_dir is not None
            ):
                raise ValueError("field_dir and stereo relaxation apply to EMC repeats")
            from radonpy.ff import gaff

            db = Path(gaff.__file__).parent / "ff_dat" / (field + ".json")
            report["hashes"] = {db.name: sha256(db.read_bytes()).hexdigest()}
            result = (
                assign_polymer_parameters(
                    structure,
                    field,
                    dp=dp,
                    atom_order="radonpy"
                    if head_end_group is None and tail_end_group is None
                    else "island",
                    charge="gasteiger",
                    head_end_group=head_end_group,
                    tail_end_group=tail_end_group,
                )
                if isinstance(structure, str)
                else assign_radonpy(structure, field, charge="gasteiger")
            )
            report.update(engine=result.version, charge=result.charge_method)
    except (ImportError, OSError) as exc:
        report.update(assignment="backend_unavailable", diagnostic=str(exc))
    except (IslandError, ValueError) as exc:
        from island.forcefields.assignment_diagnostics import assignment_failure

        report["failure"] = assignment_failure(exc, locals().get("graph"))
        msg = str(exc)
        report.update(
            assignment="failed",
            diagnostic=msg,
            source_sites=getattr(exc, "source_sites", []),
            emc_sites=re.findall(r"\{group, site\} = \{(\w+), (\d+)\}", msg),
            missing_increments=re.findall(r"increment pair \{([^}]+)\}", msg),
        )
    else:
        report["assignment"] = "passed"
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    b = sub.add_parser("baseline")
    b.add_argument("evidence", type=Path)
    q = sub.add_parser("probe")
    q.add_argument("field")
    group = q.add_mutually_exclusive_group(required=True)
    group.add_argument("--repeat")
    group.add_argument(
        "--mol", type=Path, help="Explicit molecular graph in MDL mol format"
    )
    q.add_argument("--dp", type=int)
    q.add_argument("--field-dir", type=Path)
    q.add_argument("--emc-mode", choices=("graph", "repeat"), default="graph")
    q.add_argument("--cache-dir", type=Path)
    q.add_argument("--head-end-group")
    q.add_argument("--tail-end-group")
    q.add_argument(
        "--stereo-policy", choices=("strict", "polyjarvis"), default="strict"
    )
    a = parser.parse_args()
    if a.command == "baseline":
        result = benchmark(a.evidence)
    else:
        graph = a.repeat
        if a.mol:
            from rdkit import Chem

            graph = Chem.MolFromMolFile(str(a.mol), removeHs=False)
            if graph is None:
                parser.error("Invalid mol file")
        result = probe(
            graph,
            a.field,
            dp=a.dp,
            stereo_policy=a.stereo_policy,
            field_dir=a.field_dir,
            emc_mode=a.emc_mode,
            cache_dir=a.cache_dir,
            head_end_group=a.head_end_group,
            tail_end_group=a.tail_end_group,
        )
        if a.mol:
            result["input_sha256"] = sha256(a.mol.read_bytes()).hexdigest()
    print(json.dumps(result, indent=2))
    return int(result.get("assignment", "passed") != "passed")


if __name__ == "__main__":
    raise SystemExit(main())
