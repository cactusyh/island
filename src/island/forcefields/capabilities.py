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
                for pair in set(
                    re.findall(r"increment pair \{([^}]+)\}", r["diagnostic"])
                ):
                    gaps[field][",".join(sorted(pair.split(", ")))] += 1
    fixed = json.loads((root / "phase_4j94/summary.json").read_text())
    cure = json.loads((root / "phase_4j95_j97/cure-study.json").read_text())
    return {
        "panel_sha256": sha256(panel.read_bytes()).hexdigest(),
        "linear": dict(fields),
        "published_positive": dict(published),
        "missing_increments": dict(gaps),
        "fixed_graph_md": fixed,
        "periodic_epoxy": [
            {k: c[k] for k in ("dgeba", "pacm", "seed", "statistics")}
            for c in cure["cases"]
        ],
    }


def probe(structure, field, *, dp=None, stereo_policy="strict", field_dir=None):
    """Accept a repeat string or an already built graph through its native backend."""
    report = {
        "field": field,
        "input_kind": "repeat" if isinstance(structure, str) else "graph",
        "psmiles": structure if isinstance(structure, str) else None,
        "dp": dp,
        "stereo_policy": stereo_policy,
        "hashes": {},
        "assignment": "unattempted",
        "single_point": "unmeasured",
        "dynamics": "unmeasured",
        "properties": "unmeasured",
    }
    if field not in {*FIELDS, *FAMILIES}:
        raise ValueError(f"Unknown field: {field}")
    if not isinstance(structure, str) and dp is not None:
        raise ValueError("A supplied graph already has a fixed size")
    try:
        if field in FIELDS:
            family = FIELDS[field]
            report["hashes"] = _field_sources(family, None, field_dir)[3]
            options = (
                {"stereo_policy": stereo_policy} if isinstance(structure, str) else {}
            )
            if not isinstance(structure, str) and stereo_policy != "strict":
                raise ValueError("Built graphs preserve the supplied stereo policy")
            result = assign_polymer_parameters(
                structure, field, dp=dp, field_dir=field_dir, **options
            )
            report.update(
                engine=result.emc_version,
                actual_field=result.field,
                bond_stereo_omitted=getattr(result, "bond_stereo_omitted", False),
            )
        else:
            if field_dir is not None or stereo_policy != "strict":
                raise ValueError("field_dir and stereo relaxation apply to EMC repeats")
            from radonpy.ff import gaff

            db = Path(gaff.__file__).parent / "ff_dat" / (field + ".json")
            report["hashes"] = {db.name: sha256(db.read_bytes()).hexdigest()}
            result = (
                assign_polymer_parameters(
                    structure, field, dp=dp, atom_order="radonpy", charge="gasteiger"
                )
                if isinstance(structure, str)
                else assign_radonpy(structure, field, charge="gasteiger")
            )
            report.update(engine=result.version, charge=result.charge_method)
    except (ImportError, OSError) as exc:
        report.update(assignment="backend_unavailable", diagnostic=str(exc))
    except (IslandError, ValueError) as exc:
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
        )
        if a.mol:
            result["input_sha256"] = sha256(a.mol.read_bytes()).hexdigest()
    print(json.dumps(result, indent=2))
    return int(result.get("assignment", "passed") != "passed")


if __name__ == "__main__":
    raise SystemExit(main())
