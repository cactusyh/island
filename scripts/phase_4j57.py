"""Reproduce J57 attestations and actual QM input geometries; no QM or MD run."""

import argparse
import io
import json
import lzma
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

from island.pcff_completion import (
    canonical,
    digest,
    inventory,
    pin,
    require,
    signed,
    write,
)
from island.pcff_completion.qm import ADAPTERS
from island.pcff_qm.__main__ import load_cases
from island.pcff_qm.geometry import geometries
from island.pcff_qm.inputs import emit, generate
from island.pcff_qm.readiness import closure
from island.pcff_qm.records import (
    BASE,
    BRANCH,
    PREREQUISITE,
    Record,
    confirmation_inventory,
)
from island.pcff_qm.reference import fit_readiness
from island.pcff_qm.schema import (
    AUTHOR_ARTIFACT_SCHEMA,
    FITTED_PARAMETER_SCHEMA,
    REFERENCE_SCHEMA,
    TARGET_SCHEMA,
)
from island.pcff_qm.supporting import supporting_scans
from island.pcff_qm.targets import FIXTURES
from island.pcff_qm.targets import inventory as target_inventory


def reconstruction(cases):
    """Separate-process signature binds the coordinates of every fixture."""
    return {
        "inventory": target_inventory(cases).payload,
        "fixture_geometry_checksums": {
            name: digest(
                geometries(
                    *cases[name],
                    supporting_scans(*cases[name], name)[0].payload["scan"],
                )
            )
            for name in FIXTURES
        },
    }


def archive(directory, output):
    """Deterministic regular-file-only POSIX archive, independent of host paths."""
    root = Path(directory)
    raw = io.BytesIO()
    with tarfile.open(fileobj=raw, mode="w", format=tarfile.USTAR_FORMAT) as tar:
        for path in sorted(root.rglob("*")):
            if path.is_file():
                data = path.read_bytes()
                info = tarfile.TarInfo(path.relative_to(root).as_posix())
                info.size, info.mode, info.mtime = len(data), 0o644, 0
                tar.addfile(info, io.BytesIO(data))
    Path(output).write_bytes(lzma.compress(raw.getvalue(), preset=6))


def audit(repo, distribution, output, *, fresh=True):
    repo, distribution, out = Path(repo), Path(distribution), Path(output)
    require(not out.exists(), "output_directory_already_exists")
    out.mkdir(parents=True)
    cases = load_cases(repo, distribution)
    chosen = target_inventory(cases)
    write(out / "author-attested-confirmations.json", confirmation_inventory())
    write(out / "qm-target-inventory.json", chosen.payload)
    for name, schema in (
        ("reference-data.schema.json", REFERENCE_SCHEMA),
        ("qm-target.schema.json", TARGET_SCHEMA),
        ("original-author-artifact.schema.json", AUTHOR_ARTIFACT_SCHEMA),
        ("fitted-parameter.schema.json", FITTED_PARAMETER_SCHEMA),
    ):
        write(out / name, schema)
    targets = [r["target"] for r in chosen.payload["targets"]]
    write(
        out / "fitting-readiness.json",
        signed(
            {
                "targets": [
                    fit_readiness(Record(canonical(t))).payload for t in targets
                ],
                "datasets": [],
                "fitted_parameters": [],
                "all_terms_status": "qm_fit_required",
            }
        ),
    )
    bundle_records, supporting = [], []
    with tempfile.TemporaryDirectory(prefix="island-j57-geometries-") as tmp:
        root = Path(tmp) / "bundle"
        for i, t in enumerate(targets):
            manifest = generate(
                *cases[t["molecular_fixture"]],
                Record(canonical(t)),
                root / f"targets/target-{i:02d}",
            )
            bundle_records.append(
                {
                    "query": t["query"],
                    "directory": f"targets/target-{i:02d}",
                    "target_identity": t["identity"],
                    "manifest": manifest,
                }
            )
        for fixture in FIXTURES:
            s, a = cases[fixture]
            write(out / f"{fixture}-closure.json", closure(s, a, fixture).payload)
            for plan in supporting_scans(s, a, fixture):
                p = plan.payload
                directory = f"supporting/{fixture}/{p['label']}"
                manifest = emit(s, a, p, root / directory)
                supporting.append(
                    {
                        "fixture": fixture,
                        "directory": directory,
                        "metadata": p,
                        "manifest": manifest,
                    }
                )
        files = inventory(root)
        write(
            out / "geometry-input-inventory.json",
            signed(
                {
                    "files": files,
                    "targets": bundle_records,
                    "supporting_scans": supporting,
                    "QM_outputs_present": False,
                }
            ),
        )
        archive(root, out / "geometry-input-bundle.tar.xz")
    historical_queries = json.loads(
        (
            repo / "docs/evidence/phase_4j56/audit/47-query-classifications.json"
        ).read_text()
    )["queries"]
    confirmed = {r["query"] for r in confirmation_inventory()["confirmations"]}
    require(
        confirmed
        == {
            q["query"]
            for q in historical_queries
            if q["classification"] == "author_confirmed"
        },
        "historical_confirmed_query_membership_mutation",
    )
    rows = [
        {
            **q,
            "historical_status": q["classification"],
            "classification": "author_confirmed_audited"
            if q["query"] in confirmed
            else "qm_fit_required",
            "status": "author_confirmed_audited"
            if q["query"] in confirmed
            else "qm_fit_required",
        }
        for q in historical_queries
    ]
    write(
        out / "47-query-classification.json",
        signed(
            {
                "queries": rows,
                "historical_J56_sha256": pin(
                    repo
                    / "docs/evidence/phase_4j56/audit/47-query-classifications.json"
                ),
            }
        ),
    )
    write(
        out / "provenance-and-uncertainty.json",
        signed(
            {
                "source_kind": "user-attested external authority plus independently pinned historical source files",
                "original_author_artifacts": [],
                "original_author_artifact_status": "not_yet_supplied",
                "reference_datasets": [],
                "fitted_parameters": [],
                "uncertainty_values": [],
                "original_artifact_upgrade": "additive parent-attestation-bound evidence; no overwrite",
                "target_protocol_authority": "proposed J57 preparation policy; no source authority or numerical adequacy claim",
                "future_parameter_origin": "fitted_external",
                "all_nineteen_numerical_terms": "qm_fit_required",
                "source_inputs": {
                    str(p.relative_to(repo)): pin(p)
                    for phase in (54, 55, 56)
                    for p in sorted(
                        (repo / f"docs/evidence/phase_4j{phase}/audit").rglob("*.json")
                    )
                },
            }
        ),
    )
    fresh_records = []
    if fresh:
        with tempfile.TemporaryDirectory(prefix="island-j57-relocated-") as tmp:
            relocated = Path(tmp)
            for fixture in FIXTURES:
                src = repo / "docs/evidence/phase_4j54/audit" / fixture
                shutil.copytree(
                    src, relocated / "repo/docs/evidence/phase_4j54/audit" / fixture
                )
            shutil.copytree(distribution, relocated / "sources")
            expected_reconstruction = reconstruction(cases)
            script = (
                "from pathlib import Path; import sys,json; from island.pcff_qm.__main__ import load_cases; "
                "from scripts.phase_4j57 import reconstruction; "
                "print(json.dumps(reconstruction(load_cases(Path(sys.argv[1]), Path(sys.argv[2])))))"
            )
            for seed in (1, 99397):
                run = subprocess.run(
                    [
                        sys.executable,
                        "-c",
                        script,
                        str(relocated / "repo"),
                        str(relocated / "sources"),
                    ],
                    env={**os.environ, "PYTHONHASHSEED": str(seed)},
                    capture_output=True,
                    text=True,
                    check=False,
                )
                require(
                    run.returncode == 0
                    and json.loads(run.stdout) == expected_reconstruction,
                    "relocation_or_fresh_process_mutation",
                    stderr=run.stderr,
                )
                fresh_records.append(
                    {
                        "hash_seed": seed,
                        "returncode": run.returncode,
                        "inventory_identity": chosen.payload["identity"],
                        "relocated_sources": True,
                        "all_five_fixture_graphs_reconstructed": True,
                        "fixture_geometry_checksums": expected_reconstruction[
                            "fixture_geometry_checksums"
                        ],
                    }
                )
    summary = signed(
        {
            "schema": "island.j57-audit.v1",
            "base_commit": BASE,
            "prerequisite_snapshot": PREREQUISITE,
            "branch": BRANCH,
            "base_contains_merged_J56": False,
            "author_confirmation_status": "author_confirmed_audited",
            "author_confirmed_queries": 28,
            "confirmation_authority": "user_reported_author_confirmation",
            "primary_author_artifact_status": "not_yet_supplied",
            "qm_fit_required_queries": 19,
            "directed_source_queries": 47,
            "historical_missing_site_requests": sum(q["requests"] for q in rows),
            "missing_numerical_site_requests": sum(
                q["requests"] for q in rows if q["status"] == "qm_fit_required"
            ),
            "author_attested_site_requests": sum(
                q["requests"] for q in rows if q["status"] == "author_confirmed_audited"
            ),
            "fixture_graphs": {
                name: a.payload["graph_identity"] for name, (_, a) in cases.items()
            },
            "source_hashes": cases["hmds"][1].payload["source_hashes"],
            "target_inventory_identity": chosen.payload["identity"],
            "numerical_target_samples": sum(
                r["manifest"]["sample_count"] for r in bundle_records
            ),
            "supporting_scan_count": len(supporting),
            "supporting_samples": sum(
                r["manifest"]["sample_count"] for r in supporting
            ),
            "adapters": list(ADAPTERS),
            "geometry_bundle_sha256": pin(out / "geometry-input-bundle.tar.xz"),
            "fresh_processes": fresh_records,
            "registered_profiles": [],
            "created_terms": [],
            "reference_dataset_status": "not_supplied",
            "parameter_completeness": "incomplete",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "QM_executed": False,
            "MD_executed": False,
            "NPT_executed": False,
        }
    )
    write(out / "summary.json", summary)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--distribution", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--no-fresh", action="store_true")
    a = p.parse_args()
    print(
        json.dumps(
            audit(a.repo, a.distribution, a.output, fresh=not a.no_fresh),
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
