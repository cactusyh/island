"""Explicit eight-case retained-artifact experiment; no QM or parameterization.

The declaration is written before reading scientific artifacts. A separate
--declare-only operation supports inspecting it before --execute-declared.
Exit zero means eight diagnostic projections/comparisons, not scientific acceptance.
"""

import argparse
from pathlib import Path

from island.charge_references import (
    audit_charge_projections,
    observe_retained_case,
    project_charge_observation,
    save_record,
)
from island.charge_references.projection import CONSERVATION_TOLERANCE_E, POLICY
from island.workflows import storage
from scripts.validate_oligomer_charges import MATRIX


def declaration(manifest):
    return {
        "schema": "island_conservation_experiment_v1",
        "cases": MATRIX,
        "raw_stage": "typed.mol2 post-BCC",
        "policy": POLICY,
        "conservation_tolerance_e": CONSERVATION_TOLERANCE_E,
        "target_dps": [20, 50, 100],
        "comparisons": [
            "baseline DP3/5/7 pairwise",
            "DP5 both seeds",
            "raw versus projected",
            "post-SQM geometry identity and pairwise-distance RMS",
        ],
        "source_manifest_sha256": storage.checksum(manifest.read_bytes()),
        "required_cases": 8,
        "retry_policy": "none; no QM/tool execution",
        "outputs": [
            "observations/<case>.json",
            "projections/<case>.json",
            "audit.json",
            "outcomes.json",
        ],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def execute(source, manifest_path, evidence_path, output):
    manifest = storage.read_json(manifest_path)
    evidence = storage.read_json(evidence_path)
    if evidence["original_manifest_sha256"] != storage.checksum(
        manifest_path.read_bytes()
    ):
        raise ValueError("Reviewed evidence/acceptance manifest disagreement")
    rows = manifest["cases"]
    if (
        len(rows) != 8
        or [{k: r[k] for k in ("chemistry", "psmiles", "dp", "seed")} for r in rows]
        != MATRIX
    ):
        raise ValueError("Original matrix differs from declaration")
    if evidence["cases"] and {
        r["case"]: r["source_checksums"] for r in evidence["cases"]
    } != {r["case"]: r["source_checksums"] for r in rows}:
        raise ValueError("Reviewed artifact manifests disagree")
    # Verify all original files before processing the first scientific case.
    originals = {}
    for row in rows:
        directory = storage.child(source, row["case"])
        for relative, expected in row["source_checksums"].items():
            path = storage.child(directory, relative)
            if storage.checksum(path.read_bytes()) != expected:
                raise ValueError("Original artifact checksum mismatch")
            originals[path] = expected
    (output / "observations").mkdir()
    (output / "projections").mkdir()
    projections = []
    outcomes = []
    for row in rows:
        name = row["case"]
        outcome = {
            "case": name,
            "historical_status": row["status"],
            "historical_failure": row.get("failure"),
            "status": "failed",
        }
        try:
            observation = observe_retained_case(source, manifest_path, name)
            projection = project_charge_observation(observation, policy=POLICY)
            save_record(observation, output / "observations" / f"{name}.json")
            save_record(projection, output / "projections" / f"{name}.json")
            projections.append(projection)
            outcome.update(
                status="passed",
                observation_identity=observation.identity,
                projection_identity=projection.identity,
                raw_total_e=projection.payload["projection"]["raw_total_e"],
                projected_residual_e=projection.payload["projection"][
                    "projected_residual_e"
                ],
            )
        except Exception as error:  # noqa: BLE001 -- durable diagnostic gate
            outcome["failure"] = f"{type(error).__name__}: {error}"
        outcomes.append(outcome)
    unchanged = all(storage.checksum(p.read_bytes()) == h for p, h in originals.items())
    complete = len(projections) == 8 and unchanged
    if projections:
        audit = audit_charge_projections(projections)
        save_record(audit, output / "audit.json")
        if len(audit.payload["audit"]["raw"]["comparisons"]) != 8:
            complete = False
    storage.publish(
        output / "outcomes.json",
        storage.json_bytes(
            {
                "schema": "island_conservation_outcomes_v1",
                "cases": outcomes,
                "complete_diagnostic_matrix": complete,
                "source_hashes_unchanged": unchanged,
                "checked_source_files": len(originals),
                "historical_valid_references": sum(
                    r["status"] == "passed" for r in rows
                ),
                "original_raw_conformation_acceptance": "unmet",
                "scientific_acceptance": "not_established",
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
    )
    return complete


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument(
        "--manifest", type=Path, default=Path("docs/phase_4f3_acceptance.json")
    )
    parser.add_argument(
        "--evidence", type=Path, default=Path("docs/phase_4f3_1_evidence.json")
    )
    parser.add_argument("--policy", choices=[POLICY], required=True)
    parser.add_argument("--output", type=Path, required=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--declare-only", action="store_true")
    group.add_argument("--execute-declared", action="store_true")
    args = parser.parse_args(argv)
    for source in (args.source,):
        if (
            args.output.resolve() == source.resolve()
            or args.output.resolve() in source.resolve().parents
            or source.resolve() in args.output.resolve().parents
        ):
            parser.error("Output must be separate from historical input directories")
    execution_started = False
    try:
        expected = declaration(args.manifest)
        if args.execute_declared:
            if storage.read_json(args.output / "declaration.json") != expected:
                raise ValueError("Experiment declaration changed")
            if {p.name for p in args.output.iterdir()} != {"declaration.json"}:
                raise ValueError(
                    "Execution requires a fresh declaration-only directory"
                )
        else:
            args.output.mkdir(parents=True, exist_ok=False)
            storage.publish(
                args.output / "declaration.json", storage.json_bytes(expected)
            )
        if args.declare_only:
            return 0
        with storage.writer(args.output):
            execution_started = True
            return (
                0
                if execute(
                    args.source,
                    args.manifest,
                    args.evidence,
                    args.output,
                )
                else 1
            )
    except Exception as error:  # noqa: BLE001 -- no success on missing evidence
        # Never replace an existing diagnostic, declaration or successful output.
        if (
            execution_started
            and args.output.is_dir()
            and not (args.output / "failure.json").exists()
        ):
            storage.publish(
                args.output / "failure.json",
                storage.json_bytes({"failure": f"{type(error).__name__}: {error}"}),
            )
        print(f"{type(error).__name__}: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
