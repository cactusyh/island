"""Eight fixed whole-oligomer AM1-BCC references; no retries or charge fallback.

Generate: --output NEW --amberhome PATH
Read-only reconstruction/audit: --audit-existing MATRIX_ROOT --output NEW
Partial matrices produce diagnostic audits and exit 1, never successful coverage.
"""

import argparse
import os
import time
from math import isfinite
from pathlib import Path

from island.charge_references import (
    audit_charge_references,
    create_charge_reference,
    repeat_correspondence,
    save_record,
)
from island.workflows import storage
from island.workflows.bundle import preparation_from, system_data, system_from

MATRIX = [
    {"chemistry": c, "psmiles": p, "dp": d, "seed": s}
    for c, p in [("pe", "[*:1]CC[*:2]"), ("peo", "[*:1]CCO[*:2]")]
    for d, s in [(3, 2026), (5, 2026), (7, 2026), (5, 80317)]
]
SETTINGS = {
    "matrix": MATRIX,
    "timeout_seconds_per_stage": 600,
    "retries": 0,
    "charge_tolerance_e": 0.002,
    "method": "GAFF2 whole-oligomer AM1-BCC",
    "max_atoms": 100,
}


def reconstruct(case):
    system = system_from(storage.decode(storage.read_json(case / "input-system.json")))
    p = storage.decode(storage.read_json(case / "preparation.json"))
    r = p["record"]
    artifacts = storage.child(case, Path(r["artifact_dir"]).name)
    checks = {
        **r["artifact_sha256"],
        "input.mol2": r["input_mol2_sha256"],
        "lineage.json": r["input_lineage_sha256"],
        "sqm.out": r["charge_outcome"]["sqm_out_sha256"],
    }
    for name, expected in checks.items():
        if storage.checksum(storage.child(artifacts, name).read_bytes()) != expected:
            raise ValueError(f"Artifact checksum mismatch: {name}")
    return system, preparation_from(system, p, artifacts / "result.prmtop")


def generate(root, amberhome):
    from island.builders import build_linear_polymer
    from island.forcefields import AmberToolsOptions, AmberToolsParameterizationEngine

    storage.publish(root / "declared-settings.json", storage.json_bytes(SETTINGS))
    rows = []
    original_path = os.environ.get("PATH", "")
    if amberhome:
        os.environ["PATH"] = str(amberhome / "bin") + os.pathsep + original_path
    try:
        for c in MATRIX:
            name = f"{c['chemistry']}-dp{c['dp']}-seed{c['seed']}"
            row = dict(c, case=name, status="failed")
            case = root / name
            case.mkdir()
            start = time.perf_counter()
            try:
                system = build_linear_polymer(
                    c["psmiles"],
                    dp=c["dp"],
                    coordinate_method="local_templates",
                    template_seed=c["seed"],
                    assembly_seed=c["seed"],
                )
                mapping = repeat_correspondence(system)
                if not mapping["compatible"]:
                    raise ValueError(mapping["diagnostics"])
                row["sites"] = system.number_of_sites
                storage.publish(
                    case / "input-system.json",
                    storage.json_bytes(storage.encode(system_data(system))),
                )
                prep = AmberToolsParameterizationEngine().parameterize(
                    system,
                    AmberToolsOptions(
                        "gaff2",
                        "am1bcc",
                        amberhome=amberhome,
                        work_root=case,
                        retain_success_artifacts=True,
                        timeout_seconds=600,
                        charge_tolerance=0.002,
                        max_atoms=100,
                    ),
                )
                storage.publish(
                    case / "preparation.json",
                    storage.json_bytes(
                        storage.encode(
                            {
                                "record": dict(prep.record),
                                "record_signature": prep.record_signature,
                                "import_source": prep.imported_result.source,
                                "import_provenance": dict(
                                    prep.imported_result.provenance
                                ),
                            }
                        )
                    ),
                )
                row.update(
                    status="passed",
                    record_signature=prep.record_signature,
                    import_signature=prep.imported_result.result_signature,
                    residual=prep.imported_result.charge_result.total_charge_residual,
                )
            except Exception as error:  # noqa: BLE001 -- preserve failed actual calculations
                row["failure"] = f"{type(error).__name__}: {error}"
            row["elapsed_seconds"] = time.perf_counter() - start
            rows.append(row)
            storage.publish(
                root / "matrix.json",
                storage.json_bytes(
                    {
                        "cases": rows,
                        "complete": len(rows) == 8
                        and all(r["status"] == "passed" for r in rows),
                    }
                ),
                replace=True,
            )
            print(name, row["status"], flush=True)
    finally:
        os.environ["PATH"] = original_path
    return rows


def audit_existing(source, output):
    if storage.read_json(source / "declared-settings.json") != SETTINGS:
        raise ValueError("Declaration does not match this fixed acceptance matrix")
    matrix = storage.read_json(source / "matrix.json")
    if (
        not isinstance(matrix, dict)
        or set(matrix) != {"cases", "complete"}
        or type(matrix["complete"]) is not bool
    ):
        raise ValueError("Invalid matrix summary structure")
    rows = matrix["cases"]
    if not isinstance(rows, list) or len(rows) > len(MATRIX):
        raise ValueError("Invalid matrix row container/count")
    if matrix["complete"] != (
        len(rows) == 8
        and all(isinstance(r, dict) and r.get("status") == "passed" for r in rows)
    ):
        raise ValueError("Matrix completion summary contradicts rows")
    refs = []
    outcomes = []
    seen = set()
    for row in rows:
        if not isinstance(row, dict) or not all(
            type(row.get(k)) is str for k in ("case", "chemistry", "psmiles", "status")
        ):
            raise ValueError("Invalid matrix row text fields")
        if any(type(row.get(k)) is not int for k in ("dp", "seed")) or row[
            "status"
        ] not in {"passed", "failed"}:
            raise ValueError("Invalid matrix row types/status")
        if "sites" in row and (type(row["sites"]) is not int or row["sites"] <= 0):
            raise ValueError("Invalid reported site count")
        for key in ("residual", "elapsed_seconds"):
            if key in row and (
                type(row[key]) not in (int, float)
                or not isfinite(row[key])
                or (key == "elapsed_seconds" and row[key] < 0)
            ):
                raise ValueError(f"Invalid matrix {key}")
        if (
            row["status"] == "passed"
            and not {
                "sites",
                "residual",
                "record_signature",
                "import_signature",
                "elapsed_seconds",
            }
            <= row.keys()
        ):
            raise ValueError("Successful case summary is incomplete")
        if row["status"] == "failed" and (
            type(row.get("failure")) is not str or not row["failure"]
        ):
            raise ValueError("Failed case needs a diagnostic")
        if row["status"] == "passed":
            if "failure" in row:
                raise ValueError("Passed case contradicts failure diagnostic")
            for key in ("record_signature", "import_signature"):
                value = row[key]
                if (
                    type(value) is not str
                    or len(value) != 64
                    or any(c not in "0123456789abcdef" for c in value)
                ):
                    raise ValueError("Invalid reported signature")
        name = row["case"]
        declaration = {k: row[k] for k in ("chemistry", "psmiles", "dp", "seed")}
        if (
            declaration not in MATRIX
            or name in seen
            or name != f"{row['chemistry']}-dp{row['dp']}-seed{row['seed']}"
        ):
            raise ValueError("Unexpected/duplicate matrix case")
        seen.add(name)
        result = dict(row)
        case_path = storage.child(source, name)
        result["source_checksums"] = {
            str(path.relative_to(case_path)): storage.checksum(path.read_bytes())
            for path in sorted(case_path.rglob("*"))
            if path.is_file() and not path.is_symlink()
        }
        if row["status"] == "passed":
            try:
                system, prep = reconstruct(storage.child(source, name))
                if (
                    prep.record_signature != row["record_signature"]
                    or prep.imported_result.result_signature != row["import_signature"]
                ):
                    raise ValueError("Matrix signatures disagree")
                ref = create_charge_reference(system, prep)
                if (
                    ref.payload["correspondence"]["dp"] != row["dp"]
                    or ref.payload["correspondence"]["definition"] != row["psmiles"]
                    or ref.payload["seeds"]
                    != {"template_seed": row["seed"], "assembly_seed": row["seed"]}
                ):
                    raise ValueError("Declared input differs from saved source")
                payload = ref.payload
                r = payload["preparation"]["record"]
                charge = payload["charge_result"]
                if (
                    r["charge_validation_tolerance_e"] != SETTINGS["charge_tolerance_e"]
                    or charge["tolerance"] != SETTINGS["charge_tolerance_e"]
                ):
                    raise ValueError(
                        "Actual charge tolerance differs from declared experiment"
                    )
                if (
                    r["requested_force_field"] != "gaff2"
                    or r["charge_method"] != "am1bcc"
                ):
                    raise ValueError(
                        "Actual calculation method differs from declaration"
                    )
                actual_limit = r.get("size_policy", {}).get("max_atoms", 100)
                if (
                    actual_limit != SETTINGS["max_atoms"]
                    or system.number_of_sites > actual_limit
                ):
                    raise ValueError("Actual size policy differs from declaration")
                if (
                    row["sites"] != system.number_of_sites
                    or row["residual"] != charge["total_charge_residual"]
                ):
                    raise ValueError(
                        "Reported site count/residual differs from authoritative record"
                    )
                result["sites"] = system.number_of_sites
                result["residual"] = charge["total_charge_residual"]
                result["verified_settings"] = {
                    "charge_tolerance_e": charge["tolerance"],
                    "force_field": r["requested_force_field"],
                    "charge_method": r["charge_method"],
                    "max_atoms": actual_limit,
                }
                result["declared_only_settings"] = {
                    "timeout_seconds_per_stage": SETTINGS["timeout_seconds_per_stage"],
                    "retries": SETTINGS["retries"],
                }
                save_record(ref, output / (name + ".json"))
                refs.append(ref)
                result["reference_identity"] = ref.identity
            except Exception as error:  # noqa: BLE001 -- durable validation diagnostic
                result["status"] = "reference_validation_failed"
                result["failure"] = str(error)
        outcomes.append(result)
    if refs:
        audit = audit_charge_references(refs)
        save_record(audit, output / "audit.json")
    complete = len(outcomes) == 8 and all(r["status"] == "passed" for r in outcomes)
    storage.publish(
        output / "outcomes.json",
        storage.json_bytes(
            {
                "schema": "island_charge_matrix_outcomes_v1",
                "cases": outcomes,
                "complete": complete,
                "validated_reference_count": len(refs),
                "missing_cases": 8 - len(outcomes),
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
    )
    return complete


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--amberhome", type=Path)
    p.add_argument("--audit-existing", type=Path)
    args = p.parse_args(argv)
    if args.audit_existing and (
        args.output.resolve() == args.audit_existing.resolve()
        or args.audit_existing.resolve() in args.output.resolve().parents
        or args.output.resolve() in args.audit_existing.resolve().parents
    ):
        p.error("Audit output must be separate from source")
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        if args.audit_existing:
            complete = audit_existing(args.audit_existing, args.output)
        else:
            generate(args.output, args.amberhome)
            audit_dir = args.output / "audit"
            audit_dir.mkdir()
            complete = audit_existing(args.output, audit_dir)
        return 0 if complete else 1
    except Exception as error:  # noqa: BLE001 -- command errors also fail the gate
        storage.publish(
            args.output / "failure.json",
            storage.json_bytes({"failure": f"{type(error).__name__}: {error}"}),
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
