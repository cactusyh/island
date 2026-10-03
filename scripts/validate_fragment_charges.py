"""Declared generic fragment experiment: one AM1-BCC call per repeat, no retries."""

import argparse
import os
import time
from pathlib import Path

from island.fragment_charges import (
    AM1BCCFragmentChargeBackend,
    assign_fragment_charges,
    create_fragment_template,
    load_fragment_record,
    prepare_capped_fragment,
    save_fragment_record,
)
from island.workflows import storage
from island.workflows.bundle import system_data

PLAN = {
    "schema": "island_fragment_experiment_v1",
    "matrix": [
        ("pe", "[*:1]CC[*:2]"),
        ("peo", "[*:1]CCO[*:2]"),
        ("ps", "[*:1]CC(c1ccccc1)[*:2]"),
        ("pmma", "[*:1]CC(C)(C(=O)OC)[*:2]"),
    ],
    "fragment_seed": 2026,
    "force_field": "gaff2",
    "method": "am1bcc",
    "timeout_seconds_per_stage": 600,
    "max_fragment_atoms": 100,
    "outer_retries": 0,
    "conservation_policy": "uniform_fragment_l2_v1",
    "raw_strict_tolerance_e": 0.002,
    "hydrogen_policy": "verified_parent_group_equal_v1",
    "target_tolerance_e": 1e-12,
    "mapping_dps": [1, 2, 3, 50],
    "integration": {
        "case": "pe",
        "dp": 3,
        "force_field": "gaff2",
        "method": "provided",
        "template_seed": 2026,
        "assembly_seed": 2026,
        "charge_tolerance_e": 1e-12,
        "timeout_seconds_per_stage": 600,
    },
    "production_validated": False,
    "simulation_readiness": "not_established",
}


def execute(root, amberhome):
    from island.builders import build_linear_polymer
    from island.forcefields import AmberToolsOptions, AmberToolsParameterizationEngine

    rows = []
    for name, psmiles in PLAN["matrix"]:
        case = root / name
        case.mkdir()
        row = {"case": name, "definition": psmiles, "status": "failed", "mapping": []}
        start = time.perf_counter()
        try:
            fragment = prepare_capped_fragment(psmiles, seed=2026)
            save_fragment_record(fragment, case / "fragment.json")
            row["fragment_sites"] = fragment.payload["data"]["atom_count"]
            # Exactly one charge-backend invocation, always on the small fragment.
            result = AM1BCCFragmentChargeBackend(
                force_field="gaff2",
                work_root=case / "calculation",
                amberhome=amberhome,
                timeout_seconds=600,
            ).calculate(fragment)
            save_fragment_record(result, case / "charges.json")
            row["raw_residual_e"] = result.payload["data"]["raw_residual_e"]
            row["raw_strict_gate"] = abs(row["raw_residual_e"]) <= 0.002
            row["backend_invocations"] = 1
            template = create_fragment_template(
                result,
                conservation_policy=PLAN["conservation_policy"],
                charge_tolerance=0.002,
            )
            save_fragment_record(template, case / "template.json")
            key = template.cache_key
            row["template_identity"] = template.identity
            row["cache_key"] = key
            row["force_field"] = result.payload["backend"]["force_field"]
            for dp in PLAN["mapping_dps"]:
                cached = load_fragment_record(
                    case / "template.json", expected_cache_key=key
                )
                system = build_linear_polymer(psmiles, dp=dp, generate_3d=False)
                assignment = assign_fragment_charges(
                    cached,
                    system,
                    force_field=row["force_field"],
                    charge_tolerance=1e-12,
                )
                save_fragment_record(assignment, case / f"assignment-dp{dp}.json")
                row["mapping"].append(
                    {
                        "dp": dp,
                        "sites": system.number_of_sites,
                        "total_e": assignment.payload["data"]["total_e"],
                        "assignment_identity": assignment.identity,
                    }
                )
            row["status"] = "passed"
            if name == "pe":
                integration = case / "integration"
                integration.mkdir()
                try:
                    system = build_linear_polymer(
                        psmiles,
                        dp=3,
                        coordinate_method="local_templates",
                        template_seed=2026,
                        assembly_seed=2026,
                    )
                    assignment = assign_fragment_charges(
                        template, system, force_field=row["force_field"]
                    )
                    save_fragment_record(assignment, integration / "assignment.json")
                    prep = AmberToolsParameterizationEngine().parameterize(
                        system,
                        AmberToolsOptions(
                            "gaff2",
                            "provided",
                            provided_charges=assignment.payload["data"]["assignments"],
                            charge_source="fragment-derived; AM1-BCC; uniform_fragment_l2_v1; removed_cap_to_bonded_parent_v1; verified_parent_group_equal_v1; assignment="
                            + assignment.identity,
                            charge_tolerance=1e-12,
                            timeout_seconds=600,
                            amberhome=amberhome,
                            work_root=integration,
                            retain_success_artifacts=True,
                        ),
                    )
                    prep.validate_integrity(system)
                    if (
                        prep.record["force_field_data"]["sha256"]
                        != row["force_field"]["data_sha256"]
                    ):
                        raise ValueError("Integration FF data version changed")
                    storage.publish(
                        integration / "preparation.json",
                        storage.json_bytes(
                            storage.encode(
                                {
                                    "system": system_data(system),
                                    "record": dict(prep.record),
                                    "record_signature": prep.record_signature,
                                }
                            )
                        ),
                    )
                    row["integration"] = {
                        "status": "passed",
                        "record_signature": prep.record_signature,
                        "charge_outcome": dict(prep.record["charge_outcome"]),
                    }
                except Exception as error:  # noqa: BLE001 -- independent integration gate
                    row["integration"] = {
                        "status": "failed",
                        "failure": f"{type(error).__name__}: {error}",
                    }
        except Exception as error:  # noqa: BLE001 -- no seed/method replacement
            row["failure"] = f"{type(error).__name__}: {error}"
        row["wall_seconds"] = time.perf_counter() - start
        row["file_sha256"] = {
            str(p.relative_to(case)): storage.checksum(p.read_bytes())
            for p in sorted(case.rglob("*"))
            if p.is_file()
        }
        rows.append(row)
        storage.publish(
            root / "outcomes.json", storage.json_bytes(summarize(rows)), replace=True
        )
        print(
            name,
            row["status"],
            row.get("failure", ""),
            row.get("integration"),
            flush=True,
        )
    return summarize(rows)["complete"]


def summarize(rows):
    mapping = len(rows) == 4 and all(r["status"] == "passed" for r in rows)
    integration = (
        bool(rows) and rows[0].get("integration", {}).get("status") == "passed"
    )
    return {
        "schema": "island_fragment_outcomes_v1",
        "cases": rows,
        "complete_fragment_mapping": mapping,
        "integration_passed": integration,
        "complete": mapping and integration,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--amberhome", type=Path)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--declare-only", action="store_true")
    group.add_argument("--execute-declared", action="store_true")
    args = parser.parse_args(argv)
    started = False
    old = os.environ.get("PATH", "")
    try:
        if args.execute_declared:
            if storage.decode(
                storage.read_json(args.output / "declaration.json")
            ) != PLAN or {p.name for p in args.output.iterdir()} != {
                "declaration.json"
            }:
                raise ValueError("Require fresh unchanged declaration; no retries")
        else:
            args.output.mkdir(parents=True, exist_ok=False)
            storage.publish(
                args.output / "declaration.json",
                storage.json_bytes(storage.encode(PLAN)),
            )
        if args.declare_only:
            return 0
        if args.amberhome:
            os.environ["PATH"] = (
                str(args.amberhome.resolve() / "bin") + os.pathsep + old
            )
        with storage.writer(args.output):
            started = True
            return 0 if execute(args.output, args.amberhome) else 1
    except Exception as error:  # noqa: BLE001 -- explicit unsuccessful CLI gate
        if started and not (args.output / "failure.json").exists():
            storage.publish(
                args.output / "failure.json",
                storage.json_bytes({"failure": str(error)}),
            )
        print(f"{type(error).__name__}: {error}")
        return 1
    finally:
        os.environ["PATH"] = old


if __name__ == "__main__":
    raise SystemExit(main())
