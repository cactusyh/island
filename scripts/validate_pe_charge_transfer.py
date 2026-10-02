"""Freeze PE models first; separately execute four declared held-out calculations.

Exit 0 for --execute requires four verified observation/projected targets and
prediction coverage, not four successful raw imports or scientific accuracy.
No retry, alternate seed, charge-method fallback or model refit is performed.
"""

import argparse
import os
import time
from pathlib import Path

from island.charge_references import (
    create_charge_reference,
    fit_pe_template,
    load_charge_projection,
    load_pe_record,
    observe_retained_case,
    predict_pe_charges,
    project_charge_observation,
    save_record,
    validate_pe_templates,
)
from island.charge_references.pe_template import FLAGS, HELD_OUT, TRAINING
from island.charge_references.projection import POLICY
from island.workflows import storage
from island.workflows.bundle import system_data

PLAN = {
    "schema": "island_pe_transfer_experiment_v1",
    "held_out": [
        {"dp": d, "template_seed": s, "assembly_seed": s, "expected_sites": 6 * d + 2}
        for d, s in HELD_OUT
    ],
    "psmiles": "[*:1]CC[*:2]",
    "force_field": "gaff2",
    "charge_method": "am1bcc",
    "max_atoms": 100,
    "charge_tolerance_e": 0.002,
    "timeout_seconds_per_stage": 600,
    "outer_retries": 0,
    "projection_policy": POLICY,
    "prediction_only_dps": [20, 50, 100],
    "conservation_tolerance_e": 1e-12,
    "comparisons": [
        "raw and projected separately",
        "heavy max/RMSE",
        "parent-H mean/sum",
        "roles and distances from ends",
        "repeat and molecular totals",
        "DP9 final-SQM geometry",
    ],
    "scientific_accuracy_threshold": None,
    **FLAGS,
}


def freeze(training_root, original_source, acceptance, output):
    manifest = storage.read_json(acceptance)
    training = []
    checked = {}
    for dp, seed in TRAINING:
        case = f"pe-dp{dp}-seed{seed}"
        row = next(r for r in manifest["audit"]["cases"] if r["case"] == case)
        for sub in ("observations", "projections"):
            relative = f"{sub}/{case}.json"
            path = training_root / relative
            expected = manifest["output_checksums"][relative]
            if storage.checksum(path.read_bytes()) != expected:
                raise ValueError("Training record checksum mismatch")
            checked[relative] = expected
        projection = load_charge_projection(
            training_root / "projections" / f"{case}.json"
        )
        p = projection.payload
        if (
            projection.identity != row["projection_identity"]
            or p["observation_identity"] != row["observation_identity"]
        ):
            raise ValueError("Training identity mismatch")
        # Check original files read-only as well as embedded validated evidence.
        for relative, expected in p["observation"]["data"]["artifact_sha256"].items():
            if (
                storage.checksum(
                    storage.child(original_source / case, relative).read_bytes()
                )
                != expected
            ):
                raise ValueError("Historical source checksum mismatch")
        training.append(projection)
    models = [fit_pe_template(training, mode=m) for m in ("baseline", "conserving")]
    output.mkdir(parents=True, exist_ok=False)
    for model in models:
        save_record(model, output / (model.payload["fit"]["mode"] + ".json"))
    declaration = {
        "plan": PLAN,
        "training_acceptance_sha256": storage.checksum(acceptance.read_bytes()),
        "training_file_sha256": checked,
        "model_identities": [m.identity for m in models],
        "model_files_sha256": {
            m.payload["fit"]["mode"] + ".json": storage.checksum(m.json_text.encode())
            for m in models
        },
        "fits": [m.payload["fit"] for m in models],
    }
    storage.publish(
        output / "declaration.json", storage.json_bytes(storage.encode(declaration))
    )
    print("Frozen models and held-out plan published; no new QM executed.", flush=True)


def _check_frozen(output):
    d = storage.decode(storage.read_json(output / "declaration.json"))
    if d["plan"] != PLAN:
        raise ValueError("Declared plan changed")
    for file, checksum in d["model_files_sha256"].items():
        if storage.checksum((output / file).read_bytes()) != checksum:
            raise ValueError("Frozen model file changed")
    models = [
        load_pe_record(output / (mode + ".json")) for mode in ("baseline", "conserving")
    ]
    if [m.identity for m in models] != d["model_identities"] or [
        m.payload["fit"] for m in models
    ] != d["fits"]:
        raise ValueError("Frozen model definition changed")
    return d, models


def execute(output, amberhome):
    from island.builders import build_linear_polymer
    from island.forcefields import AmberToolsOptions, AmberToolsParameterizationEngine

    declaration, models = _check_frozen(output)
    if {p.name for p in output.iterdir()} - {
        "baseline.json",
        "conserving.json",
        "declaration.json",
        ".writer.lock",
    }:
        raise ValueError(
            "Execution requires untouched frozen directory; retries are forbidden"
        )
    cases_root = output / "held-out"
    cases_root.mkdir()
    refs = {}
    rows = []
    for dp, seed in HELD_OUT:
        name = f"pe-dp{dp}-seed{seed}"
        case = cases_root / name
        case.mkdir()
        row = {
            "case": name,
            "chemistry": "pe",
            "psmiles": PLAN["psmiles"],
            "dp": dp,
            "seed": seed,
            "status": "failed",
        }
        start = time.perf_counter()
        try:
            system = build_linear_polymer(
                PLAN["psmiles"],
                dp=dp,
                coordinate_method="local_templates",
                template_seed=seed,
                assembly_seed=seed,
            )
            row["sites"] = system.number_of_sites
            if row["sites"] != 6 * dp + 2 or row["sites"] > 100:
                raise ValueError("Unexpected authoritative site count")
            storage.publish(
                case / "input-system.json",
                storage.json_bytes(storage.encode(system_data(system))),
            )
            prep = AmberToolsParameterizationEngine().parameterize(
                system,
                AmberToolsOptions(
                    "gaff2",
                    "am1bcc",
                    max_atoms=100,
                    charge_tolerance=0.002,
                    timeout_seconds=600,
                    amberhome=amberhome,
                    work_root=case,
                    retain_success_artifacts=True,
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
                            "import_provenance": dict(prep.imported_result.provenance),
                        }
                    )
                ),
            )
            ref = create_charge_reference(system, prep)
            refs[name] = ref
            row.update(
                status="passed",
                record_signature=prep.record_signature,
                import_signature=prep.imported_result.result_signature,
                reference_identity=ref.identity,
                residual=prep.imported_result.charge_result.total_charge_residual,
            )
        except Exception as error:  # noqa: BLE001 -- preserve actual failures, no fallback
            row["failure"] = f"{type(error).__name__}: {error}"
        row["elapsed_seconds"] = time.perf_counter() - start
        row["source_checksums"] = {
            str(p.relative_to(case)): storage.checksum(p.read_bytes())
            for p in sorted(case.rglob("*"))
            if p.is_file() and not p.is_symlink()
        }
        rows.append(row)
        storage.publish(
            output / "matrix.json", storage.json_bytes({"cases": rows}), replace=True
        )
        print(name, row["status"], row.get("failure", ""), flush=True)
    targets, outcomes = [], []
    for row in rows:
        name = row["case"]
        outcome = {
            "case": name,
            "raw_import_status": row["status"],
            "raw_failure": row.get("failure"),
            "projected_target_available": False,
        }
        try:
            obs = observe_retained_case(
                cases_root,
                output / "matrix.json",
                name,
                historical_reference=refs.get(name),
            )
            projection = project_charge_observation(obs, policy=POLICY)
            save_record(obs, output / (name + "-observation.json"))
            save_record(projection, output / (name + "-projection.json"))
            if name in refs:
                save_record(refs[name], output / (name + "-reference.json"))
            targets.append(projection)
            outcome.update(
                projected_target_available=True,
                observation_identity=obs.identity,
                projection_identity=projection.identity,
                raw_total_e=projection.payload["projection"]["raw_total_e"],
                projected_total_e=projection.payload["projection"]["projected_total_e"],
            )
        except Exception as error:  # noqa: BLE001 -- eligibility is an independent gate
            outcome["observation_failure"] = f"{type(error).__name__}: {error}"
        outcomes.append(outcome)
    validation = None
    if targets:
        validation = validate_pe_templates(
            models, targets, frozen_model_identities=declaration["model_identities"]
        )
        save_record(validation, output / "validation.json")
    large = []
    for dp in PLAN["prediction_only_dps"]:
        system = build_linear_polymer(PLAN["psmiles"], dp=dp, generate_3d=False)
        result = {"dp": dp, "sites": system.number_of_sites, "models": {}}
        for model in models:
            pred = predict_pe_charges(model, system)
            name = model.payload["fit"]["mode"]
            save_record(pred, output / f"prediction-dp{dp}-{name}.json")
            if pred.payload != predict_pe_charges(model, system).payload:
                raise ValueError("Nondeterministic prediction")
            result["models"][name] = {
                "identity": pred.identity,
                "total_e": pred.payload["prediction"]["total_e"],
                "coverage": len(pred.payload["prediction"]["assignments_e"]),
            }
        large.append(result)
    _check_frozen(output)
    for row in rows:
        for path, checksum in row["source_checksums"].items():
            if (
                storage.checksum(
                    storage.child(cases_root / row["case"], path).read_bytes()
                )
                != checksum
            ):
                raise ValueError("Held-out source changed")
    complete = len(targets) == 4
    storage.publish(
        output / "outcomes.json",
        storage.json_bytes(
            {
                "cases": outcomes,
                "complete_projected_validation": complete,
                "successful_raw_imports": sum(r["status"] == "passed" for r in rows),
                "frozen_models_unchanged": True,
                "larger_prediction_only": large,
                "validation_identity": validation.identity if validation else None,
                "scientific_accuracy_accepted": False,
                **FLAGS,
            }
        ),
    )
    return complete


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--freeze", action="store_true")
    group.add_argument("--execute", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--training", type=Path)
    parser.add_argument("--original-source", type=Path)
    parser.add_argument(
        "--acceptance", type=Path, default=Path("docs/phase_4f4_acceptance.json")
    )
    parser.add_argument("--amberhome", type=Path)
    args = parser.parse_args(argv)
    started = False
    original_path = os.environ.get("PATH", "")
    try:
        if args.freeze:
            if args.training is None or args.original_source is None:
                raise ValueError("--freeze requires --training and --original-source")
            freeze(args.training, args.original_source, args.acceptance, args.output)
            return 0
        if args.amberhome:
            os.environ["PATH"] = (
                str(args.amberhome.resolve() / "bin") + os.pathsep + original_path
            )
        with storage.writer(args.output):
            started = True
            return 0 if execute(args.output, args.amberhome) else 1
    except Exception as error:  # noqa: BLE001 -- durable requested-gate failure
        if started and not (args.output / "failure.json").exists():
            storage.publish(
                args.output / "failure.json",
                storage.json_bytes(
                    {"failure": f"{type(error).__name__}: {error}", **FLAGS}
                ),
            )
        print(f"{type(error).__name__}: {error}", flush=True)
        return 1
    finally:
        os.environ["PATH"] = original_path


if __name__ == "__main__":
    raise SystemExit(main())
