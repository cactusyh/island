"""Engine-free replay of native and normalized real data after relocation."""

import json
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require
from island.pcff_qm.reference import CONVERSION, intake

from .campaign import NAMES, selected
from .importer import native


def replay(repo, distribution, evidence, work, normalized):
    repo, evidence, work, normalized = map(Path, (repo, evidence, work, normalized))
    manifest = json.loads((evidence / "plan/campaign.json").read_text())
    require(
        manifest["plan_sha256"] == pin(evidence / "plan/pilot-plan.json"),
        "pilot_plan_mutation",
    )
    require(
        manifest["adapter_sha256"] == pin(repo / "src/island/pcff_qm_pilot/engine.py"),
        "executed_adapter_mutation",
    )
    require(len(manifest["jobs"]) == 62, "pilot_budget_mutation")
    for key, identity in manifest["jobs"].items():
        request = json.loads((work / "requests" / (key + ".json")).read_text())
        require(digest(request) == identity, "native_request_mutation")
        native(work / "outputs" / key, request)
    results = {}
    for name in NAMES:
        scan = json.loads((evidence / "plan" / (name + ".json")).read_text())
        system, assessment, candidates = selected(
            repo, distribution, scan["center_coordinates"]
        )
        target = next(
            t for t in candidates if t.payload["identity"] == scan["target"]["identity"]
        )
        request = json.loads((normalized / (name + "-request.json")).read_text())
        original = json.loads((normalized / (name + "-intake.json")).read_text())
        root = normalized / name
        require(
            pin(root / "importer.py")
            == pin(repo / "src/island/pcff_qm_pilot/importer.py"),
            "executed_importer_mutation",
        )
        replayed = intake(
            system,
            assessment,
            target,
            request,
            root,
            verification_sha256=pin(root / "verification.json"),
        )
        require(replayed.payload == original, "normalized_intake_mutation")
        samples = json.loads((root / "samples.json").read_text())
        for sample in samples:
            key = name + "/" + sample["id"]
            job = json.loads((work / "requests" / (key + ".json")).read_text())
            raw, arrays = native(work / "outputs" / key, job)
            require(
                abs(
                    sample["energy"]
                    - raw["energy_hartree"] * CONVERSION["energy_factor"]
                )
                < 1e-9,
                "normalized_energy_mutation",
            )
            require(
                np.allclose(
                    sample["forces"],
                    -arrays["gradient"]
                    * CONVERSION["energy_factor"]
                    / CONVERSION["distance_factor"],
                    rtol=0,
                    atol=1e-12,
                ),
                "normalized_force_or_sign_mutation",
            )
            for path in (work / "outputs" / key).iterdir():
                require(
                    pin(path) == pin(root / "native" / sample["id"] / path.name),
                    "normalized_native_copy_mutation",
                )
        results[name] = {
            "dataset_identity": original["dataset_identity"],
            "intake_identity": original["identity"],
            "sample_count": 25,
            "status": "qm_data_available",
            "term_status": "qm_fit_required",
        }
    return {
        "schema": "island.j58-replay.v1",
        "native_completed_jobs": 62,
        "datasets": results,
        "parameter_completeness": "incomplete",
        "fit_completion": "not_run",
        "fit_validation": "not_run",
        "backend_execution": "not_run",
        "cross_engine_equivalence": "not_claimed",
        "scientific_validation": "not_established",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    import argparse
    import tempfile

    from island.pcff_completion import write

    from .bundle import unpack

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--distribution", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--json", type=Path, required=True)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="j58-relocated-") as temp:
        temp = Path(temp)
        unpack(args.evidence / "native-outputs.tar.xz", temp / "native")
        unpack(args.evidence / "reference-datasets.tar.xz", temp / "normalized")
        result = replay(
            args.repo,
            args.distribution,
            args.evidence,
            temp / "native",
            temp / "normalized",
        )
    write(args.json, result)


if __name__ == "__main__":
    main()
