"""Analyze/replay J59 using existing native, source and intake infrastructure."""

import argparse
import json
import tempfile
from pathlib import Path

from island.pcff_completion import inventory, require, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.convergence_analysis import analyze, center_report
from island.pcff_qm_pilot.convergence_importer import normalize
from island.pcff_qm_pilot.residual import reconstruct
from island.pcff_qm_pilot.residual_analysis import checks


def replay(repo, evidence, distribution):
    repo, evidence = Path(repo).resolve(), Path(evidence).resolve()
    model = json.loads((evidence / "interaction-table.json").read_text())
    reconstruct(model, repo, distribution)
    plan = json.loads((evidence / "plan/convergence-plan.json").read_text())
    require(
        checks(model, plan)
        == json.loads((evidence / "fixed-evaluator-checks.json").read_text()),
        "fixed_evaluator_replay_mutation",
    )
    with tempfile.TemporaryDirectory(prefix="j59-native-replay-") as temp:
        temp = Path(temp)
        unpack(evidence / "native-outputs.tar.xz", temp / "native")
        unpack(repo / "docs/evidence/phase_4j58/native-outputs.tar.xz", temp / "j58")
        require(
            inventory(temp / "native")
            == json.loads((evidence / "native-file-inventory.json").read_text()),
            "native_recursive_inventory_mutation",
        )
        normalized = normalize(evidence, temp / "native")
        require(
            normalized
            == json.loads((evidence / "normalized-convergence.json").read_text()),
            "normalized_native_replay_mutation",
        )
        convergence = analyze(repo, evidence, normalized, temp / "j58")
        require(
            convergence == json.loads((evidence / "QM-convergence.json").read_text()),
            "convergence_replay_mutation",
        )
        center = center_report(temp / "native", repo)
        require(
            center == json.loads((evidence / "center-effects.json").read_text()),
            "center_replay_mutation",
        )
        from scripts.independent_phase_4j59 import check

        independent = check(temp / "native", evidence / "normalized-convergence.json")
        require(
            independent
            == json.loads((evidence / "independent-unit-checks.json").read_text()),
            "independent_unit_replay_mutation",
        )
    return {
        "status": "passed",
        "model_identity": model["identity"],
        "normalized_identity": normalized["identity"],
        "convergence_identity": convergence["identity"],
        "native_gradient_jobs": len(normalized["samples"]),
        "fixed_evaluator_verified": True,
        "fit_completed": False,
        "parameter_completeness": "incomplete",
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=("analyze", "replay", "residual"))
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j59"))
    p.add_argument("--distribution", type=Path, required=True)
    p.add_argument("--work", type=Path)
    p.add_argument("--j58-work", type=Path)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    if a.operation == "residual":
        from island.pcff_qm_pilot.residual import build
        from island.pcff_qm_pilot.residual_analysis import analyze as residual_analysis

        model = build(a.repo, a.distribution)
        write(a.evidence / "interaction-table.json", model)
        plan = json.loads((a.evidence / "plan/convergence-plan.json").read_text())
        write(a.evidence / "fixed-evaluator-checks.json", checks(model, plan))
        with tempfile.TemporaryDirectory(prefix="j59-source-reconstruction-") as temp:
            normalized = Path(temp) / "J58"
            unpack(
                a.repo / "docs/evidence/phase_4j58/reference-datasets.tar.xz",
                normalized,
            )
            write(
                a.evidence / "residual-identifiability.json",
                residual_analysis(model, a.repo, normalized),
            )
        result = {"model_identity": model["identity"], "complete_residual": False}
    elif a.operation == "replay":
        result = replay(a.repo, a.evidence, a.distribution)
    else:
        normalized = normalize(a.evidence, a.work)
        write(a.evidence / "normalized-convergence.json", normalized)
        write(
            a.evidence / "QM-convergence.json",
            analyze(a.repo, a.evidence, normalized, a.j58_work),
        )
        write(a.evidence / "center-effects.json", center_report(a.work, a.repo))
        result = {
            "normalized_identity": normalized["identity"],
            "gradient_jobs": len(normalized["samples"]),
        }
    if a.output:
        write(a.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
