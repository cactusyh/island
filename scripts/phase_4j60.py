"""Execute, reconstruct and analyze the bounded J60 Hamiltonian/QM campaign."""

import argparse
import json
import tempfile
from pathlib import Path

from island.pcff_completion import inventory, require, signed, write
from island.pcff_qm_pilot.bundle import unpack
from island.pcff_qm_pilot.hamiltonian import accounting, load_model
from island.pcff_qm_pilot.hamiltonian_analysis import convergence, normalize, numerical


def replay(repo, evidence):
    repo, evidence = Path(repo), Path(evidence)
    require(
        accounting(load_model(repo))
        == json.loads((evidence / "interaction-accounting.json").read_text()),
        "accounting_replay_mutation",
    )
    require(
        numerical(repo, evidence)
        == json.loads((evidence / "numerical-checks.json").read_text()),
        "numerical_replay_mutation",
    )
    from scripts.full_graph_phase_4j60 import reconstruct

    reconstruct(repo, evidence)
    with tempfile.TemporaryDirectory(prefix="j60-native-replay-") as temp:
        native = Path(temp) / "native"
        unpack(evidence / "native-outputs.tar.xz", native)
        require(
            inventory(native)
            == json.loads((evidence / "native-file-inventory.json").read_text()),
            "native_inventory_mutation",
        )
        data = normalize(evidence, native)
        require(
            data == json.loads((evidence / "normalized-reference.json").read_text()),
            "normalized_replay_mutation",
        )
        from scripts.independent_phase_4j59 import check

        require(
            check(native, evidence / "normalized-reference.json")
            == json.loads((evidence / "independent-native-checks.json").read_text()),
            "independent_native_replay_mutation",
        )
        from scripts.residuals_phase_4j60 import residuals

        require(
            residuals(repo, evidence)
            == json.loads((evidence / "real-residual-subtraction.json").read_text()),
            "fixed_subtraction_replay_mutation",
        )
        result = convergence(repo, evidence, data)
        require(
            result
            == json.loads(
                (evidence / "primary-reference-convergence.json").read_text()
            ),
            "primary_reference_replay_mutation",
        )
        stage_data = {}
        for stage in sorted((evidence / "stages").glob("*")):
            if not (stage / "native-outputs.tar.xz").exists():
                continue
            target = Path(temp) / stage.name
            unpack(stage / "native-outputs.tar.xz", target)
            require(
                inventory(target)
                == json.loads((stage / "native-file-inventory.json").read_text()),
                "stage_native_inventory_mutation",
            )
            normalized = normalize(stage, target)
            require(
                normalized
                == json.loads((stage / "normalized-reference.json").read_text()),
                "stage_native_reconstruction_mutation",
            )
            require(
                check(target, stage / "normalized-reference.json")
                == json.loads((stage / "independent-native-checks.json").read_text()),
                "stage_independent_units_mutation",
            )
            if (stage / "real-residual-subtraction.json").exists():
                require(
                    residuals(repo, stage)
                    == json.loads(
                        (stage / "real-residual-subtraction.json").read_text()
                    ),
                    "stage_fixed_subtraction_mutation",
                )
            stage_data[stage.name] = normalized
        if "tight-d-basis" in stage_data:
            corrected = stage_data["tight-d-basis"]
            if "tight-d-numerics" in stage_data:
                corrected = signed(
                    {
                        **{k: v for k, v in corrected.items() if k != "identity"},
                        "samples": {
                            **corrected["samples"],
                            **stage_data["tight-d-numerics"]["samples"],
                        },
                    }
                )
            require(
                corrected
                == json.loads(
                    (
                        evidence / "tight-d-combined-normalized-reference.json"
                    ).read_text()
                ),
                "refinement_combination_mutation",
            )
            result = convergence(repo, evidence, corrected)
            require(
                result
                == json.loads((evidence / "refinement-convergence.json").read_text()),
                "refinement_analysis_mutation",
            )
            from scripts.refine_phase_4j60 import basis_effects, replay_fit

            require(
                basis_effects(data, stage_data["tight-d-basis"])
                == json.loads((evidence / "tight-d-basis-effects.json").read_text()),
                "basis_effect_reconstruction_mutation",
            )
            replay_fit(repo, evidence, stage_data)
            from scripts.residuals_phase_4j60 import family_check

            require(
                family_check(repo, evidence)
                == json.loads(
                    (evidence / "cross-family-reference-check.json").read_text()
                ),
                "cross_family_reconstruction_mutation",
            )
        require(
            result == json.loads((evidence / "reference-convergence.json").read_text()),
            "final_reference_replay_mutation",
        )
    return {
        "status": "passed",
        "parameter_completeness": "incomplete",
        "QM_data_integrity": "passed",
        "fit_completed": json.loads((evidence / "fitting-decision.json").read_text())[
            "fit_completed"
        ],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation", choices=("numerical", "analyze", "replay"))
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j60"))
    p.add_argument("--work", type=Path)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    if a.operation == "numerical":
        result = numerical(a.repo, a.evidence)
        write(a.evidence / "numerical-checks.json", result)
    elif a.operation == "analyze":
        data = normalize(a.evidence, a.work)
        write(a.evidence / "normalized-reference.json", data)
        result = convergence(a.repo, a.evidence, data)
        write(a.evidence / "reference-convergence.json", result)
        from scripts.independent_phase_4j59 import check

        write(
            a.evidence / "independent-native-checks.json",
            check(a.work, a.evidence / "normalized-reference.json"),
        )
    else:
        result = replay(a.repo, a.evidence)
    if a.output:
        write(a.output, result)
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
