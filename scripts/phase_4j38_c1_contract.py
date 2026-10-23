"""Additive nine-stage J38 correction; historical success is explicitly inherited."""

from island.analysis.npt import FLAGS, _signed
from scripts.phase_4j38_c1_audit import require
from scripts.phase_4j38_contract import CASES, gate, validate_identity, validate_matrix

BASE = "14b1b40ae43de22bf510731859e8e6cd337d9926"
STAGES = (
    "graph_construction",
    "packing",
    "typing",
    "charge_completeness",
    "parameter_completeness",
    "executable_backend",
    "numerical_reference",
    "equilibration",
    "scientific_validation",
)
FAILED = {"oplsaa/aromatic", "gaff2/pe", "gaff2/aromatic"}


def matrix(history, audit, inventory):
    validate_matrix(history)
    cases = {}
    for key, original in history["cases"].items():
        stages = {
            s: gate(
                "not_executed",
                "Source/charge gate blocks downstream preparation and MD",
            )
            for s in STAGES
        }
        stages["graph_construction"] = gate(
            "passed",
            "Exact frozen complete final graph reconstructed; stable IDs and formal charges bound",
        )
        stages["packing"] = gate(
            "passed",
            "Frozen production packing reconstructed without changing graph, box or atom map",
        )
        if key not in FAILED:
            for s in STAGES[:7]:
                stages[s] = gate(
                    "passed",
                    "Read-only retained J38 software evidence; not a new C1 simulation",
                )
        elif key == "oplsaa/aromatic":
            stages["typing"] = gate(
                "passed",
                "Independent pinned Foyer raw-graph typing agrees exactly with native stable-ID assignments",
            )
            stages["charge_completeness"] = gate("source_incomplete", audit["cause"])
        else:
            stages["charge_completeness"] = gate(
                "asset_unavailable", inventory["requested"][key]["reason"]
            )
        stages["equilibration"] = gate(
            "insufficient_sampling" if original["execution_met"] else "not_executed",
            "No thermodynamic equilibrium established; one short retained run where executed",
        )
        stages["scientific_validation"] = gate(
            "blocked", "No experimental, material-property or production validation"
        )
        cases[key] = dict(
            case=key,
            family=original["family"],
            benchmark=original["benchmark"],
            final_graph_identity=original["input_provenance"]["final_graph_identity"],
            packing_identity=original["input_provenance"]["packing_identity"],
            historical_row_identity=original["identity"],
            stages=stages,
            execution_origin="retained_j38"
            if original["execution_met"]
            else "not_executed",
            new_dynamics_executed=False,
            acceptance_met=False,
            **FLAGS,
        )
    return _signed(
        dict(
            schema="island.j38-c1-matrix.v1",
            base_commit=BASE,
            historical_matrix_identity=history["identity"],
            cases=cases,
            acceptance_met=False,
            **FLAGS,
        )
    )


def validate(value, history, audit, inventory):
    validate_identity(value)
    require(set(value["cases"]) == set(CASES), "Changed matrix inventory")
    require(
        value == matrix(history, audit, inventory),
        "Changed stage outcomes, identities or unsupported success claim",
    )
    return value
