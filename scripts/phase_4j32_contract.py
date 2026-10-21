"""Frozen J32 case outcomes shared by execution and evidence validation."""

from island.analysis.network import _body, _digest
from scripts.phase_4j31_contract import PCFF_MESSAGE, PCFF_SHA, valid_identity
from scripts.phase_4j32_real_fixture import require

CASES = ("oplsaa/linear", "oplsaa/network", "gaff2/network", "pcff/network")
GRAPHS = {
    "linear": "412ac72f137764b52ec5bb4bf222f9097dbbe7c95faaf8e2fb3ccdb8aa4e1e79",
    "network": "01ec38d6b0177fe2f3f40e747f6c53315034be9b7902fd168b21e69a83e47407",
}


def validate_case(name, row):
    require(name in CASES, "Unknown required case")
    family, kind = name.split("/")
    require(
        row.get("family") == family and row.get("fixture") == kind,
        "Matrix/case identity mismatch",
    )
    require(row.get("final_graph_identity") == GRAPHS[kind], "Wrong frozen graph")
    require(
        row.get("production_validated") is False
        and row.get("simulation_readiness") == "not_established",
        "Readiness promotion",
    )
    require(row.get("acceptance_met") is True, "Required acceptance not met")
    require(row.get("identity") == _digest(_body(row)), "Case identity mismatch")
    if family != "pcff":
        require(
            row.get("status") == "complete" and row.get("dynamics_executed") is True,
            "Mandatory dynamics did not complete",
        )
        require(
            row.get("same_process_relocation") == "PASS"
            and row.get("fresh_process_reconstruction") == "PASS"
            and row.get("comparison", {}).get("passed") is True
            and valid_identity(row.get("analysis_identity")),
            "Missing execution/reconstruction evidence",
        )
    else:
        require(
            row.get("status") == "blocked" and row.get("dynamics_executed") is False,
            "Unexpected PCFF outcome",
        )
        require(
            not any(
                k in row
                for k in (
                    "comparison",
                    "analysis_identity",
                    "same_process_relocation",
                    "fresh_process_reconstruction",
                )
            ),
            "Contradictory PCFF success evidence",
        )
        d = row.get("diagnostic", {})
        require(
            d.get("classification") == "expected_typing_rejection"
            and d.get("stage") == "typing"
            and d.get("exception_type") == "PCFFError"
            and d.get("message") == PCFF_MESSAGE
            and d.get("origin") == "island.forcefields.pcff.automatic.bridge_data",
            "Wrong PCFF rejection stage/reason",
        )
        e = d.get("typing_evidence", {})
        require(
            e.get("source_verified") is True
            and e.get("source_sha256") == PCFF_SHA
            and e.get("final_graph_identity") == GRAPHS[kind]
            and e.get("profile") == "island_pcff_acyclic_cho_v1"
            and e.get("complete") is False
            and e.get("typing_identity")
            == "79288e9731cf8b8e9ccc8900a87cfd2429374ea7ba5078116d0b7f8de274ee16"
            and e.get("diagnostics")
            == [{"sites": list(range(1, 64)), "reason": "cyclic_component"}],
            "Wrong PCFF source/typing/graph evidence",
        )
    return row


def validate_matrix(matrix):
    require(
        matrix.get("identity") == _digest(_body(matrix)), "Matrix identity mismatch"
    )
    require(set(matrix["cases"]) == set(CASES), "Required matrix inventory mismatch")
    for name, case in matrix["cases"].items():
        require(case["exit_code"] == 0, "Case CLI failed: " + name)
        validate_case(name, case["result"])
    require(matrix.get("acceptance_met") is True, "Matrix acceptance not met")
    return matrix
