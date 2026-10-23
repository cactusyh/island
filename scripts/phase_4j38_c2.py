"""Frozen benchmark source decisions; incompatible charge sources never enable MD."""

from collections import Counter
from copy import deepcopy
from decimal import Decimal

from island.analysis.npt import FLAGS, _digest, _signed
from scripts.phase_4j38_c1_audit import require
from scripts.phase_4j38_c1_contract import FAILED
from scripts.phase_4j38_contract import validate_identity

BASE = "4e9ce54af161669577bdfda9360e64d6e1ef5361"


class SourceIncompatibleError(ValueError):
    """A native source vector disagrees with the declared neutral final graph."""


def aromatic_decision(audit, evidence):
    """Reconstruct exact decimal/component totals from independently audited rows."""
    graph = evidence["final_graph"]
    require(
        audit["final_graph_identity"] == graph["graph_identity"],
        "Stale source/graph binding",
    )
    require(audit["units"] == "elementary_charge", "Wrong charge units")
    require(
        audit["normalized"] is audit["redistributed"] is False,
        "Altered native charge vector",
    )
    sites = graph["sites"]
    rows = audit["rows"]
    require(
        [r["atom_id"] for r in rows] == [a["id"] for a in sites],
        "Changed atom-ID ordering",
    )
    summaries = []
    for index, ids in enumerate(graph["components"]):
        subset = [r for r in rows if r["atom_id"] in ids]
        atoms = [a for a in sites if a["id"] in ids]
        require(len(subset) == len(ids) == 22, "Omitted or duplicated component atoms")
        require(
            Counter(a["element"] for a in atoms) == {"C": 12, "H": 10},
            "Wrong complete biphenyl graph",
        )
        formal = sum(a["formal_charge"] for a in atoms)
        require(
            formal == 0 and all(a["formal_charge"] == 0 for a in atoms),
            "Wrong declared formal charge",
        )
        require(
            all(
                r["component"] == index
                and r["molecule"] == graph["molecule_membership"][str(r["atom_id"])]
                for r in subset
            ),
            "Changed molecule/component assignment",
        )
        require(
            all(
                r["element"] == a["element"]
                and r["formal_charge"] == a["formal_charge"]
                for r, a in zip(subset, atoms, strict=True)
            ),
            "Changed element/formal mapping",
        )
        require(
            all(
                type(r["charge_e"]) is float
                and Decimal(r["source_charge_decimal"]).is_finite()
                and float(r["source_charge_decimal"]) == r["charge_e"]
                for r in subset
            ),
            "Rounded, nonfinite or changed source vector",
        )
        source_total = sum(
            (Decimal(r["source_charge_decimal"]) for r in subset), Decimal()
        )
        measured = audit["components"][index]
        require(
            source_total
            == Decimal(measured["source_total_decimal_e"])
            == Decimal("-0.230"),
            "Changed native residual",
        )
        summaries.append(
            {
                "component": index,
                "atom_ids": ids,
                "molecule": subset[0]["molecule"],
                "formal_total_e": formal,
                "native_source_total_e": str(source_total),
                "residual_e": str(source_total - formal),
                "repeat_formula": "C6H5",
                "chain_formula": "C12H10",
                "chemically_intended_ionic_repeat": False,
            }
        )
    total = sum((Decimal(c["native_source_total_e"]) for c in summaries), Decimal())
    return _signed(
        {
            "schema": "island.j38-c2-native-source-decision.v1",
            "charge_origin": "native_oplsaa",
            "final_graph_identity": graph["graph_identity"],
            "packing_identity": evidence["packing_plan"]["identity"],
            "source_sha256": audit["source_sha256"],
            "audit_binding": audit["evidence_identity"],
            "components": summaries,
            "formal_system_total_e": 0,
            "native_system_total_e": str(total),
            "source_compatible": False,
            "classification": "source_incompatible",
            "stage": "charge_completeness",
            "diagnostic": "Pinned native OPLS-AA selected rows sum to -0.230 e per neutral C12H10 chain (-0.460 e system); incompatible with this declared neutral final graph. No ionic repeats, charged end groups, counterions or compensating-background model are declared.",
            "interpretation": "A periodic engine may define a charged-cell background convention, but that would be a different physical model; it cannot validate this neutral benchmark. No thermodynamic or experimental conclusion is made.",
            "external_charge_substitution": False,
            "normalization": False,
            "redistribution": False,
            "native_path_allowed": False,
            "dynamics_executed": False,
            **FLAGS,
        }
    )


def require_native_source(decision):
    validate_identity(decision)
    require(
        decision["charge_origin"] == "native_oplsaa",
        "External charges cannot enter the native OPLS path",
    )
    if not decision["source_compatible"]:
        raise SourceIncompatibleError(decision["diagnostic"])
    raise ValueError(
        "No new compatible aromatic source reviewed in this frozen tranche"
    )


def matrix(c1, decision, inventory):
    validate_identity(c1)
    validate_identity(decision)
    cases = {}
    for key, previous in c1["cases"].items():
        row = deepcopy(previous)
        row["c1_row_identity"] = _digest(previous)
        row["c1_matrix_identity"] = c1["identity"]
        native = key.startswith("oplsaa/")
        missing = key in {"gaff2/pe", "gaff2/aromatic"}
        status = (
            "source_incompatible"
            if key == "oplsaa/aromatic"
            else "asset_unavailable"
            if missing
            else "retained_verified"
        )
        row["selected_charge_origin"] = "native_oplsaa" if native else "RESP/Psi4"
        row["charge_source_status"] = status
        row["charge_paths"] = {
            "native_oplsaa": {
                "status": status if native else "not_applicable",
                "origin": "native_oplsaa",
                "source_sha256": decision["source_sha256"] if native else None,
            },
            "externally_provided": {
                "status": "not_selected",
                "origin": "external_provided",
                "enabled": False,
                "reason": "No reviewed external vector introduced; an external vector could never be labeled native OPLS, RESP/Psi4 or AM1-BCC",
            },
            "resp_psi4_exact_graph": {
                "status": status if not native else "not_applicable",
                "origin": "RESP/Psi4",
                "exact_graph_verified": not native and not missing,
                "review_identity": next(
                    (
                        r["independent_review_identity"]
                        for r in inventory["inventory"]
                        if r["final_graph_identity"] == row["final_graph_identity"]
                    ),
                    None,
                )
                if not native
                else None,
            },
        }
        row["source_artifacts"] = {
            "status": "asset_unavailable" if missing else "present",
            "requested_graph_identity": row["final_graph_identity"],
        }
        if key == "oplsaa/aromatic":
            row["stages"]["charge_completeness"] = {
                "status": "source_incompatible",
                "reason": decision["diagnostic"],
            }
        row["execution_origin"] = (
            "retained_j38_readonly" if key not in FAILED else "not_executed"
        )
        row["new_dynamics_executed"] = False
        cases[key] = row
    return _signed(
        {
            "schema": "island.j38-c2-matrix.v1",
            "base_commit": BASE,
            "c1_matrix_identity": c1["identity"],
            "source_decision_identity": decision["identity"],
            "cases": cases,
            "acceptance_met": False,
            "new_md_executions": 0,
            **FLAGS,
        }
    )


def validate_matrix(value, c1, decision, inventory):
    validate_identity(value)
    require(
        value == matrix(c1, decision, inventory),
        "Changed source origin, graph, provenance, stage outcome or acceptance",
    )
    return value
