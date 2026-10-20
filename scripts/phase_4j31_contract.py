"""Frozen J31 acceptance outcomes; diagnostics alone never prove acceptance."""

import math
import re
from hashlib import sha256

from island.analysis.network import FLAGS, _digest
from island.exceptions import AnalysisError

CONTRACT = "island.phase-4j31-acceptance.v1"
PCFF_SHA = "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c"
OPLS_SHA = "c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4"
PCFF_MESSAGE = "Incomplete automatic typing cannot assign charges"
CASE_BINDINGS = {
    "branched": {
        "crosslink_plan_identity": "2bd46fad3a54229640597a6a8bd139b6a209ecfbdd0b2e919b08dd4142ae5d9d",
        "final_graph_identity": "9a620a0f2158ca56d2015ac7e28c4315c20b0eddea7a186b238e522f3a4ec6e6",
        "packing_identity": "bf4e023937a18b44f39fe4632a59f3ce08780f69bd7af1ea606cb915e8149167",
    },
    "crosslinked": {
        "crosslink_plan_identity": "c07b86b54069f63687fbfe15cbfe250f74e425a5a1c92490a3f397190405c052",
        "final_graph_identity": "7ed9500b1d1c2debf3eb3a871fc9b82e322251f6138feeca25bdf147b1e4df2a",
        "packing_identity": "e21b8a782081f4afcdbace068f8d83ce7d89cd92458f1008f8d2a512c0d91a93",
    },
}
LEGACY_PCFF = {
    "pcff/branched": "f54b393888b01fa220945449d8f537db69f204e3d9cdfc0ae21f5a7d62e19677",
    "pcff/crosslinked": "91a9c5d6f466f3522ab3f6e12046dcc5a01aa39a69b968c02e4a4b1f28f03e1c",
}
KINDS = ("branched", "crosslinked")
FAMILIES = ("pcff", "oplsaa", "gaff", "gaff2")
CASES = tuple(f"{f}/{k}" for f in FAMILIES for k in KINDS)
SUCCESS_FIELDS = (
    "backend_identity",
    "source",
    "typed_graph_identity",
    "charge_identity",
    "assignment_identity",
    "independent_assignment",
    "independent_assignment_verification",
    "execution",
)


def valid_identity(value):
    return type(value) is str and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def typing_evidence(fixture, source_path):
    """Independent read-only native typing replay, never parameterization or MD."""
    from island.forcefields.pcff import load_pcff_source, type_pcff_atoms
    from island.periodic import periodic_chemical_view

    raw = source_path.read_bytes()
    if sha256(raw).hexdigest() != PCFF_SHA:
        raise ValueError("Pinned PCFF source checksum mismatch")
    source = load_pcff_source(source_path)
    typing = type_pcff_atoms(periodic_chemical_view(fixture.packed.system), source)
    return {
        "source_sha256": source.identity["sha256"],
        "source_verified": True,
        "final_graph_identity": fixture.packed.graph.identity,
        "profile": typing.payload["profile"]["name"],
        "typing_identity": typing.identity,
        "chemical_graph_identity": typing.payload["graph_identity"],
        "complete": typing.complete,
        "diagnostics": typing.payload["diagnostics"],
    }


def expected_typing(kind, evidence):
    diagnostics = [{"sites": [3], "reason": "multiple_oxygen_carbon"}]
    if kind == "crosslinked":
        diagnostics.append({"sites": list(range(1, 29)), "reason": "cyclic_component"})
    return (
        isinstance(evidence, dict)
        and evidence.get("source_verified") is True
        and evidence.get("source_sha256") == PCFF_SHA
        and evidence.get("final_graph_identity")
        == CASE_BINDINGS[kind]["final_graph_identity"]
        and evidence.get("profile") == "island_pcff_acyclic_cho_v1"
        and evidence.get("complete") is False
        and evidence.get("diagnostics") == diagnostics
        and valid_identity(evidence.get("typing_identity"))
        and valid_identity(evidence.get("chemical_graph_identity"))
    )


def case_errors(name, row, *, legacy_pcff=False):
    """Readable failed rows return errors; only an empty list meets acceptance."""
    errors = []

    def check(ok, message):
        if not ok:
            errors.append(f"{name}: {message}")

    if name not in CASES or not isinstance(row, dict):
        return [f"{name}: invalid case/record"]
    family, kind = name.split("/")
    check(
        row.get("family") == family and row.get("fixture") == kind,
        "case identity mismatch",
    )
    check(
        all(row.get(k) == v for k, v in CASE_BINDINGS[kind].items()),
        "graph/plan/packing binding mismatch",
    )
    check(
        all(row.get(k) == v and type(row.get(k)) is type(v) for k, v in FLAGS.items()),
        "readiness mismatch",
    )
    check(
        "acceptance_met" not in row or row["acceptance_met"] is True,
        "explicit acceptance flag is false",
    )
    diagnostic = row.get("diagnostic", {})
    if not isinstance(diagnostic, dict):
        diagnostic = {}
        check(False, "malformed diagnostic")
    if family == "pcff":
        check(
            row.get("status") == "blocked" and row.get("dynamics_executed") is False,
            "frozen PCFF negative control not met",
        )
        check(
            not any(k in row for k in SUCCESS_FIELDS) and "partial_evidence" not in row,
            "contradictory success fields in negative control",
        )
        # Original J31 did not serialize stage/source diagnostics. Only its exact
        # frozen rows can use compatibility; the full validator must replay typing.
        if legacy_pcff and "diagnostic" not in row:
            check(
                _digest(row) == LEGACY_PCFF[name], "unrecognized legacy PCFF rejection"
            )
        else:
            check(
                diagnostic.get("classification") == "expected_typing_rejection"
                and diagnostic.get("stage") == "typing"
                and diagnostic.get("exception_type") == "PCFFError"
                and diagnostic.get("message") == PCFF_MESSAGE
                and diagnostic.get("origin")
                == "island.forcefields.pcff.automatic.bridge_data",
                "wrong PCFF rejection stage/exception/domain",
            )
            check(
                expected_typing(kind, diagnostic.get("typing_evidence", {})),
                "PCFF verified source/domain evidence missing",
            )
            check(
                row.get("reason") == "PCFFError: " + PCFF_MESSAGE,
                "wrong PCFF diagnostic",
            )
        return errors
    check(
        row.get("status") == "complete" and row.get("dynamics_executed") is True,
        "required success not completed",
    )
    check(
        "reason" not in row and not row.get("partial_evidence"),
        "contradictory failure fields in success",
    )
    if diagnostic:
        check(
            diagnostic.get("classification") == "completed"
            and diagnostic.get("stage") == "complete",
            "success diagnostic contradicts outcome",
        )
    for key in (
        "backend_identity",
        "typed_graph_identity",
        "charge_identity",
        "assignment_identity",
    ):
        check(valid_identity(row.get(key)), "missing assignment identity: " + key)
    source = row.get("source", {})
    if not isinstance(source, dict):
        source = {}
    check(
        isinstance(source, dict)
        and source.get("family")
        == {"oplsaa": "OPLS-AA", "gaff": "GAFF", "gaff2": "GAFF2"}[family]
        and valid_identity(source.get("sha256")),
        "missing source binding",
    )
    if family == "oplsaa":
        check(
            source.get("sha256") == source.get("xml_sha256") == OPLS_SHA,
            "pinned OPLS source mismatch",
        )
    inventory = row.get("independent_assignment", {})
    if not isinstance(inventory, dict):
        inventory = {}
    check(
        all(
            inventory.get(k) == v
            for k, v in {
                "bonds": 29 if kind == "branched" else 28,
                "angles": 52,
                "proper_paths": 57 if kind == "branched" else 66,
            }.items()
        ),
        "independent final-graph inventory missing",
    )
    check(
        inventory.get("source_terms") == "all required paths present"
        and bool(inventory.get("independent_reader"))
        and bool(row.get("independent_assignment_verification")),
        "independent assignment verification missing",
    )
    if family in ("gaff", "gaff2"):
        check(
            inventory.get("parmchk2_estimation_lines") == [],
            "estimated or unverified Amber terms",
        )
        for key, bound in (
            ("energy_error_kj_mol", 1e-6),
            ("max_force_error_kj_mol_angstrom", 1e-5),
        ):
            value = inventory.get(key)
            check(
                type(value) in (int, float)
                and math.isfinite(value)
                and 0 <= value <= bound,
                "independent numerical check missing: " + key,
            )
    execution = row.get("execution", {})
    if not isinstance(execution, dict):
        execution = {}
    comparison = execution.get("comparison", {})
    if not isinstance(comparison, dict):
        comparison = {}
    check(
        valid_identity(comparison.get("manifest_identity")),
        "comparison identity missing",
    )
    check(
        all(
            comparison.get(k) is True
            for k in (
                "passed",
                "complete",
                "exact_checkpoint_and_rng_bytes",
                "exact_frames",
            )
        )
        and all(comparison.get(k) == v for k, v in FLAGS.items()),
        "split/continuous comparison missing",
    )
    check(
        execution.get("same_process_relocation")
        == execution.get("separate_process_reconstruction")
        == "PASS",
        "relocation/reconstruction missing",
    )
    check(execution.get("frame_count") == 5, "execution frame count mismatch")
    for key in (
        "analysis_record_identity",
        "selection_identity",
        "frame_content_identity",
    ):
        check(valid_identity(execution.get(key)), "execution identity missing: " + key)
    observables = execution.get("network_observables", {})
    if not isinstance(observables, dict):
        observables = {}
    required = {
        "network_connectivity",
        "degree_distribution",
        "branch_point_statistics",
        "cluster_size_distribution",
        "gel_fraction",
        "cycle_rank",
        "strand_length_distribution",
        "dangling_end_fraction",
        "crosslink_density",
        "network_radius_of_gyration",
        "periodic_spanning",
        "contact_fraction",
        "periodic_geometry",
    }
    check(
        set(observables) == required
        and all(
            isinstance(v, dict) and v.get("status") == "ok"
            for v in observables.values()
        ),
        "network observable evidence incomplete",
    )
    return errors


def require_case(name, row, *, legacy_pcff=False):
    errors = case_errors(name, row, legacy_pcff=legacy_pcff)
    if errors:
        raise AnalysisError("; ".join(errors))


def matrix_errors(matrix, *, legacy_pcff=False):
    errors = []
    if not isinstance(matrix, dict):
        return ["J31 matrix must be an object"]
    if set(matrix) != set(CASES):
        errors.append("J31 exact eight-case inventory required")
    for name, row in matrix.items():
        errors.extend(case_errors(name, row, legacy_pcff=legacy_pcff))
    return errors


def require_matrix(matrix, *, legacy_pcff=False):
    errors = matrix_errors(matrix, legacy_pcff=legacy_pcff)
    if errors:
        raise AnalysisError("; ".join(errors))
