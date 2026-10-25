"""Backend-scoped periodic workflow decisions from caller-verified evidence.

This policy checks evidence scope, bindings and completeness. It does not turn
an assertion into an independently verified measurement. Callers must verify
referenced artifacts before supplying gates. No global production attestation
or thermodynamic equilibrium is inferred from a backend decision.
"""

import hashlib
import json
import re

from island.exceptions import AnalysisError

BINDINGS = ("final_graph", "parameters", "charges", "backend", "protocol")
OPENMM = (
    "exact_graph_parameters",
    "validated_charge_origin",
    "cutoff_box_preflight",
    "internal_energy_force",
    "pressure_estimator_controls",
    "checkpoint_continuation",
    "sufficient_sampling",
)
CROSS = (
    "actual_polymer_systems",
    "numerical_policy",
    "component_energy",
    "per_atom_force",
    "virial",
    "pressure",
)
SCOPES = {
    **{"openmm." + key: "openmm" for key in OPENMM},
    **{"cross_engine." + key: "cross_engine" for key in CROSS},
    "lammps_export.exact_frame": "lammps_export",
    "lammps_kspace.native_execution": "lammps_kspace",
}
STATUSES = ("passed", "blocked", "unavailable", "not_measured")


def _require(condition, message):
    if not condition:
        raise AnalysisError(message)


def _signed(payload):
    body = {k: v for k, v in payload.items() if k != "identity"}
    encoded = json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return dict(body, identity=hashlib.sha256(encoded.encode()).hexdigest())


def _digest(value):
    return isinstance(value, str) and re.fullmatch(r"[0-9a-f]{64}", value) is not None


def assess_periodic_backend_readiness(*, bindings, checks):
    """Evaluate a single exact system/protocol; missing evidence blocks permission.

    Each supplied check needs its canonical name, backend scope, exact bindings,
    status, reason and SHA-256 artifact references. Successful checks require
    nonempty references. A LAMMPS export/runtime witness cannot satisfy OpenMM
    or cross-engine gates. All referenced evidence remains caller-verified.
    """
    _require(
        set(bindings) == set(BINDINGS) and all(_digest(v) for v in bindings.values()),
        "Invalid exact-system bindings",
    )
    _require(set(checks) <= set(SCOPES), "Unknown backend check")
    complete = {}
    for key, scope in SCOPES.items():
        gate = checks.get(
            key,
            {
                "name": key,
                "scope": scope,
                "bindings": dict(bindings),
                "status": "not_measured",
                "reason": "Required backend evidence not supplied",
                "evidence_sha256": [],
            },
        )
        _require(
            set(gate)
            == {"name", "scope", "bindings", "status", "reason", "evidence_sha256"},
            "Invalid gate fields: " + key,
        )
        _require(
            gate["name"] == key and gate["scope"] == scope,
            "Wrong backend scope: " + key,
        )
        _require(gate["bindings"] == bindings, "Stale system/protocol binding: " + key)
        _require(gate["status"] in STATUSES, "Invalid gate status: " + key)
        _require(
            isinstance(gate["reason"], str) and bool(gate["reason"].strip()),
            "Missing gate diagnostic: " + key,
        )
        refs = gate["evidence_sha256"]
        _require(
            isinstance(refs, list) and all(_digest(v) for v in refs),
            "Invalid evidence references: " + key,
        )
        _require(
            gate["status"] != "passed" or bool(refs),
            "Successful gate lacks evidence: " + key,
        )
        # JSON copy prevents a later caller mutation from changing the record.
        complete[key] = json.loads(json.dumps(gate, allow_nan=False))
    passed = lambda key: complete[key]["status"] == "passed"
    openmm_ok = all(passed("openmm." + key) for key in OPENMM)
    cross_ok = (
        all(passed("cross_engine." + key) for key in CROSS)
        and passed("lammps_kspace.native_execution")
        and passed("lammps_export.exact_frame")
    )
    return _signed(
        {
            "schema": "island.periodic-backend-readiness.v1",
            "bindings": dict(bindings),
            "checks": complete,
            "openmm_periodic_backend_status": "requirements_passed_for_bound_system"
            if openmm_ok
            else "software_tests_only",
            "lammps_periodic_export_status": "exact_frame_exported_only"
            if passed("lammps_export.exact_frame")
            else "export_not_verified",
            "lammps_kspace_runtime_status": "native_kspace_executed_only"
            if passed("lammps_kspace.native_execution")
            else "asset_unavailable_or_not_executed",
            "cross_engine_equivalence_status": "validated_for_bound_actual_polymer_system"
            if cross_ok
            else "ewald_unvalidated",
            "pressure_validation_status": "openmm_backend_pressure_validated"
            if passed("openmm.pressure_estimator_controls")
            else "pressure_unavailable",
            "scientific_equilibrium_status": "not_established",
            "openmm_periodic_workflow_permitted": openmm_ok,
            "scientific_npt_acceptance": False,
            "blocked_stages": [
                key for key, gate in complete.items() if gate["status"] != "passed"
            ],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def validate_backend_readiness(record, *, expected_bindings):
    """Recompute persisted decisions and reject replay to another system/protocol."""
    _require(
        record.get("bindings") == expected_bindings,
        "Stale expected system/protocol binding",
    )
    expected = assess_periodic_backend_readiness(
        bindings=expected_bindings, checks=record["checks"]
    )
    _require(
        record == expected and _signed(record)["identity"] == record.get("identity"),
        "Changed backend readiness decision or evidence",
    )
    return True


def require_openmm_periodic_workflow(record, *, expected_bindings):
    """Guard an opt-in qualified workflow; low-level software tests stay separate."""
    validate_backend_readiness(record, expected_bindings=expected_bindings)
    _require(
        record["openmm_periodic_workflow_permitted"],
        "OpenMM periodic workflow blocked: "
        + ", ".join(k for k in record["blocked_stages"] if k.startswith("openmm.")),
    )
