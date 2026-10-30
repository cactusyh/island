"""Strict data intake with exact scan coverage and raw-output/importer binding.

Integrity verification is a producer attestation, not proof of a QM calculation.
No imported dataset creates a fitted parameter or grants MD authority.
"""

from pathlib import Path

import numpy as np

from island.pcff_completion import (
    _finite,
    digest,
    inventory,
    read_pinned,
    require,
)

from .geometry import geometries
from .records import record
from .targets import BINDINGS, validate_target

CONVERSION = {
    "source": "CODATA 2022 hartree/bohr; SI 2019 exact Avogadro; thermochemical calorie 4.184 J",
    "energy_factor": 4.3597447222060e-18 * 6.02214076e23 / 4184,
    "distance_factor": 5.29177210544e-11 * 1e10,
}


def _text(value):
    return isinstance(value, str) and bool(value.strip())


def intake(system, assessment, target, request, root, *, verification_sha256):
    t = validate_target(system, assessment, target)
    require(
        request.get("schema") == "island.pcff-qm-reference.v2",
        "reference_schema_required",
    )
    require(request.get("family") == "pcff", "wrong_profile_family")
    require(
        request.get("target_identity") == t["identity"], "reference_target_mutation"
    )
    for k in BINDINGS:
        require(request.get(k) == t[k], k + "_mutation")
    require(request.get("graph") == t["graph"], "graph_or_bond_order_mutation")
    require(
        request.get("functional_form_identity") == t["functional_form"]["identity"],
        "incompatible_functional_form",
    )
    require(request.get("units") == t["units"], "reference_units_missing_or_changed")
    files = inventory(root)
    require(request.get("files") == files, "reference_inventory_or_hash_mutation")
    license_file = request.get("license_file")
    require(
        license_file in files and _text((Path(root) / license_file).read_text()),
        "license_evidence_required",
    )
    require(
        request.get("permission_status") in ("licensed", "permission_granted"),
        "reference_permission_required",
    )
    methods = request.get("methods")
    require(isinstance(methods, dict) and methods, "qm_method_required")
    require(
        len(methods) == 1 or _text(request.get("mixed_method_protocol")),
        "undeclared_mixed_qm_methods",
    )
    protocol = t["qm_protocol"]
    for m in methods.values():
        require(isinstance(m, dict), "qm_method_metadata_incomplete")
        require(
            all(
                _text(m.get(k))
                for k in (
                    "method",
                    "basis_set",
                    "dispersion_treatment",
                    "software",
                    "version",
                )
            ),
            "qm_method_metadata_incomplete",
        )
        for key in ("method", "basis_set", "dispersion_treatment"):
            require(m[key] == protocol[key], "qm_protocol_mismatch", field=key)
        require(
            m.get("calculation_settings") == protocol,
            "calculation_settings_missing_or_changed",
        )
    vf = request.get("verification_file")
    require(vf in files and files[vf] == verification_sha256, "unverified_output_files")
    proof = read_pinned(root, vf, files)
    importer = proof.get("importer_file")
    require(
        importer in files
        and proof.get("importer_sha256") == files[importer]
        and _text(proof.get("implementation_version")),
        "unverified_output_importer",
    )
    require(
        proof.get("normalization")
        == {
            "energy": "hartree to kcal/mol by declared factor",
            "forces": "minus gradient; hartree/bohr to kcal/mol/angstrom",
            "conversion_source": request.get("conversion_source"),
            "energy_factor": request.get("energy_factor"),
            "distance_factor": request.get("distance_factor"),
        },
        "output_conversion_evidence_required",
    )
    require(
        _text(request.get("conversion_source")), "output_conversion_source_required"
    )
    for key in ("energy_factor", "distance_factor"):
        require(
            _finite(request.get(key)).shape == () and float(request[key]) > 0,
            "invalid_output_conversion",
            field=key,
        )
        require(
            request[key] == CONVERSION[key], "output_conversion_mutation", field=key
        )
    require(
        request["conversion_source"] == CONVERSION["source"],
        "output_conversion_source_mutation",
    )
    rows = proof.get("records", [])
    require(
        isinstance(rows, list) and all(isinstance(r, dict) for r in rows),
        "output_verification_records_required",
    )
    ids = [r.get("sample_id") for r in rows]
    require(
        all(_text(i) for i in ids) and len(set(ids)) == len(ids),
        "duplicate_output_verification",
    )
    proofs = dict(zip(ids, rows, strict=True))
    samples = read_pinned(root, request.get("samples_file", ""), files)
    require(
        isinstance(samples, list)
        and samples
        and all(isinstance(s, dict) for s in samples),
        "reference_samples_required",
    )
    expected = {s["id"]: s for s in geometries(system, assessment, t["scan"])}
    require(
        len(samples) == len(expected)
        and {s.get("id") for s in samples} == set(expected),
        "incomplete_deformation_coverage",
    )
    require(set(proofs) == set(expected), "output_verification_coverage_mismatch")
    seen, states, spins = set(), set(), set()
    n = len(t["atom_order"])
    for s in samples:
        sid = s["id"]
        for key in ("graph_identity", "atom_order"):
            require(s.get(key) == t[key], "sample_" + key + "_mutation")
        require(s.get("graph") == t["graph"], "sample_bond_order_mutation")
        require(s.get("method_id") in methods, "undeclared_qm_method")
        require(
            s.get("charge_state") == t["charge_state"]
            and type(s.get("charge_state")) is int,
            "mixed_or_changed_charge_state",
        )
        require(s.get("spin_state") == t["spin_state"], "mixed_or_changed_spin_state")
        require(
            all(type(v) is int for v in s["spin_state"].values()),
            "spin_state_type_mutation",
        )
        states.add(s["charge_state"])
        spins.add(digest(s["spin_state"]))
        xyz = _finite(s.get("coordinates", []))
        require(xyz.shape == (n, 3), "reference_coordinates_incomplete")
        xyz[xyz == 0] = 0
        ghash = digest(xyz.tolist())
        require(ghash not in seen, "duplicate_geometry")
        seen.add(ghash)
        exp = expected[sid]
        require(
            np.allclose(xyz, exp["coordinates"], rtol=0, atol=1e-10),
            "scan_geometry_or_atom_order_mutation",
        )
        require(
            s.get("offsets") == exp["offsets"] and s.get("split") == exp["split"],
            "scan_coverage_or_split_mutation",
        )
        require(
            s.get("energy") is not None and _finite(s["energy"]).shape == (),
            "energy_required",
        )
        require(_finite(s.get("forces", [])).shape == (n, 3), "incomplete_force_vector")
        if all(v == 0 for v in exp["offsets"]):
            require(
                _finite(s.get("hessian", [])).shape == (3 * n, 3 * n),
                "central_hessian_required",
            )
        for key, shape in (("hessian", (3 * n, 3 * n)), ("forces", (n, 3))):
            if s.get(key) is not None:
                require(
                    _finite(s[key]).shape == shape,
                    "observable_shape_mismatch",
                    field=key,
                )
        if s.get("frequencies") is not None:
            v = _finite(s["frequencies"])
            require(v.ndim == 1 and v.size >= 3 * n - 6, "frequencies_incomplete")
        r = proofs[sid]
        output_files = r.get("output_files", {})
        require(
            isinstance(output_files, dict)
            and output_files
            and all(
                name in files and files[name] == sha
                for name, sha in output_files.items()
            )
            and r.get("normalized_sample_checksum") == digest(s)
            and r.get("status") == "verified_by_supplied_importer",
            "unverified_output_files",
        )
        raw = r.get("raw_observables", {})
        require(
            raw.get("units")
            == {
                "energy": "hartree",
                "gradient": "hartree/bohr",
                "hessian": "hartree/bohr^2",
            },
            "raw_observable_units_required",
        )
        ef, df = CONVERSION["energy_factor"], CONVERSION["distance_factor"]
        for raw_key, key, scale in (
            ("energy", "energy", ef),
            ("gradient", "forces", -ef / df),
            ("hessian", "hessian", ef / df**2),
        ):
            if s.get(key) is not None:
                values = _finite(raw.get(raw_key, []))
                require(
                    values.shape == np.asarray(s[key]).shape
                    and np.allclose(values * scale, s[key], rtol=1e-12, atol=1e-12),
                    "raw_to_normalized_conversion_mismatch",
                    channel=key,
                )
    require(len(states) == len(spins) == 1, "mixed_reference_states")
    # Actual residual isolation, bounds and uncertainty evidence remain necessary
    # even after a whole scan has passed integrity intake.
    return record(
        {
            "schema": "island.pcff-qm-reference-intake.v2",
            "family": "pcff",
            "target_identity": t["identity"],
            "dataset_identity": digest(request),
            "request": request,
            "files": files,
            "verification_sha256": verification_sha256,
            "sample_count": len(samples),
            "deformation_coverage": "complete",
            "status": "qm_data_available",
            "term_status": "qm_fit_required",
            "parameter_completeness": "incomplete",
            "scientific_validation": "not_established",
            "validation_scope": "file and supplied-importer attestation integrity; QM execution/authenticity/accuracy not independently verified",
            "execution_authority": False,
        }
    )


def fit_readiness(
    target,
    dataset=None,
    controls=None,
    *,
    system=None,
    assessment=None,
    root=None,
    verification_sha256=None,
):
    """Plan readiness only. Completed/validated fit values cannot be asserted here.

    Advancing to fit_ready replays the actual dataset and validates typed controls
    and a pinned isolation record. Integrity remains distinct from QM accuracy.
    """
    t = target.payload
    if dataset is None:
        return record(
            {
                "target_identity": t["identity"],
                "status": "qm_data_missing",
                "term_status": "qm_fit_required",
                "blockers": ["QM reference outputs not supplied"],
            }
        )
    d = dataset.payload
    require(
        d.get("target_identity") == t["identity"]
        and d.get("schema") == "island.pcff-qm-reference-intake.v2",
        "reference_target_mutation",
    )
    required = (
        "objective",
        "weights",
        "constraints",
        "bounds",
        "initial_values",
        "optimizer",
        "identifiable_source_bound_residual",
        "uncertainty_plan",
        "held_out_thresholds",
        "independent_reconstruction_plan",
    )
    missing = [k for k in required if not (controls or {}).get(k)]
    if controls and not missing:
        require(
            system is not None
            and assessment is not None
            and root is not None
            and verification_sha256,
            "dataset_replay_required_for_fit_ready",
        )
        replay = intake(
            system,
            assessment,
            target,
            d["request"],
            root,
            verification_sha256=verification_sha256,
        )
        require(replay.payload == d, "dataset_replay_mutation")
        count = len(t["fitting_variables"])
        require(
            controls["objective"]
            == "weighted_source_bound_class2_residual_least_squares",
            "fit_objective_not_supported",
        )
        bounds = _finite(controls["bounds"])
        initial = _finite(controls["initial_values"])
        require(
            bounds.shape == (count, 2)
            and initial.shape == (count,)
            and np.all(bounds[:, 0] < bounds[:, 1])
            and np.all(initial >= bounds[:, 0])
            and np.all(initial <= bounds[:, 1]),
            "fit_bounds_invalid",
        )
        require(
            isinstance(controls["constraints"], list)
            and set(t["physical_constraints"]) <= set(controls["constraints"]),
            "fit_constraints_incomplete",
        )
        samples = read_pinned(root, d["request"]["samples_file"], d["files"])
        weights = controls["weights"]
        require(
            isinstance(weights, dict)
            and set(weights) == {s["id"] for s in samples}
            and all(_finite(v).shape == () and float(v) > 0 for v in weights.values()),
            "fit_weights_invalid",
        )
        thresholds = controls["held_out_thresholds"]
        require(
            isinstance(thresholds, dict)
            and set(thresholds) == {"energy", "forces", "hessian"}
            and all(
                _finite(v).shape == () and float(v) > 0 for v in thresholds.values()
            ),
            "held_out_thresholds_invalid",
        )
        for key in ("optimizer", "uncertainty_plan", "independent_reconstruction_plan"):
            require(
                isinstance(controls[key], dict) and controls[key],
                "structured_fit_control_required",
                field=key,
            )
        isolation = controls["identifiable_source_bound_residual"]
        require(
            isinstance(isolation, dict)
            and isolation.get("file") in d["files"]
            and isolation.get("sha256") == d["files"][isolation["file"]],
            "residual_evidence_not_pinned",
        )
        proof = read_pinned(root, isolation["file"], d["files"])
        require(
            proof.get("target_identity") == t["identity"]
            and proof.get("samples_sha256") == d["files"][d["request"]["samples_file"]]
            and proof.get("status") == "independently_isolated"
            and proof.get("functional_form_identity")
            == t["functional_form"]["identity"]
            and proof.get("rank") == count
            and _text(proof.get("method"))
            and _text(proof.get("implementation_version")),
            "residual_identifiability_not_verified",
        )
    return record(
        {
            "target_identity": t["identity"],
            "dataset_identity": d["dataset_identity"],
            "status": "qm_data_available" if missing else "fit_ready",
            "term_status": "qm_fit_required",
            "blockers": missing,
            "controls_identity": digest(controls) if controls else None,
            "validation_scope": "preparation contract and supplied evidence integrity; no fitting or scientific validation",
            "execution_authority": False,
        }
    )
