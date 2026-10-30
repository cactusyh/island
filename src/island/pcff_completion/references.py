"""Strict reference intake independent of parameter fitting or QM execution.

Normalized data must be cross-bound to an independently pinned verification
record and its raw outputs/importer. Hash checks establish file integrity, not
that an engine calculation or scientific validation actually occurred.
"""

from pathlib import Path

import numpy as np

from . import (
    _finite,
    _verified,
    code_identity,
    digest,
    inventory,
    read_pinned,
    require,
    signed,
)

BINDINGS = (
    "graph_identity",
    "typing_identity",
    "charge_identity",
    "atom_order",
    "charge_vector",
    "source_hashes",
    "domain",
)
METHOD_FIELDS = (
    "method",
    "basis",
    "dispersion",
    "software",
    "version",
    "convergence",
    "multiplicity",
)


def reference_intake(system, assessment, target, request, root, *, verification_sha256):
    p = _verified(system, assessment)
    t = target.validate(system, assessment)
    require(
        request.get("schema") == "island.pcff-qm-reference.v1",
        "reference_schema_required",
    )
    require(request.get("family") == "pcff", "wrong_profile_family")
    require(
        request.get("target_identity") == t["identity"], "reference_target_mutation"
    )
    for k in BINDINGS:
        require(request.get(k) == p[k], k + "_mutation")
    require(request.get("units") == t["units"], "reference_units_missing_or_changed")
    files = inventory(root)
    require(files == request.get("files"), "reference_inventory_or_hash_mutation")
    license_file = request.get("license_file")
    require(
        license_file in files and bool((Path(root) / license_file).read_text().strip()),
        "reference_license_required",
    )
    require(
        request.get("permission") in ("licensed", "permission_granted"),
        "reference_permission_required",
    )
    methods = request.get("methods", {})
    require(isinstance(methods, dict) and bool(methods), "reference_method_required")
    for m in methods.values():
        require(
            all(m.get(k) is not None and m[k] != "" for k in METHOD_FIELDS),
            "reference_method_incomplete",
        )
    require(
        len(methods) == 1 or bool(request.get("mixed_method_protocol")),
        "undeclared_mixed_qm_methods",
    )
    verification_file = request.get("output_verification_file")
    require(
        verification_file in files and files[verification_file] == verification_sha256,
        "unverified_output_files",
    )
    verified = read_pinned(root, verification_file, files)
    importer = verified.get("importer_file")
    require(
        importer in files
        and verified.get("importer_sha256") == files[importer]
        and bool(verified.get("implementation_version")),
        "output_importer_not_pinned",
    )
    records = verified.get("records", [])
    require(
        len({r["sample_id"] for r in records}) == len(records),
        "duplicate_output_verification",
    )
    proofs = {r["sample_id"]: r for r in records}
    samples = read_pinned(root, request.get("samples_file"), files)
    require(isinstance(samples, list) and bool(samples), "reference_samples_required")
    require(len(samples) == len(proofs), "output_verification_coverage_mismatch")
    geometries, ids, states = set(), set(), set()
    n = len(p["atom_order"])
    for sample in samples:
        sid = sample.get("id")
        require(
            isinstance(sid, str) and sid and sid not in ids, "duplicate_sample_identity"
        )
        ids.add(sid)
        require(
            sample.get("graph_identity") == p["graph_identity"],
            "sample_graph_identity_required",
        )
        require(
            sample.get("atom_order") == p["atom_order"], "sample_atom_order_mutation"
        )
        require(
            sample.get("bond_order_identity") == p["graph_identity"],
            "sample_bond_order_mutation",
        )
        xyz = _finite(sample.get("coordinates", []))
        require(xyz.shape == (n, 3), "reference_coordinates_incomplete")
        # Numeric normalization makes 0 and 0.0 duplicates; signed zero too.
        xyz[xyz == 0] = 0
        geometry = digest(xyz.tolist())
        require(geometry not in geometries, "duplicate_geometry")
        geometries.add(geometry)
        state = sample.get("charge_state")
        require(
            isinstance(state, int) and not isinstance(state, bool),
            "charge_state_required",
        )
        states.add(state)
        require(sample.get("method_id") in methods, "undeclared_qm_method")
        require(
            sample.get("split") in ("train", "held_out"), "reference_split_required"
        )
        require(
            any(
                sample.get(k) is not None
                for k in ("energy", "forces", "hessian", "frequencies")
            ),
            "reference_observable_required",
        )
        for key, shape in (
            ("energy", ()),
            ("forces", (n, 3)),
            ("hessian", (3 * n, 3 * n)),
        ):
            if sample.get(key) is not None:
                require(
                    _finite(sample[key]).shape == shape,
                    "reference_observable_incomplete",
                    channel=key,
                )
        if sample.get("frequencies") is not None:
            f = _finite(sample["frequencies"])
            require(f.ndim == 1 and f.size > 0, "reference_frequencies_incomplete")
        proof = proofs.get(sid, {})
        output = proof.get("output_file")
        require(
            output in files
            and proof.get("output_sha256") == files[output]
            and proof.get("normalized_sample_checksum") == digest(sample)
            and proof.get("status") == "verified",
            "unverified_output_files",
        )
    require(len(states) == 1, "mixed_charge_states")
    charge = sum(float(v) for v in p["charge_vector"].values())
    require(
        np.isclose(charge, next(iter(states)), rtol=0, atol=1e-10),
        "charge_state_mapping_mismatch",
    )
    require(any(s["split"] == "held_out" for s in samples), "held_out_data_required")
    return signed(
        {
            "schema": "island.pcff-qm-reference-intake.v1",
            "target_identity": t["identity"],
            "request": request,
            "target": t,
            "code_identity": code_identity(),
            "dataset_identity": digest(request),
            "files": files,
            "verification_sha256": verification_sha256,
            "methods": methods,
            "license_sha256": files[license_file],
            "sample_count": len(samples),
            "graph_identity": p["graph_identity"],
            "typing_identity": p["typing_identity"],
            "charge_identity": p["charge_identity"],
            "status": "reference_integrity_passed",
            "fit_status": "qm_fit_required",
            "parameter_completeness": "incomplete",
            "validation_scope": "user-supplied output/importer attestation; file integrity only",
            "scientific_validation": "not_established",
            "execution_authority": False,
        }
    )
