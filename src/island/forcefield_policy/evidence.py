"""Optional static comparison observations and relocatable claim snapshots."""

import json
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.forcefields.opls_registry.api import require
from island.forcefields.selection.inventory import canonical, digest

from . import _signed, _update, prepare_polymer_forcefield

BOUNDS = {"energy": 3e-4, "force": 3e-4, "finite_difference": 1e-4}
# Only these observation fields may change after preparation. Parameter and
# convention completeness is always regenerated from live source files.
OBSERVATION_FIELDS = (
    "lammps_export_status",
    "lammps_runtime_status",
    "openmm_build_status",
    "openmm_runtime_status",
    "cross_engine_equivalence",
    "cross_engine_diagnostics",
    "numerical_evidence",
    "lammps_export",
    "openmm_static",
    "lammps_static",
    "openmm_build_diagnostic",
    "openmm_runtime_diagnostic",
    "lammps_runtime_diagnostic",
)


def record_static_comparison(system, prepared, path, *, expected_sha256):
    """Explicit external observation import. It cannot grant equivalence authority.

    Failure is reported independently. Passing these bounded static observations
    alone does not prove trajectory, reciprocal, pressure or scientific agreement.
    Bounds are the unchanged J49 component bounds. Sources/charges are untouched.
    """
    p = prepared.validate_integrity(system)
    path = Path(path)
    raw = path.read_bytes()
    require(
        sha256(raw).hexdigest() == expected_sha256,
        "source_incompatible",
        "diagnostic_hash",
    )
    d = json.loads(raw)
    require(
        d["assignment_identity"] == p["assignment_identity"]
        and d["input_identity"] == prepared._system_identity,
        "source_incompatible",
        "diagnostic_input_identity",
    )
    require(
        d["reference_engine"] == "OpenMM"
        and d["observed_engine"] == "LAMMPS"
        and d["scope"] == "static_components"
        and bool(d["measurements"]),
        "source_incompatible",
        "unsupported_diagnostic_scope",
    )
    results = []
    for row in d["measurements"]:
        require(
            row["kind"] in BOUNDS and row["bound"] == BOUNDS[row["kind"]],
            "source_incompatible",
            "diagnostic_bound_mutation",
        )
        ref, obs = (
            np.asarray(row["reference"], dtype=float),
            np.asarray(row["observed"], dtype=float),
        )
        require(
            ref.shape == obs.shape
            and ref.size > 0
            and np.isfinite(ref).all()
            and np.isfinite(obs).all(),
            "source_incompatible",
            "diagnostic_shape_or_nonfinite",
        )
        error = float(np.max(np.abs(ref - obs)))
        results.append(
            {
                **row,
                "max_error": error,
                "status": "passed" if error <= row["bound"] else "failed",
            }
        )
    value = {
        "evidence_sha256": expected_sha256,
        "scope": d["scope"],
        "origin": "explicit_external_observations; not independently executed by import",
        "observations": results,
        "broad_equivalence_claim": False,
        "status": "failed"
        if any(r["status"] == "failed" for r in results)
        else "passed",
    }
    status = "failed" if value["status"] == "failed" else p["cross_engine_equivalence"]
    result = _update(
        prepared,
        system,
        cross_engine_equivalence=status,
        cross_engine_diagnostics=p["cross_engine_diagnostics"] + [value],
        numerical_evidence=p["numerical_evidence"]
        + [{"sha256": expected_sha256, "record": d}],
    )
    return replace(result, _pins=prepared._pins + ((str(path), expected_sha256),))


def save(system, prepared, directory):
    from island.charge_references.records import pack
    from island.workflows.bundle import system_data

    prepared.validate_integrity(system)
    root = Path(directory)
    require(not root.exists(), "source_incompatible", "output_directory_exists")
    root.mkdir(parents=True)
    (root / "report.json").write_text(prepared.json_text)
    (root / "system.json").write_text(pack(system_data(system)))
    value = {
        "schema": "island.parameter-policy-bundle.v1",
        "files": {
            n: sha256((root / n).read_bytes()).hexdigest()
            for n in ("report.json", "system.json")
        },
    }
    (root / "manifest.json").write_text(canonical({**value, "identity": digest(value)}))
    return sha256((root / "manifest.json").read_bytes()).hexdigest()


def reconstruct(
    directory, *, manifest_sha256, assets=None, sources=None, **legacy_options
):
    from island.charge_references.records import unpack
    from island.pdms import restore_pdms_system

    root = Path(directory).resolve()
    require(
        sha256((root / "manifest.json").read_bytes()).hexdigest() == manifest_sha256,
        "source_incompatible",
        "bundle_manifest_mutation",
    )
    manifest = json.loads((root / "manifest.json").read_text())
    require(
        manifest["schema"] == "island.parameter-policy-bundle.v1"
        and set(manifest["files"]) == {"system.json", "report.json"}
        and {p.name for p in root.iterdir()}
        == {"system.json", "report.json", "manifest.json"},
        "source_incompatible",
        "bundle_membership",
    )
    for n, pin in manifest["files"].items():
        require(
            not (root / n).is_symlink()
            and sha256((root / n).read_bytes()).hexdigest() == pin,
            "source_incompatible",
            "bundle_file_mutation",
        )
    from island.workflows.bundle import system_from

    payload = unpack((root / "system.json").read_text())
    system = (
        restore_pdms_system(payload)
        if payload["box"] is not None
        else system_from(payload)
    )
    saved = json.loads((root / "report.json").read_text())
    require(
        saved["identity"]
        == digest({k: v for k, v in saved.items() if k != "identity"}),
        "source_incompatible",
        "bundle_report_checksum",
    )
    current = prepare_polymer_forcefield(
        system, saved["family"], assets=assets, sources=sources, **legacy_options
    )
    immutable = lambda p: {
        k: v for k, v in p.items() if k not in (*OBSERVATION_FIELDS, "identity")
    }
    require(
        immutable(current.payload) == immutable(saved),
        "source_incompatible",
        "reconstructed_parameter_or_convention_mutation",
    )
    for record in saved["numerical_evidence"]:
        # Raw observations are included in the independently pinned report.
        raw = record["record"]
        require(
            raw["assignment_identity"] == current.payload["assignment_identity"]
            and raw["input_identity"] == current._system_identity,
            "source_incompatible",
            "reconstructed_diagnostic_binding",
        )
    # This preserves observations only, never a live engine or MD permission.
    p = {k: v for k, v in saved.items() if k != "identity"}
    return system, _signed(
        p,
        current.assignment,
        system,
        current._audit,
        current._pins
        + tuple((str(root / n), pin) for n, pin in manifest["files"].items()),
    )
