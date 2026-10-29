"""Canonical OPLS assignment → unchanged, source-bound J47 static workflow."""

from pathlib import Path

from island.charge_references.records import unpack
from island.exceptions import IslandError
from island.forcefields.opls_registry import ForceFieldRequest, prepare_forcefield
from island.forcefields.opls_registry.api import require
from island.forcefields.pdms_opls import PROFILE
from island.pdms_lammps import LAMMPSOptions, pin, seal
from island.pdms_lammps import reconstruct as _reconstruct
from island.pdms_lammps import run_static as _run_static


def run_static(system, preparation, directory, *, engine):
    backend = preparation.validate_integrity(system)
    require(
        preparation.payload["selected_profile"] == PROFILE,
        "source_incompatible",
        "no_lammps_adapter_for_selected_profile",
    )
    result = _run_static(system, backend, directory, engine=engine)
    root = Path(directory)
    (root / "canonical-assignment.json").write_text(preparation.json_text)
    return {
        **result,
        "canonical_assignment_identity": preparation.identity,
        "manifest_sha256": seal(root),
    }


def reconstruct(directory, *, manifest_sha256, assets, engine):
    """Re-resolve canonical chemistry and rehash raw source in a fresh process."""
    from island.pdms import restore_pdms_system

    root = Path(directory)
    # Read only inert JSON before the independent manifest is checked by J47.
    try:
        require(
            pin(root / "manifest.json") == manifest_sha256,
            "source_incompatible",
            "manifest_hash_mismatch",
        )
        manifest = unpack((root / "manifest.json").read_text())
        for name in ("system.json", "canonical-assignment.json"):
            require(
                pin(root / name) == manifest["files"][name],
                "source_incompatible",
                "workflow_input_mutation",
                file=name,
            )
        system = restore_pdms_system(unpack((root / "system.json").read_text()))
        prepared = prepare_forcefield(system, ForceFieldRequest("opls"), assets=assets)
        backend = prepared.validate_integrity(system)
        restored, _, result = _reconstruct(
            root,
            manifest_sha256=manifest_sha256,
            request=backend.request,
            engine=engine,
        )
        require(
            (root / "canonical-assignment.json").read_text() == prepared.json_text,
            "source_incompatible",
            "canonical_assignment_mutation",
        )
        prepared.validate_integrity(restored)
        return restored, prepared, result
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        IslandError,
    ) as error:
        from island.forcefields.opls_registry import OPLSResolutionError

        if isinstance(error, OPLSResolutionError):
            raise
        raise OPLSResolutionError(
            "asset_unavailable"
            if isinstance(error, OSError)
            else "source_incompatible",
            [{"stage": "workflow_read_failed", "reason": str(error)}],
        ) from error


def run_md(system, preparation, directory, **kwargs):
    """No directory, engine, checkpoint or trajectory before static qualification."""
    preparation.require_md(system)
    raise AssertionError("No qualified canonical OPLS LAMMPS MD adapter registered")


__all__ = ["LAMMPSOptions", "reconstruct", "run_md", "run_static"]
