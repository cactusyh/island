"""J51 scoped resolution boundary, retaining historical numerical workflows."""

import json
from pathlib import Path

from island.exceptions import IslandError
from island.forcefields.opls_registry import OPLSResolutionError
from island.forcefields.opls_registry.api import require
from island.forcefields.selection.inventory import canonical, digest
from island.opls_execution import LAMMPSOptions
from island.opls_execution import generate_lammps as _generate
from island.opls_execution import reconstruct as _reconstruct
from island.opls_execution import run_static as _run_static
from island.opls_resolution import (
    ResolutionPreparation,
    prepare_polymer_forcefield,
)
from island.pdms_lammps import seal


def _gate(system, preparation, action="static", research_opt_in=False):
    require(
        type(preparation) is ResolutionPreparation,
        "source_incompatible",
        "scoped_convention_resolution_required",
    )
    preparation.require_execution(system, action, research_opt_in=research_opt_in)


def generate_lammps(
    system, preparation, directory, *, action="static", research_opt_in=False
):
    _gate(system, preparation, action, research_opt_in)
    return _generate(
        system,
        preparation.prior,
        directory,
        action=action,
        research_opt_in=research_opt_in,
    )


def run_static(system, preparation, directory, *, engine, research_opt_in=False):
    _gate(system, preparation, research_opt_in=research_opt_in)
    result = _run_static(
        system,
        preparation.prior,
        directory,
        engine=engine,
        research_opt_in=research_opt_in,
    )
    root = Path(directory)
    (root / "resolution.json").write_text(preparation.json_text)
    value = {
        "schema": "island.j51-static.v1",
        "resolution_identity": preparation.identity,
        "four_layer_static_identity": result["identity"],
        "diagnostic_only": True,
        "MD_authorized": False,
    }
    value["identity"] = digest(value)
    (root / "resolution-result.json").write_text(canonical(value))
    return {**result, **value, "manifest_sha256": seal(root)}


def run_md(
    system, preparation, directory, *, engine=None, action="nvt", research_opt_in=False
):
    """Reject before engine access, output, checkpoint or trajectory creation."""
    require(
        action in ("minimize", "nvt", "npt"),
        "source_incompatible",
        "MD_action_required",
    )
    _gate(system, preparation, action, research_opt_in)
    require(False, "scientifically_unqualified", "no_source_reconciled_MD_adapter")


def reconstruct(directory, *, manifest_sha256, assets, sources):
    """Authenticate every byte, rehash both PDFs and ZIP, then reconstruct in place."""
    try:
        system, prior = _reconstruct(
            directory, manifest_sha256=manifest_sha256, assets=assets
        )
        value = json.loads((Path(directory) / "resolution.json").read_text())
        prepared = prepare_polymer_forcefield(
            system,
            forcefield="opls",
            assets=assets,
            sources=sources,
            execution_mode=value["execution_mode"],
            backend_policy=value["backend_policy"],
        )
        require(
            prepared.prior.json_text == prior.json_text
            and prepared.json_text == (Path(directory) / "resolution.json").read_text(),
            "source_incompatible",
            "reconstructed_resolution_mutation",
        )
        prepared.validate_integrity(system)
        return system, prepared
    except (
        OSError,
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        IslandError,
    ) as error:
        if isinstance(error, OPLSResolutionError):
            raise
        raise OPLSResolutionError(
            "asset_unavailable"
            if isinstance(error, OSError)
            else "source_incompatible",
            [{"stage": "resolution_reconstruction_failed", "reason": str(error)}],
        ) from error


__all__ = [
    "LAMMPSOptions",
    "generate_lammps",
    "prepare_polymer_forcefield",
    "reconstruct",
    "run_md",
    "run_static",
]
