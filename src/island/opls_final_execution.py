"""J52 fail-closed execution facade with relocatable final audit identity."""

import json
from pathlib import Path

from island.forcefields.opls_registry import OPLSResolutionError
from island.forcefields.opls_registry.api import require
from island.forcefields.selection.inventory import canonical
from island.opls_conventions import signed
from island.opls_final_resolution import (
    FinalConventionPreparation,
    prepare_polymer_forcefield,
)
from island.opls_resolved_execution import LAMMPSOptions
from island.opls_resolved_execution import generate_lammps as _generate
from island.opls_resolved_execution import reconstruct as _reconstruct
from island.opls_resolved_execution import run_static as _run_static
from island.pdms_lammps import seal


def _gate(system, preparation, action="static", research_opt_in=False):
    require(
        type(preparation) is FinalConventionPreparation,
        "source_incompatible",
        "final_convention_resolution_required",
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
    (root / "final-resolution.json").write_text(preparation.json_text)
    value = signed(
        {
            "schema": "island.j52-static.v1",
            "final_resolution_identity": preparation.identity,
            "J51_static_identity": result["identity"],
            "diagnostic_only": True,
            "MD_authorized": False,
        }
    )
    (root / "final-result.json").write_text(canonical(value))
    return {**result, **value, "manifest_sha256": seal(root)}


def run_md(
    system, preparation, directory, *, engine=None, action="nvt", research_opt_in=False
):
    require(
        action in ("minimize", "nvt", "npt"),
        "source_incompatible",
        "MD_action_required",
    )
    _gate(system, preparation, action, research_opt_in)
    require(False, "scientifically_unqualified", "no_source_reconciled_MD_adapter")


def reconstruct(directory, *, manifest_sha256, assets, sources):
    system, prior = _reconstruct(
        directory, manifest_sha256=manifest_sha256, assets=assets, sources=sources
    )
    try:
        raw = (Path(directory) / "final-resolution.json").read_text()
        value = json.loads(raw)
        prepared = prepare_polymer_forcefield(
            system,
            forcefield="opls",
            assets=assets,
            sources=sources,
            execution_mode=value["execution_mode"],
            backend_policy=value["backend_policy"],
        )
        require(
            prepared.prior.json_text == prior.json_text and prepared.json_text == raw,
            "source_incompatible",
            "reconstructed_final_resolution_mutation",
        )
        prepared.validate_integrity(system)
        return system, prepared
    except (OSError, ValueError, KeyError, TypeError) as error:
        if isinstance(error, OPLSResolutionError):
            raise
        raise OPLSResolutionError(
            "asset_unavailable"
            if isinstance(error, OSError)
            else "source_incompatible",
            [{"stage": "final_resolution_reconstruction_failed", "reason": str(error)}],
        ) from error


__all__ = [
    "LAMMPSOptions",
    "generate_lammps",
    "prepare_polymer_forcefield",
    "reconstruct",
    "run_md",
    "run_static",
]
