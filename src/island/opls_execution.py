"""Explicit PDMS convention policy → checksummed static LAMMPS experiment.

No engine or output directory is created before the four-layer policy gate.
The current source audit does not authorize minimization, NVT, NPT or production.
"""

import json
import math
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.forcefields.opls_registry.api import require
from island.forcefields.pdms_opls import export_lammps
from island.forcefields.selection.inventory import canonical, digest
from island.opls_conventions import (
    ConventionPreparation,
    prepare_polymer_forcefield,
)
from island.pdms_lammps import (
    STATIC,
    LAMMPSOptions,
    _execute,
    engine_record,
    parse_static,
    pin,
    seal,
)
from island.workflows.bundle import system_data


def _validate(system, preparation, action, *, research_opt_in=False):
    require(
        type(preparation) is ConventionPreparation,
        "source_incompatible",
        "four_layer_preparation_required",
    )
    return preparation.require_execution(
        system, action, research_opt_in=research_opt_in
    )


def input_text(system, preparation, *, action="static", research_opt_in=False):
    """No fallback: every setting has a resolved, replay-validated owner."""
    backend = _validate(system, preparation, action, research_opt_in=research_opt_in)
    s = {k: r["value"] for k, r in backend["settings"].items()}
    header = (
        f"# Four-layer convention identity {preparation.identity}\n"
        f"# Execution mode {preparation.mode}; backend policy {preparation.policy}\n"
        "# Backend-only choices are not author authority. MD is blocked.\n"
        f"units {s['lammps_units']}\natom_style full\n"
        f"boundary {s['boundary_conditions']}\nnewton on\natom_modify sort 0 0\n"
        f"dielectric {s['dielectric']}\n"
        f"bond_style {s['bond_style']}\nangle_style {s['angle_style']}\n"
        f"dihedral_style {s['dihedral_style']}\nimproper_style harmonic\n"
        f"pair_style {s['pair_style']} {s['lj_inner_angstrom']} {s['pair_cutoff_angstrom']}\n"
        "pair_modify mix geometric shift no tail no "
        f"table {s['pair_table_bits']}\nspecial_bonds {s['special_bonds']}\n"
        f"kspace_style {s['kspace_style']} {s['reciprocal_tolerance']}\n"
    )
    if backend["ewald_kmax"] is not None:
        alpha = math.sqrt(-math.log(2e-7)) / 11
        header += f"kspace_modify gewald {alpha:.17g} kmax/ewald 20 20 20\n"
    header += (
        "neighbor 2 bin\nneigh_modify every 1 delay 0 check no\n"
        "read_data system.data\ninclude coefficients.in\n"
    )
    return header


def generate_lammps(
    system, preparation, directory, *, action="static", research_opt_in=False
):
    text = input_text(
        system, preparation, action=action, research_opt_in=research_opt_in
    )
    source = preparation.assignment.validate_integrity(system).validate_integrity(
        system
    )
    root = Path(directory)
    export_lammps(system, source, root)
    (root / "system.in").write_text(text)
    (root / "static.in").write_text(text + STATIC)
    return {
        n: pin(root / n)
        for n in ("system.data", "coefficients.in", "system.in", "static.in")
    }


def run_static(system, preparation, directory, *, engine, research_opt_in=False):
    backend = _validate(system, preparation, "static", research_opt_in=research_opt_in)
    require(
        type(engine) is LAMMPSOptions
        and engine.expected_sha256 == backend["engine_sha256"],
        "source_incompatible",
        "backend_policy_executable_mismatch",
        required_sha256=backend["engine_sha256"],
    )
    info = engine_record(engine)
    require(
        "KSPACE" in info["help_text"] and "double" in info["help_text"],
        "asset_unavailable",
        "pinned_KSPACE_double_engine_required",
    )
    root = Path(directory)
    require(not root.exists(), "source_incompatible", "output_directory_exists")
    # generator gates again before output creation, preventing a bypass via files.
    files = generate_lammps(
        system, preparation, root / "lammps", research_opt_in=research_opt_in
    )
    source = preparation.assignment.validate_integrity(system).validate_integrity(
        system
    )
    code = _execute(engine, root / "lammps", "static.in")
    require(code == 0, "asset_unavailable", "static_process_failed", returncode=code)
    observed = parse_static(root / "lammps", source)
    result = {
        "schema": "island.j50-static.v1",
        "convention_identity": preparation.identity,
        "family": "opls",
        "selected_profile": preparation.payload["selected_profile"],
        "charge_origin": preparation.payload["charge_origin"],
        "domain": preparation.payload["domain"],
        "execution_mode": preparation.mode,
        "backend_policy": preparation.policy,
        "deterministic_inputs": files,
        "engine": info,
        "observed": observed,
        "source_conventions_resolved": not preparation.payload["unresolved_settings"],
        "MD_executed": False,
        "NPT_executed": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    result["identity"] = digest(result)
    (root / "result.json").write_text(canonical(result))
    (root / "system.json").write_text(pack(system_data(system)))
    (root / "conventions.json").write_text(preparation.json_text)
    (root / "canonical-assignment.json").write_text(preparation.assignment.json_text)
    return {**result, "manifest_sha256": seal(root)}


def run_md(
    system, preparation, directory, *, engine=None, research_opt_in=False, action="nvt"
):
    """Opt-in is necessary, never sufficient. No current MD adapter is authorized."""
    require(
        action in ("minimize", "nvt", "npt"),
        "source_incompatible",
        "MD_action_required",
    )
    _validate(system, preparation, action, research_opt_in=research_opt_in)
    raise AssertionError("No source-reconciled PDMS MD adapter registered")


def _reconstruct(directory, *, manifest_sha256, assets):
    """Relocatable, nonexecuting replay; trust an independently pinned manifest."""
    from island.forcefields.opls_registry import use_opls_assets

    root = Path(directory).resolve()
    require(
        pin(root / "manifest.json") == manifest_sha256,
        "source_incompatible",
        "manifest_hash_mismatch",
    )
    manifest = unpack((root / "manifest.json").read_text())
    # J47 seal format; authenticate all bytes before interpreting any record.
    files = manifest["files"]
    actual = {
        str(p.relative_to(root))
        for p in root.rglob("*")
        if p.is_file() and p != root / "manifest.json"
    }
    require(actual == set(files), "source_incompatible", "workflow_membership_mutation")
    for name, expected in files.items():
        path = root / name
        require(
            not Path(name).is_absolute()
            and ".." not in Path(name).parts
            and not path.is_symlink()
            and path.resolve().is_relative_to(root)
            and pin(path) == expected,
            "source_incompatible",
            "workflow_file_mutation",
            file=name,
        )
    from island.pdms import restore_pdms_system

    system = restore_pdms_system(unpack((root / "system.json").read_text()))
    p = json.loads((root / "conventions.json").read_text())
    with use_opls_assets(assets):
        preparation = prepare_polymer_forcefield(
            system,
            forcefield="opls",
            execution_mode=p["execution_mode"],
            backend_policy=p["backend_policy"],
        )
    require(
        preparation.json_text == (root / "conventions.json").read_text()
        and preparation.assignment.json_text
        == (root / "canonical-assignment.json").read_text(),
        "source_incompatible",
        "reconstructed_convention_or_assignment_mutation",
    )
    require(
        (root / "lammps/system.in").read_text() == input_text(system, preparation),
        "source_incompatible",
        "reconstructed_input_mutation",
    )
    return system, preparation


def reconstruct(directory, *, manifest_sha256, assets):
    """Classified relocation boundary; malformed records never reach an engine."""
    from island.exceptions import IslandError
    from island.forcefields.opls_registry import OPLSResolutionError

    try:
        return _reconstruct(directory, manifest_sha256=manifest_sha256, assets=assets)
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
            [{"stage": "convention_reconstruction_failed", "reason": str(error)}],
        ) from error


__all__ = [
    "LAMMPSOptions",
    "generate_lammps",
    "input_text",
    "prepare_polymer_forcefield",
    "reconstruct",
    "run_md",
    "run_static",
]
