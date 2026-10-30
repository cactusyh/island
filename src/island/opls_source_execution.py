"""Side-effect-free source gates and relocatable, non-executable preparations."""

import json
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.forcefields.opls_registry.api import require
from island.forcefields.selection.inventory import canonical
from island.opls_conventions import signed
from island.opls_final_execution import generate_lammps as _generate
from island.opls_final_execution import run_static as _static
from island.opls_source_intake import IntakePreparation, prepare_polymer_forcefield
from island.pdms import restore_pdms_system
from island.pdms_lammps import pin
from island.workflows.bundle import system_data


def _gate(system, preparation, action="static", research_opt_in=False):
    require(
        type(preparation) is IntakePreparation,
        "source_incompatible",
        "source_intake_preparation_required",
    )
    preparation.require_execution(system, action, research_opt_in=research_opt_in)


def generate_lammps(system, preparation, directory):
    _gate(system, preparation)
    return _generate(system, preparation.prior, directory)


def run_static(system, preparation, directory, *, engine):
    _gate(system, preparation)
    return _static(system, preparation.prior, directory, engine=engine)


def run_md(
    system, preparation, directory, *, engine=None, action="nvt", research_opt_in=False
):
    _gate(system, preparation, action, research_opt_in)
    require(False, "scientifically_unqualified", "no_new_qualified_MD_adapter")


def save_preparation(system, preparation, directory):
    """Serialize a checksummed assignment report, not a checkpoint or MD bundle."""
    preparation.validate_integrity(system)
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=False)
    (root / "preparation.json").write_text(preparation.json_text)
    (root / "system.json").write_text(pack(system_data(system)))
    manifest = signed(
        {
            "schema": "island.opls-intake-preparation-bundle.v1",
            "files": {
                name: pin(root / name) for name in ("preparation.json", "system.json")
            },
            "preparation_identity": preparation.identity,
            "MD_authorized": False,
        }
    )
    (root / "manifest.json").write_text(canonical(manifest))
    return pin(root / "manifest.json")


def reconstruct(directory, *, manifest_sha256, **kwargs):
    root = Path(directory)
    require(
        not root.is_symlink()
        and not (root / "manifest.json").is_symlink()
        and pin(root / "manifest.json") == manifest_sha256,
        "source_incompatible",
        "intake_bundle_manifest_mutation",
    )
    manifest = json.loads((root / "manifest.json").read_text())
    require(
        manifest == signed({k: v for k, v in manifest.items() if k != "identity"})
        and set(manifest["files"]) == {"preparation.json", "system.json"}
        and not manifest["MD_authorized"],
        "source_incompatible",
        "invalid_intake_bundle",
    )
    for name, expected in manifest["files"].items():
        require(
            not (root / name).is_symlink() and pin(root / name) == expected,
            "source_incompatible",
            "intake_bundle_file_mutation",
            path=name,
        )
    system = restore_pdms_system(unpack((root / "system.json").read_text()))
    raw = (root / "preparation.json").read_text()
    report = json.loads(raw)
    prepared = prepare_polymer_forcefield(
        system,
        forcefield="opls",
        execution_mode=report["execution_mode"],
        backend_policy=report["backend_policy"],
        **kwargs,
    )
    require(
        prepared.json_text == raw
        and prepared.identity == manifest["preparation_identity"],
        "source_incompatible",
        "intake_reconstruction_mutation",
    )
    prepared.validate_integrity(system)
    return system, prepared
