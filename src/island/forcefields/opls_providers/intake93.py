"""J93 metadata gate for the real PE5 bundle; no graph inference or execution."""

import hashlib
from pathlib import Path

from island.forcefields.assignment_core.schema import is_hash, require
from island.forcefields.bundle_adapters.schema import loads

from .intake91 import outside_git
from .registry import _result, external_file, fence, pinned
from .schema import FLAGS

SOURCE_DIR = Path(
    "/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-opls-bundles/pe5"
)
DATA_FILE = "LIG.lammps.lmp"
SETTINGS_FILE = "LIG.lammps.in"
METADATA_FILES = (
    "request.json",
    "manifest.json",
    "provider_metadata.json",
    "independent_audit.json",
)
STRUCTURE_SUFFIXES = (".smi", ".smiles", ".mol", ".sdf", ".pdb")
PRIOR_RECEIPT = "qualification-receipt.json"
PRIOR_BLOCKERS = "missing-terms-and-blockers.json"
PRIOR_RESULT = "result.json"
PRIOR_HASHES = {
    PRIOR_RECEIPT: "92889699a7fa34eec502bdc82bba008de53aa5987dd69b7620e5054917aaac88",
    PRIOR_BLOCKERS: "acc551511b5d6260a69e927ca836db2f678044815b2fd9ac47535531d3c61f4c",
    PRIOR_RESULT: "53c80d4fa4f1cdf141e6d8d75ce930933486ff91db64f4865223ab94929ff6f7",
}
READY = {
    **FLAGS,
    "automatic_assignment": False,
    "accepted_coefficients": 0,
    "runtime_registration": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "chemical_family_transfer": "not_established",
}


def _sha256(raw):
    return hashlib.sha256(raw).hexdigest()


def _prior(prior_root):
    root = Path(prior_root)
    raw = {}
    for name, expected in PRIOR_HASHES.items():
        value = (root / name).read_bytes()
        require(_sha256(value) == expected, "J92_historical_artifact_mutation:" + name)
        raw[name] = loads(value)
    require(
        raw[PRIOR_RESULT]["decision"] == "parameter_bundle_incomplete",
        "J92_decision_must_remain_incomplete",
    )
    return raw


def _inventory(root):
    files = []
    for path in sorted(root.rglob("*")):
        require(not path.is_symlink(), "source_symlink_rejected")
        relative = path.relative_to(root)
        if path.is_file() and ".ipynb_checkpoints" not in relative.parts:
            files.append(relative.as_posix())
    return files


def intake(
    source_dir=SOURCE_DIR,
    *,
    data_sha256=None,
    settings_sha256=None,
    prior_root=Path("docs/evidence/phase_4j92_real_pe5"),
):
    fence()
    report = {
        "schema": "island.j93.pe5-metadata-gate.v1",
        "decision": "source_bundle_required",
        "source_hashes": {},
        "source_inventory": [],
        "provider": None,
        "provider_version": None,
        "exact_graph_identity": None,
        "formal_charge": None,
        "terminal_caps": None,
        "pair14_count": 0,
        "exclusion_count": 0,
        "missing_terms": [],
        "convention_blockers": [],
        "blockers": [],
        "metadata_presence": {},
        "original_structure_files": [],
        "J92_receipt_preserved": False,
        "source_files_parsed_this_step": 0,
        "scientific_execution": {
            "LAMMPS_commands": 0,
            "MD_steps": 0,
            "QM_executions": 0,
            "energy_force_evaluations": 0,
            "minimizations": 0,
            "physical_validations": 0,
        },
        **READY,
    }
    try:
        root = outside_git(source_dir)
        if not root.is_dir():
            report["blockers"] = ["explicit_external_source_directory_missing"]
            return _result(report)
        require(is_hash(data_sha256), "pinned_LAMMPS_data_sha256_required")
        require(is_hash(settings_sha256), "pinned_LAMMPS_settings_sha256_required")
        data_raw = pinned(external_file(root, DATA_FILE), data_sha256)
        settings_raw = pinned(external_file(root, SETTINGS_FILE), settings_sha256)
        del data_raw, settings_raw
        report["source_hashes"] = {
            DATA_FILE: data_sha256,
            SETTINGS_FILE: settings_sha256,
        }
        report["source_inventory"] = _inventory(root)
        prior = _prior(prior_root)
        prior_result = prior[PRIOR_RESULT]
        report.update(
            decision="parameter_bundle_incomplete",
            provider=prior_result["provider"],
            provider_version=prior_result["provider_version"],
            missing_terms=prior_result["missing_terms"],
            convention_blockers=prior_result["conventions"].get("blockers", []),
            J92_receipt_preserved=True,
            J92_artifact_hashes=dict(PRIOR_HASHES),
        )
        metadata = {name: (root / name).is_file() for name in METADATA_FILES}
        structures = [
            name
            for name in report["source_inventory"]
            if Path(name).suffix.lower() in STRUCTURE_SUFFIXES
            and ".ipynb_checkpoints" not in Path(name).parts
        ]
        report["metadata_presence"] = metadata
        report["original_structure_files"] = structures
        missing = [
            name + "_missing" for name, present in metadata.items() if not present
        ]
        if not structures:
            missing.extend(
                (
                    "original_input_structure_missing",
                    "terminal_cap_information_missing",
                    "input_PSMILES_missing",
                    "exact_graph_identity_unavailable",
                    "formal_charge_unverified",
                )
            )
        report["blockers"] = sorted(
            set(
                missing
                + [
                    "special_bonds_expansion_blocked_without_exact_graph",
                    "independent_audit_required_for_audited_exact_graph",
                ]
            )
        )
        # This gate intentionally stops before parsing or reconstructing a graph
        # whenever an authoritative structure/cap contract is absent.
        if missing:
            return _result(report)
        report["blockers"] = [
            "complete_metadata_detected_run_manifest_bound_J92_qualification"
        ]
        return _result(report)
    except (ValueError, TypeError, KeyError, OSError, UnicodeError) as error:
        report["decision"] = "parameter_bundle_incomplete"
        report["blockers"] = [
            str(error)
            if type(error) is ValueError
            else type(error).__name__ + "_invalid_input"
        ]
        return _result(report)
