"""Bounded native AM1 precision investigation, not a charge assignment.

The pinned AmberTools 24.8 CLI loses precision before antechamber reads the
charges. This additive diagnostic never publishes a vector, even for a component
whose rounded total happens to pass. Supporting a new native representation
requires independently reviewed provenance and extraction semantics.
"""

from pathlib import Path

from island.forcefields.ambertools.charge_audit import audit_component, require
from island.forcefields.ambertools.runtime import checksum, digest, signed

SCHEMA = 'island.native-am1bcc-precision-diagnostic.v1'
STUDY_SHA256 = 'f14edf758e03a41aa8a8739bca5e44a719cec5397f0fa708e5aee88828fa7374'
STUDY_ID = 'afa85536755f098b68a14f8f218e8fd4036980b047b2ed09f4c9a11419f99df0'
RUNTIME_ID = 'cdb25fdb28877c64f69f8fbc03e0cfa854b742b706de15d9f6fa5db199ee42ab'


def validate_study(study, path):
    """Anchor this specific investigation, rather than trusting re-signed claims."""
    require(checksum(path) == STUDY_SHA256, 'Unreviewed native precision study')
    require(study.get('identity') == digest({k: v for k, v in study.items() if k != 'identity'}), 'Changed study identity')
    for name, sha in study['artifacts'].items():
        target = Path(path).parent / name
        require(target.resolve().is_relative_to(Path(path).parent.resolve()), 'Unsafe study artifact path')
        require(checksum(target) == sha, 'Changed source/probe artifact: ' + name)
    require(study['runtime_identity'] == RUNTIME_ID, 'Wrong study runtime')
    return study


def inspect_component(system, directory, *, final_graph_identity, runtime, study, case):
    """Replay native outputs read-only against an exact frozen investigation.

    This is deliberately not a general fresh-preparation success path. Hashes
    bind all intermediate bytes, including zero-padding, atom order and inputs.
    Native text charges below are diagnostic observations only.
    """
    require(study.get('identity') == STUDY_ID and digest({k: v for k, v in study.items() if k != 'identity'}) == STUDY_ID, 'Unreviewed or changed native precision study')
    require(runtime.get('identity') == RUNTIME_ID, 'Wrong AmberTools executable/environment identity')
    require(digest({k: v for k, v in runtime.items() if k != 'identity'}) == RUNTIME_ID, 'Changed executable/source/runtime manifest')
    require(final_graph_identity == study['final_graph_identity'], 'Stale final graph identity')
    require(case in study['audits'], 'Unbound native component case')
    expected = study['audits'][case]
    for name, sha in expected['files'].items():
        require(checksum(Path(directory) / name) == sha, 'Changed native intermediate or stale input: ' + name)
    require({p.name for p in Path(directory).iterdir() if p.is_file()} == set(expected['files']), 'Native intermediate inventory changed')
    audit = audit_component(system, directory, runtime_identity=RUNTIME_ID)
    require(audit == expected, 'Changed graph, atom-ID mapping or native charge lineage')
    return signed({'schema': SCHEMA, 'case': case, 'status': 'blocked',
                   'acceptance_met': False, 'stage': 'native_charge_precision',
                   'classification': 'no_verified_native_high_precision_export',
                   'evidence_mode': study['evidence_mode'],
                   'final_graph_identity': final_graph_identity,
                   'study_identity': study['identity'], 'runtime_identity': RUNTIME_ID,
                   'executables': runtime['executables'], 'sources': runtime['sources'],
                   'method': 'AM1-BCC', 'atom_mapping': [
                       {k: row[k] for k in ('atom_id', 'name', 'element', 'formal_charge')}
                       for row in audit['atoms']],
                   'input_sha256': audit['input_sha256'], 'intermediate_hashes': audit['files'],
                   'observed_component_totals_e': audit['totals'],
                   'original_import_failure': audit['failure'], 'charge_tolerance_e': audit['tolerance_e'],
                   'authoritative_charge_vector': None, 'charge_record_published': False,
                   'fresh_preparation_executed': False, 'npt_executed': False,
                   'production_validated': False, 'simulation_readiness': 'not_established'})


def validate_diagnostic(record, system, directory, **bindings):
    require(record.get('identity') == digest({k: v for k, v in record.items() if k != 'identity'}), 'Changed diagnostic identity')
    require(record == inspect_component(system, directory, **bindings), 'Native precision diagnostic no longer reconstructs')
    return record


def require_authoritative_vector(record):
    """No rounded SQM/AC/MOL2/QOUT vector can authorize fresh preparation."""
    require(record.get('schema') == SCHEMA, 'Unsupported native precision record')
    raise ValueError('Blocked at native_charge_precision: pinned AmberTools24.8 has no verified supported high-precision AM1-BCC export; no charge vector published')
