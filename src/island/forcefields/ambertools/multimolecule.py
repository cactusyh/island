"""Explicit disjoint assembly of independently native-prepared GAFF2 molecules.

The original single-molecule engine and its v2/v3 records are unchanged. No
parameter or charge is inferred between disconnected components.
"""

import json
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from island.core import Coordinates
from island.exceptions import AmberToolsInputError, InvalidAmberImportResultError
from island.forcefields.ambertools.runtime import checksum, validate_runtime

SCHEMA = 'island.ambertools-multimolecule-preparation.v1'


def require(ok, message):
    if not ok:
        raise AmberToolsInputError(message)


def component_system(system, ids):
    result = system.copy()
    for site in set(result.topology.sites) - set(ids):
        result.topology.remove_site(site)
    result.coordinates = Coordinates({i: system.coordinates.get(i) for i in ids})
    result.box = None
    result.validate()
    return result


def unwrapped_component(system, ids, lengths):
    from island.analysis.network import _image_geometry

    sub = component_system(system, ids)
    adjacency = {i: sub.topology.neighbors(i) for i in ids}
    frame = SimpleNamespace(box_vectors_angstrom=np.diag(lengths), coordinates_angstrom=np.array([system.coordinates.get(i) for i in ids]))
    geometry = _image_geometry(frame, {i: n for n, i in enumerate(ids)}, adjacency, 1e-8)
    require(geometry['consistent'], 'Finite GAFF2 molecule has inconsistent periodic images')
    sub.coordinates = Coordinates(geometry['positions'])
    return sub, geometry['images']


def validate_record(system, record, signature, *, imported_signature, imported_provenance, imported_mapping, source_sha256, charge_result):
    """Offline validation of native component lineage and assembled charges."""
    from island.forcefields.ambertools.models import digest, validate_preparation_record

    def check(ok, message):
        if not ok:
            raise InvalidAmberImportResultError('Multimolecule Amber preparation: ' + message)

    check(signature == digest(record), 'record signature mismatch')
    check(record['schema'] == SCHEMA and record['engine_version'] == 'multimolecule-1', 'unsupported schema')
    check(record['requested_force_field'] == 'gaff2' and record['charge_method'] == 'am1bcc', 'unsupported family/charge mode')
    check(record['assembly_policy'] == 'parmed_disjoint_union_no_parameter_generation', 'unsupported assembly policy')
    check(record['imported_result_signature'] == imported_signature and record['artifact_sha256']['result.prmtop'] == source_sha256, 'assembled import identity changed')
    inner = {k: v for k, v in record.items() if k != 'imported_result_signature'}
    check(digest(inner) == digest(dict(imported_provenance['ambertools_preparation'])), 'contradictory signed provenance')
    check(imported_provenance['force_field'] == 'gaff2' and imported_provenance['charge_method'] == 'AM1-BCC', 'import method mismatch')
    expected = sorted([sorted(c) for c in system.topology.connected_components()])
    check(len(expected) >= 2 and [c['atom_ids'] for c in record['components']] == expected, 'component coverage/membership changed')
    mapping = {int(k): v for k, v in record['mapping'].items()}
    check(mapping == dict(imported_mapping), 'assembled atom mapping changed')
    xyz = {str(i): system.coordinates.get(i).tolist() for i in sorted(system.topology.sites)}
    check(record['input_coordinate_identity'] == digest(xyz), 'stale packed coordinates')
    check(len(record['runtime_identity']) == 64, 'missing runtime identity')
    for entry in record['components']:
        child = entry['record']
        sub, images = unwrapped_component(system, entry['atom_ids'], record['box_lengths_angstrom'])
        check(images == entry['image_shifts'], 'changed component periodic image assignment')
        validate_preparation_record(sub, child, entry['record_signature'], imported_signature=child['imported_result_signature'], imported_provenance={'force_field': 'gaff2', 'charge_method': 'AM1-BCC', 'ambertools_preparation': {k: v for k, v in child.items() if k != 'imported_result_signature'}}, imported_mapping=child['lineage']['prmtop_index_to_site_id'], source_sha256=child['artifact_sha256']['result.prmtop'], charge_result=SimpleNamespace(tolerance=charge_result.tolerance))
        check(child['parmchk2_estimation_lines'] == [], 'estimated component parameters')
        for key in ('force_field_data', 'leaprc', 'executable_sha256', 'tool_package'):
            check(child[key] == record[key], 'component source/environment mismatch: ' + key)
        check(set(entry['charges']) == {str(i) for i in entry['atom_ids']}, 'incomplete component charges')
        for site in entry['atom_ids']:
            check(abs(entry['charges'][str(site)] - charge_result.assignments[site].charge) <= 1e-12, 'assembled charge differs from native component')


def parameterize_multimolecule(system, options, output, *, runtime, runtime_identity, box_lengths):
    """Run native AM1-BCC/GAFF2 on every component of this exact input graph."""
    from island.forcefields.amber import import_amber_prmtop
    from island.forcefields.ambertools import (
        AmberToolsParameterizationEngine,
        AmberToolsPreparationResult,
    )
    from island.forcefields.ambertools.models import digest

    validate_runtime(runtime, expected_identity=runtime_identity)
    require(options.force_field == 'gaff2' and options.charge_method == 'am1bcc', 'Multimolecule preparation requires GAFF2 AM1-BCC')
    require(system.box is None, 'Prepare the explicit chemical view, then bind its periodic final graph')
    from island.forcefields.ambertools.policy import validate_size_policy

    validate_size_policy(options.max_atoms, options.charge_method, len(system.topology.sites))
    require(len(box_lengths) == 3 and all(np.isfinite(v) and v > 0 for v in box_lengths), 'Positive orthorhombic box required for component preparation')
    components = sorted([sorted(c) for c in system.topology.connected_components()])
    require(len(components) >= 2, 'Explicit multimolecule assembly requires multiple complete molecules')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    import parmed as pmd

    entries, ordered_ids, combined = [], [], None
    for index, ids in enumerate(components):
        sub, images = unwrapped_component(system, ids, box_lengths)
        directory = output / ('component-' + str(index))
        directory.mkdir()
        native = AmberToolsParameterizationEngine().parameterize(sub, replace(options, work_root=directory))
        require(native.record['parmchk2_estimation_lines'] == [], 'Estimated GAFF2 component interactions are outside J34')
        native.validate_integrity(sub)
        structure = pmd.load_file(str(Path(native.record['artifact_dir']) / 'result.prmtop'), xyz=str(Path(native.record['artifact_dir']) / 'result.rst7'))
        atom_ids = [native.imported_result.mapping[i] for i in range(len(structure.atoms))]
        ordered_ids.extend(atom_ids)
        entries.append({'atom_ids': ids, 'image_shifts': images, 'record': deepcopy(dict(native.record)), 'record_signature': native.record_signature, 'charges': {str(i): native.imported_result.charge_result.assignments[i].charge for i in ids}})
        combined = structure if combined is None else combined + structure
    require(len(ordered_ids) == len(set(ordered_ids)) == len(system.topology.sites), 'Component atom mapping is not bijective')
    combined.save(str(output / 'result.prmtop'), format='amber')
    combined.save(str(output / 'result.rst7'), format='rst7')
    mapping = dict(enumerate(ordered_ids))
    imported = import_amber_prmtop(system, output / 'result.prmtop', mapping, source='J34 ParmEd disjoint union of fresh native GAFF2 component topologies', force_field='gaff2', charge_method='AM1-BCC', charge_tolerance=options.charge_tolerance)
    imported.validate_integrity(system)
    first = entries[0]['record']
    record = {'schema': SCHEMA, 'engine_version': 'multimolecule-1', 'requested_force_field': 'gaff2', 'charge_method': 'am1bcc', 'assembly_policy': 'parmed_disjoint_union_no_parameter_generation', 'box_lengths_angstrom': list(box_lengths), 'runtime_identity': runtime_identity, 'components': entries, 'mapping': mapping, 'input_coordinate_identity': digest({str(i): system.coordinates.get(i).tolist() for i in sorted(system.topology.sites)}), 'artifact_dir': str(output), 'artifact_sha256': {n: checksum(output / n) for n in ('result.prmtop', 'result.rst7')}, 'parmchk2_estimation_lines': [], **{k: first[k] for k in ('force_field_data', 'leaprc', 'executable_sha256', 'tool_package', 'tool_versions', 'executables', 'amberhome')}}
    # JSON normalization matches the existing native signed provenance convention.
    imported = replace(imported, provenance={**dict(imported.provenance), 'ambertools_preparation': record}, result_signature='')
    imported = replace(imported, result_signature=imported.content_signature())
    record = {**record, 'imported_result_signature': imported.result_signature}
    result = AmberToolsPreparationResult(imported, record, digest(record))
    result.validate_integrity(system)
    from island.workflows import storage

    (output / 'record.json').write_text(json.dumps(storage.encode(record), indent=2, sort_keys=True) + '\n')
    return result


def validate_backend_record(payload):
    """Restore integer keys lost by the unified JSON view and verify v1 lineage."""
    from island.graph.final import _system_from_payload
    from island.periodic import periodic_chemical_view

    record = deepcopy(payload['parameter_assignment']['native_payload']['record'])
    record['mapping'] = {int(k): v for k, v in record['mapping'].items()}
    for entry in record['components']:
        child = entry['record']
        child['expected_cip_by_site'] = {int(k): v for k, v in child['expected_cip_by_site'].items()}
        for key in ('typed_mol2_index_to_site_id', 'prmtop_index_to_site_id'):
            child['lineage'][key] = {int(k): v for k, v in child['lineage'][key].items()}
    charge_result = SimpleNamespace(tolerance=record['components'][0]['record']['charge_validation_tolerance_e'], assignments={s['id']: SimpleNamespace(charge=s['charge']) for s in payload['numerical_data']['sites']})
    validate_record(periodic_chemical_view(_system_from_payload(payload['system'])), record, payload['typed_graph']['evidence'][0], imported_signature=record['imported_result_signature'], imported_provenance={'force_field': 'gaff2', 'charge_method': 'AM1-BCC', 'ambertools_preparation': {k: v for k, v in record.items() if k != 'imported_result_signature'}}, imported_mapping=record['mapping'], source_sha256=payload['config']['source']['sha256'], charge_result=charge_result)
    return record


def restore_multimolecule(system, directory, *, archive_sha256, runtime_identity):
    """Restore a declared native preparation on its exact graph, never a trajectory.

    This explicit input path is separate from fresh native parameterization. Its
    immutable archive is a current-phase source asset, not an automatic fallback.
    """
    from island.forcefields.amber import import_amber_prmtop
    from island.forcefields.ambertools.models import AmberToolsPreparationResult, digest
    from island.workflows import storage

    directory = Path(directory)
    record = storage.decode(json.loads((directory / 'record.json').read_text()))
    require(record['runtime_identity'] == runtime_identity, 'Frozen preparation runtime mismatch')
    for name, expected in record['artifact_sha256'].items():
        require(Path(name).name == name and checksum(directory / name) == expected, 'Frozen preparation artifact changed: ' + name)
    for index, entry in enumerate(record['components']):
        child = entry['record']
        native_dir = directory / ('component-' + str(index)) / Path(child['artifact_dir']).name
        for name, expected in child['artifact_sha256'].items():
            require(Path(name).name == name and checksum(native_dir / name) == expected, 'Frozen component artifact changed: ' + name)
    imported = import_amber_prmtop(system, directory / 'result.prmtop', record['mapping'], source='J34 ParmEd disjoint union of fresh native GAFF2 component topologies', force_field='gaff2', charge_method='AM1-BCC', charge_tolerance=record['components'][0]['record']['charge_validation_tolerance_e'])

    def assemble(body):
        adopted = replace(imported, provenance={**dict(imported.provenance), 'ambertools_preparation': body}, result_signature='')
        adopted = replace(adopted, result_signature=adopted.content_signature())
        final = {**body, 'imported_result_signature': adopted.result_signature}
        result = AmberToolsPreparationResult(adopted, final, digest(final))
        result.validate_integrity(system)
        return result

    original = assemble({k: v for k, v in record.items() if k != 'imported_result_signature'})
    require(original.imported_result.result_signature == record['imported_result_signature'], 'Frozen native import identity changed')
    body = {k: v for k, v in record.items() if k != 'imported_result_signature'}
    body['artifact_dir'] = str(directory)
    body['retained_native_preparation'] = {'archive_sha256': archive_sha256, 'original_record_signature': digest(record), 'original_artifact_dir': record['artifact_dir'], 'policy': 'explicit_exact_graph_native_input_no_trajectory_reuse'}
    return assemble(body)
