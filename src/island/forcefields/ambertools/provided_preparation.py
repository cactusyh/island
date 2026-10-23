"""Additive GAFF2 provided-charge preparation, preserving legacy native records."""

import json
from copy import deepcopy
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

from island.forcefields.ambertools.models import digest
from island.forcefields.ambertools.multimolecule import require, unwrapped_component
from island.forcefields.ambertools.provided import GAFF2ProvidedCharges
from island.forcefields.ambertools.runtime import checksum, validate_runtime

SCHEMA = 'island.gaff2-provided-preparation.v1'


def validate_record(system, record, signature, *, imported_signature, imported_provenance, imported_mapping, source_sha256, charge_result):
    from island.forcefields.ambertools.models import validate_preparation_record
    from island.graph.final import FinalChemicalGraph

    require(signature == digest(record) and record['schema'] == SCHEMA, 'Changed provided preparation identity')
    require(record['charge_method'] == 'provided' and record['native_am1bcc_success'] is False and record['diagnostic'] == 'provided_charge_validated', 'Contradictory provided-charge mode')
    require(record['requested_force_field'] == 'gaff2' and record['assembly_policy'] == 'native_parameters_unchanged_external_charges_no_normalization', 'Changed provided preparation policy')
    require(record['imported_result_signature'] == imported_signature and record['artifact_sha256']['result.prmtop'] == source_sha256, 'Changed assembled topology identity')
    require(digest({k: v for k, v in record.items() if k != 'imported_result_signature'}) == digest(dict(imported_provenance['ambertools_preparation'])), 'Contradictory preparation provenance')
    require(imported_provenance['force_field'] == 'gaff2' and imported_provenance['charge_method'] == 'provided', 'Changed imported charge method')
    graph = FinalChemicalGraph(record['packed_graph_raw'])
    # Compare the exact chemical graph and ancestry, allowing only the deliberate
    # removal of the box in the backend's chemical view.
    from island.core import SimulationBox
    from island.graph import final_graph

    packed = system.copy()
    packed.box = SimulationBox(*record['box_lengths_angstrom'])
    graph.validate_integrity(packed)
    require(final_graph(packed, molecule_membership={int(k): v for k, v in graph.payload['molecule_membership'].items()}, transformations=graph.payload['transformations']).identity == graph.identity, 'Changed packed final graph')
    provided = GAFF2ProvidedCharges(record['provided_charge_raw'])
    charges = provided.validate(graph, expected_reference_identity=record['reference_identity'])
    require(record['mapping'] == dict(imported_mapping), 'Changed atom mapping')
    require(record['input_coordinate_identity'] == digest({str(i): system.coordinates.get(i).tolist() for i in sorted(system.topology.sites)}), 'Stale coordinates')
    expected = sorted(sorted(c) for c in system.topology.connected_components())
    require([e['atom_ids'] for e in record['components']] == expected, 'Changed component coverage')
    for entry in record['components']:
        sub, images = unwrapped_component(system, entry['atom_ids'], record['box_lengths_angstrom'])
        require(images == entry['image_shifts'], 'Changed image mapping')
        child = entry['record']
        validate_preparation_record(sub, child, entry['record_signature'], imported_signature=child['imported_result_signature'], imported_provenance={'force_field': 'gaff2', 'charge_method': 'provided', 'ambertools_preparation': {k: v for k, v in child.items() if k != 'imported_result_signature'}}, imported_mapping=child['lineage']['prmtop_index_to_site_id'], source_sha256=child['artifact_sha256']['result.prmtop'], charge_result=SimpleNamespace(tolerance=charge_result.tolerance, assignments={int(i): SimpleNamespace(charge=q) for i, q in entry['native_bridge_charges'].items()}))
        require(child['charge_outcome']['qm_run'] is False and child['parmchk2_estimation_lines'] == [], 'Unexpected QM/estimated parameters')
        require(child['provided_charge_input']['charges'] == {str(i): charges[i] for i in entry['atom_ids']}, 'Changed component charge input')
        for key in ('force_field_data', 'leaprc', 'executable_sha256', 'tool_package'):
            require(child[key] == record[key], 'Changed component source/runtime')
    expected_mapping = dict(enumerate(site for entry in record['components'] for _, site in sorted(entry['record']['lineage']['prmtop_index_to_site_id'].items(), key=lambda item: int(item[0]))))
    require(record['mapping'] == expected_mapping, 'Assembled mapping differs from component lineage')
    for i, q in charges.items():
        require(abs(q - charge_result.assignments[i].charge) <= 1e-14, 'Provided charge precision lost in assembled topology')


def parameterize_provided(packed, graph, provided, output, *, runtime, runtime_identity, expected_reference_identity):
    """Retype each exact component; original external vector is authoritative.

    Amber's rc/MOL2 bridge supplies typing and parameters only. Charges in the
    final native topology are written from the unchanged validated input at
    16-digit scientific precision, independently checked after import.
    """
    import parmed as pmd
    from parmed.amber import AmberParm
    from parmed.amber.amberformat import FortranFormat

    from island.forcefields.amber import import_amber_prmtop
    from island.forcefields.ambertools import (
        AmberToolsOptions,
        AmberToolsParameterizationEngine,
        AmberToolsPreparationResult,
    )
    from island.periodic import periodic_chemical_view
    from island.workflows import storage

    graph.validate_integrity(packed)
    charges = provided.validate(graph, expected_reference_identity=expected_reference_identity)
    validate_runtime(runtime, expected_identity=runtime_identity)
    system = periodic_chemical_view(packed)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    entries, ordered_ids, combined = [], [], None
    for index, ids in enumerate(sorted(sorted(c) for c in system.topology.connected_components())):
        sub, images = unwrapped_component(system, ids, packed.box.lengths)
        directory = output / f'component-{index}'
        directory.mkdir()
        options = AmberToolsOptions('gaff2', 'provided', provided_charges={i: charges[i] for i in ids}, charge_source=provided.identity, amberhome=runtime['amberhome'], work_root=directory, retain_success_artifacts=True, max_atoms=100, timeout_seconds=1800)
        native = AmberToolsParameterizationEngine().parameterize(sub, options)
        require(native.record['parmchk2_estimation_lines'] == [], 'parameterization_incomplete: estimated GAFF2 terms')
        native.validate_integrity(sub)
        structure = pmd.load_file(str(Path(native.record['artifact_dir']) / 'result.prmtop'), xyz=str(Path(native.record['artifact_dir']) / 'result.rst7'))
        atom_ids = [native.imported_result.mapping[i] for i in range(len(structure.atoms))]
        for atom, site_id in zip(structure.atoms, atom_ids, strict=True):
            atom.charge = charges[site_id]
        ordered_ids.extend(atom_ids)
        entries.append({'atom_ids': ids, 'image_shifts': images, 'record': deepcopy(dict(native.record)), 'record_signature': native.record_signature, 'native_bridge_charges': {str(i): native.imported_result.charge_result.assignments[i].charge for i in ids}})
        combined = structure if combined is None else combined + structure
    require(len(ordered_ids) == len(set(ordered_ids)) == len(charges), 'Nonbijective component mapping')
    amber = AmberParm.from_structure(combined)
    # AmberParm addition can reconstruct atoms from its rounded parm_data and
    # discard component atom.charge edits. Assign the unchanged authoritative
    # vector after all native structure copies/conversions, before serialization.
    for atom, site_id in zip(amber.atoms, ordered_ids, strict=True):
        atom.charge = charges[site_id]
    amber.formats['CHARGE'] = FortranFormat('3E24.16')
    amber.write_parm(str(output / 'result.prmtop'))
    amber.save(str(output / 'result.rst7'), format='rst7')
    mapping = dict(enumerate(ordered_ids))
    imported = import_amber_prmtop(system, output / 'result.prmtop', mapping, source='J35 external provided electrostatics with native GAFF2 parameters; no AM1-BCC claim', force_field='gaff2', charge_method='provided')
    imported.validate_integrity(system)
    first = entries[0]['record']
    record = {'schema': SCHEMA, 'engine_version': 'provided-1', 'requested_force_field': 'gaff2', 'charge_method': 'provided', 'native_am1bcc_success': False, 'diagnostic': 'provided_charge_validated', 'assembly_policy': 'native_parameters_unchanged_external_charges_no_normalization', 'charge_serialization': 'Amber native CHARGE 3E24.16; original vector unchanged; absolute roundtrip <=1e-14 e', 'packed_graph_raw': graph.json_text, 'provided_charge_raw': provided.raw, 'reference_identity': expected_reference_identity, 'box_lengths_angstrom': list(packed.box.lengths), 'runtime_identity': runtime_identity, 'components': entries, 'mapping': mapping, 'input_coordinate_identity': digest({str(i): system.coordinates.get(i).tolist() for i in sorted(system.topology.sites)}), 'artifact_dir': str(output), 'artifact_sha256': {n: checksum(output / n) for n in ('result.prmtop', 'result.rst7')}, 'parmchk2_estimation_lines': [], **{k: first[k] for k in ('force_field_data', 'leaprc', 'executable_sha256', 'tool_package', 'tool_versions', 'executables', 'amberhome')}}
    imported = replace(imported, provenance={**dict(imported.provenance), 'ambertools_preparation': record}, result_signature='')
    imported = replace(imported, result_signature=imported.content_signature())
    record = {**record, 'imported_result_signature': imported.result_signature}
    result = AmberToolsPreparationResult(imported, record, digest(record))
    result.validate_integrity(system)
    (output / 'record.json').write_text(json.dumps(storage.encode(record), sort_keys=True, indent=2) + '\n')
    return result


def validate_backend_record(payload):
    """Restore JSON integer keys and validate the additive record offline."""
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
    validate_record(periodic_chemical_view(_system_from_payload(payload['system'])), record, payload['typed_graph']['evidence'][0], imported_signature=record['imported_result_signature'], imported_provenance={'force_field': 'gaff2', 'charge_method': 'provided', 'ambertools_preparation': {k: v for k, v in record.items() if k != 'imported_result_signature'}}, imported_mapping=record['mapping'], source_sha256=payload['config']['source']['sha256'], charge_result=charge_result)
    return record
