"""Additive J33-C1 checkpoint continuation; historical J33 records stay readable.

The opaque Reference checkpoint owns the complete integrator/barostat RNG state.
No trajectory, state or thermodynamic value is borrowed from the reference run.
"""

import base64
import json
import shutil
import tempfile
import xml.etree.ElementTree as ET
from hashlib import sha256
from pathlib import Path

from island.analysis.npt import (
    FLAGS,
    FRAME_SCHEMA,
    NPTConfiguration,
    _body,
    _digest,
    _signed,
    _write_json,
    analyze_npt_frames,
    load_npt_workflow,
)
from island.exceptions import AnalysisError
from island.periodic_md import _backend, _engine, _publish_directory

SCHEMA = 'island.npt-continuation-restart.v1'
COMPARISON = 'island.npt-continuation-comparison.v1'
TOLERANCES = {'absolute': 0, 'relative': 0, 'checkpoint_and_rng': 'byte_exact'}


def require(ok, reason):
    if not ok:
        raise AnalysisError(reason)


def checksum(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def parent_identity(manifest):
    return _digest({k: v for k, v in manifest.items() if k not in {'identity', 'stage_counts'}})


def restart_assets(backend, config, model, integrator, mm, initial_checkpoint):
    require(config.production_steps == 0, 'Bounded C1 continuation requires production_steps=0')
    return {
        'npt-initial-checkpoint.txt': base64.b64encode(initial_checkpoint).decode(),
        'backend.json': json.dumps(backend, sort_keys=True, allow_nan=False),
        'system.xml': mm.XmlSerializer.serialize(model),
        'integrator.xml': mm.XmlSerializer.serialize(integrator),
    }


def write_restart(root, assets, manifest, frames, boundary):
    for name, raw in assets.items():
        (root / name).write_text(raw)
    backend = json.loads(assets['backend.json'])
    record = _signed({
        'schema': SCHEMA,
        'parent_workflow_identity': parent_identity(manifest),
        'manifest_identity': manifest['identity'],
        'npt_config_identity': manifest['config']['identity'],
        'source_identity': manifest['source_identity'],
        'assignment_identity': backend['assignment_identity'],
        'typed_graph_identity': backend['typed_graph_identity'],
        'charge_identity': backend['charge_identity'],
        'final_graph_identity': backend['final_graph_identity'],
        'backend_identity': backend['identity'],
        'split_boundary_npt_step': boundary,
        'checkpoint_identity': frames[-1]['checkpoint_sha256'],
        'rng_state_sha256': frames[-1]['checkpoint_sha256'],
        'rng_evidence': 'opaque OpenMM checkpoint RNG bytes; adaptive barostat state reconstructed by independent split-prefix integration',
        'barostat_state_restoration': 'independent_prefix_replay_then_boundary_checkpoint',
        'frame_content_identity': _digest(frames),
        'files': {name: checksum(root / name) for name in (*assets, 'manifest.json', 'frames.json', 'analysis.json')},
        **FLAGS,
    })
    _write_json(root / 'restart.json', record)
    return record


def load_continuation(root, *, expected_restart=None, expected_parent=None, expected_config=None):
    """Validate nested graph, assignment, frames and restart bindings offline."""
    root = Path(root)
    workflow = load_npt_workflow(root)
    record = json.loads((root / 'restart.json').read_text())
    require(record.get('schema') == SCHEMA and record.get('identity') == _digest(_body(record)), 'Restart identity mismatch')
    require(all(record.get(k) == v and type(record.get(k)) is type(v) for k, v in FLAGS.items()), 'Changed readiness flags')
    if expected_restart is not None:
        require(record['identity'] == expected_restart, 'Stale checkpoint/restart identity')
    if expected_parent is not None:
        require(record['parent_workflow_identity'] == expected_parent, 'Stale parent workflow identity')
    if expected_config is not None:
        require(record['npt_config_identity'] == expected_config, 'Changed NPT configuration')
    require(set(record['files']) == {'npt-initial-checkpoint.txt', 'backend.json', 'system.xml', 'integrator.xml', 'manifest.json', 'frames.json', 'analysis.json'}, 'Restart asset inventory mismatch')
    require(all(checksum(root / name) == digest for name, digest in record['files'].items()), 'Restart asset hash mismatch')
    manifest, frames = workflow['manifest'], workflow['frames']
    backend = json.loads((root / 'backend.json').read_text())
    _backend(backend)
    for key in ('assignment_identity', 'typed_graph_identity', 'charge_identity', 'final_graph_identity', 'backend_identity'):
        require(record[key] == backend['identity' if key == 'backend_identity' else key], 'Stale restart ' + key)
    require(backend['identity'] == manifest['backend_identity'], 'Workflow backend mismatch')
    require(record['source_identity'] == manifest['source_identity'] and backend['config']['source']['sha256'] == record['source_identity'], 'Changed source identity')
    require(record['parent_workflow_identity'] == parent_identity(manifest) and record['manifest_identity'] == manifest['identity'], 'Stale parent workflow identity')
    require(record['npt_config_identity'] == manifest['config']['identity'], 'Changed NPT configuration')
    require(record['frame_content_identity'] == _digest(frames), 'Restart frame identity mismatch')
    require(record['checkpoint_identity'] == record['rng_state_sha256'] == frames[-1]['checkpoint_sha256'], 'Checkpoint/RNG state mismatch')
    cfg = NPTConfiguration.from_payload(manifest['config'])
    require(record.get('barostat_state_restoration') == 'independent_prefix_replay_then_boundary_checkpoint', 'Unsupported barostat-state restoration')
    require(backend['final_graph_identity'] == manifest['final_graph_identity'], 'Changed final graph binding')
    system_xml = ET.fromstring((root / 'system.xml').read_text())
    bars = [f for f in system_xml.findall('./Forces/Force') if f.get('type') == cfg.barostat_type]
    require(len(bars) == 1, 'Missing barostat')
    bar = bars[0]
    require(float(bar.get('pressure')) == cfg.pressure_bar and float(bar.get('temperature')) == cfg.temperature_kelvin and int(bar.get('frequency')) == cfg.pressure_interval and int(bar.get('randomSeed')) == cfg.barostat_seed, 'Changed pressure/barostat configuration')
    integrator_xml = ET.fromstring((root / 'integrator.xml').read_text())
    require(integrator_xml.get('type') == cfg.thermostat_type and float(integrator_xml.get('stepSize')) == cfg.timestep_fs / 1000 and float(integrator_xml.get('temperature')) == cfg.temperature_kelvin and int(integrator_xml.get('randomSeed')) == cfg.seed, 'Changed integrator configuration')
    expected = analyze_npt_frames(frames, total_mass_dalton=manifest['total_mass_dalton'], target_temperature_kelvin=cfg.temperature_kelvin, target_pressure_bar=cfg.pressure_bar)
    require(all(workflow['analysis'].get(k) == v for k, v in expected.items()), 'Copied or changed density/thermodynamic analysis')
    boundary = record['split_boundary_npt_step']
    if boundary is not None:
        require(type(boundary) is int and 0 < boundary < cfg.npt_steps and boundary % cfg.checkpoint_interval == 0, 'Invalid split boundary')
        require(frames[-1]['stage'] == 'npt' and frames[-1]['step'] == cfg.nvt_steps + boundary, 'Split checkpoint/stage mismatch')
    for frame in frames:
        if frame['stage'] == 'npt':
            require(frame['target_pressure_bar'] == cfg.pressure_bar and frame['target_temperature_kelvin'] == cfg.temperature_kelvin, 'Changed pressure/temperature')
            expected_temperature = 2 * frame['kinetic_energy_kj_mol'] / (3 * len(manifest['atom_ids']) * 0.00831446261815324)
            require(frame['measured_temperature_kelvin'] == expected_temperature, 'Copied measured temperature')
    lineage_path = root / 'continuation.json'
    if lineage_path.exists():
        lineage = json.loads(lineage_path.read_text())
        require(lineage.get('schema') == 'island.npt-continuation-lineage.v1' and lineage.get('identity') == _digest(_body(lineage)), 'Continuation lineage identity mismatch')
        parent = lineage['parent_restart']
        require(parent.get('identity') == _digest(_body(parent)), 'Parent restart signature mismatch')
        require(lineage['parent_restart_identity'] == parent['identity'] and lineage['parent_workflow_identity'] == record['parent_workflow_identity'] == parent['parent_workflow_identity'], 'Stale parent workflow identity')
        for key in ('npt_config_identity', 'source_identity', 'assignment_identity', 'backend_identity'):
            require(record[key] == parent[key], 'Continuation changed ' + key)
        prefix = [f for f in frames if f['step'] <= cfg.nvt_steps + parent['split_boundary_npt_step']]
        require(_digest(prefix) == parent['frame_content_identity'] and prefix[-1]['checkpoint_sha256'] == parent['checkpoint_identity'], 'Continuation prefix/checkpoint changed')
        require(lineage['resumed_frame_identity'] == _digest(frames), 'Resumed trajectory identity mismatch')
    return workflow, record


def relocate_continuation(source, destination, **expectations):
    source, destination = Path(source), Path(destination)
    original = load_continuation(source, **expectations)
    require(not destination.exists(), 'Relocation destination exists')
    stage = Path(tempfile.mkdtemp(prefix='.npt-relocate-', dir=destination.parent))
    try:
        shutil.copytree(source, stage, dirs_exist_ok=True)
        require(load_continuation(stage, **expectations) == original, 'Relocation reconstruction mismatch')
        _publish_directory(stage, destination)
    finally:
        if stage.exists():
            shutil.rmtree(stage)
    return original[1]


def capture(ctx, config, atom_count, step, mm, unit):
    state = ctx.getState(getPositions=True, getVelocities=True, getEnergy=True)
    raw = ctx.createCheckpoint()
    vectors = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer).tolist()
    import numpy as np

    kinetic = float(state.getKineticEnergy().value_in_unit(unit.kilojoule_per_mole))
    return {
        'schema': FRAME_SCHEMA, 'stage': 'npt', 'step': step,
        'time_ps': float(state.getTime().value_in_unit(unit.picosecond)),
        'positions_nm': state.getPositions(asNumpy=True).value_in_unit(unit.nanometer).tolist(),
        'velocities_nm_ps': state.getVelocities(asNumpy=True).value_in_unit(unit.nanometer / unit.picosecond).tolist(),
        'box_vectors_nm': vectors, 'volume_nm3': float(np.linalg.det(np.asarray(vectors))),
        'potential_energy_kj_mol': float(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)),
        'kinetic_energy_kj_mol': kinetic,
        'target_temperature_kelvin': config.temperature_kelvin, 'target_pressure_bar': config.pressure_bar,
        'measured_temperature_kelvin': 2 * kinetic / (3 * atom_count * 0.00831446261815324),
        'measured_temperature_source': 'kinetic_energy', 'measured_pressure_bar': None,
        'pressure_observation_source': 'unavailable', 'minimum_distance_angstrom': None,
        'checkpoint_sha256': sha256(raw).hexdigest(), 'checkpoint_base64': base64.b64encode(raw).decode(),
    }


def resume_continuation(source, destination, *, expected_restart, expected_parent, expected_config):
    """Resume solely from the relocated bundle; publish a separate complete bundle."""
    source, destination = Path(source), Path(destination)
    workflow, restart = load_continuation(source, expected_restart=expected_restart, expected_parent=expected_parent, expected_config=expected_config)
    require(not destination.exists(), 'Resume destination exists')
    manifest, frames = workflow['manifest'], workflow['frames']
    cfg = NPTConfiguration.from_payload(manifest['config'])
    boundary = restart['split_boundary_npt_step']
    require(boundary is not None and cfg.production_steps == 0, 'No resumable NPT boundary')
    import openmm as mm
    from openmm import unit

    require(_engine(mm) == manifest['engine'], 'Unsupported engine for exact checkpoint restoration')
    model = mm.XmlSerializer.deserialize((source / 'system.xml').read_text())
    integrator = mm.XmlSerializer.deserialize((source / 'integrator.xml').read_text())
    bars = [f for f in model.getForces() if isinstance(f, mm.MonteCarloBarostat)]
    require(len(bars) == 1, 'Missing barostat')
    bar = bars[0]
    require(bar.getDefaultPressure().value_in_unit(unit.bar) == cfg.pressure_bar and bar.getDefaultTemperature().value_in_unit(unit.kelvin) == cfg.temperature_kelvin and bar.getFrequency() == cfg.pressure_interval and bar.getRandomNumberSeed() == cfg.barostat_seed, 'Changed pressure/barostat configuration')
    require(integrator.getRandomNumberSeed() == cfg.seed and integrator.getStepSize().value_in_unit(unit.picosecond) == cfg.timestep_fs / 1000 and integrator.getTemperature().value_in_unit(unit.kelvin) == cfg.temperature_kelvin, 'Changed integrator configuration')
    ctx = mm.Context(model, integrator, mm.Platform.getPlatformByName('Reference'))
    try:
        # OpenMM does not checkpoint the adaptive barostat's proposal scale and
        # acceptance counters. Reconstruct those by independent prefix integration,
        # using only this split bundle, before restoring its boundary checkpoint.
        ctx.loadCheckpoint(base64.b64decode((source / 'npt-initial-checkpoint.txt').read_text(), validate=True))
        for _ in range(boundary):
            integrator.step(1)
        replayed = capture(ctx, cfg, len(manifest['atom_ids']), frames[-1]['step'], mm, unit)
        require(replayed == frames[-1], 'Independent barostat-state replay differs from split boundary')
        ctx.loadCheckpoint(base64.b64decode(frames[-1]['checkpoint_base64'], validate=True))
        restored = capture(ctx, cfg, len(manifest['atom_ids']), frames[-1]['step'], mm, unit)
        require(restored == frames[-1], 'Restored checkpoint state/RNG differs from retained boundary')
        for i in range(boundary + 1, cfg.npt_steps + 1):
            integrator.step(1)
            if i % cfg.checkpoint_interval == 0 or i == cfg.npt_steps:
                frames.append(capture(ctx, cfg, len(manifest['atom_ids']), cfg.nvt_steps + i, mm, unit))
    finally:
        del ctx, integrator
    manifest = _body(manifest)
    manifest['stage_counts'] = {stage: sum(f['stage'] == stage for f in frames) for stage in {f['stage'] for f in frames}}
    manifest = _signed(manifest)
    analysis = analyze_npt_frames(frames, total_mass_dalton=manifest['total_mass_dalton'], target_temperature_kelvin=cfg.temperature_kelvin, target_pressure_bar=cfg.pressure_bar)
    analysis.update(total_mass_dalton=manifest['total_mass_dalton'], bindings=workflow['analysis']['bindings'], manifest_identity=manifest['identity'])
    stage = Path(tempfile.mkdtemp(prefix='.npt-resume-', dir=destination.parent))
    try:
        _write_json(stage / 'manifest.json', manifest)
        _write_json(stage / 'frames.json', _signed({'schema': 'island.periodic-npt-frames.v1', 'frames': frames}))
        _write_json(stage / 'analysis.json', _signed(analysis))
        assets = {name: (source / name).read_text() for name in ('npt-initial-checkpoint.txt', 'backend.json', 'system.xml', 'integrator.xml')}
        write_restart(stage, assets, manifest, frames, None)
        _write_json(stage / 'continuation.json', _signed({'schema': 'island.npt-continuation-lineage.v1', 'parent_restart': restart, 'parent_restart_identity': restart['identity'], 'parent_workflow_identity': restart['parent_workflow_identity'], 'split_boundary_npt_step': boundary, 'resumed_frame_identity': _digest(frames), **FLAGS}))
        load_continuation(stage)
        _publish_directory(stage, destination)
    except BaseException as error:
        (stage / 'failure.txt').write_text(str(error))
        error.add_note('Failed resume candidate preserved: ' + str(stage))
        raise
    return load_continuation(destination)[1]


def compare_continuations(first, resumed, split):
    a, ar = load_continuation(first)
    b, br = load_continuation(resumed)
    s, sr = load_continuation(split)
    lineage = json.loads((Path(resumed) / 'continuation.json').read_text())
    require(lineage.get('identity') == _digest(_body(lineage)), 'Continuation lineage identity mismatch')
    require(lineage['parent_restart_identity'] == sr['identity'] and lineage['parent_workflow_identity'] == sr['parent_workflow_identity'] == ar['parent_workflow_identity'] == br['parent_workflow_identity'], 'Stale parent workflow identity')
    require(lineage['resumed_frame_identity'] == _digest(b['frames']), 'Resumed lineage frame mismatch')
    require(s['frames'] == b['frames'][:len(s['frames'])], 'Split prefix changed during continuation')
    require(a['manifest'] == b['manifest'], 'Uninterrupted/resumed workflow identities differ')
    require(a['frames'][-1]['step'] == b['frames'][-1]['step'] == a['manifest']['config']['nvt_steps'] + a['manifest']['config']['npt_steps'], 'Incomplete NPT history')
    mismatches = []
    if len(a['frames']) != len(b['frames']):
        mismatches.append({'field': 'frame_count'})
    for i, (af, bf) in enumerate(zip(a['frames'], b['frames'])):
        for key in sorted(set(af) | set(bf)):
            if af.get(key) != bf.get(key):
                mismatches.append({'frame_index': i, 'step': af['step'], 'field': key})
    if a['analysis'] != b['analysis']:
        mismatches.append({'field': 'analysis_volume_density_observables'})
    result = _signed({
        'schema': COMPARISON, 'passed': not mismatches, 'mismatches': mismatches,
        'parent_workflow_identity': ar['parent_workflow_identity'],
        'npt_config_identity': ar['npt_config_identity'],
        'source_identity': ar['source_identity'], 'assignment_identity': ar['assignment_identity'],
        'final_graph_identity': ar['final_graph_identity'], 'split_restart_identity': sr['identity'],
        'split_boundary_npt_step': sr['split_boundary_npt_step'],
        'uninterrupted_trajectory_identity': _digest(a['frames']), 'resumed_trajectory_identity': _digest(b['frames']),
        'checkpoint_hashes': [f['checkpoint_sha256'] for f in a['frames']],
        'resumed_checkpoint_hashes': [f['checkpoint_sha256'] for f in b['frames']],
        'rng_state_hashes': [f['checkpoint_sha256'] for f in a['frames']],
        'resumed_rng_state_hashes': [f['checkpoint_sha256'] for f in b['frames']],
        'density_frame_identity': a['analysis']['frame_content_identity'],
        'tolerances': TOLERANCES, 'rng_evidence': ar['rng_evidence'], **FLAGS,
    })
    return result
