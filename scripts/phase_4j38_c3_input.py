"""Exact J38 PE/aromatic quantum inputs; no charge bundle inputs."""

import json
from collections import deque

import numpy as np

from island.forcefields.ambertools.provided import reference_digest, require
from scripts.phase_4j36_input import POLICY
from scripts.phase_4j38_fixture import build


def build_input(kind):
    require(kind in ('pe', 'aromatic'), 'Only exact blocked J38 graphs allowed')
    fixture, _ = build(kind)
    system, graph = fixture.packed.system, fixture.packed.graph
    graph.validate_integrity(system)
    lengths = np.array(system.box.lengths)
    coords = {i: np.array(system.coordinates.get(i)) for i in system.topology.sites}
    components = []
    for ids in sorted(sorted(c) for c in graph.payload['components']):
        assigned = {ids[0]: coords[ids[0]]}
        queue = deque([ids[0]])
        while queue:
            i = queue.popleft()
            for j in sorted(system.topology.neighbors(i)):
                delta = coords[j] - coords[i]
                delta -= np.rint(delta / lengths) * lengths
                pos = assigned[i] + delta
                if j in assigned:
                    require(np.max(abs(pos - assigned[j])) < 1e-6, 'Inconsistent input component periodic image')
                else:
                    assigned[j] = pos
                    queue.append(j)
        atoms = [{'atom_id': i, 'element': system.topology.sites[i].element,
                  'position_angstrom': assigned[i].tolist(), 'formal_charge': system.topology.sites[i].formal_charge} for i in ids]
        components.append({'atoms': atoms, 'formal_charge': sum(a['formal_charge'] for a in atoms), 'multiplicity': 1})
    policy = json.loads(json.dumps(POLICY))
    return {'schema': 'island.j38-c3-quantum-input.v1', 'fixture': kind,
            'final_graph_identity': graph.identity, 'graph': graph.payload, 'components': components,
            **{k: v for k, v in policy.items() if k != 'schema'}}


def validate_input(data):
    expected = build_input(data['fixture'])
    require(reference_digest(data) == reference_digest(expected), 'Stale quantum graph/geometry/mapping, changed formal charge, method, basis or fitting policy')
    return data
