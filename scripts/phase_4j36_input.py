"""Fresh quantum input from the retained complete graphs, never charge bundles."""

import json
from collections import deque

import numpy as np

from island.forcefields.ambertools.provided import reference_digest, require
from scripts.phase_4j34_real_fixture import build_fixture

POLICY = {
    'schema': 'island.psi4-resp-policy.v1', 'method': 'HF', 'basis': '6-31G*',
    'optimization_method': 'HF', 'optimization_basis': 'STO-3G',
    'geometry_optimization': {'geom_maxiter': 100, 'g_convergence': 'GAU_LOOSE', 'scf_type': 'df', 'reference': 'rhf'},
    'resp_options': {'METHOD_ESP': 'hf', 'BASIS_ESP': '6-31g*',
                     'PSI4_OPTIONS': {'scf_type': 'df', 'e_convergence': 1e-8, 'd_convergence': 1e-8},
                     'VDW_SCALE_FACTORS': [1.4, 1.6, 1.8, 2.0], 'VDW_POINT_DENSITY': 1.0,
                     'RESP_A': 0.0005, 'RESP_B': 0.1, 'IHFREE': True, 'TOLER': 1e-8, 'MAX_IT': 200},
}


def build_input(kind):
    fixture = build_fixture(kind, 'GAFF2').fixture
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
    return {'schema': 'island.j36-quantum-input.v1', 'fixture': kind,
            'final_graph_identity': graph.identity, 'graph': graph.payload, 'components': components,
            **{k: v for k, v in policy.items() if k != 'schema'}}


def validate_input(data):
    expected = build_input(data['fixture'])
    require(reference_digest(data) == reference_digest(expected), 'Stale quantum graph/geometry/mapping, changed formal charge, method, basis or fitting policy')
    return data
