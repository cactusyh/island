"""Load independently validated RESP/Psi4 evidence; never infer quantum provenance.

This adapter is a loader, not a quantum or RESP implementation. The caller must
obtain the manifest identity from an independent review of the actual quantum
calculation and fit. Self-signed metadata is not a trust anchor. Runtime and
intermediate assets remain necessary when loading offline after relocation.
"""

import json
from dataclasses import dataclass
from pathlib import Path

from island.forcefields.ambertools.provided import (
    GAFF2ProvidedCharges,
    reference_digest,
    require,
)

SCHEMA = 'island.resp-psi4-result.v1'


def file_hash(path):
    import hashlib

    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def bound_file(root, name, expected):
    root = Path(root).resolve()
    path = (root / name).resolve()
    require(path.is_relative_to(root) and path.is_file(), 'Missing or unsafe RESP asset: ' + name)
    require(file_hash(path) == expected, 'Mutated RESP asset: ' + name)
    return path


@dataclass(frozen=True)
class RESPCharges:
    """Immutable validated result; an explicit RESP-only entry to provided mode."""

    raw: str

    @property
    def payload(self):
        return json.loads(self.raw)

    @classmethod
    def load(cls, directory, graph, *, expected_identity):
        p = json.loads((Path(directory) / 'resp-result.json').read_text())
        require(p.get('schema') == SCHEMA, 'Unsupported RESP result schema')
        require(p.get('charge_origin') == 'resp_psi4' and p.get('charge_method') == 'RESP'
                and p.get('quantum_software') == 'Psi4', 'Only actual RESP/Psi4 evidence is accepted; external OPLS-AA is excluded')
        require(isinstance(expected_identity, str) and len(expected_identity) == 64
                and p.get('identity') == expected_identity == reference_digest({k: v for k, v in p.items() if k != 'identity'}), 'Independent RESP review identity mismatch')
        for key in ('psi4_version', 'method', 'basis', 'resp_implementation', 'resp_version', 'provenance'):
            require(isinstance(p.get(key), str) and bool(p[key].strip()), 'Missing quantum provenance: ' + key)
        opt = p.get('geometry_optimization', {})
        require(type(opt.get('executed')) is bool and bool(opt.get('settings')) and bool(opt.get('input_geometry_identity')), 'Missing geometry optimization settings')
        sites = graph.payload['sites']
        require(type(p.get('multiplicity')) is int and p['multiplicity'] > 0, 'Invalid multiplicity')
        require(type(p.get('formal_charge')) is int and p['formal_charge'] == sum(s['formal_charge'] for s in sites), 'Wrong quantum formal charge')
        require(p.get('final_graph_identity') == graph.identity, 'Stale quantum final graph')
        mapping = [{'atom_id': s['id'], 'element': s['element']} for s in sites]
        require(p.get('atom_mapping') == mapping and all(type(r['atom_id']) is int for r in p['atom_mapping']), 'Stale/reordered quantum atom mapping')
        require(p.get('normalization') is False and p.get('redistribution') is False, 'Charge repair is prohibited')
        # Every role is mandatory and separately fingerprinted. No inference from
        # a filename, a "RESP" label, or the total charge establishes provenance.
        roles = {'psi4_input', 'psi4_output', 'esp_grid', 'resp_input', 'resp_output',
                 'charge_vector', 'independent_verification', 'psi4_executable',
                 'resp_implementation_source', 'environment_manifest'}
        require(set(p.get('assets', {})) == roles, 'Incomplete quantum/fit evidence')
        paths = {}
        for role, asset in p['assets'].items():
            require(set(asset) == {'path', 'sha256'}, 'Malformed RESP asset')
            paths[role] = bound_file(directory, asset['path'], asset['sha256'])
        require(len(set(paths.values())) == len(roles), 'Evidence roles must refer to distinct assets')
        review = json.loads(paths['independent_verification'].read_text())
        require(review.get('schema') == 'island.resp-independent-verification.v1'
                and review.get('actual_quantum_and_fit_executed') is True
                and review.get('source_kind') == 'Psi4_ESP_RESP_fit'
                and bool(review.get('reviewer')) and bool(review.get('evidence_references')), 'Missing independent quantum/RESP verification')
        settings = {k: p[k] for k in ('psi4_version', 'method', 'basis', 'geometry_optimization', 'formal_charge', 'multiplicity', 'atom_mapping', 'resp_implementation', 'resp_version')}
        require(review.get('settings') == settings, 'Quantum/RESP settings differ from independent review')
        require(review.get('asset_hashes') == {k: v['sha256'] for k, v in p['assets'].items() if k != 'independent_verification'}, 'Independent source/intermediate verification mismatch')
        reference = p['reference']
        require(reference.get('method') == 'RESP/Psi4' and reference.get('provenance') == p['provenance'], 'Wrong RESP source/method provenance')
        require(reference.get('source_sha256') == p['assets']['psi4_output']['sha256']
                and reference.get('generator_sha256') == p['assets']['resp_implementation_source']['sha256'], 'Wrong quantum source/fit identity')
        require(p.get('evidence_references') == review['evidence_references'], 'Tampered evidence references')
        # Authoritative full-precision exported rows must match exactly; this
        # loader neither reformats nor repairs the independent fit's vector.
        vector = json.loads(paths['charge_vector'].read_text())
        require(vector == reference.get('rows'), 'Rounded, truncated or reordered RESP charge file')
        provided = GAFF2ProvidedCharges.from_reference(graph, reference, expected_reference_identity=reference['identity'])
        require(review.get('reference_identity') == reference['identity'], 'Stale independently reviewed charge vector')
        require(p.get('native_am1bcc_success') is False, 'RESP is not native AM1-BCC')
        result = {**p, 'provided_charge': provided.payload}
        return cls(json.dumps(result, sort_keys=True, allow_nan=False))

    def provided(self, directory, graph, *, expected_identity):
        """Only a validated load creates this bridge; no external OPLS fallback."""
        checked = type(self).load(directory, graph, expected_identity=expected_identity)
        require(checked.payload == self.payload, 'Loaded RESP result changed')
        return GAFF2ProvidedCharges(json.dumps(checked.payload['provided_charge'], sort_keys=True, allow_nan=False))
