"""Explicit, source-bound provided electrostatics for GAFF2.

No native charge generation, normalization, redistribution or decimal rounding.
A caller must supply an independently verified reference identity; a self-signed
charge vector alone cannot establish its scientific source.
"""

import json
import math
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal

from island.forcefields.ambertools.models import digest

SCHEMA = 'island.gaff2-provided-charge.v1'


class ProvidedChargeError(ValueError):
    """A failed provided-charge gate with a machine-readable classification."""

    def __init__(self, message, classification='provided_charge_source_unverified'):
        super().__init__(message)
        self.classification = classification


def require(ok, message):
    if not ok:
        raise ProvidedChargeError(message)


def reference_digest(value):
    import hashlib

    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def validate_reference(reference, graph, expected_identity):
    """Validate against the caller's independently established trust anchor."""
    require(isinstance(expected_identity, str) and len(expected_identity) == 64, 'Independent reference identity required')
    require(reference.get('identity') == expected_identity == reference_digest({k: v for k, v in reference.items() if k != 'identity'}), 'Changed independent charge reference')
    require(reference.get('final_graph_identity') == graph.identity, 'Stale final graph')
    require(reference.get('units') == 'elementary_charge' and reference.get('normalization') is False, 'Wrong units or normalization policy')
    for key in ('method', 'provenance', 'source_sha256', 'generator_sha256'):
        require(isinstance(reference.get(key), str) and bool(reference[key]), 'Missing source/method/provenance: ' + key)
    require(len(reference['source_sha256']) == len(reference['generator_sha256']) == 64, 'Missing source fingerprints')
    require(bool(reference.get('independent_checks')), 'Missing independent evidence')
    sites = graph.payload['sites']
    rows = reference.get('rows', [])
    require([r.get('atom_id') for r in rows] == [s['id'] for s in sites] and all(type(r['atom_id']) is int for r in rows), 'Missing, extra, duplicate or reordered stable atom IDs')
    values = {}
    for row, site in zip(rows, sites, strict=True):
        value = row.get('charge_e')
        require(type(value) in (float, int) and math.isfinite(value), 'One finite non-boolean charge per atom required')
        require(row.get('element') == site['element'], 'Atom element mapping changed')
        try:
            decimal = Decimal(row['decimal_e'])
        except Exception as error:
            raise ProvidedChargeError('Malformed authoritative charge decimal') from error
        require(decimal.is_finite() and float(decimal) == value, 'Rounded/truncated charge differs from source decimal')
        values[site['id']] = decimal
    components = []
    for ids in sorted(sorted(c) for c in graph.payload['components']):
        formal = sum(s['formal_charge'] for s in sites if s['id'] in ids)
        total = sum((values[i] for i in ids), Decimal(0))
        require(total == formal, 'Provided component charge differs from formal total')
        components.append({'atom_ids': ids, 'formal_total_e': formal, 'provided_total_e': str(total)})
    require(reference.get('components') == components, 'Changed component totals/membership')
    require(sum(values.values(), Decimal(0)) == sum(s['formal_charge'] for s in sites), 'Total charge differs from graph')


@dataclass(frozen=True)
class GAFF2ProvidedCharges:
    """JSON-owned immutable record; access returns a fresh, mutable copy."""

    raw: str

    @property
    def payload(self):
        return json.loads(self.raw)

    @property
    def identity(self):
        return self.payload['identity']

    def validate(self, graph, *, expected_reference_identity):
        p = self.payload
        if p.get('schema') == 'island.resp-finite-precision-charge.v1':
            from .resp_policy import FiniteRESPCharges

            return FiniteRESPCharges(self.raw).validate(graph, expected_reference_identity=expected_reference_identity)
        require(p.get('schema') == SCHEMA and p.get('identity') == digest({k: v for k, v in p.items() if k != 'identity'}), 'Provided record identity changed')
        validate_reference(p['reference'], graph, expected_reference_identity)
        require(p['final_graph_identity'] == graph.identity and p['reference_identity'] == expected_reference_identity, 'Changed graph/reference binding')
        require(all(type(r.get('atom_id')) is int and type(r.get('charge_e')) in (int, float) and math.isfinite(r['charge_e']) for r in p['rows']), 'Invalid provided atom ID/charge')
        require(digest(p['rows']) == digest(p['reference']['rows']), 'Provided vector differs from independent source; no rounding or normalization allowed')
        require(p['source_identity'] == p['reference']['source_sha256'] and p['method_identity'] == digest(p['reference']['method']), 'Source/method mutation')
        require(p['diagnostic'] == 'provided_charge_validated' and p['native_am1bcc_success'] is False and p['units'] == 'elementary_charge', 'Contradictory charge mode')
        require(p['evidence_references'] == [expected_reference_identity] and p['provenance'] == p['reference']['provenance'], 'Changed evidence/provenance')
        return {r['atom_id']: r['charge_e'] for r in p['rows']}

    @classmethod
    def from_reference(cls, graph, reference, *, expected_reference_identity):
        validate_reference(reference, graph, expected_reference_identity)
        body = {'schema': SCHEMA, 'final_graph_identity': graph.identity, 'reference_identity': expected_reference_identity, 'reference': deepcopy(reference), 'rows': deepcopy(reference['rows']), 'units': 'elementary_charge', 'source_identity': reference['source_sha256'], 'method_identity': digest(reference['method']), 'provenance': reference['provenance'], 'evidence_references': [expected_reference_identity], 'diagnostic': 'provided_charge_validated', 'native_am1bcc_success': False}
        result = cls(json.dumps({**body, 'identity': digest(body)}, sort_keys=True, allow_nan=False))
        result.validate(graph, expected_reference_identity=expected_reference_identity)
        return result
