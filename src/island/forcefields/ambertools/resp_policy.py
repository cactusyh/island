"""Versioned finite-binary64 RESP conservation, never charge repair.

The absolute allowance is the sum of half-ULPs of the source values. The
relative allowance is Higham's gamma_n reduction bound, scaled by sum(abs(q)).
This is a floating-point consistency check, not a chemical accuracy claim.
"""

import json
import math
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

from island.forcefields.ambertools.provided import reference_digest, require
from island.forcefields.ambertools.resp_precision import exact_decimal, signed

SCHEMA = 'island.resp-finite-precision-charge.v1'
POLICY = 'island.resp-binary64-conservation-policy.v1'
UNIT_ROUNDOFF = Fraction(1, 2**53)


def measure(values, formal):
    require(values and all(type(q) is float and math.isfinite(q) for q in values), 'Finite nonboolean native binary64 charges required')
    require(type(formal) is int, 'Integer formal charge required')
    n = len(values)
    require(n < 1000000, 'Unsupported reduction size')
    exact = sum((Fraction.from_float(q) for q in values), Fraction())
    l1 = sum((abs(Fraction.from_float(q)) for q in values), Fraction())
    absolute = sum((Fraction.from_float(math.ulp(q)) / 2 for q in values), Fraction())
    relative = n * UNIT_ROUNDOFF / (1 - n * UNIT_ROUNDOFF)
    bound = absolute + relative * l1
    residual = exact - formal
    with localcontext() as ctx:
        ctx.prec = 2200
        decimal = sum((Decimal.from_float(q) for q in values), Decimal())
    return {'atom_count': n, 'formal_total_e': formal, 'native_exact_total_e': exact_decimal(exact),
            'residual_e': exact_decimal(residual), 'absolute_residual_e': exact_decimal(abs(residual)),
            'relative_residual_l1': float(abs(residual) / l1) if l1 else (0.0 if not residual else None),
            'relative_residual_formal': float(abs(residual) / abs(formal)) if formal else None,
            'naive_sum_e': sum(values), 'compensated_fsum_e': math.fsum(values),
            'decimal_from_binary_total_e': str(decimal), 'decimal_source_text': False,
            'l1_charge_e': exact_decimal(l1), 'absolute_tolerance_e': exact_decimal(absolute),
            'relative_tolerance': float(relative), 'relative_tolerance_fraction': [relative.numerator, relative.denominator],
            'allowed_residual_e': float(bound), 'residual_fraction': [residual.numerator, residual.denominator],
            'bound_fraction': [bound.numerator, bound.denominator],
            'within_source_half_ulp_envelope': abs(residual) <= absolute,
            'accepted': abs(residual) <= bound}


def make_policy(audit, runtime):
    require(audit['charge_origin'] == 'resp_psi4' and audit['quantum_software'] == 'Psi4' and audit['charge_method'] == 'RESP', 'RESP/Psi4 origin required')
    require(audit['source_dtype'] == 'float64' and audit['source_significand_bits'] == 53, 'Unsupported source precision')
    require(runtime['identity'] == audit['runtime_identity'] == reference_digest({k: v for k, v in runtime.items() if k != 'identity'}), 'Wrong executable/source environment')
    rows = audit['rows']
    components = [measure([r['charge_e'] for r in rows if r['atom_id'] in c['atom_ids']], c['formal_total_e']) for c in audit['components']]
    system = measure([r['charge_e'] for r in rows], audit['formal_system_total_e'])
    return signed({'schema': POLICY, 'charge_origin': 'RESP/Psi4', 'source_precision': 'IEEE-754 binary64, 53 significand bits',
                   'output_precision': 'lossless native NPY binary64; JSON/hex roundtrip; no rounded table authority',
                   'source_and_executable_hashes': runtime['files'], 'runtime_identity': runtime['identity'],
                   'summation_method': 'exact binary Fraction decision; math.fsum compensated reduction; Decimal.from_float cross-check',
                   'absolute_tolerance_rule': 'sum(0.5*ulp(q_i)) e, round-to-nearest source representation envelope',
                   'relative_tolerance_rule': 'gamma_n=n*u/(1-n*u), u=2^-53; scale=sum(abs(q_i)); bound=atol+gamma_n*scale',
                   'justification': 'Binary64 representation and standard forward error bound for finite floating-point summation (Higham, Accuracy and Stability of Numerical Algorithms, 2nd ed., section 3.1). Residual is in source values, not solely accumulation error. Native constrained float64 fit plus rounding is numerically consistent when within this declared envelope; no chemical accuracy inferred.',
                   'evidence_references': [audit['identity'], runtime['identity']], 'components': components, 'system': system,
                   'accepted': system['accepted'] and all(c['accepted'] for c in components),
                   'normalization': False, 'redistribution': False})


def create(audit, graph, runtime, independent):
    require(audit['identity'] == reference_digest({k: v for k, v in audit.items() if k != 'identity'}), 'Changed source audit identity')
    policy = make_policy(audit, runtime)
    rows = audit['rows']
    require([r['atom_id'] for r in rows] == [s['id'] for s in graph.payload['sites']], 'Wrong source atom order')
    for row, site in zip(rows, graph.payload['sites'], strict=True):
        require(type(row['atom_id']) is int and row['element'] == site['element'] and row['ieee754_hex'] == row['charge_e'].hex(), 'Atom mapping/value changed')
    require(audit['final_graph_identity'] == graph.identity, 'Stale graph')
    expected_components = sorted(sorted(c) for c in graph.payload['components'])
    require([c['atom_ids'] for c in audit['components']] == expected_components, 'Changed component membership')
    require([c['formal_total_e'] for c in audit['components']] == [sum(s['formal_charge'] for s in graph.payload['sites'] if s['id'] in ids) for ids in expected_components], 'Changed component formal totals')
    require(audit['formal_system_total_e'] == sum(s['formal_charge'] for s in graph.payload['sites']), 'Changed system formal total')
    require(independent['source_audit_identity'] == audit['identity'] and independent['rows'] == rows and independent['mapping_verified'] is True, 'Independent raw parsing mismatch')
    require(independent['component_residual_fractions'] == [c['residual_fraction'] for c in policy['components']] and independent['system_residual_fraction'] == policy['system']['residual_fraction'], 'Independent component/total mismatch')
    require(independent['identity'] == reference_digest({k: v for k, v in independent.items() if k != 'identity'}), 'Changed independent evidence')
    require(policy['accepted'], 'RESP conservation exceeds declared finite-precision bound')
    return FiniteRESPCharges(json.dumps(signed({'schema': SCHEMA, 'final_graph_identity': graph.identity,
        'charge_origin': 'resp_psi4', 'charge_method': 'RESP/Psi4', 'native_am1bcc_success': False,
        'units': 'elementary_charge', 'source_audit': audit, 'runtime': runtime,
        'independent': independent, 'policy': policy, 'rows': rows,
        'normalization': False, 'redistribution': False}), sort_keys=True, allow_nan=False))


@dataclass(frozen=True)
class FiniteRESPCharges:
    raw: str

    @property
    def payload(self):
        return json.loads(self.raw)

    @property
    def identity(self):
        return self.payload['identity']

    def validate(self, graph, *, expected_reference_identity):
        p = self.payload
        require(p['schema'] == SCHEMA and p['identity'] == expected_reference_identity == reference_digest({k: v for k, v in p.items() if k != 'identity'}), 'Changed RESP policy/result trust anchor')
        require(p['charge_origin'] == 'resp_psi4' and p['charge_method'] == 'RESP/Psi4' and p['native_am1bcc_success'] is False, 'Wrong charge origin')
        require(p['units'] == 'elementary_charge' and p['normalization'] is p['redistribution'] is False, 'Wrong units/charge repair')
        audit = p['source_audit']
        require(audit['identity'] == reference_digest({k: v for k, v in audit.items() if k != 'identity'}), 'Changed native audit')
        require(p['policy'] == make_policy(audit, p['runtime']), 'Changed policy, tolerances, source or executable')
        rebuilt = create(audit, graph, p['runtime'], p['independent'])
        require(rebuilt.payload == p, 'Changed RESP charge vector/provenance')
        return {r['atom_id']: r['charge_e'] for r in p['rows']}

    def save(self, path):
        with Path(path).open('x') as stream:
            stream.write(self.raw + '\n')

    @classmethod
    def load(cls, path, graph, *, expected_identity):
        result = cls(Path(path).read_text())
        result.validate(graph, expected_reference_identity=expected_identity)
        return result
