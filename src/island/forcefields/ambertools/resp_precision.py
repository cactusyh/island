"""Lossless native RESP audit and reviewed-charge loader; no charge repair.

Exact rational arithmetic distinguishes a native float64 constraint residual
from serialization loss. A failed audit cannot publish provided charges.
Historical RESP/provided schemas and tolerances remain unchanged.
"""

import json
from dataclasses import dataclass
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np

from island.forcefields.ambertools.provided import (
    ProvidedChargeError,
    reference_digest,
    require,
)
from island.forcefields.ambertools.psi4_execution import load_execution
from island.forcefields.ambertools.psi4_runtime import validate_runtime
from island.forcefields.ambertools.resp import RESPCharges, file_hash

SCHEMA = 'island.resp-native-precision.v1'


def signed(body):
    return {**body, 'identity': reference_digest(body)}


def exact_decimal(value):
    """Exact terminating decimal for a binary rational, independent of context."""
    denominator = value.denominator
    require(denominator & (denominator - 1) == 0, 'Expected binary source rational')
    exponent = denominator.bit_length() - 1
    coefficient = value.numerator * 5 ** exponent
    return str(Decimal((int(coefficient < 0), tuple(int(c) for c in str(abs(coefficient))), -exponent)))


def audit_native(directory, graph, *, expected_identity, expected_input, runtime_identity, worker_sha256):
    directory = Path(directory)
    # Every finite binary64 value has a terminating decimal spanning fewer than
    # 1,400 digits. This context prevents ambient Decimal rounding; it changes
    # neither source values nor the zero charge tolerance.
    with localcontext() as context:
        context.prec = 2000 + len(str(len(graph.payload['sites'])))
        execution = load_execution(directory, graph, expected_identity=expected_identity, expected_input=expected_input,
                                   runtime_identity=runtime_identity, worker_sha256=worker_sha256)
    require(execution['charge_generation_executed'] is True, 'Native RESP generation did not execute')
    rows, components, source_total = [], [], Fraction(0)
    for index, component in enumerate(expected_input['components']):
        root = directory / f'component-{index}'
        atoms = component['atoms']
        # The independent worker's geometry text uses round-trip 17 significant
        # digits. Inspect it directly, rather than trusting copied mapping labels.
        lines = (root / 'geometry-input.dat').read_text().splitlines()
        require(lines[0] == f"{component['formal_charge']} {component['multiplicity']}", 'Changed formal charge/multiplicity input')
        xyz = np.array([a['position_angstrom'] for a in atoms])
        centered = xyz - xyz.mean(axis=0)
        geometry_rows = [line.split() for line in lines[1:len(atoms) + 1]]
        require([r[0] for r in geometry_rows] == [a['element'] for a in atoms], 'Native geometry element ordering changed')
        require(all(len(r) == 4 for r in geometry_rows) and np.array_equal(np.array([[float(v) for v in r[1:]] for r in geometry_rows]), centered), 'Native input coordinate mapping changed')
        require(lines[len(atoms) + 1:] == ['units angstrom', 'symmetry c1', 'no_reorient', 'no_com'], 'Native geometry orientation/units policy changed')
        values = np.load(root / 'native-resp-charges.npy', allow_pickle=False)
        saved = json.loads((root / 'charge-vector.json').read_text())
        require(all(type(r['charge_e']) is float and type(r['atom_id']) is int for r in saved), 'Boolean/nonfloat charge or atom ID rejected')
        with localcontext() as context:
            context.prec = 2000 + len(str(len(atoms)))
            repr_decimal_total = str(sum((Decimal(r['decimal_e']) for r in saved), Decimal(0)))
        # load_execution checks dtype, shape, finiteness, IDs, elements, native
        # optimized geometry, graph, source files, and all artifact hashes.
        require(np.array_equal(values.view(np.uint64), np.array([r['charge_e'] for r in saved], dtype=np.float64).view(np.uint64)), 'Charge import lost native binary bits')
        native_sum = sum((Fraction.from_float(float(q)) for q in values), Fraction(0))
        repr_sum = sum((Fraction(Decimal(r['decimal_e'])) for r in saved), Fraction(0))
        formal = component['formal_charge']
        source_total += native_sum
        components.append({'atom_ids': [a['atom_id'] for a in atoms], 'formal_total_e': formal,
                           'native_exact_total_e': exact_decimal(native_sum),
                           'native_total_fraction': [native_sum.numerator, native_sum.denominator],
                           'repr_decimal_total_e': repr_decimal_total,
                           'repr_decimal_total_fraction': [repr_sum.numerator, repr_sum.denominator],
                           'exact_conservation': native_sum == formal, 'json_roundtrip_bit_exact': True,
                           'input_coordinate_mapping_exact': True,
                           'input_geometry_sha256': file_hash(root / 'geometry-input.dat'),
                           'optimized_geometry_sha256': file_hash(root / 'optimized-geometry-bohr.npy'),
                           'native_array_sha256': file_hash(root / 'native-resp-charges.npy')})
        rows.extend({'atom_id': a['atom_id'], 'element': a['element'], 'charge_e': float(q),
                     'decimal_e': str(Decimal.from_float(float(q))), 'ieee754_hex': float(q).hex()} for a, q in zip(atoms, values, strict=True))
    rows.sort(key=lambda row: row['atom_id'])
    formal_total = sum(c['formal_charge'] for c in expected_input['components'])
    exact = source_total == formal_total and all(c['exact_conservation'] for c in components)
    # Existing provided-charge import remains a second, unchanged gate. An
    # exact-native vector alone does not establish reviewed RESP provenance.
    return signed({'schema': SCHEMA, 'execution_identity': execution['identity'], 'final_graph_identity': graph.identity,
                   'input_identity': execution['input_identity'], 'runtime_identity': runtime_identity,
                   'worker_sha256': worker_sha256, 'optimized_geometry_identity': execution['optimized_geometry_identity'],
                   'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4',
                   'precision_source': 'native resp.resp API float64 ndarray; native-resp-charges.npy',
                   'source_dtype': 'float64', 'source_significand_bits': 53, 'units': 'elementary_charge',
                   'rows': rows, 'components': components, 'formal_system_total_e': formal_total,
                   'native_exact_system_total_e': exact_decimal(source_total),
                   'native_system_total_fraction': [source_total.numerator, source_total.denominator],
                   'charge_conservation_passed': exact, 'acceptance_met': False,
                   'classification': 'native_charge_totals_exact' if exact else 'authoritative_native_precision_insufficient',
                   'stage': 'native_component_and_system_total_validation', 'serialization_loss_detected': False,
                   'normalization': False, 'redistribution': False, 'charge_tolerance_e': '0',
                   'native_am1bcc_success': False, 'dynamics_executed': False,
                   'production_validated': False, 'simulation_readiness': 'not_established'})


def load_audit(directory, graph, *, expected_audit_identity, **kwargs):
    saved = json.loads((Path(directory) / 'precision-audit.json').read_text())
    rebuilt = audit_native(directory, graph, **kwargs)
    require(saved == rebuilt and rebuilt['identity'] == expected_audit_identity, 'Changed native precision/provenance audit')
    return rebuilt


@dataclass(frozen=True)
class PrecisionRESPCharges:
    """Accepted native bridge additionally requires the existing independent review."""

    raw: str

    @property
    def payload(self):
        return json.loads(self.raw)

    @classmethod
    def load(cls, directory, graph, *, runtime, expected_runtime_identity,
             expected_audit_identity, review_directory, expected_review_identity, **kwargs):
        loader = load_in_memory_audit if 'native_directory' in kwargs else load_audit
        audit = loader(directory, graph, expected_audit_identity=expected_audit_identity, **kwargs)
        require(audit['runtime_identity'] == expected_runtime_identity, 'Native runtime binding changed')
        validate_runtime(runtime, expected_identity=expected_runtime_identity)
        if not audit['charge_conservation_passed']:
            raise ProvidedChargeError('Authoritative float64 RESP vector has nonzero exact component/system residuals; lossless import cannot repair them. See precision-audit.json.', 'authoritative_native_precision_insufficient')
        with localcontext() as context:
            context.prec = 2000 + len(str(len(audit['rows'])))
            reviewed = RESPCharges.load(review_directory, graph, expected_identity=expected_review_identity)
        expected_rows = [{k: r[k] for k in ('atom_id', 'element', 'charge_e', 'decimal_e')} for r in audit['rows']]
        require(reviewed.payload['reference']['rows'] == expected_rows, 'Reviewed RESP vector differs from exact native source')
        require(reviewed.payload['psi4_version'] == runtime['psi4_version'] and reviewed.payload['method'] == kwargs['expected_input']['method']
                and reviewed.payload['basis'] == kwargs['expected_input']['basis'], 'Reviewed native settings differ')
        body = {'schema': 'island.precision-resp-charge.v1', 'audit_identity': audit['identity'],
                'review_identity': expected_review_identity, 'runtime_identity': expected_runtime_identity,
                'charge_origin': 'resp_psi4', 'native_am1bcc_success': False, 'reviewed_result': reviewed.payload}
        return cls(json.dumps(signed(body), sort_keys=True, allow_nan=False))

    def provided(self, directory, graph, **kwargs):
        checked = type(self).load(directory, graph, **kwargs)
        require(checked.payload == self.payload, 'Accepted native RESP record changed')
        with localcontext() as context:
            context.prec = 2000 + len(str(len(checked.payload['reviewed_result']['reference']['rows'])))
            return RESPCharges(json.dumps(checked.payload['reviewed_result'], sort_keys=True, allow_nan=False)).provided(
                kwargs['review_directory'], graph, expected_identity=kwargs['expected_review_identity'])


def audit_in_memory(directory, graph, *, native_directory, in_memory_worker_sha256, **kwargs):
    """Validate a fresh native ESP/RESP fit anchored to retained optimized geometry."""
    from island.forcefields.ambertools.resp import bound_file
    directory = Path(directory)
    original = audit_native(native_directory, graph, **kwargs)
    data = kwargs['expected_input']
    raw = json.loads((directory / 'result.json').read_text())
    require(raw['schema'] == 'island.native-esp-resp.v1' and raw['charge_origin'] == 'resp_psi4' and raw['charge_method'] == 'RESP', 'Wrong native ESP charge origin')
    require(raw['final_graph_identity'] == graph.identity and raw['method'] == data['method'] and raw['basis'] == data['basis'], 'Changed native ESP graph/method/basis')
    require(raw['psi4_version'] == '1.9.1' and raw['resp_module_version'] == '0.8', 'Wrong native ESP runtime version')
    require(raw['normalization'] is False and raw['redistribution'] is False and raw['dynamics_executed'] is False and raw['acceptance_met'] is False, 'Unproven native ESP outcome')
    require(raw['input_sha256'] == file_hash(Path(native_directory) / 'input.json') == file_hash(directory / 'input.json'), 'Changed native ESP input')
    require(file_hash(directory / 'worker-source.py') == in_memory_worker_sha256, 'Changed native ESP worker')
    required = {'input.json'}
    for i in range(len(data['components'])):
        required.update(f'component-{i}/{name}' for name in ('geometry-bohr.npy', 'geometry-input.dat', 'grid-bohr.npy', 'esp-native.npy', 'fit-options.json', 'fit-input.npz', 'charges-native.npy', 'charges-native.json', 'psi4.out', 'grid.dat', 'grid_esp.dat'))
    require(required <= set(raw['artifacts']), 'Missing native ESP evidence')
    for name, expected in raw['artifacts'].items():
        bound_file(directory, name, expected)
    require(len(raw['components']) == len(data['components']), 'Incomplete native ESP component inventory')
    rows, components, total = [], [], Fraction(0)
    for index, component in enumerate(data['components']):
        root = directory / f'component-{index}'
        atoms = component['atoms']
        entry = raw['components'][index]
        geometry = np.load(root / 'geometry-bohr.npy', allow_pickle=False)
        require(file_hash(root / 'geometry-bohr.npy') == file_hash(Path(native_directory) / f'component-{index}/optimized-geometry-bohr.npy') == entry['source_geometry_sha256'] == entry['geometry_sha256'], 'Stale optimized geometry source/mapping')
        require(geometry.shape == (len(atoms), 3) and np.isfinite(geometry).all(), 'Invalid native ESP geometry')
        require(entry['atom_ids'] == [a['atom_id'] for a in atoms] and entry['elements'] == [a['element'] for a in atoms]
                and type(entry['formal_charge']) is int and entry['formal_charge'] == component['formal_charge'] and entry['multiplicity'] == component['multiplicity'], 'Changed native ESP atom/formal-charge mapping')
        expected_geometry = f"{component['formal_charge']} {component['multiplicity']}\n" + '\n'.join(f"{a['element']} {r[0]:.17g} {r[1]:.17g} {r[2]:.17g}" for a, r in zip(atoms, geometry, strict=True)) + '\nunits bohr\nsymmetry c1\nno_com\nno_reorient\n'
        require((root / 'geometry-input.dat').read_text() == expected_geometry, 'Native ESP coordinate mapping changed')
        charges = np.load(root / 'charges-native.npy', allow_pickle=False)
        require(charges.dtype == np.float64 and charges.shape == (len(atoms),) and np.isfinite(charges).all(), 'Invalid native ESP/RESP charges')
        current = [{'atom_id': a['atom_id'], 'element': a['element'], 'charge_e': float(q), 'ieee754_hex': float(q).hex()} for a, q in zip(atoms, charges, strict=True)]
        saved = json.loads((root / 'charges-native.json').read_text())
        require(saved == current and all(type(r['charge_e']) is float and type(r['atom_id']) is int for r in saved), 'Reordered/rounded native ESP charge vector')
        require(entry['charge_array_sha256'] == file_hash(root / 'charges-native.npy'), 'Changed native ESP charge array')
        points = np.load(root / 'grid-bohr.npy', allow_pickle=False)
        esp = np.load(root / 'esp-native.npy', allow_pickle=False)
        require(points.dtype == esp.dtype == np.float64 and points.ndim == 2 and points.shape[1] == 3 and esp.shape == (len(points),) and np.isfinite(points).all() and np.isfinite(esp).all(), 'Invalid native grid/ESP array')
        require(np.array_equal(np.loadtxt(root / 'grid.dat'), points), 'Native ESP grid text lost coordinates')
        with np.load(root / 'fit-input.npz', allow_pickle=False) as fit:
            require(np.array_equal(fit['esp'].view(np.uint64), esp.view(np.uint64)) and fit['invr'].shape == (len(points), len(atoms)) and np.isfinite(fit['invr']).all(), 'Native ESP fitting import changed source values')
        options = {**data['resp_options'], 'CONSTRAINT_CHARGE': [], 'CONSTRAINT_GROUP': [], 'WEIGHT': [1], 'RESTRAINT': True}
        require(json.loads((root / 'fit-options.json').read_text()) == options, 'Changed native fitting options')
        difference = float(np.max(abs(esp - np.loadtxt(root / 'grid_esp.dat'))))
        require(difference == entry['native_esp_vs_official_text_max_difference_au'] and difference <= 5.1e-11, 'Official ESP API units/semantics check failed')
        require(entry['fit_labels'] == ['ESP', 'RESP'] and entry['fit_note'] == '', 'Native fit failed to converge')
        value = sum((Fraction.from_float(float(q)) for q in charges), Fraction(0))
        require(entry['total_fraction'] == [value.numerator, value.denominator] and entry['exact_conservation'] is (value == component['formal_charge']), 'Native charge-total result changed')
        total += value
        components.append({'atom_ids': entry['atom_ids'], 'formal_total_e': entry['formal_charge'], 'native_total_fraction': entry['total_fraction'],
                           'native_exact_total_e': exact_decimal(value), 'exact_conservation': value == component['formal_charge'],
                           'native_array_sha256': entry['charge_array_sha256'], 'optimized_geometry_sha256': entry['geometry_sha256']})
        rows.extend({'atom_id': a['atom_id'], 'element': a['element'], 'charge_e': float(q), 'decimal_e': str(Decimal.from_float(float(q))), 'ieee754_hex': float(q).hex()} for a, q in zip(atoms, charges, strict=True))
    rows.sort(key=lambda r: r['atom_id'])
    formal = sum(c['formal_charge'] for c in data['components'])
    exact = total == formal and all(c['exact_conservation'] for c in components)
    require(raw['charge_conservation_passed'] is exact and raw['classification'] == ('native_charge_totals_exact' if exact else 'authoritative_native_precision_insufficient'), 'Contradictory native ESP charge acceptance')
    return signed({'schema': 'island.resp-native-memory-precision.v1', 'original_audit_identity': original['identity'],
                   'final_graph_identity': graph.identity, 'runtime_identity': kwargs['runtime_identity'],
                   'optimized_geometry_identity': original['optimized_geometry_identity'],
                   'native_result_sha256': file_hash(directory / 'result.json'), 'worker_sha256': in_memory_worker_sha256,
                   'charge_origin': 'resp_psi4', 'charge_method': 'RESP', 'quantum_software': 'Psi4',
                   'method': data['method'], 'basis': data['basis'], 'precision_source': 'Psi4 ESPPropCalc native Vector -> CDS espfit.fit ndarray',
                   'source_dtype': 'float64', 'source_significand_bits': 53, 'rows': rows, 'components': components,
                   'native_exact_system_total_e': exact_decimal(total), 'native_system_total_fraction': [total.numerator, total.denominator],
                   'formal_system_total_e': formal, 'charge_conservation_passed': exact,
                   'classification': raw['classification'], 'stage': 'native_component_and_system_total_validation',
                   'serialization_loss_detected': False, 'normalization': False, 'redistribution': False,
                   'acceptance_met': False, 'charge_tolerance_e': '0', 'native_am1bcc_success': False, 'dynamics_executed': False,
                   'production_validated': False, 'simulation_readiness': 'not_established'})


def load_in_memory_audit(directory, graph, *, expected_audit_identity, **kwargs):
    saved = json.loads((Path(directory) / 'precision-audit.json').read_text())
    rebuilt = audit_in_memory(directory, graph, **kwargs)
    require(saved == rebuilt and rebuilt['identity'] == expected_audit_identity, 'Changed native in-memory precision audit')
    return rebuilt
