"""Independent raw RESP provenance review using only the Python standard library.

No ISLAND/importer/policy imports. Trust is rooted in reviewed external archive
bytes, not in self-signed charge metadata. Offline review does not claim that a
hash can authenticate an untrusted external calculation or reproduce quantum MD.
"""

import argparse
import ast
import hashlib
import json
import math
import re
import struct
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path

# C3 trust roots are supplied explicitly and pinned by the additive receipt.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def read(path):
    return json.loads(Path(path).read_text())


def array(path):
    """Parse native NPY independently of NumPy, including coordinate matrices."""
    raw = Path(path).read_bytes()
    require(raw[:8] == b'\x93NUMPY\x01\x00', 'NPY v1 required')
    size = struct.unpack('<H', raw[8:10])[0]
    header = ast.literal_eval(raw[10:10 + size].decode('latin1').strip())
    shape = header['shape']
    require(header['descr'] == '<f8' and header['fortran_order'] is False and 1 <= len(shape) <= 2, 'Ordered binary64 required')
    count = math.prod(shape)
    require(len(raw) == 10 + size + 8 * count, 'Truncated/extended native array')
    values = list(struct.unpack('<' + 'd' * count, raw[10 + size:]))
    require(all(math.isfinite(v) for v in values), 'Nonfinite native values')
    return values, shape


def bits(value):
    require(type(value) is float and math.isfinite(value), 'Finite nonboolean float required')
    return struct.pack('<d', value).hex()


def residual(values, formal):
    """Independent rational reduction and C2 bound reconstruction."""
    require(type(formal) is int, 'Integer formal total required')
    exact = sum(map(Fraction.from_float, values), Fraction())
    l1 = sum((abs(Fraction(v)) for v in values), Fraction())
    half_ulp = sum((Fraction(math.ulp(v)) / 2 for v in values), Fraction())
    nu = Fraction(len(values), 2**53)
    gamma = nu / (1 - nu)
    bound = half_ulp + gamma * l1
    delta = exact - formal
    with localcontext() as context:
        context.prec = 2200
        decimal_total = sum(map(Decimal.from_float, values), Decimal())
    require(Fraction(decimal_total) == exact, 'Independent Decimal/rational disagreement')
    return {'atom_count': len(values), 'formal_total_e': formal,
            'residual_fraction': [delta.numerator, delta.denominator],
            'bound_fraction': [bound.numerator, bound.denominator],
            'absolute_tolerance_e': str(float(half_ulp)),
            'absolute_residual_e': float(abs(delta)), 'relative_residual_l1': float(abs(delta) / l1),
            'compensated_fsum_e': math.fsum(values), 'decimal_binary_total_e': str(decimal_total),
            'accepted': abs(delta) <= bound}


def compare_policy(values, formal, policy):
    measured = residual(values, formal)
    require(measured['accepted'], 'Raw charge conservation exceeds C2 bound')
    for field in ('atom_count', 'formal_total_e', 'residual_fraction', 'bound_fraction', 'compensated_fsum_e', 'accepted'):
        require(policy[field] == measured[field], 'C2 policy mismatch: ' + field)
    require(Fraction(Decimal(policy['absolute_tolerance_e'])) == sum((Fraction(math.ulp(v)) / 2 for v in values), Fraction()), 'Changed absolute tolerance')
    require(Fraction(Decimal(policy['absolute_residual_e'])) == abs(Fraction(*measured['residual_fraction'])), 'Changed residual')
    gamma = Fraction(len(values), 2**53 - len(values))
    require(policy['relative_tolerance_fraction'] == [gamma.numerator, gamma.denominator], 'Changed relative tolerance')
    return measured


def audit(root, kind, charges, runtime, anchors):
    root = Path(root)
    anchor = anchors['cases'][kind]
    for name, expected in anchor['files'].items():
        require(sha(root / name) == expected, 'Unreviewed raw source: ' + name)
    require(digest({k: v for k, v in runtime.items() if k != 'identity'}) == runtime['identity'], 'Changed runtime identity')
    # Match full runtime contents independently, including executable/source hashes.
    pinned = read(Path('docs/evidence/phase_4j36/environment.json'))
    require(sha('docs/evidence/phase_4j36/environment.json') == anchors['runtime_sha256'] and runtime == pinned, 'Changed executable/source provenance')
    require(charges['charge_origin'] == 'resp_psi4' and charges['charge_method'] == 'RESP/Psi4' and charges['native_am1bcc_success'] is False, 'Wrong charge origin')
    require(charges['runtime'] == runtime and charges['normalization'] is charges['redistribution'] is False and charges['units'] == 'elementary_charge', 'Changed runtime/units/charge repair')
    require(charges['identity'] == digest({k: v for k, v in charges.items() if k != 'identity'}), 'Changed charge identity')
    original = root / ('quantum-' + kind)
    memory = root / ('native-esp-' + kind)
    data, result = read(memory / 'input.json'), read(memory / 'result.json')
    execution = read(original / 'execution-record.json')
    graph = data['graph']
    chemical = {k: v for k, v in graph.items() if k not in ('graph_identity', 'history_identity')}
    chemical['transformations'] = []
    graph_id = digest({'payload': chemical, 'sha256': digest(chemical)})
    require(graph_id == graph['graph_identity'] == data['final_graph_identity'] == charges['final_graph_identity'], 'Changed final graph')
    require(execution['runtime_identity'] == runtime['identity'] and execution['charge_generation_executed'] is True, 'Wrong external execution provenance')
    require(result['psi4_version'] == runtime['psi4_version'] == '1.9.1' and result['resp_module_version'] == runtime['resp_module_version'] == '0.8', 'Wrong quantum versions')
    require(result['charge_origin'] == 'resp_psi4' and result['charge_method'] == 'RESP' and result['method'] == data['method'] == 'HF' and result['basis'] == data['basis'] == '6-31G*', 'Wrong quantum method/origin')
    for directory, record, key in ((memory, result, 'artifacts'), (original, execution, 'artifact_sha256')):
        for name, expected in record[key].items():
            path = (directory / name).resolve()
            require(path.is_relative_to(directory.resolve()) and sha(path) == expected, 'Mutated raw intermediate: ' + name)
    adj = {a['id']: set() for a in graph['sites']}
    sites = {a['id']: a for a in graph['sites']}
    for bond in graph['bonds']:
        a, b = bond['sites'];adj[a].add(b);adj[b].add(a)
    components, remaining = [], set(adj)
    while remaining:
        pending, seen = [min(remaining)], set()
        while pending:
            a = pending.pop()
            if a not in seen:
                seen.add(a);pending.extend(adj[a] - seen)
        remaining -= seen;components.append(sorted(seen))
    require(components == graph['components'] == [[a['atom_id'] for a in c['atoms']] for c in data['components']], 'Component membership mismatch')
    policy = charges['policy']
    require(policy['schema'] == 'island.resp-binary64-conservation-policy.v1' and policy['charge_origin'] == 'RESP/Psi4' and policy['source_and_executable_hashes'] == runtime['files'], 'Wrong finite precision policy provenance')
    allvalues, rows, summaries, coordinates = [], [], [], []
    for index, component in enumerate(data['components']):
        atoms = component['atoms'];ids = components[index];directory = memory / f'component-{index}'
        require([a['element'] for a in atoms] == [sites[i]['element'] for i in ids], 'Element mapping mismatch')
        formal = sum(sites[i]['formal_charge'] for i in ids)
        require(component['formal_charge'] == formal == sum(a['formal_charge'] for a in atoms), 'Formal/component total mismatch')
        initial = (original / f'component-{index}/geometry-input.dat').read_text().splitlines()
        require(initial[0] == f"{formal} {component['multiplicity']}", 'Initial formal charge changed')
        centers = [math.fsum(a['position_angstrom'][j] for a in atoms) / len(atoms) for j in range(3)]
        for atom, line in zip(atoms, initial[1:len(ids) + 1], strict=True):
            fields = line.split()
            require(fields[0] == atom['element'] and max(abs(float(fields[j + 1]) - (atom['position_angstrom'][j] - centers[j])) for j in range(3)) <= 5e-14, 'Initial input coordinate mapping changed')
        geometry, shape = array(directory / 'geometry-bohr.npy')
        require(shape == (len(ids), 3) and sha(directory / 'geometry-bohr.npy') == sha(original / f'component-{index}/optimized-geometry-bohr.npy'), 'Stale optimized geometry')
        xyz = [geometry[i:i + 3] for i in range(0, len(geometry), 3)]
        lines = (directory / 'geometry-input.dat').read_text().splitlines()
        require(lines[0] == f"{formal} {component['multiplicity']}" and lines[len(ids) + 1:] == ['units bohr', 'symmetry c1', 'no_com', 'no_reorient'], 'Geometry units/formal charge policy changed')
        for atom, point, line in zip(atoms, xyz, lines[1:len(ids) + 1], strict=True):
            fields = line.split()
            require(fields[0] == atom['element'] and [float(v) for v in fields[1:]] == point, 'Geometry atom/coordinate order mismatch')
        # Independent native Psi4 log cross-check: this printed geometry is a
        # mapping witness, never authoritative charge/geometry precision.
        log = (directory / 'psi4.out').read_text()
        marker = f"Geometry (in Bohr), charge = {formal}, multiplicity = {component['multiplicity']}:"
        require(marker in log and 'Energy and wave function converged.' in log and 'Name: 6-31G*' in log and 'RHF Reference' in log, 'Quantum log/settings incomplete')
        printed = re.findall(r'^\s+([CHO])\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+[\d.]+\s*$', log.split(marker, 1)[1], re.MULTILINE)[:len(ids)]
        require(len(printed) == len(ids), 'Missing Psi4 coordinate witness')
        require(all(r[0] == a['element'] and max(abs(float(v) - q) for v, q in zip(r[1:], point, strict=True)) <= 5.1e-12 for r, a, point in zip(printed, atoms, xyz, strict=True)), 'Psi4 log geometry mismatch')
        values, shape = array(directory / 'charges-native.npy')
        require(shape == (len(ids),), 'Wrong vector size')
        saved = read(directory / 'charges-native.json')
        require([r['atom_id'] for r in saved] == ids and [r['element'] for r in saved] == [a['element'] for a in atoms], 'Reordered raw charges')
        require([bits(r['charge_e']) for r in saved] == [bits(v) for v in values], 'Rounded raw charges')
        summaries.append(compare_policy(values, formal, policy['components'][index]))
        allvalues.extend(values)
        rows.extend({'atom_id': a['atom_id'], 'element': a['element'], 'ieee754_hex': q.hex(), 'charge_e': q} for a, q in zip(atoms, values, strict=True))
        coordinates.append({'atom_ids': ids, 'elements': [a['element'] for a in atoms], 'geometry_sha256': sha(directory / 'geometry-bohr.npy'), 'formal_charge': formal})
    rows.sort(key=lambda r: r['atom_id'])
    for field in ('rows',):
        imported = charges[field]
        require([(r['atom_id'], r['element'], bits(r['charge_e'])) for r in imported] == [(r['atom_id'], r['element'], bits(r['charge_e'])) for r in rows], 'Imported vector/order differs from raw bits')
    total = sum(c['formal_charge'] for c in data['components'])
    system = compare_policy(allvalues, total, policy['system'])
    require([c['atom_ids'] for c in charges['source_audit']['components']] == components and [c['formal_total_e'] for c in charges['source_audit']['components']] == [c['formal_charge'] for c in data['components']], 'Imported component totals changed')
    # Frozen reviewed policy and charge record forbid changing even descriptive
    # fields or re-signing an altered origin/tolerance/source into acceptance.
    pinned_charges = read(root / f'{kind}-charges.json')
    require(sha(root / f'{kind}-charges.json') == anchor['charges_sha256'] and charges == pinned_charges, 'Unreviewed policy/charge provenance')
    body = {'schema': 'island.j38-c3-independent-provenance.v1', 'fixture': kind,
            'final_graph_identity': graph_id, 'charge_identity': charges['identity'], 'runtime_identity': runtime['identity'],
            'anchor_identity': digest(anchors), 'auditor_sha256': sha(__file__), 'raw_files': anchor['files'],
            'charge_origin': 'RESP/Psi4', 'method': data['method'], 'basis': data['basis'],
            'components': summaries, 'system': system, 'mapping': coordinates,
            'initial_centering_crosscheck_tolerance_angstrom': 5e-14,
            'printed_geometry_witness_tolerance_bohr': 5.1e-12,
            'import_bit_exact': True, 'normalization': False, 'redistribution': False,
            'raw_provenance_verified': True, 'policy_matches': True,
            'limitation': 'Reviewed raw artifact provenance; hash integrity is not independent experimental validation.',
            'production_validated': False, 'simulation_readiness': 'not_established'}
    return {**body, 'identity': digest(body)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--fixture', choices=('pe', 'aromatic'), required=True)
    parser.add_argument('--charges', type=Path, required=True)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--anchors', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(args.root, args.fixture, read(args.charges), read(args.runtime), read(args.anchors))
        code = 0
    except (ValueError, KeyError, OSError, TypeError) as error:
        result = {'status': 'blocked', 'raw_provenance_verified': False, 'diagnostic': str(error)}
        code = 1
    with args.output.open('x') as stream:
        json.dump(result, stream, indent=2, allow_nan=False);stream.write('\n')
    return code


if __name__ == '__main__':
    raise SystemExit(main())
