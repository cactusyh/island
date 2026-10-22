"""Independent stdlib NPY/JSON parser; no production charge/policy imports."""

import argparse
import ast
import hashlib
import json
import math
import struct
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def parse_npy(path):
    raw = Path(path).read_bytes()
    if raw[:8] != b'\x93NUMPY\x01\x00':
        raise ValueError('Expected retained NPY v1 binary64 representation')
    size = struct.unpack('<H', raw[8:10])[0]
    header = ast.literal_eval(raw[10:10 + size].decode('latin1').strip())
    if header['descr'] != '<f8' or header['fortran_order'] or len(header['shape']) != 1:
        raise ValueError('Expected ordered little-endian native float64 vector')
    count = header['shape'][0]
    if len(raw) != 10 + size + 8 * count:
        raise ValueError('Truncated/extended native array')
    return list(struct.unpack('<' + 'd' * count, raw[10 + size:])), hashlib.sha256(raw).hexdigest()


def check(directory, audit):
    directory = Path(directory)
    data = json.loads((directory / 'input.json').read_text())
    memory = audit['schema'] == 'island.resp-native-memory-precision.v1'
    allrows, components, hashes, formal_total = [], [], {}, 0
    for index, component in enumerate(data['components']):
        path = directory / f'component-{index}'
        filename = 'charges-native.npy' if memory else 'native-resp-charges.npy'
        values, sha = parse_npy(path / filename)
        hashes[f'component-{index}/{filename}'] = sha
        if sha != audit['components'][index]['native_array_sha256']:
            raise ValueError('Wrong authoritative native array')
        saved = json.loads((path / ('charges-native.json' if memory else 'charge-vector.json')).read_text())
        ids = [a['atom_id'] for a in component['atoms']]
        if [r['atom_id'] for r in saved] != ids or [r['element'] for r in saved] != [a['element'] for a in component['atoms']]:
            raise ValueError('Wrong ID/element ordering')
        if [struct.pack('<d', v) for v in values] != [struct.pack('<d', r['charge_e']) for r in saved] or not all(math.isfinite(v) for v in values):
            raise ValueError('Raw JSON/binary value mismatch')
        formal = component['formal_charge'];formal_total += formal
        residual = sum((Fraction.from_float(v) for v in values), Fraction()) - formal
        with localcontext() as ctx:
            ctx.prec = 2200
            decimal_total = sum((Decimal.from_float(v) for v in values), Decimal())
        components.append({'residual_fraction': [residual.numerator, residual.denominator], 'decimal_binary_total_e': str(decimal_total), 'math_fsum_e': math.fsum(values), 'atom_ids': ids})
        allrows.extend({'atom_id': a['atom_id'], 'element': a['element'], 'charge_e': v, 'decimal_e': str(Decimal.from_float(v)), 'ieee754_hex': v.hex()} for a, v in zip(component['atoms'], values, strict=True))
    allrows.sort(key=lambda r: r['atom_id'])
    if allrows != audit['rows']:
        raise ValueError('Independent source rows differ from native audit')
    residual = sum((Fraction.from_float(r['charge_e']) for r in allrows), Fraction()) - formal_total
    body = {'schema': 'island.resp-finite-independent.v1', 'source_audit_identity': audit['identity'], 'mapping_verified': True,
            'rows': allrows, 'raw_array_hashes': hashes, 'components': components,
            'component_residual_fractions': [c['residual_fraction'] for c in components],
            'system_residual_fraction': [residual.numerator, residual.denominator],
            'decimal_text_authority': False, 'precision_note': 'Decimal.from_float reconstructs binary values; printed decimal tables are not authority.',
            'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    return {**body, 'identity': digest(body)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    with args.output.open('x') as stream:
        json.dump(check(args.directory, json.loads(args.audit.read_text())), stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
