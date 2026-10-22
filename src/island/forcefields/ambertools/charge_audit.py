"""Read-only native AM1-BCC charge-path diagnostics; never repairs charges.

These records explain an unsuccessful preparation. They are not a charge
assignment and cannot authorize dynamics or replace the strict importer.
"""

import json
import math
import re
from decimal import Decimal
from pathlib import Path

from island.exceptions import IncompleteChargeAssignmentError
from island.forcefields.amber import import_amber_prmtop
from island.forcefields.ambertools.lineage import parse_and_validate_mol2, prepare_input
from island.forcefields.ambertools.runtime import checksum, digest, signed

SCHEMA = 'island.am1bcc-charge-audit.v1'


def require(value, message):
    if not value:
        raise ValueError(message)


def _ac(path, names):
    lines = [line for line in path.read_text().splitlines() if line.startswith('ATOM ')]
    rows = [line.split() for line in lines]
    require(len(rows) == len(names), 'Incomplete native AC charge table')
    # AC is fixed width: a four-character name touches the residue name.
    require([line[13:17].strip() for line in lines] == list(names), 'Native AC atom-name/order mismatch')
    return [row[-2] for row in rows]


def _sqm(path, elements):
    text = path.read_text()
    tables = re.findall(r'Atom\s+Element\s+Mulliken Charge\s*\n(.*?)Total Mulliken Charge\s*=\s*([-+\d.]+)', text, re.DOTALL)
    require(tables, 'Missing native SQM Mulliken charge table')
    rows = [line.split() for line in tables[-1][0].splitlines() if line.strip()]
    require(len(rows) == len(elements), 'Incomplete SQM atom table')
    require([int(row[0]) for row in rows] == list(range(1, len(elements) + 1)), 'Wrong SQM atom order')
    require([row[1] for row in rows] == elements, 'Wrong SQM element association')
    require(all(len(row) == 3 for row in rows), 'Unsupported SQM charge table')
    return [row[2] for row in rows], tables[-1][1]


def audit_component(system, directory, *, runtime_identity, tolerance=1e-4):
    """Reconstruct the exact input, native lineage and strict import failure.

    Decimal totals below retain the native text precision. The actual decision
    is made independently by the existing importer using its unchanged tolerance.
    Missing or mismatched artifacts raise, rather than suggesting SQM precision
    loss without a complete trace.
    """
    import parmed

    directory = Path(directory)
    prepared = prepare_input(system, None, max_atoms=100)
    require((directory / 'input.mol2').read_text() == prepared.mol2_text, 'Native input differs from exact final component')
    names = prepared.names
    ids = list(names.values())
    elements = [system.topology.sites[i].element for i in ids]
    formal = [system.topology.sites[i].formal_charge for i in ids]
    sqm_input = (directory / 'sqm.in').read_text()
    qmcharge = re.search(r'qmcharge\s*=\s*([-+\d.]+)', sqm_input, re.IGNORECASE)
    require(qmcharge is not None and Decimal(qmcharge[1]) == sum(formal), 'SQM formal charge was lost')
    # SQM input carries the same stable names. Native antechamber first writes
    # coordinates through fixed-width AC (three decimal places), then SQM input.
    # Compare the actual native texts exactly, without changing coordinates or
    # introducing a geometry tolerance.
    from parmed.periodic_table import AtomicNum

    sqm_rows = [line.split() for line in sqm_input.splitlines() if re.match(r'^\s*\d+\s+\S+\s+[-+\d.]', line)]
    ac_rows = [line for line in (directory / 'ANTECHAMBER_AC.AC').read_text().splitlines() if line.startswith('ATOM ')]
    require(len(sqm_rows) == len(ac_rows) == len(names), 'Incomplete SQM input lineage')
    require([row[1] for row in sqm_rows] == [line[13:17].strip() for line in ac_rows] == list(names), 'SQM input atom-name/order mismatch')
    require([int(row[0]) for row in sqm_rows] == [AtomicNum[e] for e in elements], 'SQM input element mismatch')
    require(all([Decimal(x) for x in row[2:5]] == [Decimal(x) for x in line.split()[4:7]] for row, line in zip(sqm_rows, ac_rows, strict=True)), 'SQM input/native AC coordinate mismatch')
    sqm, native_total = _sqm(directory / 'sqm.out', elements)
    pre = _ac(directory / 'ANTECHAMBER_AM1BCC_PRE.AC', names)
    bcc = _ac(directory / 'ANTECHAMBER_AM1BCC.AC', names)
    mapping, parsed, _ = parse_and_validate_mol2(directory / 'typed.mol2', system, names, prepared.expected_cip, stage='charge audit')
    require(set(mapping.values()) == set(ids), 'Incomplete stable-ID mapping')
    section = (directory / 'typed.mol2').read_text().split('@<TRIPOS>ATOM\n')[1].split('@<TRIPOS>BOND')[0]
    mol2 = {line.split()[1]: line.split()[8] for line in section.splitlines() if line.strip()}
    require(set(mol2) == set(names), 'Wrong MOL2 charge names')
    parm = parmed.load_file(str(directory / 'result.prmtop'))
    require(len(parm.atoms) == len(ids) and {a.name for a in parm.atoms} == set(names), 'Wrong prmtop atom names')
    prmtop_map = {a.idx: names[a.name] for a in parm.atoms}
    prmtop = {a.name: a.charge for a in parm.atoms}
    # Exact graph/element/name validation is also performed by the original
    # importer below. No partial charge result is ever published by this audit.
    failure = None
    try:
        import_amber_prmtop(system, directory / 'result.prmtop', prmtop_map,
                            source='Read-only native AM1-BCC charge-path audit',
                            force_field='gaff2', charge_method='AM1-BCC',
                            charge_tolerance=tolerance)
    except IncompleteChargeAssignmentError as error:
        failure = {'exception': type(error).__name__, 'message': str(error),
                   'stage': 'import/component_total_validation'}
    rows = [{'atom_id': i, 'name': name, 'element': element, 'formal_charge': charge,
             'sqm_mulliken_e': s, 'pre_bcc_e': p, 'bcc_e': b,
             'mol2_e': mol2[name], 'prmtop_e': prmtop[name]}
            for (name, i), element, charge, s, p, b in zip(names.items(), elements, formal, sqm, pre, bcc, strict=True)]
    def total(values):
        return str(sum((Decimal(v) for v in values), Decimal(0)))
    totals = {'expected_e': sum(formal), 'sqm_reported_total_e': native_total,
              'sqm_printed_atom_sum_e': total(sqm), 'pre_bcc_e': total(pre),
              'bcc_e': total(bcc), 'mol2_e': total(mol2.values()),
              'prmtop_e': math.fsum(prmtop.values())}
    native_precision_loss = (
        failure is not None and Decimal(native_total) == sum(formal)
        and abs(Decimal(totals['sqm_printed_atom_sum_e']) - sum(formal)) > Decimal(str(tolerance))
        # Native antechamber averages equivalent atoms before PRE.AC; that
        # established family convention can change individual values while
        # preserving the already incorrect printed SQM total.
        and Decimal(totals['sqm_printed_atom_sum_e']) == Decimal(totals['pre_bcc_e'])
        and Decimal(totals['pre_bcc_e']) == Decimal(totals['bcc_e']) == Decimal(totals['mol2_e'])
        and all(Decimal(b) == Decimal(mol2[name]) for name, b in zip(names, bcc, strict=True))
        and max(abs(parsed[names[name]] - prmtop[name]) for name in names) < 1e-7
    )
    return signed({'schema': SCHEMA, 'runtime_identity': runtime_identity,
                   'input_sha256': checksum(directory / 'input.mol2'),
                   'classification': 'native_sqm_printed_charge_precision_loss' if native_precision_loss else 'charge_path_observed',
                   'failure': failure, 'tolerance_e': tolerance,
                   'sqm_input_lineage_validated': True,
                   'native_pre_bcc_changed_atom_count': sum(Decimal(s) != Decimal(p) for s, p in zip(sqm, pre, strict=True)),
                   'totals': totals, 'atoms': rows,
                   'published_success_marker': (directory / 'ISLAND_SUCCESS').exists(),
                   'files': {p.name: checksum(p) for p in sorted(directory.iterdir()) if p.is_file()},
                   'production_validated': False, 'simulation_readiness': 'not_established'})


def validate_audit(record, system, directory, *, runtime_identity):
    require(record.get('identity') == digest({k: v for k, v in record.items() if k != 'identity'}), 'Charge audit identity changed')
    require(record == audit_component(system, directory, runtime_identity=runtime_identity), 'Native charge audit no longer reconstructs')
    return record


def write_audit(record, path):
    Path(path).write_text(json.dumps(record, indent=2, sort_keys=True, allow_nan=False) + '\n')
