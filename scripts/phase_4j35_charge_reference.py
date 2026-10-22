"""Independent external electrostatic reference; never reads frozen Amber inputs.

This script uses upstream Foyer SMARTS and source XML decimals, with a separate
local-neighbourhood oracle for the bounded PEO fixtures. No ISLAND typing,
parameter assignment or provided-charge implementation is imported.
"""

import argparse
import hashlib
import json
from decimal import Decimal
from importlib.metadata import version
from pathlib import Path
from xml.etree import ElementTree

from scripts.phase_4j32_real_fixture import XML
from scripts.phase_4j34_real_fixture import build_fixture

SOURCE_SHA = 'c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4'
METHOD = 'external_oplsaa_source_decimals_foyer_and_neighbourhood_v1'


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def neighbourhood_type(system, i):
    """Independent bounded chemical rule oracle, not a general atom typer."""
    sites = system.topology.sites
    neighbours = system.topology.neighbors(i)
    element = sites[i].element
    elements = sorted(sites[n].element for n in neighbours)
    if element == 'O':
        return {('C', 'C'): 'opls_180', ('C', 'H'): 'opls_154'}[tuple(elements)]
    if element == 'C':
        oxygens = [n for n in neighbours if sites[n].element == 'O']
        if not oxygens:
            return {3: 'opls_135', 2: 'opls_136'}[elements.count('H')]
        assert len(oxygens) == 1
        alcohol = any(sites[n].element == 'H' for n in system.topology.neighbors(oxygens[0]))
        if alcohol:
            assert elements.count('H') == 2
            return 'opls_157'
        return {3: 'opls_181', 2: 'opls_182'}[elements.count('H')]
    assert element == 'H' and len(neighbours) == 1
    parent = next(iter(neighbours))
    if sites[parent].element == 'O':
        return 'opls_155'
    return 'opls_185' if any(sites[n].element == 'O' and neighbourhood_type(system, n) == 'opls_180' for n in system.topology.neighbors(parent)) else 'opls_140'


def reconstruct(kind):
    import foyer
    from foyer.atomtyper import AtomTypingRulesProvider, find_atomtypes
    from foyer.topology_graph import TopologyGraph

    pin = json.loads(Path('src/island/forcefields/oplsaa/pin.json').read_text())
    for relative, expected in pin['implementation_files'].items():
        assert hashlib.sha256((Path(foyer.__file__).parent / relative).read_bytes()).hexdigest() == expected
    assert hashlib.sha256(XML.read_bytes()).hexdigest() == SOURCE_SHA
    fixture = build_fixture(kind, 'GAFF2').fixture
    system, graph = fixture.packed.system, fixture.packed.graph
    tree = ElementTree.parse(XML).getroot()
    types = {e.attrib['name']: e.attrib for e in tree.find('AtomTypes')}
    charges = {e.attrib['type']: e.attrib['charge'] for e in tree.find('NonbondedForce')}
    rules = AtomTypingRulesProvider({k: v['def'] for k, v in types.items() if v.get('def')}, {k: {s.strip() for s in v.get('overrides', '').split(',') if s.strip()} for k, v in types.items()}, set())
    upstream = TopologyGraph()
    for i, site in system.topology.sites.items():
        upstream.add_atom(i, site.element, atomic_number=site.atomic_number, symbol=site.element)
    for bond in system.topology.bonds.values():
        upstream.add_bond(*bond.key, bond_order=bond.order)
    matches = find_atomtypes(upstream, rules, max_iter=len(types) * len(system.topology.sites) + 1)
    rows = []
    for i in sorted(system.topology.sites):
        atomtype = matches[i]['atomtype']
        assert atomtype == neighbourhood_type(system, i), (i, atomtype, neighbourhood_type(system, i))
        rows.append({'atom_id': i, 'element': system.topology.sites[i].element, 'atom_type': atomtype, 'decimal_e': charges[atomtype], 'charge_e': float(Decimal(charges[atomtype]))})
    components = []
    for ids in sorted(sorted(c) for c in system.topology.connected_components()):
        total = sum((Decimal(r['decimal_e']) for r in rows if r['atom_id'] in ids), Decimal(0))
        formal = sum(system.topology.sites[i].formal_charge for i in ids)
        assert total == formal, (ids, total, formal)
        components.append({'atom_ids': ids, 'formal_total_e': formal, 'provided_total_e': str(total)})
    value = {'schema': 'island.j35-independent-charge-reference.v1', 'fixture': kind, 'final_graph_identity': graph.identity, 'rows': rows, 'components': components, 'units': 'elementary_charge', 'method': METHOD, 'source_sha256': SOURCE_SHA, 'source_revision': pin['revision'], 'implementation_hashes': pin['implementation_files'], 'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'versions': {n: version(n) for n in ('foyer', 'parmed', 'numpy')}, 'provenance': 'External OPLS-AA XML source charges with GAFF2 bonded/LJ parameters; not AM1-BCC; no physical hybrid-model validation', 'independent_checks': ['upstream_foyer_smarts', 'bounded_local_chemical_neighbourhood', 'exact_decimal_component_totals'], 'normalization': False}
    return {**value, 'identity': digest(value)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture', choices=('linear', 'network'), required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    value = reconstruct(args.fixture)
    with args.output.open('x') as stream:
        stream.write(json.dumps(value, sort_keys=True, indent=2) + '\n')
    print(value['identity'])


if __name__ == '__main__':
    main()
