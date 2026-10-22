"""Fresh complete J34 linear fixture; retain the reviewed J33 connected network."""

from dataclasses import dataclass

from island.analysis.npt import NPTConfiguration
from island.builders import build_linear_polymer
from island.graph import final_graph
from island.packing import PeriodicPackingConfig, pack_multichain_periodic
from scripts.phase_4j32_real_fixture import Fixture, validate_complete_chain
from scripts.phase_4j33_real_fixture import J33Fixture, build_j33_fixture


@dataclass
class LinearFixture(Fixture):
    def evidence(self):
        value = super().evidence()
        value.update(dp=4, chain_count=2, psmiles="[*:1]CCO[*:2]", namespace='J34-linear-peo-dp4-two-chain')
        return value


def build_fixture(kind, family):
    if kind == 'network':
        return build_j33_fixture(kind, family)
    if kind != 'linear':
        raise ValueError('Unsupported J34 fixture')
    chains = tuple(build_linear_polymer("[*:1]CCO[*:2]", dp=4, random_seed=2032, chain_id=f'J34-chain-{i}') for i in range(2))
    for chain in chains:
        validate_complete_chain(chain)
    units = [(f'J34-unit-{i}', chain, final_graph(chain)) for i, chain in enumerate(chains)]
    packed = pack_multichain_periodic(units, config=PeriodicPackingConfig(target_density=0.3, seed=2032, minimum_interunit_distance=1.8), provenance='J34 complete explicit-H PEO DP4 linear chains; setup density only', evidence=tuple(g.identity for _, _, g in units))
    fixture = LinearFixture('linear', chains, units, packed)
    config = NPTConfiguration(family=family, source_identity='0' * 64, final_graph_identity=packed.graph.identity, packing_identity=packed.plan.identity, box_lengths_nm=tuple(v / 10 for v in packed.system.box.lengths))
    return J33Fixture('linear', fixture, config)
