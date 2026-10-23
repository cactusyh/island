"""Frozen PSMILES benchmark inventory; end hydrogens complete finite oligomers."""

from dataclasses import asdict, dataclass

from island.builders import build_linear_polymer
from island.graph import final_graph
from island.packing import PeriodicPackingConfig, pack_multichain_periodic
from scripts.phase_4j32_real_fixture import Fixture, validate_complete_chain
from scripts.phase_4j34_real_fixture import build_fixture as retained


@dataclass(frozen=True)
class Benchmark:
    key: str
    name: str
    psmiles: str
    dp: int
    chains: int
    seed: int
    construction: str


BENCHMARKS = (
    Benchmark("pe", "polyethylene oligomer", "[*:1]CC[*:2]", 4, 2, 2038, "fresh"),
    Benchmark(
        "peo",
        "poly(ethylene oxide) oligomer",
        "[*:1]CCO[*:2]",
        4,
        2,
        2032,
        "retained_linear_recipe",
    ),
    Benchmark(
        "aromatic",
        "poly(para-phenylene) oligomer",
        "[*:1]c1ccc([*:2])cc1",
        2,
        2,
        2038,
        "fresh",
    ),
    Benchmark(
        "network",
        "controlled terminal-linked PEO network",
        "[*:1]CCO[*:2]",
        3,
        3,
        2032,
        "retained_network_recipe",
    ),
)
SPECS = {s.key: s for s in BENCHMARKS}


def build(key):
    spec = SPECS[key]
    if key in {"peo", "network"}:
        fixture = retained("linear" if key == "peo" else "network", "OPLS-AA").fixture
    else:
        chains = tuple(
            build_linear_polymer(
                spec.psmiles,
                dp=spec.dp,
                random_seed=spec.seed + i,
                chain_id=f"J38-{key}-{i}",
            )
            for i in range(spec.chains)
        )
        for chain in chains:
            validate_complete_chain(chain)
        units = [
            (f"J38-{key}-{i}", chain, final_graph(chain))
            for i, chain in enumerate(chains)
        ]
        packed = pack_multichain_periodic(
            units,
            config=PeriodicPackingConfig(
                target_density=0.3, seed=spec.seed, minimum_interunit_distance=1.8
            ),
            provenance="J38 complete explicit-H finite PSMILES oligomer packing",
            evidence=tuple(g.identity for _, _, g in units),
        )
        fixture = Fixture(key, chains, units, packed)
    evidence = fixture.evidence()
    evidence.update(
        benchmark=asdict(spec),
        kind=key,
        psmiles=spec.psmiles,
        dp=spec.dp,
        chain_count=spec.chains,
        namespace="J38-" + key,
        end_group_policy="hydrogen valence completion with RDKit AddHs; network terminal H removed before documented J20 C/O closure",
        stereochemistry="no controlled stereocenters in selected repeat units",
        replica_count=1,
        provenance_scope="fresh construction from frozen recipe; retained recipes deliberately preserve J37 graph/charge bindings",
    )
    evidence["construction_seeds"] = (
        [spec.seed, spec.seed]
        if key == "peo"
        else [spec.seed + i for i in range(spec.chains)]
    )
    evidence["endpoints"] = []
    if key != "network":
        for chain, mapping in zip(
            fixture.chains,
            fixture.packed.plan.payload["site_id_remapping"],
            strict=True,
        ):
            p = chain.metadata["polymer"]
            evidence["endpoints"].append(
                [mapping[str(p["head_site_id"])], mapping[str(p["tail_site_id"])]]
            )
    return fixture, evidence
