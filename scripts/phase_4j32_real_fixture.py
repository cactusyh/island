"""Fresh complete PEO chains and optional J20 terminal C/O network for J32."""

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from island.builders import build_linear_polymer
from island.chemistry.rdkit_graph import system_to_rdkit_graph
from island.core import Coordinates, MolecularSystem, Topology
from island.crosslinking import ReactiveSiteRule, generate_crosslink_network
from island.graph import final_graph
from island.packing import PeriodicPackingConfig, pack_multichain_periodic

PSMILES = "[*:1]CCO[*:2]"
CHAIN_COUNT = 3
DP = 4
NETWORK_DP = 3
DENSITY = 0.4
SEED = 2032
ROOT = Path(__file__).resolve().parents[2]
FRC = ROOT / "lammps/lammps/tools/msi2lmp/frc_files/pcff.frc"
XML = ROOT / "island-validation/foyer-4g1-upstream/foyer/forcefields/xml/oplsaa.xml"
AMBER = ROOT / "miniforge3/envs/md_rg"
FOYER = ROOT / "island-validation/phase4g1-env/bin/python"
RULE = ReactiveSiteRule(
    "J32 terminal C/O closure",
    "j32_reactive",
    True,
    ("C", "O"),
    1,
    1.0,
    "Explicit terminal H removal before J20",
    ("docs/phase_4j32.md",),
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def validate_complete_chain(system):
    """No implicit H, radicals, disconnected fragments or missing coordinates."""
    system.validate()
    require(
        len(system.topology.connected_components()) == 1, "Disconnected chain input"
    )
    converted = system_to_rdkit_graph(system)
    for atom in converted.mol.GetAtoms():
        require(
            atom.GetNumRadicalElectrons() == 0
            and atom.GetNumImplicitHs() == 0
            and atom.GetNumExplicitHs() == 0,
            "Incomplete explicit-H chemical graph",
        )
    return converted


@dataclass
class Fixture:
    kind: str
    chains: tuple
    units: list
    packed: object
    input_system: object = None
    input_graph: object = None
    crosslinked: object = None

    def evidence(self):
        return {
            "kind": self.kind,
            "psmiles": PSMILES,
            "dp": NETWORK_DP if self.kind == "network" else DP,
            "chain_count": CHAIN_COUNT,
            "complete_chains": [s.to_dict() for s in self.chains],
            "packing_units": [
                {"label": label, "system": s.to_dict(), "graph": g.payload}
                for label, s, g in self.units
            ],
            "packing_plan": self.packed.plan.payload,
            "packed_system": self.packed.system.to_dict(),
            "final_graph": self.packed.graph.payload,
            "crosslink_plan": self.crosslinked.plan.payload
            if self.crosslinked
            else None,
            "crosslink_input": self.input_system.to_dict()
            if self.input_system
            else None,
            "crosslink_input_graph": self.input_graph.payload
            if self.input_graph
            else None,
            "setup": setup_metadata(self.packed),
        }


def setup_metadata(packed):
    packed.graph.validate_integrity(packed.system)
    p = packed.graph.payload
    members = p["molecule_membership"]
    ids = sorted(packed.system.topology.sites)
    lengths = np.array(packed.system.box.lengths)
    require(
        packed.system.box.periodic == (True, True, True), "Missing fully periodic box"
    )
    minimum = None
    for n, a in enumerate(ids):
        for b in ids[n + 1 :]:
            if members[str(a)] != members[str(b)]:
                delta = packed.system.coordinates.get(
                    a
                ) - packed.system.coordinates.get(b)
                delta -= np.rint(delta / lengths) * lengths
                d = float(np.linalg.norm(delta))
                minimum = d if minimum is None else min(minimum, d)
    mass = sum(s.mass for s in packed.system.topology.sites.values())
    return {
        "atom_count": len(ids),
        "chain_count": len(set(members.values())),
        "component_count": len(p["components"]),
        "crosslink_count": len(p["crosslink_provenance"]),
        "box_lengths_angstrom": lengths.tolist(),
        "volume_angstrom3": float(np.prod(lengths)),
        "mass_dalton": mass,
        "density_g_cm3": mass / float(np.prod(lengths)) * 1.66053906660,
        "minimum_interchain_distance_angstrom": minimum,
        "minimum_interunit_distance_angstrom": packed.config.minimum_interunit_distance,
        "density_meaning": "derived setup quantity, not an equilibrium or experimental property",
    }


def build_fixture(kind="linear", *, link_seed=0):
    require(kind in ("linear", "network"), "Unknown fixture kind")
    dp = NETWORK_DP if kind == "network" else DP
    chains = tuple(
        build_linear_polymer(
            PSMILES, dp=dp, random_seed=SEED + i, chain_id=f"J32-chain-{i}"
        )
        for i in range(CHAIN_COUNT)
    )
    for s in chains:
        validate_complete_chain(s)
    units = [(f"J32-unit-{i}", s, final_graph(s)) for i, s in enumerate(chains)]
    input_system = input_graph = linked = None
    if kind == "network":
        topology, xyz = Topology(), {}
        for i, complete in enumerate(chains):
            s = complete.copy()
            for atom_id in (1, dp * 3):
                hid = next(
                    h
                    for h in s.topology.neighbors(atom_id)
                    if s.topology.sites[h].element == "H"
                )
                s.topology.remove_site(hid)
                s.coordinates._positions.pop(hid)
                s.topology.sites[atom_id].metadata["j32_reactive"] = True
            mapping = {a: i * 1000 + a for a in s.topology.sites}
            for a in s.topology.sites.values():
                copied = deepcopy(a)
                copied.id = mapping[a.id]
                topology.add_site(copied)
                xyz[copied.id] = s.coordinates.get(a.id) + np.array([i * 20, 0, 0])
            for b in s.topology.bonds.values():
                topology.add_bond(mapping[b.site1], mapping[b.site2], order=b.order)
        input_system = MolecularSystem(topology, Coordinates(xyz))
        input_graph = final_graph(input_system)
        linked = generate_crosslink_network(
            input_system,
            input_graph,
            RULE,
            target_crosslinks=CHAIN_COUNT,
            seed=link_seed,
            provenance="J32 terminal closure, final graph re-embedded before J21",
            evidence=("docs/phase_4j32.md",),
        )
        converted = validate_complete_chain(linked.system)
        require(
            all(
                {linked.system.topology.sites[a].element for a in row["sites"]}
                == {"C", "O"}
                for row in linked.graph.payload["crosslink_provenance"]
            ),
            "Closure must contain only C/O bonds",
        )
        from rdkit.Chem import AllChem

        params = AllChem.ETKDGv3()
        params.randomSeed = SEED
        require(
            AllChem.EmbedMolecule(converted.mol, params) == 0,
            "Final graph embedding failed",
        )
        conf = converted.mol.GetConformer()
        linked.system.coordinates = Coordinates(
            {
                a: list(conf.GetAtomPosition(j))
                for a, j in converted.site_id_to_rdkit_index.items()
            }
        )
        units = [("J32-network", linked.system, linked.graph)]
    packed = pack_multichain_periodic(
        units,
        config=PeriodicPackingConfig(
            target_density=DENSITY,
            seed=650 if kind == "network" else SEED,
            minimum_interunit_distance=1.8 if kind == "linear" else 0,
        ),
        provenance="J32 complete rigid units packed with J21",
        evidence=tuple(g.identity for _, _, g in units),
    )
    return Fixture(kind, chains, units, packed, input_system, input_graph, linked)
