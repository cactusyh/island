"""Construct a two-atom molecular system without external chemistry toolkits."""

from island import AtomSite, Coordinates, MolecularSystem, Topology

topology = Topology()
topology.add_site(AtomSite(id=10, name="C1", mass=12.011, element="C"))
topology.add_site(AtomSite(id=20, name="C2", mass=12.011, element="C"))
topology.add_bond(10, 20, order=1)

system = MolecularSystem(
    topology=topology,
    coordinates=Coordinates({10: [0.0, 0.0, 0.0], 20: [1.54, 0.0, 0.0]}),
)
system.validate()
print(system.to_dict())
