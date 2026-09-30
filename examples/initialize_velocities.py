"""Explicit thermal velocities with only core ISLAND and NumPy dependencies."""

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.dynamics import initialize_velocities

graph = Topology()
graph.add_site(AtomSite(17, "C", 12.01, element="C", atomic_number=6))
graph.add_site(AtomSite(91, "H", 1.008, element="H", atomic_number=1))
graph.add_bond(17, 91)
system = MolecularSystem(graph, Coordinates({17: (0.0, 0.0, 0.0), 91: (1.1, 0.0, 0.0)}))
sample = initialize_velocities(system, temperature_kelvin=300, seed=78123)
print("velocities (angstrom/ps):", dict(sample.velocities))
print(
    "target/instantaneous temperature (K):",
    sample.temperature_kelvin,
    sample.instantaneous_temperature_kelvin,
)
print("kinetic energy (kJ/mol), DOF:", sample.kinetic_energy, sample.degrees_of_freedom)
print(
    "velocity seed/RNG/NumPy:",
    sample.velocity_seed,
    sample.rng_algorithm,
    sample.numpy_version,
)
print("initialization identity:", sample.initialization_fingerprint)
print(
    "No COM/rotation removal or rescaling; this sample is not an equilibrated system."
)
