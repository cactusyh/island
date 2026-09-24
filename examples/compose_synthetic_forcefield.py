"""Compose synthetic parameters, charges, and an explicit nonbonded policy."""

from island.builders import build_linear_polymer
from island.exceptions import IncompleteChargeAssignmentError
from island.forcefields import (
    AtomTypeChargeEngine,
    NonbondedPolicy,
    ParameterAssignmentEngine,
    ParameterizedSystem,
    ProvidedChargeEngine,
    RDKitSmartsAtomTypingEngine,
    island_demo_charges_v1,
    island_demo_parameters_v1,
    island_demo_v1_ruleset,
)

system = build_linear_polymer("[*]CC[*]", dp=3, generate_3d=False)
ruleset = island_demo_v1_ruleset()
typing = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
parameters = ParameterAssignmentEngine().assign(
    system, typing, ruleset, island_demo_parameters_v1()
)
charges = AtomTypeChargeEngine().assign(
    system, typing, ruleset, island_demo_charges_v1()
)
policy = NonbondedPolicy(
    name="island_synthetic_nonbonded_v1",
    version="1.0.0",
    mixing_rule="lorentz_berthelot",
    lj_scale_12=0.0,
    coulomb_scale_12=0.0,
    lj_scale_13=0.0,
    coulomb_scale_13=0.0,
    lj_scale_14=0.5,
    coulomb_scale_14=0.833333333333,
    source="SYNTHETIC SOFTWARE-TEST POLICY — NOT FOR SCIENTIFIC SIMULATION.",
)
snapshot = ParameterizedSystem.from_components(system, parameters, charges, policy)

print("parameter coverage:", parameters.complete_supported_scope)
print("charge coverage:", charges.complete)
print("charge method:", charges.method, charges.source)
print("nonbonded policy:", policy.name, policy.mixing_formula)
print("local exceptions:", len(policy.local_pair_scalings(system.topology)))
print("aggregate signature:", snapshot.aggregate_signature)
print("production validated:", snapshot.metadata["aggregate"]["production_validated"])
print("simulation readiness:", snapshot.metadata["aggregate"]["simulation_readiness"])

bad_values = {site_id: 0.0 for site_id in system.topology.sites}
bad_values[min(bad_values)] = 0.25
try:
    ProvidedChargeEngine().assign(
        system, bad_values, source="INTENTIONALLY INVALID EXAMPLE"
    )
except IncompleteChargeAssignmentError as error:
    print("intentional charge failure:", error.result.diagnostics[-1])

print("Coordinates remain angstrom; parameter lengths remain nm; no conversion ran.")
