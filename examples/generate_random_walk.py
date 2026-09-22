"""Generate non-destructive self-avoiding random-walk polymer coordinates."""

from island.builders import build_linear_polymer
from island.conformations import generate_polymer_conformation

chemical_system = build_linear_polymer("[*]CC[*]", dp=20, generate_3d=False)
result = generate_polymer_conformation(chemical_system, seed=2026)
conformed_system = result.apply_to(chemical_system)

print(result.attempts, result.rejected_trials, result.rollback_count)
print(result.minimum_nonbonded_distance)
print(conformed_system.number_of_sites)
