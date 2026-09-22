"""Build a finite linear polyethylene chain from a labeled PSMILES repeat unit."""

from island.builders import build_linear_polymer

system = build_linear_polymer("[*:1]CC[*:2]", dp=5, random_seed=2026)

print(system.number_of_sites)
print(system.metadata["polymer"])
