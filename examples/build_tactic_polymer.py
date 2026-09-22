"""Build a stereochemically defined finite linear polymer."""

from island.builders import build_linear_polymer

system = build_linear_polymer(
    "[*:1]N[C@H](F)C[*:2]",
    dp=6,
    tacticity="syndiotactic",
    random_seed=2026,
)

print(system.metadata["polymer"]["stereochemical_sequence"])
