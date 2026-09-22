"""Build a DP=50 polymer from short local 3D templates."""

from island.builders import build_linear_polymer

system = build_linear_polymer(
    "[*:1]CC[*:2]",
    dp=50,
    coordinate_method="local_templates",
    template_seed=2026,
    assembly_seed=2026,
)

print(system.number_of_sites)
print(system.metadata["polymer"]["coordinates"])
print(system.metadata["polymer"]["coordinate_generation"])
