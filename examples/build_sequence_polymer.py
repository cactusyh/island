"""Build a sequence-defined finite linear copolymer."""

from island.builders import build_polymer_from_sequence

system = build_polymer_from_sequence(
    repeat_units={
        "A": "[*:1]CC[*:2]",
        "B": "[*:1]CO[*:2]",
    },
    sequence=["A", "A", "B", "A", "B"],
    random_seed=2026,
)

print(system.metadata["polymer"]["sequence"])
print(system.metadata["polymer"]["composition_counts"])
