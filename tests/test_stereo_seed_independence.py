from island.builders import PolymerSequence, build_polymer_from_sequence


def test_stereo_seed_is_independent_of_sequence_seed_metadata() -> None:
    repeat_units = {"A": "[*:1]N[C@H](F)C[*:2]"}
    first = build_polymer_from_sequence(
        repeat_units,
        PolymerSequence(("A",) * 10, {"sequence_seed": 1}),
        tacticity="atactic",
        stereo_seed=77,
        random_seed=11,
        generate_3d=False,
    )
    second = build_polymer_from_sequence(
        repeat_units,
        PolymerSequence(("A",) * 10, {"sequence_seed": 999}),
        tacticity="atactic",
        stereo_seed=77,
        random_seed=22,
        generate_3d=False,
    )
    assert (
        first.metadata["polymer"]["stereochemical_sequence"]
        == second.metadata["polymer"]["stereochemical_sequence"]
    )
    assert (
        first.metadata["polymer"]["sequence_generation"]
        != second.metadata["polymer"]["sequence_generation"]
    )
