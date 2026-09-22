from island.forcefields import ForceFieldBackend, ParameterizedSystem


def test_forcefield_abstraction_is_importable() -> None:
    assert issubclass(ForceFieldBackend, object)
    assert ParameterizedSystem.__name__ == "ParameterizedSystem"
