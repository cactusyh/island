from island.builders import build_linear_polymer


def test_polyethylene_terminal_backbone_atoms_are_hydrogen_capped() -> None:
    system = build_linear_polymer("[*]CC[*]", dp=3, random_seed=2026)
    metadata = system.metadata["polymer"]

    for terminal_key in ("head_site_id", "tail_site_id"):
        terminal_id = metadata[terminal_key]
        neighbors = [
            system.topology.get_site(site_id)
            for site_id in system.topology.neighbors(terminal_id)
        ]
        assert sum(site.element == "H" for site in neighbors) == 3
        assert sum(site.element == "C" for site in neighbors) == 1
