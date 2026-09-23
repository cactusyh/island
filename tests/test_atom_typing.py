from copy import deepcopy

import pytest

pytest.importorskip("rdkit")
from rdkit.Chem import AllChem

from island import AtomSite, Topology
from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.exceptions import (
    IncompleteAtomTypingError,
    UnsupportedAtomTypingError,
)
from island.forcefields.typing import (
    AtomTypingRule,
    AtomTypingRuleSet,
    RDKitSmartsAtomTypingEngine,
    island_demo_v1_ruleset,
)


def assigned_types(result: object) -> dict[int, str]:
    return {
        site_id: assignment.atom_type
        for site_id, assignment in result.assignments.items()
    }


def one_rule_set(
    *rules: AtomTypingRule, hydrogen_policy: str = "implicit_allowed"
) -> AtomTypingRuleSet:
    return AtomTypingRuleSet(
        "test_rules",
        "1",
        rules,
        hydrogen_policy=hydrogen_policy,  # type: ignore[arg-type]
    )


@pytest.mark.parametrize(
    ("smiles", "expected"),
    [
        (
            "CCO",
            {
                "demo_c_aliphatic",
                "demo_o_hydroxyl",
                "demo_h_on_carbon",
                "demo_h_on_hetero",
            },
        ),
        (
            "CC(=O)C",
            {
                "demo_c_aliphatic",
                "demo_c_carbonyl",
                "demo_o_carbonyl",
                "demo_h_on_carbon",
            },
        ),
        ("c1ccccc1", {"demo_c_aromatic", "demo_h_on_carbon"}),
        (
            "CN",
            {
                "demo_c_aliphatic",
                "demo_n_neutral",
                "demo_h_on_carbon",
                "demo_h_on_hetero",
            },
        ),
        (
            "C[NH3+]",
            {
                "demo_c_aliphatic",
                "demo_n_positive",
                "demo_h_on_carbon",
                "demo_h_on_hetero",
            },
        ),
    ],
)
def test_demonstration_types_representative_molecules(
    smiles: str, expected: set[str]
) -> None:
    system = from_smiles(smiles, random_seed=7)
    result = RDKitSmartsAtomTypingEngine().type_system(system, island_demo_v1_ruleset())
    assert result.complete
    assert set(assigned_types(result).values()) == expected
    assert result.metadata["parameterized"] is False


def test_explicit_override_resolution_records_diagnostics() -> None:
    system = from_smiles("CC(=O)C", random_seed=7)
    result = RDKitSmartsAtomTypingEngine().type_system(system, island_demo_v1_ruleset())
    carbonyl_id = next(
        site_id
        for site_id, assignment in result.assignments.items()
        if assignment.atom_type == "demo_c_carbonyl"
    )
    diagnostic = result.diagnostics[carbonyl_id]
    assert diagnostic.matched_rule_ids == ("demo_c_carbonyl", "demo_c_generic")
    assert diagnostic.eliminated_rule_ids == ("demo_c_generic",)
    assert diagnostic.surviving_rule_ids == ("demo_c_carbonyl",)
    assert diagnostic.eliminated_by == {"demo_c_generic": ("demo_c_carbonyl",)}


def test_transitive_override_resolution() -> None:
    topology = Topology()
    topology.add_site(AtomSite(10, "C", 12.011, element="C", atomic_number=6))
    rules = one_rule_set(
        AtomTypingRule("broad", "broad", "[#6:1]"),
        AtomTypingRule("middle", "middle", "[#6:1]", overrides=("broad",)),
        AtomTypingRule("specific", "specific", "[#6:1]", overrides=("middle",)),
    )
    result = RDKitSmartsAtomTypingEngine().type_topology(topology, rules)
    assert result.assignments[10].atom_type == "specific"
    assert result.diagnostics[10].eliminated_rule_ids == ("broad", "middle")


def test_untyped_and_ambiguous_diagnostics_never_fabricate_assignments() -> None:
    topology = Topology()
    topology.add_site(AtomSite(10, "C", 12.011, element="C", atomic_number=6))
    topology.add_site(AtomSite(20, "O", 15.999, element="O", atomic_number=8))
    topology.add_bond(10, 20, order=2)
    rules = one_rule_set(
        AtomTypingRule("carbon_a", "a", "[#6:1]"),
        AtomTypingRule("carbon_b", "b", "[#6:1]"),
    )
    engine = RDKitSmartsAtomTypingEngine()
    result = engine.type_topology(topology, rules, strict=False)
    assert not result.complete
    assert result.ambiguous_site_ids == (10,)
    assert result.untyped_site_ids == (20,)
    assert result.assignments == {}
    assert result.diagnostics[10].surviving_atom_types == ("a", "b")
    assert result.diagnostics[20].matched_rule_ids == ()
    with pytest.raises(IncompleteAtomTypingError) as captured:
        engine.type_topology(topology, rules)
    assert captured.value.result.diagnostics == result.diagnostics


def test_surviving_rules_that_agree_on_type_are_not_ambiguous() -> None:
    topology = Topology()
    topology.add_site(AtomSite(10, "C", 12.011, element="C", atomic_number=6))
    rules = one_rule_set(
        AtomTypingRule("carbon_a", "same", "[#6:1]"),
        AtomTypingRule("carbon_b", "same", "[#6:1]"),
    )
    result = RDKitSmartsAtomTypingEngine().type_topology(topology, rules)
    assert result.assignments[10].selected_rule_ids == ("carbon_a", "carbon_b")


def test_noncontiguous_stable_ids_are_assignment_keys() -> None:
    ids = [10, 20, 35, 47, 59, 61, 73, 89, 101]
    system = from_smiles("CCO", random_seed=7, site_ids=ids)
    result = RDKitSmartsAtomTypingEngine().type_system(system, island_demo_v1_ruleset())
    assert set(result.assignments) == set(ids)
    assert set(result.diagnostics) == set(ids)


def test_graph_only_typing_does_not_need_coordinates_or_embed(
    monkeypatch: object,
) -> None:
    topology = Topology()
    topology.add_site(AtomSite(10, "C", 12.011, element="C", atomic_number=6))
    rules = one_rule_set(AtomTypingRule("carbon", "carbon", "[#6:1]"))

    def fail_embed(*args: object, **kwargs: object) -> None:
        raise AssertionError("typing must never embed coordinates")

    monkeypatch.setattr(AllChem, "EmbedMolecule", fail_embed)
    result = RDKitSmartsAtomTypingEngine().type_topology(topology, rules)
    assert result.assignments[10].atom_type == "carbon"


def test_typing_does_not_mutate_input_on_success_or_strict_failure() -> None:
    system = from_smiles("CCO", random_seed=7)
    before = system.to_dict()
    engine = RDKitSmartsAtomTypingEngine()
    engine.type_system(system, island_demo_v1_ruleset())
    assert system.to_dict() == before
    incomplete = one_rule_set(AtomTypingRule("carbon", "carbon", "[#6:1]"))
    with pytest.raises(IncompleteAtomTypingError):
        engine.type_system(system, incomplete)
    assert system.to_dict() == before
    assert all(
        "atom_type" not in site.metadata for site in system.topology.sites.values()
    )


def test_explicit_hydrogen_policy_rejects_omitted_hydrogens() -> None:
    system = from_smiles("CCO", add_hydrogens=False, random_seed=7)
    with pytest.raises(UnsupportedAtomTypingError, match="explicit hydrogen"):
        RDKitSmartsAtomTypingEngine().type_system(system, island_demo_v1_ruleset())


def test_implicit_hydrogen_policy_allows_implicit_hydrogens() -> None:
    system = from_smiles("CC", add_hydrogens=False, random_seed=7)
    rules = one_rule_set(AtomTypingRule("carbon", "carbon", "[#6:1]"))
    result = RDKitSmartsAtomTypingEngine().type_system(system, rules)
    assert set(result.assignments) == set(system.topology.sites)


def test_explicit_policy_accepts_chemically_hydrogen_free_structure() -> None:
    system = from_smiles("O=C=O", add_hydrogens=False, random_seed=7)
    result = RDKitSmartsAtomTypingEngine().type_system(system, island_demo_v1_ruleset())
    assert result.complete
    assert sorted(assigned_types(result).values()) == [
        "demo_c_carbonyl",
        "demo_o_carbonyl",
        "demo_o_carbonyl",
    ]


def reordered_topology(source: Topology) -> Topology:
    reordered = Topology()
    for site_id in reversed(list(source.sites)):
        reordered.add_site(deepcopy(source.sites[site_id]))
    for key in reversed(list(source.bonds)):
        bond = source.bonds[key]
        reordered.add_bond(
            bond.site2, bond.site1, order=bond.order, aromatic=bond.aromatic
        )
    reordered.rebuild_derived_interactions()
    return reordered


def test_rule_and_topology_insertion_order_do_not_affect_results() -> None:
    system = from_smiles("CC(=O)C", random_seed=7)
    original_rules = island_demo_v1_ruleset()
    reversed_rules = AtomTypingRuleSet(
        original_rules.name,
        original_rules.version,
        tuple(reversed(original_rules.rules)),
        hydrogen_policy=original_rules.hydrogen_policy,
        description=original_rules.description,
    )
    engine = RDKitSmartsAtomTypingEngine()
    first = engine.type_topology(system.topology, original_rules)
    second = engine.type_topology(reordered_topology(system.topology), reversed_rules)
    assert assigned_types(first) == assigned_types(second)
    assert first.diagnostics == second.diagnostics
    assert first.graph_signature == second.graph_signature
    assert first.ruleset_signature == second.ruleset_signature
    assert first.typing_signature == second.typing_signature


def test_matching_is_exhaustive_beyond_rdkit_default_limit() -> None:
    topology = Topology()
    for site_id in range(2000, 3105):
        topology.add_site(
            AtomSite(
                site_id,
                f"C{site_id}",
                12.011,
                element="C",
                atomic_number=6,
                metadata={"no_implicit_hydrogens": True},
            )
        )
    rules = one_rule_set(AtomTypingRule("carbon", "carbon", "[#6:1]"))
    result = RDKitSmartsAtomTypingEngine().type_topology(topology, rules)
    assert result.complete
    assert len(result.assignments) == 1105


@pytest.mark.parametrize("dp", [50, 100])
def test_demonstration_rules_type_long_polyethylene(dp: int) -> None:
    system = build_linear_polymer("[*]CC[*]", dp=dp, generate_3d=False)
    result = RDKitSmartsAtomTypingEngine().type_system(system, island_demo_v1_ruleset())
    types = assigned_types(result)
    assert result.complete
    assert sum(value == "demo_c_aliphatic" for value in types.values()) == 2 * dp
    assert sum(value == "demo_h_on_carbon" for value in types.values()) == 4 * dp + 2


def test_graph_only_adapter_has_explicit_stable_id_mappings() -> None:
    from island.chemistry import topology_to_rdkit_graph

    topology = Topology()
    topology.add_site(AtomSite(35, "O", 15.999, element="O", atomic_number=8))
    topology.add_site(AtomSite(10, "C", 12.011, element="C", atomic_number=6))
    topology.add_bond(35, 10, order=2)
    converted = topology_to_rdkit_graph(topology)
    assert converted.site_id_to_rdkit_index == {10: 0, 35: 1}
    assert converted.rdkit_index_to_site_id == {0: 10, 1: 35}
    assert converted.mol.GetNumConformers() == 0
    assert converted.mol.GetBondWithIdx(0).GetBondTypeAsDouble() == 2.0


def test_coarse_grained_topology_is_rejected_by_typing_engine() -> None:
    from island import BeadSite

    topology = Topology()
    topology.add_site(BeadSite(10, "B", 28.0, bead_type="demo"))
    rules = one_rule_set(AtomTypingRule("carbon", "carbon", "[#6:1]"))
    with pytest.raises(UnsupportedAtomTypingError, match="AtomSite"):
        RDKitSmartsAtomTypingEngine().type_topology(topology, rules)
