"""Illustrative SMARTS rules; not a production force field."""

from island.forcefields.typing.models import AtomTypingRule, AtomTypingRuleSet


def island_demo_v1_ruleset() -> AtomTypingRuleSet:
    """Return a small explicit-hydrogen demonstration ruleset.

    These ``demo_*`` labels exercise matching and override behavior only. They do
    not provide charges, parameters, or compatibility with an established force
    field.
    """
    rules = (
        AtomTypingRule(
            "demo_c_generic",
            "demo_carbon",
            "[#6:1]",
            description="Broad carbon rule used to demonstrate overrides.",
        ),
        AtomTypingRule(
            "demo_c_aliphatic",
            "demo_c_aliphatic",
            "[C;!$(C=O):1]",
            overrides=("demo_c_generic",),
            description="Non-aromatic, non-carbonyl carbon.",
        ),
        AtomTypingRule(
            "demo_c_aromatic",
            "demo_c_aromatic",
            "[c:1]",
            overrides=("demo_c_generic",),
            description="Aromatic carbon.",
        ),
        AtomTypingRule(
            "demo_c_carbonyl",
            "demo_c_carbonyl",
            "[C:1](=[O])",
            overrides=("demo_c_generic",),
            description="Carbonyl carbon.",
        ),
        AtomTypingRule(
            "demo_o_generic",
            "demo_oxygen",
            "[#8:1]",
            description="Broad oxygen rule used to demonstrate overrides.",
        ),
        AtomTypingRule(
            "demo_o_carbonyl",
            "demo_o_carbonyl",
            "[O:1]=[C]",
            overrides=("demo_o_generic",),
            description="Carbonyl oxygen.",
        ),
        AtomTypingRule(
            "demo_o_hydroxyl",
            "demo_o_hydroxyl",
            "[O;X2:1]-[H]",
            overrides=("demo_o_generic",),
            description="Oxygen bonded to an explicit hydrogen.",
        ),
        AtomTypingRule(
            "demo_o_ether",
            "demo_o_ether",
            "[O;X2:1]([#6])[#6]",
            overrides=("demo_o_generic",),
            description="Oxygen bonded to two carbons.",
        ),
        AtomTypingRule(
            "demo_n_generic",
            "demo_nitrogen",
            "[#7:1]",
            description="Broad nitrogen rule used to demonstrate charge rules.",
        ),
        AtomTypingRule(
            "demo_n_neutral",
            "demo_n_neutral",
            "[#7+0:1]",
            overrides=("demo_n_generic",),
            description="Neutral nitrogen.",
        ),
        AtomTypingRule(
            "demo_n_positive",
            "demo_n_positive",
            "[#7+1:1]",
            overrides=("demo_n_generic",),
            description="Positively charged nitrogen.",
        ),
        AtomTypingRule(
            "demo_h_generic",
            "demo_hydrogen",
            "[#1:1]",
            description="Broad explicit-hydrogen rule.",
        ),
        AtomTypingRule(
            "demo_h_on_carbon",
            "demo_h_on_carbon",
            "[#1:1]-[#6]",
            overrides=("demo_h_generic",),
            description="Explicit hydrogen bonded to carbon.",
        ),
        AtomTypingRule(
            "demo_h_on_hetero",
            "demo_h_on_hetero",
            "[#1:1]-[#7,#8]",
            overrides=("demo_h_generic",),
            description="Explicit hydrogen bonded to nitrogen or oxygen.",
        ),
    )
    return AtomTypingRuleSet(
        name="island_demo_v1",
        version="1.0.0",
        rules=rules,
        supported_representation="atomistic",
        hydrogen_policy="explicit_required",
        description=(
            "Illustrative Phase 4A rules for deterministic typing diagnostics; "
            "not suitable for production parameterization."
        ),
    )
