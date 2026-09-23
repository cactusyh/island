import pytest

pytest.importorskip("rdkit")

from island.exceptions import InvalidAtomTypingRuleError
from island.forcefields.typing import (
    AtomTypingRule,
    AtomTypingRuleSet,
    RDKitSmartsAtomTypingEngine,
)


def rule(rule_id: str, atom_type: str = "type", **kwargs: object) -> AtomTypingRule:
    return AtomTypingRule(rule_id, atom_type, "[#6:1]", **kwargs)


def ruleset(*rules: AtomTypingRule) -> AtomTypingRuleSet:
    return AtomTypingRuleSet("test", "1", rules, hydrogen_policy="implicit_allowed")


@pytest.mark.parametrize("field", ["rule_id", "atom_type", "smarts"])
def test_rule_requires_nonempty_identity_fields(field: str) -> None:
    values = {"rule_id": "r", "atom_type": "t", "smarts": "[#6:1]"}
    values[field] = ""
    with pytest.raises(InvalidAtomTypingRuleError, match=field):
        AtomTypingRule(**values)


def test_duplicate_rule_ids_are_invalid() -> None:
    with pytest.raises(InvalidAtomTypingRuleError, match="Duplicate"):
        ruleset(rule("same"), rule("same"))


def test_unknown_and_self_overrides_are_invalid() -> None:
    with pytest.raises(InvalidAtomTypingRuleError, match="unknown override"):
        ruleset(rule("a", overrides=("missing",)))
    with pytest.raises(InvalidAtomTypingRuleError, match="cannot override itself"):
        rule("a", overrides=("a",))


def test_override_cycle_is_invalid() -> None:
    with pytest.raises(InvalidAtomTypingRuleError, match="Override cycle"):
        ruleset(
            rule("a", overrides=("b",)),
            rule("b", overrides=("c",)),
            rule("c", overrides=("a",)),
        )


@pytest.mark.parametrize(
    "smarts",
    ["not valid smarts", "[#6]", "[#6:1]-[#8:1]"],
)
def test_smarts_must_compile_with_exactly_one_target(smarts: str) -> None:
    candidate = AtomTypingRule("bad", "bad", smarts)
    with pytest.raises(InvalidAtomTypingRuleError):
        RDKitSmartsAtomTypingEngine().validate_ruleset(ruleset(candidate))


def test_ruleset_identity_and_hydrogen_policy_are_validated() -> None:
    with pytest.raises(InvalidAtomTypingRuleError, match="name"):
        AtomTypingRuleSet("", "1", (rule("a"),))
    with pytest.raises(InvalidAtomTypingRuleError, match="version"):
        AtomTypingRuleSet("test", "", (rule("a"),))
    with pytest.raises(InvalidAtomTypingRuleError, match="hydrogen_policy"):
        AtomTypingRuleSet(
            "test",
            "1",
            (rule("a"),),
            hydrogen_policy="sometimes",  # type: ignore[arg-type]
        )
