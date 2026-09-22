"""RDKit-free public data models for graph-based atom typing."""

from dataclasses import dataclass, field
from typing import Any, Literal

from island.exceptions import InvalidAtomTypingRuleError

HydrogenPolicy = Literal["explicit_required", "implicit_allowed"]
TypingStatus = Literal["assigned", "untyped", "ambiguous"]


@dataclass(frozen=True)
class AtomTypingRule:
    """One SMARTS rule with an explicitly marked target atom ``:1``."""

    rule_id: str
    atom_type: str
    smarts: str
    overrides: tuple[str, ...] = ()
    description: str = ""
    source: str = ""

    def __post_init__(self) -> None:
        for field_name in ("rule_id", "atom_type", "smarts"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise InvalidAtomTypingRuleError(
                    f"Atom typing rule {field_name} must be a non-empty string"
                )
        for field_name in ("description", "source"):
            if not isinstance(getattr(self, field_name), str):
                raise InvalidAtomTypingRuleError(
                    f"Atom typing rule {field_name} must be a string"
                )
        if isinstance(self.overrides, str) or any(
            not isinstance(rule_id, str) or not rule_id.strip()
            for rule_id in self.overrides
        ):
            raise InvalidAtomTypingRuleError(
                "Atom typing rule overrides must contain non-empty rule IDs"
            )
        normalized = tuple(sorted(set(self.overrides)))
        if self.rule_id in normalized:
            raise InvalidAtomTypingRuleError(
                f"Rule {self.rule_id!r} cannot override itself"
            )
        object.__setattr__(self, "overrides", normalized)


@dataclass(frozen=True)
class AtomTypingRuleSet:
    """Versioned rule collection and its supported chemical representation."""

    name: str
    version: str
    rules: tuple[AtomTypingRule, ...]
    supported_representation: str = "atomistic"
    hydrogen_policy: HydrogenPolicy = "explicit_required"
    description: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise InvalidAtomTypingRuleError("Ruleset name must be a non-empty string")
        if not isinstance(self.version, str) or not self.version.strip():
            raise InvalidAtomTypingRuleError(
                "Ruleset version must be a non-empty string"
            )
        rules = tuple(self.rules)
        if any(not isinstance(rule, AtomTypingRule) for rule in rules):
            raise InvalidAtomTypingRuleError(
                "Ruleset rules must be AtomTypingRule instances"
            )
        if not rules:
            raise InvalidAtomTypingRuleError("An atom-typing ruleset cannot be empty")
        object.__setattr__(self, "rules", rules)
        if self.supported_representation != "atomistic":
            raise InvalidAtomTypingRuleError(
                "Phase 4A rulesets support atomistic representation only"
            )
        if self.hydrogen_policy not in {
            "explicit_required",
            "implicit_allowed",
        }:
            raise InvalidAtomTypingRuleError(
                "hydrogen_policy must be 'explicit_required' or 'implicit_allowed'"
            )
        by_id = {rule.rule_id: rule for rule in rules}
        if len(by_id) != len(rules):
            duplicates = sorted(
                rule_id
                for rule_id in by_id
                if sum(rule.rule_id == rule_id for rule in rules) > 1
            )
            raise InvalidAtomTypingRuleError(
                f"Duplicate atom-typing rule IDs: {duplicates}"
            )
        unknown = {
            referenced
            for rule in rules
            for referenced in rule.overrides
            if referenced not in by_id
        }
        if unknown:
            raise InvalidAtomTypingRuleError(
                f"Rules reference unknown override IDs: {sorted(unknown)}"
            )
        _validate_acyclic_overrides(by_id)

    @property
    def rules_by_id(self) -> dict[str, AtomTypingRule]:
        return {rule.rule_id: rule for rule in self.rules}


@dataclass(frozen=True)
class AtomTypeAssignment:
    """Resolved external atom-type assignment for one stable site."""

    site_id: int
    atom_type: str
    selected_rule_ids: tuple[str, ...]
    matched_rule_ids: tuple[str, ...]


@dataclass(frozen=True)
class SiteTypingDiagnostic:
    """Deterministic rule-resolution trace for one stable site."""

    site_id: int
    status: TypingStatus
    matched_rule_ids: tuple[str, ...]
    eliminated_rule_ids: tuple[str, ...]
    surviving_rule_ids: tuple[str, ...]
    surviving_atom_types: tuple[str, ...]
    eliminated_by: dict[str, tuple[str, ...]] = field(default_factory=dict)


@dataclass(frozen=True)
class AtomTypingResult:
    """Stable-site assignments, diagnostics, and compatibility signatures."""

    ruleset_name: str
    ruleset_version: str
    assignments: dict[int, AtomTypeAssignment]
    diagnostics: dict[int, SiteTypingDiagnostic]
    complete: bool
    untyped_site_ids: tuple[int, ...]
    ambiguous_site_ids: tuple[int, ...]
    graph_signature: str
    ruleset_signature: str
    typing_signature: str
    engine_name: str
    engine_version: str
    dependency_versions: dict[str, str]
    metadata: dict[str, Any] = field(default_factory=dict)

    def is_compatible_with(
        self,
        topology_or_system: object,
        ruleset: AtomTypingRuleSet,
    ) -> bool:
        """Check reusable typing inputs without consulting coordinates."""
        from island.forcefields.typing.signatures import (
            graph_signature,
            ruleset_signature,
        )

        topology = getattr(topology_or_system, "topology", topology_or_system)
        representation = getattr(topology_or_system, "representation", None)
        if representation is not None and (
            representation != ruleset.supported_representation
        ):
            return False
        return self.graph_signature == graph_signature(topology) and (
            self.ruleset_signature == ruleset_signature(ruleset)
        )


def _validate_acyclic_overrides(by_id: dict[str, AtomTypingRule]) -> None:
    visiting: set[str] = set()
    visited: set[str] = set()

    def visit(rule_id: str, path: tuple[str, ...]) -> None:
        if rule_id in visiting:
            cycle_start = path.index(rule_id)
            cycle = (*path[cycle_start:], rule_id)
            raise InvalidAtomTypingRuleError(
                f"Override cycle detected: {' -> '.join(cycle)}"
            )
        if rule_id in visited:
            return
        visiting.add(rule_id)
        for overridden in by_id[rule_id].overrides:
            visit(overridden, (*path, rule_id))
        visiting.remove(rule_id)
        visited.add(rule_id)

    for rule_id in sorted(by_id):
        visit(rule_id, ())
