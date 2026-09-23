"""RDKit SMARTS atom typing with deterministic override resolution."""

from importlib.metadata import PackageNotFoundError, version

import rdkit
from rdkit import Chem

from island.chemistry.rdkit_graph import topology_to_rdkit_graph
from island.core import Topology
from island.exceptions import (
    IncompleteAtomTypingError,
    InvalidAtomTypingRuleError,
    RDKitConversionError,
    UnsupportedAtomTypingError,
)
from island.forcefields.typing.base import AtomTypingEngine
from island.forcefields.typing.models import (
    AtomTypeAssignment,
    AtomTypingResult,
    AtomTypingRule,
    AtomTypingRuleSet,
    SiteTypingDiagnostic,
)
from island.forcefields.typing.signatures import (
    graph_signature,
    ruleset_signature,
    typing_signature,
)


class RDKitSmartsAtomTypingEngine(AtomTypingEngine):
    """Apply ordinary RDKit SMARTS rules to an atomistic chemical graph."""

    engine_name = "rdkit_smarts"
    engine_version = "1"

    def validate_ruleset(
        self, ruleset: AtomTypingRuleSet
    ) -> dict[str, tuple[Chem.Mol, int]]:
        """Compile patterns and validate exactly one target atom mapped ``:1``."""
        compiled: dict[str, tuple[Chem.Mol, int]] = {}
        for rule in ruleset.rules:
            query = Chem.MolFromSmarts(rule.smarts)
            if query is None:
                raise InvalidAtomTypingRuleError(
                    f"Rule {rule.rule_id!r} contains invalid SMARTS: {rule.smarts!r}"
                )
            targets = [
                atom.GetIdx() for atom in query.GetAtoms() if atom.GetAtomMapNum() == 1
            ]
            if len(targets) != 1:
                raise InvalidAtomTypingRuleError(
                    f"Rule {rule.rule_id!r} must contain exactly one target atom "
                    "marked :1"
                )
            compiled[rule.rule_id] = (query, targets[0])
        return compiled

    def type_topology(
        self,
        topology: Topology,
        ruleset: AtomTypingRuleSet,
        *,
        strict: bool = True,
    ) -> AtomTypingResult:
        """Type every authoritative site with exhaustive SMARTS matching."""
        compiled = self.validate_ruleset(ruleset)
        try:
            converted = topology_to_rdkit_graph(
                topology, representation=ruleset.supported_representation
            )
        except RDKitConversionError as error:
            raise UnsupportedAtomTypingError(
                f"Cannot type the supplied topology: {error}"
            ) from error
        self._validate_hydrogen_policy(
            converted.mol, converted.rdkit_index_to_site_id, ruleset
        )
        matches_by_site = {site_id: set() for site_id in topology.sites}
        for rule in sorted(ruleset.rules, key=lambda item: item.rule_id):
            query, target_query_index = compiled[rule.rule_id]
            matches = converted.mol.GetSubstructMatches(
                query,
                uniquify=False,
                useChirality=False,
                maxMatches=2_147_483_647,
            )
            target_indices = {match[target_query_index] for match in matches}
            for rdkit_index in target_indices:
                matches_by_site[converted.rdkit_index_to_site_id[rdkit_index]].add(
                    rule.rule_id
                )

        rules_by_id = ruleset.rules_by_id
        override_closure = _override_closure(rules_by_id)
        assignments: dict[int, AtomTypeAssignment] = {}
        diagnostics: dict[int, SiteTypingDiagnostic] = {}
        untyped: list[int] = []
        ambiguous: list[int] = []
        for site_id in sorted(topology.sites):
            matched = set(matches_by_site[site_id])
            eliminated_by = {
                candidate: tuple(
                    sorted(
                        winner
                        for winner in matched
                        if candidate in override_closure[winner]
                    )
                )
                for candidate in matched
            }
            eliminated_by = {
                candidate: winners
                for candidate, winners in eliminated_by.items()
                if winners
            }
            eliminated = set(eliminated_by)
            surviving = matched - eliminated
            surviving_types = {rules_by_id[item].atom_type for item in surviving}
            if not matched:
                status = "untyped"
                untyped.append(site_id)
            elif len(surviving_types) == 1:
                status = "assigned"
                atom_type = next(iter(surviving_types))
                assignments[site_id] = AtomTypeAssignment(
                    site_id=site_id,
                    atom_type=atom_type,
                    selected_rule_ids=tuple(sorted(surviving)),
                    matched_rule_ids=tuple(sorted(matched)),
                )
            else:
                status = "ambiguous"
                ambiguous.append(site_id)
            diagnostics[site_id] = SiteTypingDiagnostic(
                site_id=site_id,
                status=status,
                matched_rule_ids=tuple(sorted(matched)),
                eliminated_rule_ids=tuple(sorted(eliminated)),
                surviving_rule_ids=tuple(sorted(surviving)),
                surviving_atom_types=tuple(sorted(surviving_types)),
                eliminated_by=dict(sorted(eliminated_by.items())),
            )

        graph_digest = graph_signature(topology)
        rules_digest = ruleset_signature(ruleset)
        result = AtomTypingResult(
            ruleset_name=ruleset.name,
            ruleset_version=ruleset.version,
            assignments=assignments,
            diagnostics=diagnostics,
            complete=not untyped and not ambiguous,
            untyped_site_ids=tuple(untyped),
            ambiguous_site_ids=tuple(ambiguous),
            graph_signature=graph_digest,
            ruleset_signature=rules_digest,
            typing_signature=typing_signature(
                graph_digest,
                rules_digest,
                engine_name=self.engine_name,
                engine_version=self.engine_version,
            ),
            engine_name=self.engine_name,
            engine_version=self.engine_version,
            dependency_versions={
                "island": _package_version("island"),
                "rdkit": rdkit.__version__,
            },
            metadata={
                "matching": "ordinary_rdkit_smarts",
                "resolution": "transitive_explicit_overrides",
                "substructure_matching": "exhaustive_no_default_limit",
                "use_chirality": False,
                "parameterized": False,
            },
        )
        if strict and not result.complete:
            raise IncompleteAtomTypingError(
                "Atom typing is incomplete: "
                f"untyped={list(result.untyped_site_ids)}, "
                f"ambiguous={list(result.ambiguous_site_ids)}",
                result=result,
            )
        return result

    @staticmethod
    def _validate_hydrogen_policy(
        mol: Chem.Mol,
        rdkit_index_to_site_id: dict[int, int],
        ruleset: AtomTypingRuleSet,
    ) -> None:
        if ruleset.hydrogen_policy == "implicit_allowed":
            return
        omitted = []
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                continue
            count = atom.GetNumImplicitHs() + atom.GetNumExplicitHs()
            if count:
                omitted.append((rdkit_index_to_site_id[atom.GetIdx()], count))
        if omitted:
            raise UnsupportedAtomTypingError(
                f"Ruleset {ruleset.name!r} requires explicit hydrogen sites; "
                f"omitted counts by site: {omitted[:12]}"
            )


def _override_closure(
    rules_by_id: dict[str, AtomTypingRule],
) -> dict[str, frozenset[str]]:
    memo: dict[str, frozenset[str]] = {}

    def closure(rule_id: str) -> frozenset[str]:
        if rule_id not in memo:
            direct = rules_by_id[rule_id].overrides
            memo[rule_id] = frozenset(
                set(direct).union(*(closure(item) for item in direct))
            )
        return memo[rule_id]

    return {rule_id: closure(rule_id) for rule_id in sorted(rules_by_id)}


def _package_version(package: str) -> str:
    try:
        return version(package)
    except PackageNotFoundError:
        return "unknown"
