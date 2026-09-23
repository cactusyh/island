"""Deterministic signatures for atom-typing input compatibility."""

import hashlib
import json

from island.core import AtomSite, Topology
from island.forcefields.typing.models import AtomTypingRuleSet


def graph_signature(topology: Topology) -> str:
    """Hash all graph attributes consulted by the Phase 4A SMARTS matcher."""
    topology.validate()
    sites = []
    for site_id in sorted(topology.sites):
        site = topology.sites[site_id]
        if not isinstance(site, AtomSite):
            kind = type(site).__name__
            sites.append({"id": site_id, "kind": kind})
            continue
        sites.append(
            {
                "id": site_id,
                "kind": "AtomSite",
                "atomic_number": site.atomic_number,
                "element": site.element if site.atomic_number is None else None,
                "formal_charge": site.formal_charge,
                "aromatic": bool(site.metadata.get("aromatic", False)),
                "no_implicit_hydrogens": bool(
                    site.metadata.get("no_implicit_hydrogens", False)
                ),
                "explicit_hydrogen_count": site.metadata.get(
                    "explicit_hydrogen_count", 0
                ),
            }
        )
    bonds = [
        {
            "sites": list(bond_key),
            "order": topology.bonds[bond_key].order,
            "aromatic": topology.bonds[bond_key].aromatic,
        }
        for bond_key in sorted(topology.bonds)
    ]
    return _digest({"sites": sites, "bonds": bonds})


def ruleset_signature(ruleset: AtomTypingRuleSet) -> str:
    """Hash ruleset identity and normalized rule content independent of order."""
    rules = [
        {
            "rule_id": rule.rule_id,
            "atom_type": rule.atom_type,
            "smarts": rule.smarts,
            "overrides": list(rule.overrides),
            "description": rule.description,
            "source": rule.source,
        }
        for rule in sorted(ruleset.rules, key=lambda item: item.rule_id)
    ]
    return _digest(
        {
            "name": ruleset.name,
            "version": ruleset.version,
            "supported_representation": ruleset.supported_representation,
            "hydrogen_policy": ruleset.hydrogen_policy,
            "description": ruleset.description,
            "rules": rules,
        }
    )


def typing_signature(
    graph_digest: str,
    ruleset_digest: str,
    *,
    engine_name: str,
    engine_version: str,
) -> str:
    return _digest(
        {
            "graph_signature": graph_digest,
            "ruleset_signature": ruleset_digest,
            "engine_name": engine_name,
            "engine_version": engine_version,
        }
    )


def _digest(value: object) -> str:
    serialized = json.dumps(value, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(serialized.encode()).hexdigest()
