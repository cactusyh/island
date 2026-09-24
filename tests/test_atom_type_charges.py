import os
import subprocess
import sys
from dataclasses import replace

import pytest
from test_parameter_assignment_graph_only import graph_typing

from island.forcefields import (
    AtomTypeChargeEngine,
    AtomTypeChargeEntry,
    AtomTypeChargeTable,
)
from island.forcefields.charges.signatures import charge_table_signature
from island.forcefields.typing.signatures import ruleset_signature

SOURCE = "SYNTHETIC SOFTWARE-TEST CHARGES"


def charge_table(ruleset: object) -> AtomTypeChargeTable:
    common = {"source": SOURCE, "table_name": "test_charges", "table_version": "1"}
    return AtomTypeChargeTable(
        "test_charges",
        "1",
        ruleset.name,
        ruleset.version,
        ruleset_signature(ruleset),
        (
            AtomTypeChargeEntry("charge_A", "A", 0.2, **common),
            AtomTypeChargeEntry("charge_B", "B", -0.2, **common),
            AtomTypeChargeEntry("charge_C", "C", 0.3, **common),
            AtomTypeChargeEntry("charge_D", "D", -0.3, **common),
        ),
        SOURCE,
    )


def test_atom_type_charge_lookup_records_entries_and_values() -> None:
    topology, typing, ruleset = graph_typing()
    table = charge_table(ruleset)
    result = AtomTypeChargeEngine().assign(topology, typing, ruleset, table)
    assert result.complete
    assert {site_id: item.charge for site_id, item in result.assignments.items()} == {
        10: 0.2,
        30: -0.2,
        70: 0.3,
        120: -0.3,
    }
    assert result.assignments[10].entry_id == "charge_A"
    assert result.assignments[10].atom_type == "A"
    assert result.table_signature == charge_table_signature(table)
    assert result.total_assigned_charge == pytest.approx(0.0)


def test_explicit_zero_entry_is_assignment_while_absence_is_missing() -> None:
    topology, typing, ruleset = graph_typing()
    table = charge_table(ruleset)
    zero_entries = tuple(replace(entry, charge=0.0) for entry in table.entries)
    zero_table = replace(table, entries=zero_entries)
    zero_result = AtomTypeChargeEngine().assign(topology, typing, ruleset, zero_table)
    assert zero_result.complete
    assert zero_result.assignments[10].charge == 0.0

    missing_table = replace(
        zero_table,
        entries=tuple(entry for entry in zero_entries if entry.atom_type != "A"),
    )
    partial = AtomTypeChargeEngine().assign(
        topology, typing, ruleset, missing_table, strict=False
    )
    assert not partial.complete
    assert 10 not in partial.assignments
    assert any(
        item.reason == "missing" and item.site_ids == (10,)
        for item in partial.diagnostics
    )


def test_conflicting_type_entries_are_ambiguous_not_order_resolved() -> None:
    topology, typing, ruleset = graph_typing()
    table = charge_table(ruleset)
    conflict = replace(table.entries[0], entry_id="charge_A_conflict", charge=0.25)
    conflicting = replace(table, entries=(*table.entries, conflict))
    result = AtomTypeChargeEngine().assign(
        topology, typing, ruleset, conflicting, strict=False
    )
    diagnostic = next(item for item in result.diagnostics if item.site_ids == (10,))
    assert diagnostic.reason == "ambiguous"
    assert diagnostic.candidate_entry_ids == ("charge_A", "charge_A_conflict")
    assert 10 not in result.assignments


def test_charge_table_changes_invalidate_compatibility_and_signatures() -> None:
    topology, typing, ruleset = graph_typing()
    table = charge_table(ruleset)
    result = AtomTypeChargeEngine().assign(topology, typing, ruleset, table)
    changed_entry = replace(table.entries[0], charge=0.25)
    changed_table = replace(table, entries=(changed_entry, *table.entries[1:]))
    assert not result.is_input_compatible_with(
        topology, typing_result=typing, table=changed_table
    )
    changed_result = AtomTypeChargeEngine().assign(
        topology, typing, ruleset, changed_table, strict=False
    )
    assert changed_result.table_signature != result.table_signature
    assert changed_result.result_signature != result.result_signature

    reordered = replace(table, entries=tuple(reversed(table.entries)))
    reordered_result = AtomTypeChargeEngine().assign(
        topology, typing, ruleset, reordered
    )
    assert reordered_result.table_signature == result.table_signature
    assert reordered_result.result_signature == result.result_signature


def test_graph_only_charge_lookup_operates_without_rdkit() -> None:
    code = r"""
import importlib.abc
import sys

class BlockRDKit(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "rdkit" or fullname.startswith("rdkit."):
            raise RuntimeError("charge assignment imported RDKit")
        return None

sys.meta_path.insert(0, BlockRDKit())
from island import AtomSite, Topology
from island.forcefields import AtomTypeChargeEngine, AtomTypeChargeEntry, AtomTypeChargeTable
from island.forcefields.typing import (
    AtomTypeAssignment, AtomTypingResult, AtomTypingRule, AtomTypingRuleSet,
    SiteTypingDiagnostic,
)
from island.forcefields.typing.signatures import graph_signature, ruleset_signature, typing_signature

topology = Topology()
topology.add_site(AtomSite(10, "C", 12.0, element="C", atomic_number=6))
ruleset = AtomTypingRuleSet("r", "1", (AtomTypingRule("rule", "A", "[#6:1]"),), hydrogen_policy="implicit_allowed")
graph = graph_signature(topology)
rules = ruleset_signature(ruleset)
typing = AtomTypingResult(
    "r", "1", {10: AtomTypeAssignment(10, "A", ("rule",), ("rule",))},
    {10: SiteTypingDiagnostic(10, "assigned", ("rule",), (), ("rule",), ("A",))},
    True, (), (), graph, rules,
    typing_signature(graph, rules, engine_name="precomputed", engine_version="1"),
    "precomputed", "1", {},
)
entry = AtomTypeChargeEntry("a", "A", 0.0, "synthetic", "t", "1")
table = AtomTypeChargeTable("t", "1", "r", "1", rules, (entry,), "synthetic")
result = AtomTypeChargeEngine().assign(topology, typing, ruleset, table)
assert result.complete and result.assignments[10].charge == 0.0
assert "rdkit" not in sys.modules
"""
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        check=False,
        env=environment,
        text=True,
    )
    assert completed.returncode == 0, completed.stderr
