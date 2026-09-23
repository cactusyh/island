import os
import subprocess
import sys
from dataclasses import replace

from island import AtomSite, Topology
from island.forcefields.parameters import (
    HarmonicAngleParameter,
    HarmonicBondParameter,
    LennardJonesParameter,
    ParameterAssignmentEngine,
    ParameterLibrary,
    PeriodicTorsionTerm,
    ProperTorsionParameter,
    parameter_library_signature,
)
from island.forcefields.parameters.signatures import typing_assignment_content_signature
from island.forcefields.typing import (
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

SOURCE = "synthetic graph-only test parameters"
COMMON = {"source": SOURCE, "library_name": "graph_test", "library_version": "1"}


def graph_typing(
    *, reversed_insertion: bool = False
) -> tuple[Topology, AtomTypingResult, AtomTypingRuleSet]:
    logical = ((10, "A"), (30, "B"), (70, "C"), (120, "D"))
    topology = Topology()
    for site_id, atom_type in reversed(logical) if reversed_insertion else logical:
        topology.add_site(
            AtomSite(site_id, atom_type, 12.0, element="C", atomic_number=6)
        )
    edges = ((10, 30), (30, 70), (70, 120))
    for site1, site2 in reversed(edges) if reversed_insertion else edges:
        topology.add_bond(site2, site1, order=1)
    ruleset = AtomTypingRuleSet(
        "graph_rules",
        "1",
        tuple(
            AtomTypingRule(f"rule_{atom_type}", atom_type, "[#6:1]")
            for atom_type in ("A", "B", "C", "D")
        ),
        hydrogen_policy="implicit_allowed",
    )
    assignments = {}
    diagnostics = {}
    for site_id, atom_type in logical:
        rule_id = f"rule_{atom_type}"
        assignments[site_id] = AtomTypeAssignment(
            site_id, atom_type, (rule_id,), (rule_id,)
        )
        diagnostics[site_id] = SiteTypingDiagnostic(
            site_id,
            "assigned",
            (rule_id,),
            (),
            (rule_id,),
            (atom_type,),
        )
    graph_digest = graph_signature(topology)
    rule_digest = ruleset_signature(ruleset)
    result = AtomTypingResult(
        ruleset_name=ruleset.name,
        ruleset_version=ruleset.version,
        assignments=assignments,
        diagnostics=diagnostics,
        complete=True,
        untyped_site_ids=(),
        ambiguous_site_ids=(),
        graph_signature=graph_digest,
        ruleset_signature=rule_digest,
        typing_signature=typing_signature(
            graph_digest,
            rule_digest,
            engine_name="precomputed_test",
            engine_version="1",
        ),
        engine_name="precomputed_test",
        engine_version="1",
        dependency_versions={},
    )
    return topology, result, ruleset


def graph_library(ruleset: AtomTypingRuleSet) -> ParameterLibrary:
    records = tuple(
        LennardJonesParameter(f"lj_{kind}", kind, 0.2, 0.3, **COMMON)
        for kind in ("A", "B", "C", "D")
    ) + (
        HarmonicBondParameter("bond_ba", ("B", "A"), 10.0, 0.15, **COMMON),
        HarmonicBondParameter("bond_cb", ("C", "B"), 10.0, 0.15, **COMMON),
        HarmonicBondParameter("bond_dc", ("D", "C"), 10.0, 0.15, **COMMON),
        HarmonicAngleParameter("angle_cba", ("C", "B", "A"), 20.0, 110.0, **COMMON),
        HarmonicAngleParameter("angle_dcb", ("D", "C", "B"), 20.0, 110.0, **COMMON),
        ProperTorsionParameter(
            "torsion_dcba",
            ("D", "C", "B", "A"),
            (
                PeriodicTorsionTerm(1.0, 1, 0.0),
                PeriodicTorsionTerm(0.5, 2, 180.0),
            ),
            **COMMON,
        ),
    )
    return ParameterLibrary(
        "graph_test",
        "1",
        ruleset.name,
        ruleset.version,
        ruleset_signature(ruleset),
        records,
        description="synthetic exact reversal tests",
        source=SOURCE,
    )


def test_exact_patterns_have_bond_angle_and_torsion_reversal_symmetry() -> None:
    topology, typing, ruleset = graph_typing()
    result = ParameterAssignmentEngine().assign(
        topology, typing, ruleset, graph_library(ruleset)
    )
    assert result.bond_assignments[(10, 30)].parameter_id == "bond_ba"
    assert result.angle_assignments[(10, 30, 70)].parameter_id == "angle_cba"
    assert (
        result.proper_torsion_assignments[(10, 30, 70, 120)].parameter_id
        == "torsion_dcba"
    )
    assert (
        len(result.proper_torsion_assignments[(10, 30, 70, 120)].parameter.terms) == 2
    )


def test_noncontiguous_ids_and_insertion_order_do_not_change_assignment() -> None:
    first_topology, first_typing, ruleset = graph_typing()
    second_topology, second_typing, _ = graph_typing(reversed_insertion=True)
    library = graph_library(ruleset)
    engine = ParameterAssignmentEngine()
    first = engine.assign(first_topology, first_typing, ruleset, library)
    second = engine.assign(second_topology, second_typing, ruleset, library)
    assert first.assignment_signature == second.assignment_signature
    assert first.site_assignments == second.site_assignments
    assert first.bond_assignments == second.bond_assignments
    assert first.angle_assignments == second.angle_assignments
    assert first.proper_torsion_assignments == second.proper_torsion_assignments


def test_assignment_ignores_empty_and_invalid_stale_derived_caches() -> None:
    from island import Angle, Dihedral

    topology, typing, ruleset = graph_typing()
    topology.angles.clear()
    topology.dihedrals.clear()
    topology.angles[(10, 30, 999)] = Angle(10, 30, 999)
    topology.dihedrals[(10, 30, 999, 120)] = Dihedral(10, 30, 999, 120)
    result = ParameterAssignmentEngine().assign(
        topology, typing, ruleset, graph_library(ruleset)
    )
    assert result.complete_supported_scope
    assert set(result.angle_assignments) == {(10, 30, 70), (30, 70, 120)}
    assert set(result.proper_torsion_assignments) == {(10, 30, 70, 120)}


def test_reuse_check_includes_library_and_mutable_typing_content() -> None:
    topology, typing, ruleset = graph_typing()
    library = graph_library(ruleset)
    result = ParameterAssignmentEngine().assign(topology, typing, ruleset, library)
    assert result.is_compatible_with(topology, typing, library)

    changed_record = replace(library.records[0], epsilon=0.21)
    changed_library = replace(library, records=(changed_record, *library.records[1:]))
    assert parameter_library_signature(changed_library) != result.library_signature
    assert not result.is_compatible_with(topology, typing, changed_library)

    typing.assignments[10] = replace(typing.assignments[10], atom_type="changed")
    assert typing_assignment_content_signature(typing) != (
        result.typing_assignment_signature
    )
    assert not result.is_compatible_with(topology, typing, library)


def test_library_record_order_does_not_change_assignment() -> None:
    topology, typing, ruleset = graph_typing()
    library = graph_library(ruleset)
    reversed_library = replace(library, records=tuple(reversed(library.records)))
    engine = ParameterAssignmentEngine()
    first = engine.assign(topology, typing, ruleset, library)
    second = engine.assign(topology, typing, ruleset, reversed_library)
    assert first.library_signature == second.library_signature
    assert first.assignment_signature == second.assignment_signature
    assert first.site_assignments == second.site_assignments


def test_explicit_zero_records_are_assigned_and_signed_not_treated_as_missing() -> None:
    topology, typing, ruleset = graph_typing()
    library = graph_library(ruleset)
    zero_records = []
    for record in library.records:
        if record.parameter_id == "lj_A":
            record = replace(record, epsilon=0.0)
        elif record.parameter_id == "torsion_dcba":
            record = replace(
                record,
                terms=(replace(record.terms[0], force_constant=0.0),),
            )
        zero_records.append(record)
    zero_library = replace(library, records=tuple(zero_records))
    zero_result = ParameterAssignmentEngine().assign(
        topology, typing, ruleset, zero_library
    )
    assert zero_result.complete_supported_scope
    assert zero_result.site_assignments[10].parameter.epsilon == 0.0
    torsion = zero_result.proper_torsion_assignments[(10, 30, 70, 120)].parameter
    assert torsion.terms[0].force_constant == 0.0

    nonzero_result = ParameterAssignmentEngine().assign(
        topology, typing, ruleset, library
    )
    assert zero_result.library_signature != nonzero_result.library_signature
    assert zero_result.assignment_signature != nonzero_result.assignment_signature

    missing_library = replace(
        zero_library,
        records=tuple(
            record for record in zero_library.records if record.parameter_id != "lj_A"
        ),
    )
    missing_result = ParameterAssignmentEngine().assign(
        topology, typing, ruleset, missing_library, strict=False
    )
    assert not missing_result.complete_supported_scope
    assert any(
        diagnostic.family == "site"
        and diagnostic.site_ids == (10,)
        and diagnostic.reason == "missing"
        for diagnostic in missing_result.diagnostics
    )


def test_precomputed_typing_assignment_operates_without_rdkit() -> None:
    code = r"""
import importlib.abc
import sys

class BlockRDKit(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "rdkit" or fullname.startswith("rdkit."):
            raise RuntimeError("parameter assignment imported RDKit")
        return None

sys.meta_path.insert(0, BlockRDKit())
from island import AtomSite, Topology
from island.forcefields.parameters import (
    LennardJonesParameter, ParameterAssignmentEngine, ParameterLibrary,
)
from island.forcefields.typing import (
    AtomTypeAssignment, AtomTypingResult, AtomTypingRule, AtomTypingRuleSet,
    SiteTypingDiagnostic,
)
from island.forcefields.typing.signatures import (
    graph_signature, ruleset_signature, typing_signature,
)

topology = Topology()
topology.add_site(AtomSite(10, "C", 12.0, element="C", atomic_number=6))
ruleset = AtomTypingRuleSet(
    "r", "1", (AtomTypingRule("rule", "A", "[#6:1]"),),
    hydrogen_policy="implicit_allowed",
)
graph = graph_signature(topology)
rules = ruleset_signature(ruleset)
typing = AtomTypingResult(
    "r", "1",
    {10: AtomTypeAssignment(10, "A", ("rule",), ("rule",))},
    {10: SiteTypingDiagnostic(10, "assigned", ("rule",), (), ("rule",), ("A",))},
    True, (), (), graph, rules,
    typing_signature(graph, rules, engine_name="precomputed", engine_version="1"),
    "precomputed", "1", {},
)
record = LennardJonesParameter("lj", "A", 0.2, 0.3, "synthetic", "p", "1")
library = ParameterLibrary(
    "p", "1", "r", "1", rules, (record,),
    description="synthetic test", source="synthetic",
)
result = ParameterAssignmentEngine().assign(topology, typing, ruleset, library)
assert result.complete_supported_scope
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
