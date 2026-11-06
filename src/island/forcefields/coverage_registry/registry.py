"""Evidence-aware, graph-scoped coverage registry. No runtime registration API."""

from dataclasses import dataclass, field

from island.forcefields.assignment_core import builtin_contracts, resolve, type_graph
from island.forcefields.assignment_core.inventory import requests
from island.forcefields.assignment_core.schema import (
    CROSS_TERMS,
    READINESS,
    Coverage,
    EvidenceRecord,
    ForceFieldProfile,
    data,
    freeze_sequences,
    identity,
    require,
    unique,
)

from .catalog import ChemicalMotifCatalog, MotifMatch

DIMENSIONS = (
    ("bonded_terms", ("bond", "angle")),
    ("torsion", ("torsion",)),
    ("improper", ("improper",)),
    ("nonbonded", ("nonbonded",)),
    ("charge", ("charge",)),
    ("cross_terms", CROSS_TERMS),
)


@dataclass(frozen=True)
class CapabilityCell:
    profile_id: str
    profile_identity: str
    family: str
    functional_class: str
    motif: str
    topology_class: str
    sites: tuple[int, ...]
    coverage: tuple[tuple[str, Coverage], ...]
    request_ids: tuple[str, ...]
    unresolved_term_ids: tuple[str, ...]
    evidence_statuses: tuple[str, ...]
    evidence_ids: tuple[str, ...]
    coverage_decision: str
    reasons: tuple[str, ...]


@dataclass(frozen=True)
class ProfileCapabilityMatrix:
    cells: tuple[CapabilityCell, ...]
    applicability: str = "exact_final_graph_only_not_general_motif_support"

    def __post_init__(self):
        freeze_sequences(self)
        unique(
            [(c.profile_id, c.motif, c.topology_class) for c in self.cells],
            "duplicate_capability_cell_identity",
        )

    def lookup(self, profile, motif, topology):
        return next(
            (
                c
                for c in self.cells
                if (c.profile_id, c.motif, c.topology_class)
                == (profile, motif, topology)
            ),
            None,
        )


@dataclass(frozen=True)
class CoverageReport:
    graph_identity: str
    registry_identity: str
    catalog_identity: str
    topology_class: str
    motif_matches: tuple[MotifMatch, ...]
    capability_matrix: ProfileCapabilityMatrix
    profile_reports: tuple
    unresolved_motifs: tuple[tuple[str, tuple[str, ...]], ...]
    unresolved_terms: tuple[tuple[str, tuple], ...]
    evidence_records: tuple[tuple[str, tuple[EvidenceRecord, ...]], ...]
    source_hashes: tuple[str, ...]
    ranked_profiles: tuple[str, ...]
    blocking_reasons: tuple[str, ...]
    replay_identity: str
    automatic_assignment_enabled: bool = field(default=False, init=False)
    production_validated: bool = field(default=False, init=False)
    simulation_readiness: str = field(default="not_established", init=False)
    accepted_coefficients: int = field(default=0, init=False)
    runtime_registration: bool = field(default=False, init=False)

    def __post_init__(self):
        freeze_sequences(self)


@dataclass(frozen=True)
class ForceFieldProfileRegistry:
    profiles: tuple[ForceFieldProfile, ...] = ()
    trusted_evidence: tuple[EvidenceRecord, ...] = ()
    catalog: ChemicalMotifCatalog = field(default_factory=ChemicalMotifCatalog)

    def __post_init__(self):
        freeze_sequences(self)
        unique([p.id for p in self.profiles], "duplicate_profile_identity")
        unique(
            [e.id for e in self.trusted_evidence], "duplicate_trusted_evidence_identity"
        )
        object.__setattr__(
            self, "profiles", tuple(sorted(self.profiles, key=lambda p: p.id))
        )
        object.__setattr__(
            self,
            "trusted_evidence",
            tuple(sorted(self.trusted_evidence, key=lambda e: e.id)),
        )

    @classmethod
    def contracts(cls):
        """Five empty scientific contracts, not installed force-field parameter sets."""
        return cls(builtin_contracts())

    @property
    def identity(self):
        return identity(self)

    def coverage(self, graph):
        matches = self.catalog.detect(graph)
        decision = resolve(graph, self.profiles, self.trusted_evidence)
        reports = {r.profile_id: r for r in decision.reports}
        by_motif = {}
        for m in matches:
            by_motif.setdefault(m.motif, set()).update(m.sites)
        cells = []
        for p in self.profiles:
            report = reports[p.id]
            typed = type_graph(graph, p, self.trusted_evidence)
            rows = requests(graph, typed, p)
            matched = dict(report.matched_records)
            missing = {m.request_id: m for m in report.missing_terms}
            params = {v.id: v for v in p.parameters}
            rules = {v.id: v for v in p.rules}
            trace = dict(typed.rule_trace)
            for motif, sites in sorted(by_motif.items()):
                # Include boundary interactions, then the full dependency closure.
                selected = {r.id for r in rows if set(r.sites) & sites}
                by_id = {r.id: r for r in rows}
                pending = list(selected)
                while pending:
                    for dep in by_id[pending.pop()].dependencies:
                        if dep not in selected:
                            selected.add(dep)
                            pending.append(dep)
                coverage = []
                for name, families in DIMENSIONS:
                    required = {
                        rid for rid in selected if by_id[rid].family in families
                    }
                    covered = len(required & matched.keys())
                    status = (
                        "not_required_by_class_i_contract"
                        if name == "cross_terms" and p.functional_class == "I"
                        else "complete"
                        if covered == len(required)
                        else "incomplete"
                    )
                    coverage.append((name, Coverage(len(required), covered, status)))
                count = len(sites & dict(typed.assignments).keys())
                coverage.append(
                    (
                        "atom_typing",
                        Coverage(
                            len(sites),
                            count,
                            "complete" if count == len(sites) else "incomplete",
                        ),
                    )
                )
                coverage.append(("topology", dict(report.coverage)["topology"]))
                evidence = {c.evidence_id for c in p.conventions}
                for site in sites:
                    evidence.update(rules[r].evidence_id for r in trace.get(site, ()))
                for rid in selected:
                    evidence.update(params[k].evidence_id for k in matched.get(rid, ()))
                    if rid in missing:
                        evidence.update(
                            params[k].evidence_id for k in missing[rid].candidate_ids
                        )
                # Preserve blocked typing evidence as well, even if no rule survived.
                evidence.update(r.evidence_id for r in p.rules)
                statuses = {e.status for e in p.evidence if e.id in evidence}
                unresolved = tuple(sorted(selected & missing.keys()))
                reasons = set(report.blocking_reasons)
                if motif == "unclassified":
                    reasons.add("motif_not_in_catalog")
                if report.status == "synthetic_complete" and motif != "unclassified":
                    outcome = "provisional"
                    reasons.add("synthetic_software_control_not_scientific_evidence")
                elif (
                    report.status == "coverage_complete_not_activated"
                    and motif != "unclassified"
                ):
                    outcome = "accepted"
                    reasons.add("audited_exact_graph_coverage_only_runtime_disabled")
                else:
                    outcome = "blocked"
                if unresolved or count != len(sites):
                    statuses.add("missing")
                if outcome == "blocked":
                    statuses.add("blocked")
                cells.append(
                    CapabilityCell(
                        p.id,
                        p.identity,
                        p.family,
                        p.functional_class,
                        motif,
                        graph.architecture,
                        tuple(sorted(sites)),
                        tuple(sorted(coverage)),
                        tuple(sorted(selected)),
                        unresolved,
                        tuple(sorted(statuses)),
                        tuple(sorted(evidence)),
                        outcome,
                        tuple(sorted(reasons)),
                    )
                )
        matrix = ProfileCapabilityMatrix(tuple(cells))
        unresolved_motifs = tuple(
            (
                p.id,
                tuple(
                    sorted(
                        c.motif
                        for c in cells
                        if c.profile_id == p.id and c.coverage_decision != "accepted"
                    )
                ),
            )
            for p in self.profiles
        )
        payload = {
            "graph_identity": graph.identity,
            "registry_identity": self.identity,
            "catalog_identity": self.catalog.identity,
            "topology_class": graph.architecture,
            "motif_matches": matches,
            "capability_matrix": matrix,
            "profile_reports": tuple(reports[p.id] for p in self.profiles),
            "unresolved_motifs": unresolved_motifs,
            "unresolved_terms": tuple(
                (p.id, reports[p.id].missing_terms) for p in self.profiles
            ),
            "evidence_records": tuple((p.id, p.evidence) for p in self.profiles),
            "source_hashes": tuple(
                sorted({e.source_sha256 for p in self.profiles for e in p.evidence})
            ),
            "ranked_profiles": decision.ranked_profiles,
            "blocking_reasons": (
                "J83_runtime_assignment_disabled",
                *(("no_candidate_profiles",) if not self.profiles else ()),
            ),
        }
        result = CoverageReport(**payload, replay_identity=identity(payload))
        require(
            all(data(result)[k] == v for k, v in READINESS.items()),
            "readiness_mutation",
        )
        return result
