"""Scoped convention resolution without turning engine defaults into source data.

The canonical OPLS assignment and all four J50 records remain byte-identical.
Resolution of an engine implementation is independent of source authorization.
"""

import json
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path

from island.forcefields.opls_registry.api import require
from island.forcefields.opls_registry.registry import PAPER_SHA256
from island.forcefields.selection.inventory import canonical, digest
from island.opls_conventions import LAYERS, signed
from island.opls_conventions import prepare_polymer_forcefield as legacy_prepare

from .backend_evidence import CITATIONS

TEN_CONVENTIONS = (
    "dielectric",
    "coulomb_prefactor",
    "periodic_image_policy",
    "precision",
    "neighbor_list_policy",
    "pressure_convention",
    "virial_convention",
    "pair_table_bits",
    "erfc_kernel",
    "lj_tail_correction",
)
REMAINING_SOURCE_GAPS = tuple(n for n in TEN_CONVENTIONS if n != "lj_tail_correction")
RESOLUTIONS = (
    "source_resolved",
    "converter_resolved",
    "backend_resolved",
    "derived_resolved",
    "unresolved",
    "conflicting",
)
EVIDENCE_CLASSES = (
    "explicit_source_statement",
    "source_equation_or_table",
    "source_code_or_input",
    "converter_defined",
    "backend_defined",
    "independently_derived",
    "unresolved",
)
STATUS_FOR = dict(
    zip(
        EVIDENCE_CLASSES,
        (
            "source_resolved",
            "source_resolved",
            "source_resolved",
            "converter_resolved",
            "backend_resolved",
            "derived_resolved",
            "unresolved",
        ),
    )
)
SUPPLEMENT_SHA256 = "6e86aa64085466e286b5c59878c2288f5f16da1b7a5298d37ba1833a88cd602f"
_SOURCES = ContextVar("opls_convention_sources", default=None)


@dataclass(frozen=True)
class ConventionSources:
    """Explicit local evidence installation; never discovers or downloads assets."""

    paper: object
    supplement: object

    def validate(self):
        records = {}
        for name, path, expected in (
            ("jp4c08471.pdf", self.paper, PAPER_SHA256),
            ("jp4c08471_si_001.pdf", self.supplement, SUPPLEMENT_SHA256),
        ):
            try:
                observed = sha256(Path(path).read_bytes()).hexdigest()
            except (OSError, TypeError) as error:
                require(
                    False,
                    "asset_unavailable",
                    "convention_source_unavailable",
                    file=name,
                    reason=str(error),
                )
            require(
                observed == expected,
                "source_incompatible",
                "convention_source_replacement",
                file=name,
                expected_sha256=expected,
                observed_sha256=observed,
            )
            records[name] = expected
        return records


@contextmanager
def use_convention_sources(sources):
    require(
        type(sources) is ConventionSources,
        "source_incompatible",
        "explicit_convention_sources_required",
    )
    sources.validate()
    token = _SOURCES.set(sources)
    try:
        yield sources
    finally:
        _SOURCES.reset(token)


def reconcile(setting, candidates):
    """Compare claims for ONE scope. Disagreement is retained, never prioritized away.

    This is an audit primitive, not an API for registering arbitrary conventions.
    Backend diagnostic experiments live in a separate record from source claims.
    """
    claims = sorted(candidates, key=canonical)
    require(
        all(c.get("evidence_class") in EVIDENCE_CLASSES for c in claims),
        "source_incompatible",
        "unknown_convention_evidence_class",
    )
    for c in claims:
        if c["value"] is not None:
            require(
                c["evidence_class"] != "unresolved" and c.get("evidence"),
                "source_incompatible",
                "uncited_convention_claim",
                setting=setting,
            )
    values = {canonical(c["value"]) for c in claims if c["value"] is not None}
    status = "unresolved"
    value = None
    if len(values) > 1:
        status = "conflicting"
    elif values:
        value = json.loads(next(iter(values)))
        # Equal values may have several authorities; retain every original claim.
        status = next(
            s
            for s in RESOLUTIONS[:4]
            if any(STATUS_FOR[c["evidence_class"]] == s for c in claims)
        )
    return signed(
        {
            "setting": setting,
            "resolution": status,
            "value": value,
            "candidates": claims,
            "execution_rejected": status in ("unresolved", "conflicting"),
        }
    )


def reject_conflicts(rows):
    conflicts = [r for r in rows if r["resolution"] == "conflicting"]
    require(
        not conflicts,
        "source_incompatible",
        "conflicting_conventions",
        conflicts=conflicts,
    )


def _source_evidence(evidence):
    result = []
    for e in evidence:
        if "source_file" in e:
            result.append(
                {
                    "doi": "10.1021/acs.jpcb.4c08471",
                    "citation": "Xiang et al., JPCB 129 (2025), 1864-1873",
                    **e,
                }
            )
        elif "file" in e:
            result.extend(
                {
                    "source_file": e["file"],
                    "source_sha256": e["sha256"],
                    "doi": None,
                    "citation": "Pinned LAMMPS definition",
                    **c,
                }
                for c in e["citations"]
            )
    return result


def _table(prior):
    p = prior.payload
    source, backend = p[LAYERS[0]]["settings"], p[LAYERS[2]]["settings"]
    rows = []
    for name, s in sorted(source.items()):
        missing = s["value"] is None
        kind = (
            "unresolved"
            if missing
            else (
                "converter_defined"
                if s["classification"] == "derived_by_documented_formula"
                else "explicit_source_statement"
                if name == "boundary_conditions"
                else "source_code_or_input"
            )
        )
        evidence = _source_evidence(s["evidence"])
        value = s["value"]
        if name == "lj_tail_correction":
            # A documented restriction of the explicitly selected functional
            # form, NOT inference from pair_modify's default tail flag.
            missing, kind, value = False, "converter_defined", "no"
            evidence = (
                _source_evidence(source["pair_style"]["evidence"]) + CITATIONS[name]
            )
        claim = {"value": value, "evidence_class": kind, "evidence": evidence}
        r = reconcile(name, [claim])
        b = backend[name]
        backend_cites = CITATIONS.get(name, evidence)
        implementation = {
            "version": "equation_preserving_real_to_OpenMM_v1",
            "hashes": p[LAYERS[1]]["implementation_hashes"],
        }
        if name == "lj_tail_correction":
            implementation = {
                "version": "island.opls-resolution.charmm-no-tail.v1",
                "code": "island.opls_resolution._table",
                "engine_definition_sha256": CITATIONS[name][-1]["source_sha256"],
                "source_command": "pair_style lj/charmm/coul/long 9.0 11.0",
            }
        # A backend implementation can be reconstructed without proving what
        # software/settings the author actually ran. These statuses never merge.
        rows.append(
            signed(
                {
                    **{k: v for k, v in r.items() if k != "identity"},
                    "source_authoritative": kind
                    in ("explicit_source_statement", "source_code_or_input"),
                    "source_status": "unresolved" if missing else r["resolution"],
                    "units": _units(name),
                    "conversion_formula": (
                        "lj/charmm/coul/long cannot use pair_modify tail: switched LJ is zero beyond outer cutoff; no extra long-range tail term"
                        if name == "lj_tail_correction"
                        else s.get("formula", "identity; no physical parameter change")
                    ),
                    "conversion_implementation": implementation,
                    "independent_reconstruction": {
                        "method": "raw ZIP/PDF/engine lines outside production parser",
                        "result": "source_absence_retained"
                        if missing
                        else (
                            "J51_independent_style_restriction_reconstruction"
                            if name == "lj_tail_correction"
                            else "verified_in_pinned_J50_audit"
                        ),
                        "receipt_sha256": "f61eda5b5b423b2501d29d665a196da9be6cc1374ab4678458b0160e8d31d7b5",
                    },
                    "missing_evidence": s.get("missing_evidence") if missing else None,
                    "searched_include_chain": s.get("searched_include_chain", []),
                    "absence_inventory_sha256": s.get("absence_evidence_sha256"),
                    "backend": {
                        "resolution": "backend_resolved",
                        "evidence_class": "backend_defined",
                        "value": b["value"],
                        "source_authoritative": False,
                        "evidence": backend_cites,
                        "policy": prior.policy,
                        "executable_sha256": p[LAYERS[2]]["engine_sha256"],
                        "reconstruction": "independent engine-source/build audit required; not author proof",
                    },
                    "evidence_checksum": digest(
                        {"source": claim, "backend": backend_cites}
                    ),
                }
            )
        )
    return rows


def _units(name):
    if name == "coulomb_prefactor":
        return "kcal angstrom/(mol e^2); multiply by 4.184 for kJ"
    if "angstrom" in name:
        return "angstrom"
    if name == "pressure_convention":
        return "atm (configurational proxy only)"
    if name == "virial_convention":
        return "kJ/mol; P_atm*V_angstrom^3*4.184/68568.415"
    return "explicit in interpreted value; otherwise dimensionless"


def _document(prior, sources):
    pins = sources.validate()
    rows = _table(prior)
    reject_conflicts(rows)
    unresolved = [r for r in rows if r["resolution"] == "unresolved"]
    require(
        {r["setting"] for r in unresolved} == set(REMAINING_SOURCE_GAPS),
        "source_incompatible",
        "unexpected_source_gap_inventory",
    )
    return signed(
        {
            "schema": "island.opls-convention-resolution.v1",
            "family": "opls",
            "selected_profile": prior.payload["selected_profile"],
            "charge_origin": prior.payload["charge_origin"],
            "domain": prior.payload["domain"],
            "four_layer_identity": prior.identity,
            **{k: prior.payload[k] for k in LAYERS},
            "source_assets": pins,
            "rows": rows,
            "ten_conventions": list(TEN_CONVENTIONS),
            "ten_convention_resolutions": [
                next(r for r in rows if r["setting"] == n) for n in TEN_CONVENTIONS
            ],
            "unresolved_conventions": unresolved,
            "conflicts": [],
            "source_scope_review": {
                "status": "unresolved",
                "archive_pair_style": "lj/charmm/coul/long 9.0 11.0",
                "SI_Table_S4_pair_style": "lj/cut/coul/long 11 11",
                "reason": "SI OPLS-AA table is not explicitly mapped to serialized DFT/NPA branch; do not combine or choose by a default",
                "supplement_sha256": SUPPLEMENT_SHA256,
                "evidence": {
                    "source_file": "jp4c08471_si_001.pdf",
                    "source_sha256": SUPPLEMENT_SHA256,
                    "doi": "10.1021/acs.jpcb.4c08471",
                    "page": 6,
                    "table": "S4",
                    "raw_text": "pair_style lj/cut/coul/long 11 11",
                },
                "required": "author DFT/NPA-scoped reconciliation of Table S4 and archive include chain",
            },
            "execution_mode": prior.mode,
            "backend_policy": prior.policy,
            "research_MD_available": False,
            "default_MD_enabled": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "limitations": prior.payload["limitations"],
        }
    )


@dataclass(frozen=True)
class ResolutionPreparation:
    json_text: str
    prior: object = field(repr=False)
    sources: ConventionSources = field(repr=False)

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    @property
    def assignment(self):
        return self.prior.assignment

    def validate_integrity(self, system):
        self.prior.validate_integrity(system)
        expected = _document(self.prior, self.sources)
        require(
            self.json_text == canonical(expected),
            "source_incompatible",
            "mutated_convention_resolution",
            expected_identity=expected["identity"],
        )
        return self.assignment.validate_integrity(system)

    def require_execution(self, system, action="static", *, research_opt_in=False):
        self.validate_integrity(system)
        p = self.payload
        reject_conflicts(p["rows"])
        if self.prior.mode == "research_opt_in" and research_opt_in is True:
            require(
                not p["unresolved_conventions"],
                "parameter_incomplete",
                "source_conventions_unresolved",
                missing=p["unresolved_conventions"],
                source_scope_review=p["source_scope_review"],
            )
        return self.prior.require_execution(
            system, action, research_opt_in=research_opt_in
        )


def prepare_polymer_forcefield(system, forcefield="opls", *, sources=None, **kwargs):
    """Canonical family only; sources are installed once, with no profile/charge input."""
    sources = sources if sources is not None else _SOURCES.get()
    require(
        type(sources) is ConventionSources,
        "asset_unavailable",
        "explicit_convention_sources_required",
    )
    prior = legacy_prepare(system, forcefield=forcefield, **kwargs)
    return ResolutionPreparation(canonical(_document(prior, sources)), prior, sources)


__all__ = [
    "RESOLUTIONS",
    "TEN_CONVENTIONS",
    "ConventionSources",
    "ResolutionPreparation",
    "prepare_polymer_forcefield",
    "use_convention_sources",
]
