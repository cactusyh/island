"""Automatic OPLS selection, explicit asset installation, immutable assignment.

No source discovery or numerical engine construction occurs in this module.
Profile and charge defaults come from the curated registry, never family guesses.
"""

import json
import os
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace

from island.exceptions import IslandError
from island.forcefields.external_pdms import LIMITS
from island.forcefields.pdms_opls import PROFILE, PDMSOPLSError
from island.forcefields.selection.inventory import canonical, digest, requirements

from .registry import (
    ARCHIVE_SHA256,
    NATIVE,
    canonical_registry,
    validate_registry,
)


class OPLSResolutionError(IslandError):
    def __init__(self, status, diagnostics, *, result=None):
        self.status, self.diagnostics, self.result = status, diagnostics, result
        super().__init__(status + ": " + canonical(diagnostics))


def require(condition, status, stage, **details):
    if not condition:
        raise OPLSResolutionError(status, [{"stage": stage, **details}])


@dataclass(frozen=True)
class OPLSAssets:
    """Explicit one-time local installation; no profile or charge selection.

    Paths never participate in assignment identities, so source relocation is safe.
    License scope stays user-authorized local research, not redistribution.
    """

    archive_path: object = None
    archive_sha256: str | None = None
    native_source_path: object = None
    native_source_sha256: str | None = None
    use_scope: str = "local_research"


_ASSETS = ContextVar("island_opls_assets", default=None)


@contextmanager
def use_opls_assets(assets):
    require(
        type(assets) is OPLSAssets,
        "asset_unavailable",
        "explicit_asset_configuration_required",
    )
    token = _ASSETS.set(assets)
    try:
        yield
    finally:
        _ASSETS.reset(token)


def configured_assets():
    assets = _ASSETS.get()
    if assets is not None:
        return assets
    path = os.environ.get("ISLAND_OPLS_ASSETS")
    if not path:
        return OPLSAssets()
    expected = os.environ.get("ISLAND_OPLS_ASSETS_SHA256")
    raw = Path(path).read_bytes()
    require(
        expected and sha256(raw).hexdigest() == expected,
        "asset_unavailable",
        "independent_asset_configuration_hash_required",
    )
    return OPLSAssets(**json.loads(raw))


@dataclass(frozen=True)
class ForceFieldRequest:
    family: str = "opls"
    final_graph: object = None
    typed_graph: object = None
    graph_charges: object = None
    # Expert-only constraints; never needed for automatic selection.
    profile: str | None = None
    charge_method: str | None = None


def graph_binding(system):
    """Chemical identity excludes coordinates, names and box; retains stable order.

    Execution uses a separate coordinate/box-bound J47 preparation.
    """
    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    system.validate()
    system.topology.validate_bond_graph()
    system_to_rdkit_graph(system)
    sites = []
    for i, atom in system.topology.sites.items():
        require(
            i == atom.id and type(i) is int,
            "typing_ambiguous",
            "stable_atom_ids_required",
        )
        sites.append(
            {
                "id": i,
                "element": atom.element,
                "atomic_number": atom.atomic_number,
                "mass": atom.mass,
                "formal_charge": atom.formal_charge,
                "radical_electrons": atom.metadata.get("radical_electrons", 0),
                "isotope": atom.metadata.get("isotope", 0),
                "molecule": atom.metadata.get("chain_id"),
            }
        )
    return {
        "atom_order": list(system.topology.sites),
        "sites": sites,
        "bonds": [
            [*k, b.order, b.aromatic] for k, b in sorted(system.topology.bonds.items())
        ],
        "components": sorted(sorted(c) for c in system.topology.connected_components()),
        "impropers": [list(k) for k in system.topology.impropers],
        "polymer_charge_partition": {
            k: system.metadata.get("polymer", {}).get(k)
            for k in (
                "architecture",
                "head_site_id",
                "tail_site_id",
                "chain_id",
                "number_of_repeat_units",
                "end_groups",
                "source_psmiles",
            )
        },
    }


def rank_profiles(entries):
    """Pure deterministic conflict detection; cannot grant registry authority."""
    require(bool(entries), "parameter_incomplete", "no_compatible_opls_profile")
    entries = sorted(
        entries,
        key=lambda e: (-e["specificity"], -e["priority"], e["profile_identity"]),
    )
    best = entries[0]
    tied = [
        e
        for e in entries
        if (e["specificity"], e["priority"]) == (best["specificity"], best["priority"])
    ]
    # Compare all physical/provenance declarations, not just profile names.
    require(
        len({canonical(e) for e in tied}) == 1,
        "profile_conflict",
        "equal_priority_opls_profiles_disagree",
        profiles=sorted(e["profile_identity"] for e in tied),
    )
    return best


def _anchors(system, graph):
    """Recognize actual source domain before authorizing a charge partition."""
    from island.forcefields.oplsaa.pdms_conventions import _partition
    from island.forcefields.pdms_opls import final_graph

    source_graph = final_graph(system)
    elements = {a["id"]: a["element"] for a in graph["sites"]}
    anchors = []
    for component in source_graph["components"]:
        oxygens = [i for i in component if elements[i] == "O"]
        anchor = None
        if len(oxygens) > 1:
            polymer = system.metadata.get("polymer", {})
            # Only a builder-recorded head site identifies the source E fragment.
            # No arbitrary tie-break on symmetric ends, geometry or molecule name.
            anchor = polymer.get("head_site_id")
            require(
                not polymer
                or (
                    polymer.get("architecture") == "linear"
                    and polymer.get("number_of_repeat_units") == len(oxygens)
                    and polymer.get("end_groups", {}).get("policy")
                    == "explicit_single_dummy_single_bond_v1"
                ),
                "typing_ambiguous",
                "unverified_builder_charge_partition",
            )
        subset = {
            **source_graph,
            "components": [component],
            "sites": [a for a in source_graph["sites"] if a["id"] in component],
        }
        view = SimpleNamespace(
            topology=SimpleNamespace(
                sites={i: system.topology.sites[i] for i in component},
                neighbors=system.topology.neighbors,
            )
        )
        partition = _partition(view, subset, anchor)
        anchors.append(partition["e_oxygen_site_id"])
    # Historical HMDS assignment used the unique-O rule with an empty tuple.
    if all(sum(elements[i] == "O" for i in c) == 1 for c in source_graph["components"]):
        return ()
    return tuple(anchors)


def _record(value):
    value = value.payload if hasattr(value, "payload") else value
    require(type(value) is dict, "source_incompatible", "graph_record_required")
    return value


def _compare_record(supplied, expected, stage):
    if supplied is not None:
        require(
            canonical(_record(supplied)) == canonical(expected),
            "charge_incomplete"
            if stage == "charge_record_mismatch"
            else "source_incompatible",
            stage,
            expected_identity=digest(expected),
            supplied_identity=digest(_record(supplied)),
        )


def _pin_asset(path, expected, pin):
    require(
        path is not None and expected is not None,
        "asset_unavailable",
        "explicit_source_path_and_independent_hash_required",
    )
    require(
        expected == pin,
        "source_incompatible",
        "unregistered_source_hash",
        expected=pin,
        supplied=expected,
    )
    raw = Path(path).read_bytes()
    require(
        sha256(raw).hexdigest() == pin, "source_incompatible", "source_hash_mismatch"
    )


@dataclass(frozen=True)
class Preparation:
    json_text: str
    request: ForceFieldRequest
    assets: OPLSAssets = field(repr=False)
    _backend: object = field(default=None, repr=False, compare=False)

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    def validate_integrity(self, system):
        current = prepare_forcefield(system, self.request, assets=self.assets)
        require(
            current.json_text == self.json_text,
            "source_incompatible",
            "stale_or_tampered_opls_assignment",
        )
        p = current.payload
        if p["status"] not in ("prepared_static_only", "prepared"):
            raise OPLSResolutionError(p["status"], p["diagnostics"], result=self)
        return current._backend

    def require_md(self, system):
        backend = self.validate_integrity(system)
        require(
            self.payload["stages"]["lammps_MD_verified"],
            "scientifically_unqualified",
            "lammps_static_gates_not_closed",
            missing=self.payload["missing_conventions"],
        )
        return backend


def prepare_forcefield(system, request=None, *, assets=None, registry=None):
    """Inspect exact chemistry, choose a source and charge method, fail closed.

    Optional graph records are constraints, never parameter/charge replacements.
    Returned source assignments and charges are owned immutable JSON snapshots.
    """
    request = ForceFieldRequest() if request is None else request
    registry = canonical_registry() if registry is None else registry
    p = {
        "schema": "island_canonical_opls_assignment_v1",
        "family": "opls",
        "registry_identity": canonical_registry()["identity"],
        "selected_profile": None,
        "source_kind": None,
        "source_identity": None,
        "graph_identity": None,
        "final_graph": None,
        "typed_graph": None,
        "graph_charges": None,
        "charge_origin": None,
        "charge_method": None,
        "domain": None,
        "attempted_profiles": [],
        "assigned_terms": [],
        "missing_terms": [],
        "missing_atom_environments": [],
        "missing_charges": [],
        "missing_conventions": [],
        "diagnostics": [],
        "backend_status": {},
        "stages": {
            k: False
            for k in (
                "typing_complete",
                "charge_complete",
                "parameter_complete",
                "executable",
                "numerically_verified",
                "scientifically_validated",
                "lammps_MD_verified",
            )
        },
        "production_validated": False,
        "simulation_readiness": "not_established",
        "limitations": LIMITS,
    }
    backend = None
    selected_assets = OPLSAssets()
    stage = "request"
    try:
        require(
            type(request) is ForceFieldRequest and request.family == "opls",
            "source_incompatible",
            "opls_public_family_required",
        )
        stage = "final_graph"
        graph = graph_binding(system)
        p.update(final_graph=graph, graph_identity=digest(graph))
        if request.final_graph is not None:
            if hasattr(request.final_graph, "validate_integrity"):
                request.final_graph.validate_integrity(system)
            else:
                _compare_record(request.final_graph, graph, "final_graph_mismatch")
        p["missing_terms"] = requirements(system, "opls")
        elements = {r["element"] for r in graph["sites"]}
        candidates = []
        anchors = ()
        for entry in sorted(registry["entries"], key=lambda e: e["profile_identity"]):
            attempt = {
                "profile": entry["profile_identity"],
                "status": "domain_mismatch",
            }
            p["attempted_profiles"].append(attempt)
            unsupported = [
                r
                for r in graph["sites"]
                if r["element"] not in entry["supported_elements"]
            ]
            if unsupported:
                attempt["missing_atom_environments"] = unsupported
                continue
            if entry["profile_identity"] == PROFILE:
                try:
                    anchors = _anchors(system, graph)
                except Exception as error:  # domain predicates are always classified
                    attempt.update(
                        status=getattr(error, "status", "typing_ambiguous"),
                        diagnostics=getattr(
                            error,
                            "diagnostics",
                            [{"stage": "domain_validation", "reason": str(error)}],
                        ),
                    )
                    if "Si" in elements:
                        raise OPLSResolutionError(
                            attempt["status"], attempt["diagnostics"]
                        ) from error
                    continue
            attempt["status"] = "domain_matched"
            candidates.append(entry)
        require(
            bool(candidates),
            "typing_ambiguous",
            "unsupported_opls_atom_environments",
            environments=graph["sites"],
        )
        stage = "profile_resolution"
        selected = rank_profiles(candidates)
        validate_registry(registry)
        name = selected["profile_identity"]
        require(
            request.profile is None or request.profile == name,
            "source_incompatible",
            "expert_profile_conflicts_with_exact_domain",
            requested=request.profile,
            selected=name,
        )
        require(
            request.charge_method is None
            or request.charge_method == selected["charge_method"],
            "charge_incomplete",
            "expert_charge_method_conflicts_with_profile",
        )
        p.update(
            selected_profile=name,
            source_kind=selected["source_kind"],
            source_identity=digest(selected["source_hashes"]),
            charge_method=selected["charge_method"],
            charge_origin=selected["charge_origin"],
            domain=selected["domain_restrictions"],
            backend_status=selected["backend_support"],
        )
        stage = "source"
        selected_assets = configured_assets() if assets is None else assets
        require(
            type(selected_assets) is OPLSAssets
            and selected_assets.use_scope == "local_research",
            "asset_unavailable",
            "local_research_assets_required",
        )
        if name == PROFILE:
            from island.forcefields.pdms_lammps import ForceFieldRequest as LRequest
            from island.forcefields.pdms_lammps import (
                OPLSOptions,
            )
            from island.forcefields.pdms_lammps import (
                prepare_forcefield as prepare_lammps,
            )

            _pin_asset(
                selected_assets.archive_path,
                selected_assets.archive_sha256,
                ARCHIVE_SHA256,
            )
            backend = prepare_lammps(
                system,
                LRequest(
                    "opls",
                    OPLSOptions(
                        name,
                        selected_assets.archive_path,
                        selected_assets.archive_sha256,
                        anchors,
                    ),
                ),
            )
            data = backend.payload
            require(
                data["status"] == "prepared_static_only",
                "charge_incomplete"
                if "source_charge_literal_missing"
                in canonical(data.get("diagnostics", []))
                else data["status"],
                "source_bound_preparation_failed",
                diagnostics=data["diagnostics"],
            )
            raw = backend.validate_integrity(system).payload
            p["typed_graph"] = {
                "graph_identity": p["graph_identity"],
                "profile_identity": name,
                "source_identity": p["source_identity"],
                "atom_types": [
                    {
                        k: r[k]
                        for k in (
                            "site_id",
                            "author_type_id",
                            "fragment",
                            "molecule_id",
                        )
                    }
                    for r in raw["mapping"]
                ],
                "e_oxygen_site_ids": list(anchors),
            }
            charges = {r["site_id"]: r["charge_e"] for r in raw["mapping"]}
            p["graph_charges"] = {
                "graph_identity": p["graph_identity"],
                "profile_identity": name,
                "source_identity": p["source_identity"],
                "charge_origin": raw["charge_policy"],
                "charge_generation_method": selected["charge_method"],
                "unit": "elementary_charge",
                "vector": [[i, charges[i]] for i in graph["atom_order"]],
                "component_totals": [
                    str(sum((Fraction(charges[i]) for i in c), Fraction()))
                    for c in graph["components"]
                ],
                "evidence": raw["source"]["source_hashes"],
                "normalization": False,
                "redistribution": False,
            }
            p.update(
                assigned_terms=raw["interactions"],
                source_evidence=raw["source"],
                source_assignment_identity=raw["identity"],
                execution_input_identity=data["input_identity"],
                lammps_preparation_identity=data["identity"],
                missing_terms=[],
                status="prepared_static_only",
                missing_conventions=[
                    "identical_PPPM_Ewald_solver_splitting_mesh_contract",
                    "qualified_virial_acceptance_contract",
                    "exact_LAMMPS_checkpoint_and_complete_RNG_continuation",
                ],
            )
            p["stages"].update(
                typing_complete=True,
                charge_complete=True,
                parameter_complete=True,
                executable=True,
            )
        else:
            from island.forcefields.preparation import OPLSOptions
            from island.forcefields.selection import (
                prepare_polymer_forcefield as legacy,
            )

            _pin_asset(
                selected_assets.native_source_path,
                selected_assets.native_source_sha256,
                selected["source_hashes"]["oplsaa.xml"],
            )
            backend = legacy(
                system,
                "opls",
                source_profile=NATIVE,
                charge_method="native",
                options=OPLSOptions(selected_assets.native_source_path),
            )
            data = backend.payload
            p.update(
                status=data["status"],
                diagnostics=data["diagnostics"],
                assigned_terms=[
                    r for r in data["requirements"] if r["status"] == "assigned"
                ],
                missing_terms=[
                    r for r in data["requirements"] if r["status"] != "assigned"
                ],
                typed_graph=data.get("typing_evidence")
                or {"identity": data["typed_graph_identity"]},
                graph_charges=data.get("charge_evidence"),
            )
            p["stages"].update(data["stages"])
        _compare_record(
            request.typed_graph, p["typed_graph"], "typed_graph_record_mismatch"
        )
        _compare_record(
            request.graph_charges, p["graph_charges"], "charge_record_mismatch"
        )
    except Exception as error:  # noqa: BLE001 -- classified public boundary
        if isinstance(error, (OPLSResolutionError, PDMSOPLSError)) or hasattr(
            error, "diagnostics"
        ):
            status = getattr(error, "status", "source_incompatible")
            diagnostics = error.diagnostics
        else:
            status = (
                "typing_ambiguous"
                if stage == "final_graph"
                else "asset_unavailable"
                if isinstance(error, OSError)
                else "source_incompatible"
            )
            diagnostics = [{"stage": stage, "reason": str(error)}]
        p.update(status=status, diagnostics=diagnostics)
        p["stages"]["executable"] = False
        if status in ("charge_unavailable", "charge_incomplete"):
            p["stages"]["charge_complete"] = False
            if not p["missing_terms"] and p["final_graph"] is not None:
                p["missing_terms"] = [
                    {"family": "charge", "sites": [i], "reason": diagnostics}
                    for i in p["final_graph"]["atom_order"]
                ]
        backend = None
        p["missing_atom_environments"] = [
            a
            for r in p["attempted_profiles"]
            for a in r.get("missing_atom_environments", [])
        ]
        p["missing_atom_environments"].extend(
            {"profile": r["profile"], "diagnostics": r["diagnostics"]}
            for r in p["attempted_profiles"]
            if r.get("status") in ("typing_blocked", "typing_ambiguous")
            and r.get("diagnostics")
        )
        p["missing_charges"] = [
            r for r in p["missing_terms"] if r.get("family") == "charge"
        ]
    p["identity"] = digest(p)
    return Preparation(canonical(p), request, selected_assets, backend)


def prepare_polymer_forcefield(system, forcefield="opls", **kwargs):
    if forcefield == "opls":
        return prepare_forcefield(system, ForceFieldRequest(), **kwargs)
    from island.forcefields.selection import prepare_polymer_forcefield as legacy

    return legacy(system, forcefield, **kwargs)


def prepare_polymer(
    *,
    psmiles,
    forcefield="opls",
    coordinates=None,
    box=None,
    dp=None,
    head_end_group=None,
    tail_end_group=None,
    **kwargs,
):
    """Build explicit-H chain with explicit terminal groups, then resolve chemistry.

    Returns a PolymerPreparation; dp/end groups define chemistry, not a charge model.
    """
    from island.builders import build_linear_polymer
    from island.core import Coordinates

    system = build_linear_polymer(
        psmiles,
        dp=dp,
        head_end_group=head_end_group,
        tail_end_group=tail_end_group,
        add_hydrogens=True,
        random_seed=2026,
    )
    if coordinates is not None:
        system.coordinates = (
            coordinates
            if isinstance(coordinates, Coordinates)
            else Coordinates(dict(zip(system.topology.sites, coordinates, strict=True)))
        )
    system.box = box
    return PolymerPreparation(
        system, prepare_polymer_forcefield(system, forcefield, **kwargs)
    )


@dataclass(frozen=True)
class PolymerPreparation:
    system: object
    assignment: object

    @property
    def payload(self):
        return self.assignment.payload

    def validate_integrity(self):
        return self.assignment.validate_integrity(self.system)


def create_evaluator(system, preparation):
    require(
        type(preparation) is Preparation,
        "source_incompatible",
        "canonical_preparation_required",
    )
    backend = preparation.validate_integrity(system)
    if preparation.payload["selected_profile"] == PROFILE:
        from island.forcefields.research import create_evaluator as existing

        return existing(system, backend.validate_integrity(system))
    from island.forcefields.selection import create_evaluator as existing

    return existing(system, backend)
