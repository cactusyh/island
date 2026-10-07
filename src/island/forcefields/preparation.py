"""Explicit adapters over native force-field contracts; no numerical engine."""

import json
from copy import deepcopy
from dataclasses import dataclass, replace
from hashlib import sha256
from pathlib import Path

from island.exceptions import ForceFieldError


class ForceFieldRequestError(ForceFieldError):
    """Invalid selector or backend-specific options, before scientific execution."""


class PreparedForceFieldError(ForceFieldError):
    """Malformed or incompatible facade/native preparation."""


def _require(ok, message, error=PreparedForceFieldError):
    if not ok:
        raise error(message)


@dataclass(frozen=True)
class OPLSOptions:
    """Explicit installed, pinned Foyer XML; native charges only."""

    source_path: Path
    typed_graph: object | None = None
    graph_charges: object | None = None

    def __post_init__(self):
        if self.typed_graph is not None or self.graph_charges is not None:
            from island.graph import UnifiedGraphCharges, UnifiedTypedGraph

            _require(
                type(self.typed_graph) is UnifiedTypedGraph
                and type(self.graph_charges) is UnifiedGraphCharges,
                "OPLS unified typed_graph and graph_charges are required",
                ForceFieldRequestError,
            )
            _require(
                self.typed_graph.payload["force_field"] == "OPLS-AA"
                and self.graph_charges.payload["force_field"] == "OPLS-AA",
                "OPLS records cannot borrow another force field",
                ForceFieldRequestError,
            )
        _require(
            isinstance(self.source_path, (str, Path)),
            "source_path required",
            ForceFieldRequestError,
        )
        object.__setattr__(self, "source_path", Path(self.source_path))


@dataclass(frozen=True)
class PCFFOptions:
    """Pinned FRC with mandatory caller-selected special-pair weights."""

    source_path: Path
    lj: tuple[float, float, float]
    coulomb: tuple[float, float, float]
    typing_profile: str = "island_pcff_acyclic_cho_v1"
    resolution_policy: str | None = None
    source_profile: object | None = None
    typed_graph: object | None = None
    graph_charges: object | None = None

    def __post_init__(self):
        from .pcff import special_pair_policy
        from .pcff.fallbacks import (
            COMPATIBILITY_POLICY,
            DOMAIN_POLICY,
            MSI_POLICY,
            POLICY,
        )

        if self.typed_graph is not None or self.graph_charges is not None:
            from island.graph import UnifiedGraphCharges, UnifiedTypedGraph

            from .pcff.typed_graph import PCFFGraphCharges, PCFFTypedGraph

            valid_old = (
                type(self.typed_graph) is PCFFTypedGraph
                and type(self.graph_charges) is PCFFGraphCharges
            )
            valid_unified = (
                type(self.typed_graph) is UnifiedTypedGraph
                and type(self.graph_charges) is UnifiedGraphCharges
            )
            _require(
                valid_old or valid_unified,
                "Explicit PCFF requires matching typed_graph and graph_charges records",
                ForceFieldRequestError,
            )
            if valid_unified:
                _require(
                    self.typed_graph.payload["force_field"] == "PCFF"
                    and self.graph_charges.payload["force_field"] == "PCFF",
                    "PCFF records cannot borrow another force field",
                    ForceFieldRequestError,
                )
            _require(
                self.source_profile is None
                and self.typing_profile == "island_pcff_acyclic_cho_v1",
                "External records cannot select an automatic or operational profile",
                ForceFieldRequestError,
            )
            self.typed_graph.validate_integrity()
            self.graph_charges.validate_integrity()
            charge = self.graph_charges.payload
            if valid_old:
                _require(
                    charge["typing_identity"] == self.typed_graph.identity
                    and charge["resolution_policy"] == self.resolution_policy,
                    "External typing/charge/policy mismatch",
                    ForceFieldRequestError,
                )
            else:
                _require(
                    charge["typed_graph_identity"] == self.typed_graph.identity,
                    "Unified typing/charge identity mismatch",
                    ForceFieldRequestError,
                )

        if self.source_profile is not None:
            from .pcff.operational_profile import PCFFOperationalSelection

            _require(
                type(self.source_profile) is PCFFOperationalSelection,
                "Typed operational source selection required",
                ForceFieldRequestError,
            )
            profile = self.source_profile.profile().payload
            _require(
                self.typing_profile == profile["typing_profile"]
                and self.resolution_policy == profile["resolution_policy"],
                "Operational typing/resolution options mismatch",
                ForceFieldRequestError,
            )

        _require(
            self.resolution_policy is None
            or (
                self.resolution_policy
                in (POLICY, DOMAIN_POLICY, MSI_POLICY, COMPATIBILITY_POLICY)
                and (
                    self.typed_graph is not None
                    or self.typing_profile != "island_pcff_acyclic_cho_v1"
                )
            ),
            "Invalid PCFF fallback policy/profile",
            ForceFieldRequestError,
        )
        _require(
            self.typing_profile
            not in (
                "island_pcff_source_graph_v2",
                "island_pcff_source_graph_v3",
                "island_pcff_source_graph_v4",
                "island_pcff_source_graph_v5",
                "island_pcff_source_graph_v6",
            )
            or self.resolution_policy
            in (POLICY, DOMAIN_POLICY, MSI_POLICY, COMPATIBILITY_POLICY),
            "Graph v2 requires explicit fallback policy",
            ForceFieldRequestError,
        )

        _require(
            self.typing_profile
            in (
                "island_pcff_acyclic_cho_v1",
                "island_pcff_source_graph_v1",
                "island_pcff_source_graph_v2",
                "island_pcff_source_graph_v3",
                "island_pcff_source_graph_v4",
                "island_pcff_source_graph_v5",
                "island_pcff_source_graph_v6",
            ),
            "Unsupported PCFF typing profile",
            ForceFieldRequestError,
        )

        _require(
            isinstance(self.source_path, (str, Path)),
            "source_path required",
            ForceFieldRequestError,
        )
        try:
            policy = special_pair_policy(lj=self.lj, coulomb=self.coulomb)
        except Exception as exc:
            raise ForceFieldRequestError(
                f"Invalid PCFF special-pair policy: {exc}"
            ) from exc
        object.__setattr__(self, "source_path", Path(self.source_path))
        object.__setattr__(self, "lj", tuple(policy["lj"]))
        object.__setattr__(self, "coulomb", tuple(policy["coulomb"]))


@dataclass(frozen=True)
class ForceFieldRequest:
    """Discriminated selection: gaff/gaff2, oplsaa, or pcff; no defaults."""

    family: str
    options: object
    final_graph: object | None = None
    typed_graph: object | None = None
    graph_charges: object | None = None

    def __post_init__(self):
        from .ambertools import AmberToolsOptions

        expected = {
            "gaff": AmberToolsOptions,
            "gaff2": AmberToolsOptions,
            "oplsaa": OPLSOptions,
            "pcff": PCFFOptions,
        }
        _require(
            type(self.family) is str and self.family in expected,
            "Unknown force-field selector",
            ForceFieldRequestError,
        )
        _require(
            type(self.options) is expected[self.family],
            "Options do not belong to selected force field",
            ForceFieldRequestError,
        )
        if self.family in ("gaff", "gaff2"):
            _require(
                self.options.force_field == self.family,
                "Amber selector/options force-field mismatch",
                ForceFieldRequestError,
            )
        try:
            # Re-run option validation even for reconstructed/tampered instances.
            object.__setattr__(self, "options", replace(self.options))
        except Exception as exc:
            raise ForceFieldRequestError(f"Invalid backend options: {exc}") from exc
        if self.final_graph is not None:
            from island.graph import FinalChemicalGraph

            _require(
                type(self.final_graph) is FinalChemicalGraph,
                "FinalChemicalGraph required",
                ForceFieldRequestError,
            )
            self.final_graph.validate_integrity()
        if self.typed_graph is not None or self.graph_charges is not None:
            from island.graph import UnifiedGraphCharges, UnifiedTypedGraph

            _require(
                type(self.typed_graph) is UnifiedTypedGraph
                and type(self.graph_charges) is UnifiedGraphCharges,
                "Unified typed_graph and graph_charges must be supplied together",
                ForceFieldRequestError,
            )
            _require(
                self.final_graph is not None,
                "Unified records require final_graph",
                ForceFieldRequestError,
            )
            self.typed_graph.validate_integrity(self.final_graph)
            self.graph_charges.validate_integrity(self.final_graph, self.typed_graph)
            family = {
                "pcff": "PCFF",
                "oplsaa": "OPLS-AA",
                "gaff": "GAFF",
                "gaff2": "GAFF2",
            }[self.family]
            _require(
                self.typed_graph.payload["force_field"] == family
                and self.graph_charges.payload["force_field"] == family,
                "Unified records cannot borrow another force field",
                ForceFieldRequestError,
            )
            _require(
                self.graph_charges.payload["typed_graph_identity"]
                == self.typed_graph.identity,
                "Unified graph charge/typed graph identity mismatch",
                ForceFieldRequestError,
            )
            if self.final_graph is not None:
                _require(
                    self.typed_graph.payload["graph_identity"]
                    == self.final_graph.identity,
                    "Unified typed graph/final graph identity mismatch",
                    ForceFieldRequestError,
                )


def _binding(system):
    from .typing.signatures import graph_signature

    system.validate()
    return {
        "graph": graph_signature(system.topology),
        "representation": system.representation,
        "sites": [
            (
                i,
                s.element,
                s.atomic_number,
                s.mass,
                s.metadata.get("cip_label"),
                s.metadata.get("chiral_tag"),
            )
            for i, s in sorted(system.topology.sites.items())
        ],
    }


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _description(system, family, native, source, pcff_profile=None):
    """Delegate complete scientific validation; derive all advertised metadata."""
    if family in ("gaff", "gaff2"):
        from .ambertools import AmberToolsPreparationResult

        _require(
            type(native) is AmberToolsPreparationResult and source is None,
            "Amber requires an AmberToolsPreparationResult and no separate source",
        )
        native.validate_integrity(system)
        r = native.record
        _require(r["requested_force_field"] == family, "Native Amber family mismatch")
        profile = "island_ambertools_" + family + "_adapter_v1"
        native_id = native.record_signature
        sources = {k: deepcopy(r[k]) for k in ("force_field_data", "leaprc")}
        charge = r["charge_method"]
        policy = deepcopy(native.imported_result.nonbonded_policy.__dict__)
    elif family == "oplsaa":
        from .oplsaa import FoyerOPLSSource, OPLSParameterizationResult

        _require(
            type(native) is OPLSParameterizationResult
            and type(source) is FoyerOPLSSource,
            "OPLS requires native parameters and pinned Foyer source",
        )
        native.validate_integrity(system, source)
        profile = "pinned_foyer_oplsaa_parameters_v1"
        native_id, sources = native.identity, source.identity
        charge = sources["charge_policy"]
        policy = {
            "mixing": "geometric",
            "lj14": 0.5,
            "coulomb14": 0.5,
            "scope": "finite_nonperiodic_NoCutoff",
        }
    elif family == "pcff":
        from .pcff import PCFFModelSpecification

        _require(
            type(native) is PCFFModelSpecification and source is None,
            "PCFF requires a native model specification with its owned source",
        )
        native.validate_integrity(system)
        p = native.payload
        _require(
            p["model_definition_complete"], f"Incomplete PCFF model: {p['diagnostics']}"
        )
        profile = p["compatibility_profile"]
        native_id, sources = native.identity, p["source"]
        charge = "source_native_bond_increments_automatic_pcff_v1"
        if p["schema"] == "island_pcff_source_model_v1":
            charge = "source_native_bond_increments_source_graph_v1"
        if p["schema"] == "island_pcff_source_model_v2":
            charge = (
                "source_native_bond_increments_positional_fallbacks_v"
                + profile["resolution_policy"]["name"][-1]
            )
        if p["schema"] == "island_pcff_typed_graph_model_v1":
            record = native.assignment.payload["charge_record"]
            charge = "external_types_" + record["origin"] + "_v1"
        policy = p["special_pairs"]
    else:
        raise PreparedForceFieldError("Unknown prepared family")
    description = {
        "schema": "island_prepared_forcefield_v1",
        "adapter": "island_native_dispatch_v1",
        "family": family,
        "profile": profile,
        "native_identity": native_id,
        "source": sources,
        "charge_method": charge,
        "nonbonded_policy": policy,
        "binding": _binding(system),
        "status": "prepared",
        "limitations": [
            "Native chemistry/source scope only",
            "Finite nonperiodic evaluation",
            "Not an Amber workflow bundle",
            "No scientific suitability attestation",
        ],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    if pcff_profile is not None:
        from .pcff.operational_profile import PCFFOperationalProfile

        _require(
            family == "pcff" and type(pcff_profile) is PCFFOperationalProfile,
            "Operational profile belongs to PCFF only",
        )
        pcff_profile.validate_native(system, native)
        _require(
            description["charge_method"]
            == pcff_profile.payload["authorized_charge_policy"],
            "Operational charge policy mismatch",
        )
        description.update(
            schema="island_prepared_forcefield_v2",
            pcff_operational_profile=pcff_profile.payload,
        )
    return description


@dataclass(frozen=True)
class PreparedForceField:
    """Owned facade. Native access returns copies; identity is not a native hash.

    The original preparation system is retained because Amber preparation evidence
    binds its original coordinates. Evaluation permits coordinate replacement via
    the native imported parameter contract; historical preparation is never edited.
    """

    _system: object
    _family: str
    _native: object
    _source: object
    _json_text: str
    _pcff_profile: object | None = None

    def __post_init__(self):
        for key in ("_system", "_native", "_source", "_pcff_profile"):
            object.__setattr__(self, key, deepcopy(getattr(self, key)))
        self.validate_integrity()

    @property
    def metadata(self):
        self.validate_integrity()
        return json.loads(self._json_text)

    @property
    def native_result(self):
        self.validate_integrity()
        return deepcopy(self._native)

    @property
    def source(self):
        self.validate_integrity()
        return deepcopy(self._source)

    @property
    def operational_profile(self):
        self.validate_integrity()
        return deepcopy(self._pcff_profile)

    @property
    def identity(self):
        self.validate_integrity()
        return sha256(self._json_text.encode()).hexdigest()

    def validate_integrity(self, system=None):
        try:
            expected = _description(
                self._system,
                self._family,
                self._native,
                self._source,
                self._pcff_profile,
            )
            _require(
                self._json_text == _json(expected), "Contradictory prepared metadata"
            )
            if system is not None:
                _require(
                    _binding(system) == _binding(self._system),
                    "Prepared chemical binding mismatch",
                )
                if self._family in ("gaff", "gaff2"):
                    self._native.imported_result.validate_integrity(system)
                elif self._family == "oplsaa":
                    self._native.validate_integrity(system, self._source)
                else:
                    self._native.validate_integrity(system)
        except PreparedForceFieldError:
            raise
        except Exception as exc:
            raise PreparedForceFieldError(
                f"Native preparation validation failed: {exc}"
            ) from exc


def adopt_forcefield(system, family, native_result, *, source=None, pcff_profile=None):
    """Adopt validated native results without running tools or parameterization.

    Amber: full preparation record and original input system. OPLS: resolved result
    and pinned source. PCFF: model specification retaining its H3 assignment/source.
    """
    try:
        metadata = _description(system, family, native_result, source, pcff_profile)
        return PreparedForceField(
            system, family, native_result, source, _json(metadata), pcff_profile
        )
    except PreparedForceFieldError:
        raise
    except Exception as exc:
        raise PreparedForceFieldError(f"Cannot adopt native result: {exc}") from exc


def prepare_forcefield(system, request):
    """Execute exactly the selected native backend; native failure types propagate."""
    _require(
        type(request) is ForceFieldRequest,
        "ForceFieldRequest required",
        ForceFieldRequestError,
    )
    request = ForceFieldRequest(
        request.family,
        request.options,
        request.final_graph,
        request.typed_graph,
        request.graph_charges,
    )
    if request.final_graph is not None:
        request.final_graph.validate_integrity(system)
        if request.final_graph.periodic:
            raise PreparedForceFieldError(
                f"{request.family} preparation does not support periodic final graphs; "
                "retain the neutral graph bundle and use a periodic-capable backend"
            )
    options = request.options
    if request.typed_graph is not None and request.family != "pcff":
        raise PreparedForceFieldError(
            f"{request.family} unified assignment adapter is validation-only; "
            "its existing backend does not accept externally assigned types"
        )
    if request.family in ("gaff", "gaff2"):
        from .ambertools import AmberToolsParameterizationEngine

        result = AmberToolsParameterizationEngine().parameterize(system, options)
        _require(
            result.record["charge_method"] == options.charge_method,
            "Native charge method differs from request",
        )
        return adopt_forcefield(system, request.family, result)
    if request.family == "oplsaa":
        from .oplsaa import load_oplsaa_source, parameterize_oplsaa

        source = load_oplsaa_source(options.source_path)
        return adopt_forcefield(
            system, "oplsaa", parameterize_oplsaa(system, source), source=source
        )
    from .pcff import (
        assign_automatic_pcff_charges,
        assign_pcff_parameters,
        define_pcff_model,
        load_pcff_source,
        special_pair_policy,
        type_pcff_atoms,
    )

    profile = None
    if options.source_profile is not None:
        profile = options.source_profile.profile()
        profile.validate_system(system)
    source = load_pcff_source(options.source_path)
    if request.typed_graph is not None:
        from .pcff import (
            assign_typed_pcff_charges,
            bind_pcff_types,
            provide_pcff_charges,
        )

        unified_typed = request.typed_graph
        unified_charges = request.graph_charges
        unified_typed.validate_integrity(request.final_graph)
        unified_charges.validate_integrity(request.final_graph, unified_typed)
        _require(
            unified_typed.payload["source"]["sha256"] == source.identity["sha256"]
            and unified_charges.payload["source"]["sha256"] == source.identity["sha256"],
            "Unified PCFF source hash differs from loaded FRC",
        )
        typing = bind_pcff_types(
            system,
            source,
            unified_typed.atom_types,
            provenance=unified_typed.payload["provenance"],
            evidence_references=unified_typed.payload["evidence"],
        )
        charge_payload = unified_charges.payload
        if charge_payload["charge_origin"] == "provided":
            charges = provide_pcff_charges(
                system,
                typing,
                unified_charges.charges,
                unit=charge_payload["unit"],
                provenance=charge_payload["provenance"],
                evidence_references=charge_payload["evidence"],
                component_totals={
                    int(k): v for k, v in charge_payload["component_totals"].items()
                },
                total_charge=charge_payload["total_charge"],
                resolution_policy=options.resolution_policy,
            )
        else:
            charges = assign_typed_pcff_charges(
                system, typing, resolution_policy=options.resolution_policy
            )
    elif options.typed_graph is not None:
        typing, charges = options.typed_graph, options.graph_charges
        typing.validate_integrity(system)
        charges.validate_integrity(system)
        _require(
            typing.source.identity == source.identity,
            "External record/request source mismatch",
        )
    else:
        typing = type_pcff_atoms(system, source, profile=options.typing_profile)
        charges = assign_automatic_pcff_charges(
            system, typing, resolution_policy=options.resolution_policy
        )
    parameters = assign_pcff_parameters(
        system, typing, charges, resolution_policy=options.resolution_policy
    )
    model = define_pcff_model(
        parameters,
        special_pairs=special_pair_policy(lj=options.lj, coulomb=options.coulomb),
    )
    return adopt_forcefield(system, "pcff", model, pcff_profile=profile)


def create_evaluator(system, prepared):
    """Return the native Reference bound evaluator, including its session API."""
    _require(type(prepared) is PreparedForceField, "PreparedForceField required")
    prepared.validate_integrity(system)
    if prepared._family in ("gaff", "gaff2"):
        from island.evaluation import OpenMMSinglePointEvaluator

        return OpenMMSinglePointEvaluator(system, prepared._native.imported_result)
    if prepared._family == "oplsaa":
        from island.evaluation.oplsaa import OPLSSinglePointEvaluator

        return OPLSSinglePointEvaluator(system, prepared._native, prepared._source)
    from island.evaluation import PCFFSinglePointEvaluator

    return PCFFSinglePointEvaluator(system, prepared._native)
