"""One explicit selection boundary. No numerical engines or source discovery."""

import json
from copy import deepcopy
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path

from island.exceptions import IslandError
from island.forcefields.external_pdms import LIMITS

from .completion import MODES, CompletionRecord
from .inventory import canonical, digest, final_graph, key, requirements
from .registry import PDMS, USER_OPLS, inspect_profile


def _resolver_identity():
    root = Path(__file__).parent
    return digest(
        {p.name: sha256(p.read_bytes()).hexdigest() for p in sorted(root.glob("*.py"))}
    )


class UnifiedForceFieldError(IslandError):
    def __init__(self, status, diagnostics, *, result=None):
        self.status, self.diagnostics, self.result = (
            status,
            deepcopy(diagnostics),
            result,
        )
        super().__init__(status + ": " + str(diagnostics))


class _Stop(Exception):
    def __init__(self, status, stage, **details):
        self.status, self.diagnostic = status, {"stage": stage, **details}


def stop(condition, status, stage, **details):
    if not condition:
        raise _Stop(status, stage, **details)


@dataclass(frozen=True)
class SourceSelection:
    """An explicit profile and its existing backend options; never a file search."""

    profile_identity: str
    options: object = None


@dataclass(frozen=True)
class RESPSource:
    """Explicit independently reviewed RESP/Psi4 evidence, never a renamed vector."""

    directory: object
    expected_identity: str


@dataclass(frozen=True)
class UnifiedForceFieldPreparation:
    """Owned inspection report and a separately validated backend capability.

    JSON alone never reconstructs execution authority. Repeat preparation with
    explicit source assets after relocation; compare report identities.
    """

    json_text: str
    _original: str = field(repr=False)
    _prepared: object = field(default=None, repr=False)
    _files: tuple = field(default=(), repr=False)

    @property
    def payload(self):
        self.validate_integrity()
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]

    def validate_integrity(self, system=None):
        try:
            p = json.loads(self.json_text)
            stop(
                self.json_text == self._original,
                "source_incompatible",
                "report_mutation",
            )
            stop(
                p["identity"]
                == digest({k: v for k, v in p.items() if k != "identity"}),
                "source_incompatible",
                "report_checksum",
            )
            stop(
                p["resolver_identity"] == _resolver_identity(),
                "source_incompatible",
                "resolver_code_mutation",
            )
            for name, pin in self._files:
                stop(
                    sha256(Path(name).read_bytes()).hexdigest() == pin,
                    "source_incompatible",
                    "source_mutation",
                    file=name,
                )
            if system is not None:
                current = final_graph(system) if p["final_graph"] is not None else None
                stop(
                    canonical(current) == canonical(p["final_graph"])
                    and (
                        current is not None
                        or _raw_binding(system) == p["input_binding"]
                    ),
                    "source_incompatible",
                    "stale_final_graph",
                )
                if self._prepared is not None:
                    self._prepared.validate_integrity(system)
            stop(
                not p["stages"]["executable"] or self._prepared is not None,
                "source_incompatible",
                "missing_backend_authority",
            )
            if p["stages"]["executable"]:
                stop(
                    p["status"] == "prepared"
                    and all(
                        p["stages"][s]
                        for s in (
                            "typing_complete",
                            "charge_complete",
                            "parameter_complete",
                        )
                    )
                    and not p["blocked_terms"]
                    and all(
                        r["status"] in ("assigned", "not_applicable")
                        for r in p["requirements"]
                    )
                    and p["backend_identity"] == _backend_identity(self._prepared),
                    "source_incompatible",
                    "inconsistent_execution_authority",
                )
        except _Stop as error:
            raise UnifiedForceFieldError(
                error.status, [error.diagnostic], result=self
            ) from error
        except UnifiedForceFieldError:
            raise
        except Exception as error:
            raise UnifiedForceFieldError(
                "source_incompatible",
                [{"stage": "integrity_validation", "reason": str(error)}],
                result=self,
            ) from error

    def require_executable(self, system):
        self.validate_integrity(system)
        p = self.payload
        if not p["stages"]["executable"]:
            raise UnifiedForceFieldError(p["status"], p["diagnostics"], result=self)
        if hasattr(self._prepared, "require_executable"):
            self._prepared.require_executable(system)
        return self._prepared


def _backend_identity(backend):
    if hasattr(backend, "identity"):
        return backend.identity
    return backend.payload["identity"]


def _raw_binding(system):
    try:
        data = system.to_dict()
        data.pop("coordinates", None)
        return digest(data)
    except Exception:  # noqa: BLE001 -- even an invalid graph needs a classified report
        return None


def _failure_status(error, stage):
    from island.exceptions import AtomTypingError, ChargeAssignmentError

    causes, cursor = [], error
    while cursor is not None and cursor not in causes:
        causes.append(cursor)
        cursor = cursor.__cause__ or cursor.__context__
    if stage == "final_graph" or any(isinstance(e, AtomTypingError) for e in causes):
        return "typing_ambiguous"
    if any(isinstance(e, ChargeAssignmentError) for e in causes):
        return "charge_incomplete"
    if stage == "request":
        return "source_incompatible"
    if any(
        isinstance(e, (OSError, ImportError))
        or "Unavailable" in type(e).__name__
        or "Dependency" in type(e).__name__
        for e in causes
    ):
        return "asset_unavailable"
    return "parameter_incomplete"


def _set(p, kind, sites, value, **extra):
    ident = key(kind, sites)
    row = next((r for r in p["requirements"] if r["id"] == ident), None)
    if row is None:
        row = {"id": ident, "family": kind, "sites": list(sites), "dependencies": []}
        p["requirements"].append(row)
    row.update(status="assigned", value=value, **extra)
    row.pop("reason", None)


def _pin(p, files, path, expected=None):
    stop(path is not None, "asset_unavailable", "explicit_source_path_required")
    path = Path(path)
    stop(path.is_file(), "asset_unavailable", "source_file_unavailable", file=str(path))
    pin = sha256(path.read_bytes()).hexdigest()
    if expected is not None:
        stop(
            pin == expected,
            "source_incompatible",
            "source_hash_mismatch",
            actual=pin,
            expected=expected,
        )
    files.append((str(path), pin))
    p.setdefault("source_hashes", []).append(pin)
    return pin


def _atoms(p, system, types, charges=None, masses=None):
    for i, label in sorted(types.items()):
        _set(p, "atom_type", [i], label)
        _set(
            p,
            "mass",
            [i],
            masses[i] if masses else system.topology.sites[i].mass,
            origin="source mass" if masses else "validated final-graph mass",
        )
        if charges is not None:
            _set(
                p,
                "charge",
                [i],
                charges[i],
                unit="elementary_charge",
                origin=p["charge_method"],
            )


def _policies(p, policy):
    for kind in ("mixing_rule", "scaling_14", "exclusions", "electrostatics"):
        _set(p, kind, [], policy, source="selected backend source policy")


def _pcff(system, p, options, files):
    from island.forcefields.pcff import (
        assign_automatic_pcff_charges,
        assign_pcff_parameters,
        define_pcff_model,
        load_pcff_source,
        special_pair_policy,
        type_pcff_atoms,
    )
    from island.forcefields.preparation import PCFFOptions, adopt_forcefield

    stop(
        type(options) is PCFFOptions,
        "source_incompatible",
        "explicit_pcff_options_required",
    )
    stop(
        options.typing_profile == p["profile_identity"],
        "source_incompatible",
        "typing_profile_mismatch",
    )
    stop(
        options.resolution_policy is None and options.source_profile is None,
        "source_incompatible",
        "strict_source_rejects_legacy_fallback",
    )
    stop(
        options.typed_graph is None and options.graph_charges is None,
        "source_incompatible",
        "explicit_graph_records_require_separate_adapter",
    )
    _pin(p, files, options.source_path)
    source = load_pcff_source(options.source_path)
    typing = type_pcff_atoms(system, source, profile=options.typing_profile)
    t = typing.payload
    p["typing_evidence"] = t
    p["typed_graph_identity"] = typing.identity
    p["typing_evidence"] = typing.payload
    _atoms(p, system, t["assignments"])
    p["stages"]["typing_complete"] = typing.complete
    stop(
        typing.complete,
        "typing_ambiguous",
        "native_pcff_typing",
        decisions=t["diagnostics"],
    )
    charges = assign_automatic_pcff_charges(system, typing)
    p["charge_evidence"] = charges.payload
    p["charge_identity"] = charges.identity
    p["stages"]["charge_complete"] = charges.complete
    stop(
        charges.complete,
        "charge_incomplete",
        "native_pcff_charge_increments",
        record=charges.payload["native_charge_record"],
    )
    _atoms(p, system, typing.assignments, charges.charges)
    parameters = assign_pcff_parameters(system, typing, charges)
    for a in parameters.payload["assignments"]:
        ident = key(a["family"], a["sites"])
        row = next((r for r in p["requirements"] if r["id"] == ident), None)
        if row is None:
            row = {"id": ident, "family": a["family"], "sites": a["sites"]}
            p["requirements"].append(row)
        row.update(deepcopy(a))
        if row["status"] in ("assigned", "not_applicable"):
            row.pop("reason", None)
    p["parameter_identity"] = parameters.identity
    _policies(
        p,
        {
            "mixing": "Class-II sixth-power",
            "special_pairs": special_pair_policy(
                lj=options.lj, coulomb=options.coulomb
            ),
            "electrostatics": "finite NoCutoff",
        },
    )
    stop(
        parameters.payload["parameter_coverage_complete"],
        "parameter_incomplete",
        "class2_missing_rows_or_equilibrium_dependencies",
    )
    model = define_pcff_model(
        parameters,
        special_pairs=special_pair_policy(lj=options.lj, coulomb=options.coulomb),
    )
    return adopt_forcefield(system, "pcff", model)


def _opls(system, p, options, files):
    from xml.etree import ElementTree as ET

    from island.forcefields.oplsaa import (
        assign_native_charges,
        load_oplsaa_source,
        parameterize_oplsaa,
        type_atoms,
    )
    from island.forcefields.oplsaa.parameters import POLICY, inventory, selected_record
    from island.forcefields.preparation import OPLSOptions, adopt_forcefield

    stop(
        type(options) is OPLSOptions,
        "source_incompatible",
        "explicit_opls_options_required",
    )
    stop(
        options.typed_graph is None and options.graph_charges is None,
        "source_incompatible",
        "native_opls_external_records_not_supported",
    )
    _pin(p, files, options.source_path)
    source = load_oplsaa_source(options.source_path)
    try:
        typing = type_atoms(system, source)
    except Exception as error:
        raise _Stop(
            "typing_ambiguous", "native_opls_typing", reason=str(error)
        ) from error
    p["typed_graph_identity"] = typing.identity
    types = typing.payload["data"]["assignments"]
    _atoms(p, system, types)
    p["stages"]["typing_complete"] = True
    try:
        charge = assign_native_charges(system, typing, source)
    except Exception as error:
        raise _Stop(
            "charge_incomplete", "native_opls_charges", reason=str(error)
        ) from error
    p["charge_identity"] = charge.identity
    p["charge_evidence"] = charge.payload
    p["stages"]["charge_complete"] = True
    _atoms(
        p, system, types, {i: r["charge"] for i, r in charge.payload["charges"].items()}
    )
    paths, pairs = inventory(system)
    root = ET.fromstring(source.xml)
    definitions, _ = source.entries
    # Inspect every required path, rather than reporting only the first miss.
    for kind, paths_for_kind in paths.items():
        for sites in paths_for_kind:
            family = {"bonds": "bond", "angles": "angle", "rb": "proper"}[kind]
            labels = tuple(types[i] for i in sites)
            try:
                index, attrs, candidates = selected_record(
                    root, definitions, labels, kind
                )
                _set(
                    p,
                    family,
                    sites,
                    attrs,
                    labels=labels,
                    source_record_index=index,
                    matching_records=candidates,
                )
            except Exception as error:  # noqa: BLE001 -- classified public boundary
                row = next(
                    r for r in p["requirements"] if r["id"] == key(family, sites)
                )
                row.update(labels=labels, reason=str(error))
    for i, label in sorted(types.items()):
        rows = [
            r.attrib
            for r in root.find("NonbondedForce")
            if r.attrib.get("type") == label
        ]
        if len(rows) == 1:
            _set(p, "nonbonded", [i], rows[0])
    for pair in pairs:
        _set(p, "special_pair", pair["sites"], pair)
    _policies(p, POLICY)
    stop(
        all(r["status"] == "assigned" for r in p["requirements"]),
        "parameter_incomplete",
        "opls_missing_rows_or_improper_policy",
    )
    native = parameterize_oplsaa(system, source)
    p["parameter_identity"] = native.identity
    p["parameter_evidence"] = native.payload
    return adopt_forcefield(system, "oplsaa", native, source=source)


def _external(system, p, options, files):
    if p["profile_identity"] == PDMS:
        from island.forcefields.pdms_opls import PDMSOPLSOptions, prepare_pdms_opls

        stop(
            type(options) is PDMSOPLSOptions,
            "asset_unavailable",
            "explicit_pdms_archive_options_required",
        )
        stop(
            options.expected_sha256 is not None,
            "asset_unavailable",
            "independent_archive_sha256_required",
        )
        _pin(p, files, options.archive_path, options.expected_sha256)
        result = prepare_pdms_opls(system, options)
        data = result.payload
        stop(
            data["status"] == "prepared",
            {"typing_blocked": "typing_ambiguous"}.get(data["status"], data["status"]),
            "external_pdms_preparation",
            diagnostics=data["diagnostics"],
        )
        result.require_executable(system)
        mapping = data["mapping"]
        p["typed_graph_identity"] = digest(mapping)
        p["charge_identity"] = digest(data["charge_policy"] | {"mapping": mapping})
        p["parameter_identity"] = data["identity"]
        p["source_evidence"] = data["source"]
        p["charge_evidence"] = data["charge_policy"]
        _atoms(
            p,
            system,
            {r["site_id"]: r["author_type_id"] for r in mapping},
            {r["site_id"]: r["charge_e"] for r in mapping},
            {r["site_id"]: r["mass_dalton"] for r in mapping},
        )
        for family, rows in data["interactions"].items():
            for row in rows:
                _set(
                    p,
                    {"proper_torsion": "proper"}.get(family, family),
                    row["sites"],
                    row,
                )
        from island.forcefields.oplsaa.parameters import inventory

        for pair in inventory(system)[1]:
            _set(
                p,
                "special_pair",
                pair["sites"],
                {"source_policy": data["model"], "sites": pair["sites"]},
            )
        _policies(p, data["model"])
    else:
        from island.forcefields.opls_external import (
            ExternalOPLSOptions,
            prepare_external_opls,
        )

        stop(
            type(options) is ExternalOPLSOptions,
            "asset_unavailable",
            "explicit_external_opls_options_required",
        )
        stop(
            options.expected_sha256 is not None,
            "asset_unavailable",
            "independent_manifest_hash_required",
        )
        _pin(p, files, options.source_path, options.expected_sha256)
        result = prepare_external_opls(system, options)
        document = result.payload
        stop(
            document["status"] == "prepared",
            {"typing_blocked": "typing_ambiguous"}.get(
                document["status"], document["status"]
            ),
            "external_opls_preparation",
            diagnostics=document["diagnostics"],
        )
        result.require_executable(system)
        p["external_assignment"] = document
        p["parameter_identity"] = result.identity
        data = document["assignment"]
        # Pin every companion asset too; a manifest alone is not a live source.
        from island.forcefields.opls_external import load_source

        source = load_source(options.source_path, options.expected_sha256)
        d = source.document
        _pin(p, files, source.root_path / d["xml_file"], d["xml_sha256"])
        proof = d["provenance"]["license_evidence"]
        if proof:
            _pin(p, files, source.root_path / proof["file"], proof["sha256"])
        sites = data["sites"]
        _atoms(
            p,
            system,
            {i: r["type"] for i, r in sites.items()},
            {i: r["exact_charge_e"] for i, r in sites.items()},
            {i: r["mass_dalton"] for i, r in sites.items()},
        )
        for i, r in sites.items():
            _set(p, "nonbonded", [i], r)
        for name, family in (("bonds", "bond"), ("angles", "angle"), ("rb", "proper")):
            for row in data[name]:
                _set(p, family, row["sites"], row)
        for row in data["pairs"]:
            _set(p, "special_pair", row["sites"], row)
        _policies(
            p,
            {
                "model": data["model"],
                "special_pairs": data["special_pairs"],
                "mixing": "geometric_sigma_and_epsilon",
                "electrostatics": "NoCutoff",
            },
        )
        p["typed_graph_identity"] = digest(data["typing_matches"])
        p["charge_identity"] = digest(
            {i: r["exact_charge_e"] for i, r in sites.items()}
        )
    p["stages"].update(typing_complete=True, charge_complete=True)
    return result


def _amber(system, p, options, files, resp_source):
    from island.forcefields.ambertools import AmberToolsOptions
    from island.forcefields.preparation import ForceFieldRequest, prepare_forcefield

    stop(
        type(options) is AmberToolsOptions and options.force_field == p["forcefield"],
        "source_incompatible",
        "amber_family_options_mismatch",
    )
    expected = (
        "provided"
        if p["charge_method"] in ("provided", "external_provided", "resp_psi4")
        else p["charge_method"]
    )
    stop(
        expected in ("provided", "am1bcc") and expected == options.charge_method,
        "charge_incomplete",
        "explicit_verified_charge_bridge_required",
        requested=p["charge_method"],
    )
    if expected == "provided":
        stop(
            bool(options.charge_source),
            "charge_incomplete",
            "provided_charge_origin_required",
        )
        p["charge_origin"] = {
            "method": p["charge_method"],
            "source": options.charge_source,
            "verified_method": False,
        }
        if p["charge_method"] == "resp_psi4":
            from island.forcefields.ambertools.resp import RESPCharges
            from island.graph.final import final_graph as amber_graph

            stop(
                type(resp_source) is RESPSource,
                "charge_incomplete",
                "independent_resp_evidence_required",
            )
            graph = amber_graph(system)
            record = RESPCharges.load(
                resp_source.directory,
                graph,
                expected_identity=resp_source.expected_identity,
            )
            bridge = record.provided(
                resp_source.directory,
                graph,
                expected_identity=resp_source.expected_identity,
            )
            exact = bridge.validate(
                graph,
                expected_reference_identity=record.payload["reference"]["identity"],
            )
            stop(
                dict(options.provided_charges) == dict(exact),
                "charge_incomplete",
                "resp_charge_vector_mismatch",
            )
            for asset in record.payload["assets"].values():
                _pin(
                    p,
                    files,
                    Path(resp_source.directory) / asset["path"],
                    asset["sha256"],
                )
            _pin(p, files, Path(resp_source.directory) / "resp-result.json")
            p["charge_origin"].update(
                verified_method=True, evidence_identity=record.payload["identity"]
            )
        else:
            stop(
                "resp" not in options.charge_source.lower(),
                "charge_incomplete",
                "resp_label_requires_verified_resp_bridge",
            )
    else:
        p["charge_origin"] = {
            "method": "AM1-BCC",
            "source": "native AmberTools SQM/BCC",
        }
    stop(
        options.amberhome is not None,
        "asset_unavailable",
        "explicit_amberhome_required",
    )
    home = Path(options.amberhome)
    _pin(p, files, home / "dat/leap/parm" / (p["forcefield"] + ".dat"))
    _pin(p, files, home / "dat/leap/cmd" / ("leaprc." + p["forcefield"]))
    native = prepare_forcefield(system, ForceFieldRequest(p["forcefield"], options))
    result = native.native_result.imported_result
    p["typed_graph_identity"] = digest(dict(result.atom_types))
    p["charge_identity"] = digest(
        {i: r.charge for i, r in result.charge_result.assignments.items()}
    )
    p["parameter_identity"] = result.result_signature
    _atoms(
        p,
        system,
        dict(result.atom_types),
        {i: r.charge for i, r in result.charge_result.assignments.items()},
    )
    for kind, attr in (
        ("bond", "bond_assignments"),
        ("angle", "angle_assignments"),
        ("proper", "proper_torsion_assignments"),
        ("improper", "improper_assignments"),
        ("nonbonded", "site_assignments"),
    ):
        from dataclasses import asdict

        for sites, selection in getattr(result, attr).items():
            _set(p, kind, [sites] if type(sites) is int else sites, asdict(selection))
    for r in p["requirements"]:
        if r["family"] == "special_pair":
            _set(
                p,
                "special_pair",
                r["sites"],
                {
                    "source_exclusions": result.source_exclusions,
                    "source_14_pairs": result.source_14_pairs,
                },
            )
    _policies(
        p, {"origin": "exact Amber import", "record_signature": result.result_signature}
    )
    p["stages"].update(typing_complete=True, charge_complete=True)
    return native


def prepare_polymer_forcefield(
    system,
    forcefield,
    *,
    source_profile=None,
    charge_method=None,
    resolution_policy="strict_source",
    options=None,
    completion_records=(),
    completion_root=None,
    resp_source=None,
):
    """Prepare a selected final graph; return classified reports for all failures.

    Existing typed options carry explicit sources and conventions. Completion
    evidence never substitutes for missing native rows or confers MD authority.
    """
    p = {
        "schema": "island_unified_forcefield_report_v1",
        "resolver": "island_unified_forcefield_selection_v1",
        "resolver_identity": _resolver_identity(),
        "forcefield": forcefield,
        "profile_identity": None,
        "charge_method": charge_method,
        "resolution_policy": resolution_policy,
        "final_graph": None,
        "final_graph_identity": None,
        "typed_graph_identity": None,
        "charge_identity": None,
        "parameter_identity": None,
        "input_binding": _raw_binding(system),
        "requirements": [],
        "diagnostics": [],
        "completion_records": [],
        "stages": {
            k: False
            for k in (
                "typing_complete",
                "charge_complete",
                "parameter_complete",
                "executable",
                "numerically_verified",
                "scientifically_validated",
            )
        },
        "production_validated": False,
        "simulation_readiness": "not_established",
        "limitations": deepcopy(LIMITS),
    }
    prepared, files = None, []
    stage = "request"
    try:
        stop(
            forcefield in ("pcff", "opls", "gaff", "gaff2"),
            "source_incompatible",
            "unknown_forcefield",
        )
        stop(
            resolution_policy in MODES,
            "source_incompatible",
            "unknown_resolution_policy",
        )
        if type(source_profile) is SourceSelection:
            stop(options is None, "source_incompatible", "duplicate_source_options")
            options, source_profile = (
                source_profile.options,
                source_profile.profile_identity,
            )
        defaults = {
            "pcff": "island_pcff_acyclic_cho_v1",
            "opls": "island_oplsaa_foyer_v1",
            "gaff": "island_gaff_ambertools_v1",
            "gaff2": "island_gaff2_ambertools_v1",
        }
        name = source_profile if source_profile is not None else defaults[forcefield]
        definition = inspect_profile(name)
        p.update(profile_identity=name, source_kind=definition["source_kind"])
        stop(
            definition["family"] == forcefield,
            "source_incompatible",
            "profile_family_mismatch",
        )
        stop(
            definition.get("selection_adapter", True),
            "source_incompatible",
            "legacy_profile_requires_separate_policy_review",
        )
        stage = "final_graph"
        graph = final_graph(system)
        p.update(
            final_graph=graph, final_graph_identity=digest(graph), input_binding=None
        )
        p["requirements"] = [
            {
                **r,
                "status": "missing",
                "reason": "unresolved source/type/charge dependency",
            }
            for r in requirements(system, forcefield)
        ]
        allowed = definition["supported_elements"]
        unsupported = [
            a["id"] for a in graph["sites"] if allowed and a["element"] not in allowed
        ]
        stop(
            not unsupported,
            "typing_ambiguous",
            "unsupported_atom_environments",
            atom_ids=unsupported,
        )
        stop(
            charge_method in definition["charge_methods"],
            "charge_incomplete",
            "explicit_compatible_charge_method_required",
            allowed=definition["charge_methods"],
        )
        if completion_records:
            stop(
                resolution_policy != "strict_source" and completion_root is not None,
                "source_incompatible",
                "explicit_completion_policy_and_root_required",
            )
            for record in completion_records:
                stop(
                    type(record) is CompletionRecord,
                    "source_incompatible",
                    "completion_record_type",
                )
                data = record.validate_integrity(completion_root)
                stop(
                    data["source_identity"] == name,
                    "source_incompatible",
                    "completion_source_mismatch",
                )
                kind = (
                    "derived"
                    if resolution_policy == "derived_from_declared_formula"
                    else "fitted_external"
                )
                stop(
                    data["kind"] == kind,
                    "source_incompatible",
                    "completion_policy_mismatch",
                )
                p["completion_records"].append(data)
                entries = data.get("inputs", [data.get("reference")])
                for entry in entries:
                    _pin(
                        p, files, Path(completion_root) / entry["file"], entry["sha256"]
                    )
        elif resolution_policy != "strict_source":
            raise _Stop(
                "parameter_incomplete",
                "completion_evidence_unavailable",
                policy=resolution_policy,
            )
        stage = "source"
        stop(
            options is not None, "asset_unavailable", "explicit_source_options_required"
        )
        if name in (PDMS, USER_OPLS):
            prepared = _external(system, p, options, files)
        elif forcefield == "pcff":
            prepared = _pcff(system, p, options, files)
        elif forcefield == "opls":
            prepared = _opls(system, p, options, files)
        else:
            prepared = _amber(system, p, options, files, resp_source)
        missing = [
            r["id"]
            for r in p["requirements"]
            if r["status"] not in ("assigned", "not_applicable")
        ]
        stop(
            not missing,
            "parameter_incomplete",
            "unresolved_required_interactions",
            missing=missing,
        )
        p["stages"].update(parameter_complete=True, executable=True)
        p["status"] = "prepared"
        p["backend_identity"] = _backend_identity(prepared)
    except _Stop as error:
        p.update(status=error.status, diagnostics=[error.diagnostic])
        prepared = None
    except Exception as error:  # noqa: BLE001 -- classified public boundary
        status = _failure_status(error, stage)
        p.update(
            status=status,
            diagnostics=[
                {
                    "stage": stage,
                    "exception": type(error).__name__,
                    "reason": str(error),
                }
            ],
        )
        prepared = None
    p["requirements"].sort(key=lambda r: r["id"])
    p["source_identity"] = digest(
        {
            "profile_identity": p["profile_identity"],
            "sha256": sorted(set(p.get("source_hashes", []))),
        }
    )
    p["assigned_terms"] = [
        r["id"] for r in p["requirements"] if r["status"] == "assigned"
    ]
    p["blocked_terms"] = [
        r["id"]
        for r in p["requirements"]
        if r["status"] not in ("assigned", "not_applicable")
    ]
    p["identity"] = digest(p)
    raw = canonical(p)
    return UnifiedForceFieldPreparation(raw, raw, prepared, tuple(files))


report_polymer_forcefield = prepare_polymer_forcefield


def create_evaluator(system, preparation):
    """Check the unified result before importing or constructing an engine."""
    if type(preparation) is not UnifiedForceFieldPreparation:
        raise UnifiedForceFieldError(
            "source_incompatible", [{"stage": "unified_preparation_required"}]
        )
    backend = preparation.require_executable(system)
    from island.forcefields.research import create_evaluator as existing

    return existing(system, backend)
