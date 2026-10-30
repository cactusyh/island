"""Parameter completeness, backend preparation and scientific claims are separate.

An additive J55 view: historical facades, hashes and decisions are not rewritten.
Ordinary family selection remains forcefield="opls". Numerical evidence is opt-in.
"""

import json
from dataclasses import dataclass, field
from hashlib import sha256
from pathlib import Path

from island.forcefields.opls_registry import (
    ForceFieldRequest,
    OPLSAssets,
    use_opls_assets,
)
from island.forcefields.opls_registry.api import require
from island.forcefields.opls_registry.registry import NATIVE
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import canonical, digest, requirements

from .claims import PDMS_MODEL_GAPS, SCHEMA, ClaimRequirements, assess


def _proof(value, file, **extra):
    return {
        "value": value,
        "evidence": [
            {
                "file": str(file),
                "sha256": sha256(Path(file).read_bytes()).hexdigest(),
                **extra,
            }
        ],
    }


def policy_identity():
    return digest(
        {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(__file__).parent.glob("*.py"))
        }
    )


def _amber_conventions(assignment):
    from dataclasses import asdict

    from island.forcefields.amber import importer, models

    prepared = assignment._prepared
    result = prepared.native_result.imported_result
    a = assignment.payload
    values = {
        "chemical_domain": {
            "graph": a["final_graph_identity"],
            "typing": result.atom_types,
        },
        "charge_origin": a.get("charge_origin", a.get("charge_evidence")),
        "units": "imported explicit parameter units: nm, kJ/mol, radian, e",
        "functional_forms": "harmonic bonds/angles; source periodic proper/improper; LJ12-6; Coulomb",
        "torsion_convention": "source periodic K*(1+cos(n*phi-phase)); radians; SCEE/SCNB retained",
        "exclusions": {
            "source_exclusions": result.source_exclusions,
            "pairs_14": result.source_14_pairs,
        },
        "coulomb_14": {
            "value": result.nonbonded_policy.coulomb_scale_14,
            "source": result.source_sha256,
        },
        "lj_14": {
            "value": result.nonbonded_policy.lj_scale_14,
            "source": result.source_sha256,
        },
        "lj_mixing": asdict(result.nonbonded_policy),
        "cutoff": "none: finite NoCutoff",
        "long_range": "none: finite imported evaluator; periodic execution not claimed",
        "provenance": result.provenance,
        "source_hashes": result.source_sha256,
        "evidence": {
            "result_signature": result.result_signature,
            "native_record": prepared.native_result.record,
        },
    }
    evidence = [
        {
            "file": "island/forcefields/amber/" + Path(m.__file__).name,
            "sha256": sha256(Path(m.__file__).read_bytes()).hexdigest(),
            "source_sha256": result.source_sha256,
        }
        for m in (models, importer)
    ]
    return {k: {"value": v, "evidence": evidence} for k, v in values.items()}


def _external_conventions(assignment):
    from island.forcefields.opls_external import source

    a = assignment.payload
    external = a["external_assignment"]
    data = external["assignment"]
    values = {
        "chemical_domain": data["model"],
        "charge_origin": a["charge_method"],
        "units": source.UNITS,
        "functional_forms": {"model": data["model"], "torsion": data["torsion_form"]},
        "torsion_convention": {
            "form": data["torsion_form"],
            "source_rows": data["source_row_evidence"],
        },
        "exclusions": data["special_pairs"],
        "coulomb_14": data["special_pairs"],
        "lj_14": data["special_pairs"],
        "lj_mixing": "geometric_sigma_epsilon",
        "cutoff": "none: finite NoCutoff",
        "long_range": "none: finite NoCutoff",
        "provenance": external,
        "source_hashes": a["source_hashes"],
        "evidence": data["source_row_evidence"],
    }
    evidence = {
        "file": "island/forcefields/opls_external/source.py",
        "sha256": sha256(Path(source.__file__).read_bytes()).hexdigest(),
        "external_assignment_identity": external["identity"],
    }
    return {k: {"value": v, "evidence": [evidence]} for k, v in values.items()}


def _pcff_conventions(assignment):
    from island.forcefields.pcff import model

    a = assignment.payload
    m = assignment._prepared.native_result.payload
    profile = m["compatibility_profile"]
    special = m["special_pairs"]
    values = {
        "chemical_domain": {
            "typed_graph": a["typed_graph_identity"],
            "model": profile["name"],
        },
        "charge_origin": a["charge_evidence"],
        "units": profile["units"],
        "functional_forms": profile["equations"],
        "torsion_convention": {
            "equation": profile["equations"]["torsion_3"],
            "geometry": profile["torsion_geometry"],
            "source_discrepancies": profile["frc_discrepancies"],
        },
        "exclusions": special,
        "coulomb_14": {"value": special["coulomb"][2], "origin": special["origin"]},
        "lj_14": {"value": special["lj"][2], "origin": special["origin"]},
        "lj_mixing": profile["mixing"],
        "cutoff": profile["boundary"],
        "long_range": profile["boundary"],
        "provenance": m["source"],
        "source_hashes": a["source_hashes"],
        "evidence": a["requirements"],
    }
    proof = {
        "file": "island/forcefields/pcff/model.py",
        "sha256": sha256(Path(model.__file__).read_bytes()).hexdigest(),
        "model_identity": m["identity"],
        "authority": "existing explicit converter model; not recovered source facts",
    }
    return {k: {"value": v, "evidence": [proof]} for k, v in values.items()}


def _pdms(system, assignment, sources):
    from island.opls_resolution import prepare_polymer_forcefield as resolve

    audit = resolve(
        system, forcefield="opls", assets=assignment.assets, sources=sources
    )
    return _pdms_conventions(assignment, audit), audit


def _pdms_conventions(assignment, audit):
    p = audit.payload
    settings = {r["setting"]: r for r in p["rows"]}
    conventions = {}
    names = {
        "units": "energy_units",
        "functional_forms": "lj_functional_form",
        "torsion_convention": "torsion_convention",
        "exclusions": "special_bonds",
        "coulomb_14": "coulomb_14_scaling",
        "lj_14": "lj_14_scaling",
        "lj_mixing": "lj_mixing",
        "cutoff": "pair_cutoff_angstrom",
        "long_range": "long_range_solver",
    }
    for name, setting in names.items():
        row = settings[setting]
        conventions[name] = {"value": row["value"], "evidence": [row]}
    a = assignment.payload
    for name, value in (
        ("chemical_domain", a["domain"]),
        ("charge_origin", a["graph_charges"]),
        ("provenance", a["source_evidence"]),
        ("source_hashes", a["source_evidence"]["source_hashes"]),
        ("evidence", a["source_evidence"]["rows"]),
    ):
        conventions[name] = {"value": value, "evidence": [a["source_identity"]]}
    for name in PDMS_MODEL_GAPS:
        row = settings[name]
        conventions[name] = {
            "value": row["value"],
            "required": True,
            "evidence": [row],
            "authority": row["resolution"],
        }
    # Preserve all old gaps and their scope; never turn backend choices into
    # source authority. Seven numerical/observable gaps do not gate row coverage.
    return conventions


def _native_conventions(assignment):
    from island.forcefields.oplsaa import parameters

    source = Path(assignment.assets.native_source_path)
    code = Path(parameters.__file__)
    p = assignment.payload
    values = {
        "chemical_domain": p["domain"],
        "charge_origin": p["graph_charges"],
        "units": parameters.UNITS,
        "functional_forms": "harmonic bonds/angles; RB torsion; LJ12-6; Coulomb",
        "torsion_convention": "sum C_n*cos(phi-pi)^n, n=0..5; kJ/mol",
        "exclusions": [1, 2],
        "coulomb_14": "0.5",
        "lj_14": "0.5",
        "lj_mixing": parameters.POLICY["mixing"],
        "cutoff": "none: finite NoCutoff",
        "long_range": "none: finite NoCutoff; no periodic claim",
        "provenance": p["source_identity"],
        "source_hashes": sha256(source.read_bytes()).hexdigest(),
        "evidence": p["assigned_terms"],
    }
    # Paths are evidence labels only. They are made relocation-independent.
    return {
        k: {
            **_proof(v, code),
            "evidence": [
                {
                    "file": "island/forcefields/oplsaa/parameters.py",
                    "sha256": sha256(code.read_bytes()).hexdigest(),
                    "source_xml_sha256": sha256(source.read_bytes()).hexdigest(),
                    "authority": "pinned source and documented converter",
                }
            ],
        }
        for k, v in values.items()
    }


def _pdms_rows(system, raw):
    mapping = {r["site_id"]: r for r in raw["mapping"]}
    interactions = {}
    for family, rows in raw["interactions"].items():
        for row in rows:
            family_name = {"proper_torsion": "proper"}.get(family, family)
            interactions[(family_name, tuple(row["sites"]))] = row
    result = []
    for req in requirements(system, "opls"):
        kind, sites = req["family"], req["sites"]
        value = None
        if kind in ("atom_type", "mass", "charge"):
            value = mapping[sites[0]]
        elif kind in (
            "mixing_rule",
            "scaling_14",
            "exclusions",
            "electrostatics",
            "special_pair",
        ):
            value = raw["model"]
        else:
            value = interactions.get((kind, tuple(sites)))
        result.append(
            {
                **req,
                "status": "assigned" if value is not None else "missing",
                "value": value,
                "evidence": raw["source"]["source_hashes"],
            }
        )
    return result


@dataclass(frozen=True)
class Preparation:
    json_text: str
    assignment: object = field(repr=False)
    _original: str = field(repr=False)
    _system_identity: str = field(repr=False)
    _audit: object = field(default=None, repr=False)
    _pins: tuple = field(default=(), repr=False)

    @property
    def payload(self):
        require(
            self.json_text == self._original,
            "source_incompatible",
            "policy_report_mutation",
        )
        p = json.loads(self.json_text)
        require(
            p["identity"] == digest({k: v for k, v in p.items() if k != "identity"}),
            "source_incompatible",
            "policy_report_checksum",
        )
        for path, expected in self._pins:
            require(
                Path(path).is_file()
                and sha256(Path(path).read_bytes()).hexdigest() == expected,
                "source_incompatible",
                "policy_source_or_evidence_mutation",
            )
        return p

    @property
    def identity(self):
        return self.payload["identity"]

    def validate_integrity(self, system):
        from island.workflows.bundle import system_data

        p = self.payload
        a = self.assignment.payload
        bound = {
            "selected_profile": a.get("selected_profile", a.get("profile_identity")),
            "source_identity": a.get("source_identity"),
            "source_kind": a.get("source_kind"),
            "graph_identity": a.get("graph_identity", a.get("final_graph_identity")),
            "charge_record": a.get("graph_charges", a.get("charge_evidence")),
        }
        require(
            all(p[k] == v for k, v in bound.items()),
            "source_incompatible",
            "parameter_source_profile_charge_binding",
        )
        require(
            p["execution_input_identity"] == self._system_identity,
            "source_incompatible",
            "execution_input_identity_mutation",
        )
        require(
            p["policy_identity"] == policy_identity(),
            "source_incompatible",
            "policy_code_mutation",
        )
        require(
            self.assignment.payload == p["source_assignment"]
            and self.assignment.payload["identity"] == p["assignment_identity"],
            "source_incompatible",
            "assignment_identity_mutation",
        )
        require(
            digest(system_data(system)) == self._system_identity,
            "source_incompatible",
            "policy_graph_charge_or_coordinate_mutation",
        )
        if p["source_integrity"] == "passed":
            backend = self.assignment.validate_integrity(system)
            if self._audit is not None:
                self._audit.validate_integrity(system)
                conventions = _pdms_conventions(self.assignment, self._audit)
                rows = _pdms_rows(system, backend.validate_integrity(system).payload)
            elif p["selected_profile"] == NATIVE:
                conventions = _native_conventions(self.assignment)
                rows = backend.payload["requirements"]
            elif p["family"] in ("gaff", "gaff2"):
                conventions = _amber_conventions(self.assignment)
                rows = self.assignment.payload["requirements"]
            elif p["family"] == "pcff":
                conventions = _pcff_conventions(self.assignment)
                rows = self.assignment.payload["requirements"]
            else:
                conventions = _external_conventions(self.assignment)
                rows = self.assignment.payload["requirements"]
            expected = assess(
                requirements(system, p["family"]),
                rows,
                conventions,
                typing_complete=p["source_assignment"]["stages"]["typing_complete"],
                charge_complete=p["source_assignment"]["stages"]["charge_complete"],
                source_valid=True,
            )
            require(
                conventions == p["conventions"]
                and all(expected[k] == p[k] for k in expected),
                "source_incompatible",
                "parameter_completeness_reconstruction",
            )
        return p

    def require(self, claims=None):
        (claims or ClaimRequirements()).enforce(self.payload)
        return self


def _signed(p, assignment, system, audit=None, pins=()):
    from island.workflows.bundle import system_data

    p = {**p, "identity": digest(p)}
    text = canonical(p)
    return Preparation(text, assignment, text, digest(system_data(system)), audit, pins)


def prepare_polymer_forcefield(
    system,
    forcefield="opls",
    *,
    request=None,
    assets=None,
    sources=None,
    claims=None,
    **legacy_options,
):
    """Resolve and audit parameters without running numerical diagnostics.

    Asset/source installation is an administrator task; profile/charge overrides
    are never needed for ordinary canonical OPLS selection. Legacy family options
    retain the original PCFF and AmberTools provenance and rejection boundaries.
    """
    from island.forcefields.opls_registry import prepare_forcefield as canonical_prepare
    from island.forcefields.selection import prepare_polymer_forcefield as legacy

    if forcefield == "opls" and not legacy_options:
        assignment = canonical_prepare(
            system, request or ForceFieldRequest(), assets=assets
        )
        a = assignment.payload
        name = a["selected_profile"]
        rows, conventions, audit = [], {}, None
        valid = a["status"] in ("prepared", "prepared_static_only")
        diagnostics = list(a["diagnostics"])
        if valid:
            try:
                backend = assignment.validate_integrity(system)
                if name == PROFILE:
                    rows = _pdms_rows(
                        system, backend.validate_integrity(system).payload
                    )
                    conventions, audit = _pdms(system, assignment, sources)
                elif name == NATIVE:
                    rows = backend.payload["requirements"]
                    conventions = _native_conventions(assignment)
            except Exception as error:  # noqa: BLE001 -- inspection keeps exact blockers
                diagnostics.append({"stage": "convention_audit", "reason": str(error)})
                valid = False
        stages = a["stages"]
    else:
        assignment = legacy(system, forcefield, **legacy_options)
        a = assignment.payload
        name, stages = a["profile_identity"], a["stages"]
        rows, conventions, audit = a["requirements"], {}, None
        valid = a["status"] == "prepared"
        diagnostics = list(a["diagnostics"])
        # No broad conversion adapter grants a new completeness claim. Preserve
        # the full existing row/charge evidence and exact coupling gaps.
        if forcefield == "opls" and "external_assignment" in a and valid:
            conventions = _external_conventions(assignment)
        elif forcefield in ("gaff", "gaff2") and valid:
            conventions = _amber_conventions(assignment)
        elif forcefield == "pcff" and valid:
            conventions = _pcff_conventions(assignment)
        else:
            diagnostics.append(
                {
                    "stage": "convention_audit",
                    "reason": "Source requirements are incomplete; no cross-engine prerequisite",
                }
            )
    p = assess(
        requirements(system, forcefield),
        rows,
        conventions,
        typing_complete=stages["typing_complete"],
        charge_complete=stages["charge_complete"],
        source_valid=valid,
    )
    from island.workflows.bundle import system_data

    p.update(
        execution_input_identity=digest(system_data(system)),
        schema=SCHEMA,
        policy_identity=policy_identity(),
        family=forcefield,
        selected_profile=name,
        source_kind=a.get("source_kind"),
        source_identity=a.get("source_identity"),
        charge_origin=a.get("charge_origin", a.get("charge_method")),
        graph_identity=a.get("graph_identity", a.get("final_graph_identity")),
        assignment_identity=a["identity"],
        charge_record=a.get("graph_charges", a.get("charge_evidence")),
        domain=a.get("domain"),
        source_assignment=a,
        conventions=conventions,
        source_convention_gaps=[]
        if audit is None
        else audit.payload["unresolved_conventions"],
        backend_settings={}
        if audit is None
        else audit.payload["backend_execution_settings"],
        diagnostics=diagnostics,
        lammps_export_status="not_run",
        lammps_runtime_status="not_run",
        openmm_build_status="not_run",
        openmm_runtime_status="not_run",
        cross_engine_equivalence="not_claimed",
        cross_engine_diagnostics=[],
        scientific_validation="not_established",
        production_readiness="not_established",
        production_validated=False,
        simulation_readiness="not_established",
        MD_authorized=False,
        numerical_evidence=[],
    )
    pins = []
    if hasattr(assignment, "_files"):
        pins.extend(assignment._files)
    if forcefield == "opls" and hasattr(assignment, "assets"):
        selected = assignment.assets
        for file in (selected.archive_path, selected.native_source_path):
            if file is not None and Path(file).is_file():
                pins.append((str(file), sha256(Path(file).read_bytes()).hexdigest()))
    if audit is not None:
        for file in (audit.sources.paper, audit.sources.supplement):
            pins.append((str(file), sha256(Path(file).read_bytes()).hexdigest()))
    result = _signed(p, assignment, system, audit, tuple(pins))
    return result.require(claims)


def _update(prepared, system, **changes):
    p = prepared.validate_integrity(system)
    p.pop("identity")
    p.update(changes)
    return _signed(p, prepared.assignment, system, prepared._audit, prepared._pins)


def export_lammps(system, prepared, directory):
    """Independent export status, preserving diagnostic-only backend authority."""
    p = prepared.validate_integrity(system)
    try:
        require(
            p["selected_profile"] == PROFILE and prepared._audit is not None,
            "asset_unavailable",
            "no_lammps_export_adapter",
        )
        from island.opls_resolved_execution import generate_lammps

        files = generate_lammps(system, prepared._audit, directory)
        result = _update(
            prepared,
            system,
            lammps_export_status="passed",
            lammps_export={"files": files, "diagnostic_only": True},
        )
    except Exception as error:  # noqa: BLE001 -- backend failure is a separate claim
        result = _update(
            prepared,
            system,
            lammps_export_status="blocked",
            lammps_export={"diagnostic": str(error)},
        )
    return result


def build_openmm(system, prepared):
    """Return (new report, evaluator); construction never claims equivalence."""
    prepared.validate_integrity(system)
    from island.forcefields.opls_registry import create_evaluator

    try:
        if prepared.payload["family"] == "opls" and hasattr(
            prepared.assignment, "assets"
        ):
            model = create_evaluator(system, prepared.assignment)
        else:
            from island.forcefields.selection import create_evaluator as existing

            model = existing(system, prepared.assignment)
        return _update(prepared, system, openmm_build_status="passed"), model
    except Exception as error:  # noqa: BLE001 -- preserve other claims
        return _update(
            prepared,
            system,
            openmm_build_status="failed",
            openmm_build_diagnostic=str(error),
        ), None


def check_openmm(system, prepared):
    """Backend-specific static check; no MD or scientific promotion."""
    result, model = build_openmm(system, prepared)
    if model is None:
        return _update(result, system, openmm_runtime_status="blocked")
    try:
        observed = model.evaluate()
        energy = getattr(observed, "potential_energy", None)
        # Actual evaluator result is immutable; preserve the full public record.
        raw = {
            "potential_energy": observed.potential_energy,
            "energy_components": dict(observed.energy_components),
            "forces": dict(observed.forces),
            "backend": observed.backend_name,
            "version": observed.backend_version,
            "settings": dict(observed.settings),
        }
        return _update(
            result,
            system,
            openmm_runtime_status="passed",
            openmm_static={"result": raw, "energy": energy, "scope": "static_only"},
        )
    except Exception as error:  # noqa: BLE001
        return _update(
            result,
            system,
            openmm_runtime_status="failed",
            openmm_runtime_diagnostic=str(error),
        )
    finally:
        if hasattr(model, "close"):
            model.close()


def check_lammps(system, prepared, directory, *, engine):
    prepared.validate_integrity(system)
    try:
        require(
            prepared._audit is not None, "asset_unavailable", "no_lammps_static_adapter"
        )
        from island.opls_resolved_execution import run_static

        observed = run_static(system, prepared._audit, directory, engine=engine)
        return _update(
            prepared,
            system,
            lammps_runtime_status="passed",
            lammps_export_status="passed",
            lammps_static=observed,
        )
    except Exception as error:  # noqa: BLE001
        return _update(
            prepared,
            system,
            lammps_runtime_status="failed",
            lammps_runtime_diagnostic=str(error),
        )


__all__ = [
    "ClaimRequirements",
    "ForceFieldRequest",
    "OPLSAssets",
    "Preparation",
    "build_openmm",
    "check_lammps",
    "check_openmm",
    "export_lammps",
    "prepare_polymer_forcefield",
    "use_opls_assets",
]


def prepare_forcefield(system, request=None, **kwargs):
    request = ForceFieldRequest() if request is None else request
    return prepare_polymer_forcefield(system, request.family, request=request, **kwargs)


def prepare_polymer(
    *,
    psmiles,
    forcefield="opls",
    dp=None,
    head_end_group=None,
    tail_end_group=None,
    coordinates=None,
    box=None,
    **kwargs,
):
    from island.builders import build_linear_polymer
    from island.core import Coordinates
    from island.forcefields.opls_registry.api import PolymerPreparation

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


__all__ += ["prepare_forcefield", "prepare_polymer"]
