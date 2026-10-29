"""Additive, opt-in PDMS LAMMPS boundary; historical facades stay unchanged.

Only the author DFT/NPA external OPLS-compatible model is accepted. Static
execution authority does not grant MD, pressure, native OPLS or PCFF authority.
"""

from dataclasses import dataclass

from island.charge_references.records import pack, unpack
from island.forcefields.external_pdms import LIMITS, digest
from island.forcefields.pdms_opls import (
    PROFILE,
    PDMSOPLSError,
    PDMSOPLSOptions,
    check,
    prepare_pdms_opls,
)
from island.pdms import box_lengths
from island.workflows.bundle import system_data


@dataclass(frozen=True)
class OPLSOptions:
    profile: str
    archive_path: object = None
    expected_sha256: str | None = None
    e_oxygen_site_ids: tuple[int, ...] = ()
    use_scope: str = "local_research"
    lj_mixing: str = "geometric_sigma_epsilon"
    special_lj: tuple = (0, 0, 0.5)
    special_coulomb: tuple = (0, 0, 0.5)
    lj_inner_angstrom: float = 9.0
    cutoff_angstrom: float = 11.0
    reciprocal_method: str = "pppm"
    reciprocal_tolerance: float = 1e-4

    def source_options(self):
        check(self.profile == PROFILE, "unsupported_external_profile", "typing_blocked")
        expected = {
            "lj_mixing": "geometric_sigma_epsilon",
            "special_lj": (0, 0, 0.5),
            "special_coulomb": (0, 0, 0.5),
            "lj_inner_angstrom": 9.0,
            "cutoff_angstrom": 11.0,
            "reciprocal_method": "pppm",
            "reciprocal_tolerance": 1e-4,
        }
        for name, value in expected.items():
            check(getattr(self, name) == value, "source_convention_changed", field=name)
        return PDMSOPLSOptions(
            self.archive_path,
            self.expected_sha256,
            self.e_oxygen_site_ids,
            self.use_scope,
        )


@dataclass(frozen=True)
class ForceFieldRequest:
    family: str
    options: OPLSOptions

    def validate(self):
        check(
            self.family == "opls",
            "external_opls_profile_family_mismatch",
            "typing_blocked",
        )
        check(
            type(self.options) is OPLSOptions,
            "explicit_opls_options_required",
            "typing_blocked",
        )
        return self.options.source_options()


@dataclass(frozen=True)
class Preparation:
    json_text: str
    request: ForceFieldRequest

    @property
    def payload(self):
        return unpack(self.json_text)

    def validate_integrity(self, system):
        expected = prepare_forcefield(system, self.request)
        check(
            expected.json_text == self.json_text,
            "stale_or_tampered_lammps_assignment",
            "typing_blocked",
        )
        if self.payload["status"] != "prepared_static_only":
            error = PDMSOPLSError(
                self.payload["status"], "external_preparation_blocked", result=self
            )
            error.diagnostics = self.payload["diagnostics"]
            raise error
        return prepare_pdms_opls(system, self.request.validate())

    def require_md(self, system):
        self.validate_integrity(system)
        raise PDMSOPLSError(
            "scientifically_unqualified",
            "lammps_md_acceptance_incomplete",
            result=self,
            missing=[
                "identical_long_range_solver_contract",
                "exact_component_energy_parity",
                "qualified_force_and_virial_parity",
                "qualified_minimization_and_NVT",
                "exact_fresh_process_checkpoint_and_complete_RNG_equality",
            ],
        )


def prepare_forcefield(system, request):
    """Validate source/graph before any numerical engine is constructed."""
    check(
        type(request) is ForceFieldRequest,
        "explicit_forcefield_request_required",
        "typing_blocked",
    )
    try:
        source_options = request.validate()
        source = prepare_pdms_opls(system, source_options)
        source.require_executable(system)
        p = source.payload
        lengths = box_lengths(system)
        report = {
            "schema": "island_pdms_lammps_preparation_v1",
            "status": "prepared_static_only",
            "profile_identity": PROFILE,
            "profile_kind": "external_literature_opls_compatible_research_only",
            "preparation_identity": p["identity"],
            "source_identity": digest(p["source"]),
            "archive_sha256": p["archive_sha256"],
            "graph_identity": p["graph_identity"],
            "input_identity": digest(system_data(system)),
            "box_angstrom": list(lengths),
            "atom_order": [r["site_id"] for r in p["mapping"]],
            "ordered_charge_literals": [r["charge_e"] for r in p["mapping"]],
            "charge_origin": p["charge_policy"],
            "supported_domain": "explicit-H neutral linear trimethylsilyl-capped siloxanes; explicit E-oxygen anchors",
            "assigned_interactions": p["interactions"],
            "missing_interactions": [],
            "conventions": {
                **p["model"],
                "units": "real",
                "boundary": "p p p",
                "reciprocal": {"method": "pppm", "tolerance": 1e-4},
                "numerical_pair_table": "pinned executable default; no override",
            },
            "backend_availability": {
                "LAMMPS": "requires_explicit_executable_and_SHA256",
                "OpenMM": "existing_Ewald_research_reference",
            },
            "stages": {
                "source_complete": True,
                "typing_complete": True,
                "charge_complete": True,
                "parameter_complete": True,
                "static_export_complete": True,
                "executable_MD": False,
                "numerically_verified": False,
                "workflow_verified": False,
                "scientifically_validated": False,
            },
            "diagnostics": [],
            "limitations": LIMITS,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    except PDMSOPLSError as error:
        report = {
            "schema": "island_pdms_lammps_preparation_v1",
            "status": error.status,
            "profile_identity": PROFILE,
            "diagnostics": error.diagnostics,
            "execution_authority": False,
            "limitations": LIMITS,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    return Preparation(pack({**report, "identity": digest(report)}), request)


__all__ = ["ForceFieldRequest", "OPLSOptions", "Preparation", "prepare_forcefield"]
