"""Four-layer, replay-validated PDMS convention audit and execution policy.

Backend choices describe a diagnostic experiment. They cannot supply missing
author authority, change the physical profile, or enable MD through re-signing.
"""

import json
from dataclasses import dataclass, field
from hashlib import sha256
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

from island.forcefields.external_pdms import LIMITS
from island.forcefields.opls_registry import prepare_polymer_forcefield as resolve
from island.forcefields.opls_registry.api import require
from island.forcefields.opls_registry.registry import ARCHIVE_SHA256, PAPER_SHA256
from island.forcefields.oplsaa.pdms_author import PINS
from island.forcefields.pdms_opls import PROFILE
from island.forcefields.selection.inventory import canonical, digest

from .evidence import CONVERTER_PINS, ENGINE_EVIDENCE, PAPER_PBC

SCHEMA = "island.opls-pdms-convention-reconciliation.v1"
POLICIES = ("source_replay_v1", "exact_erfc_diagnostic_v1")
MODES = ("diagnostic", "research_opt_in", "production_candidate")
LAYERS = (
    "source_authoritative_settings",
    "converter_settings",
    "backend_execution_settings",
    "scientific_verification_settings",
)
CLASSIFICATIONS = (
    "explicitly_stated_by_source",
    "explicitly_stated_by_converter",
    "inherited_from_forcefield_family",
    "derived_by_documented_formula",
    "selected_by_backend_policy",
    "unresolved",
)
AUTHOR_ENGINE = "6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191"
EXACT_ENGINE = "3a989d8e36345bd4508d932832a78e0eca3313b134d5c414c280065d615029ea"

# No implicit source defaults: None means absence of source evidence, even if
# LAMMPS would fill the field at runtime. All are explicitly set for diagnostics.
SOURCE_VALUES = {
    "lammps_units": "real",
    "distance_units": "angstrom",
    "energy_units": "kcal/mol",
    "force_units": "kcal/(mol angstrom)",
    "charge_units": "elementary_charge",
    "dielectric": None,
    "coulomb_prefactor": None,
    "lj_functional_form": "4*epsilon*((sigma/r)^12-(sigma/r)^6)",
    "lj_epsilon_convention": "well_depth; kcal/mol",
    "lj_sigma_convention": "zero_crossing; r_min=2^(1/6)*sigma; angstrom",
    "pair_style": "lj/charmm/coul/long",
    "bond_style": "harmonic",
    "angle_style": "harmonic",
    "dihedral_style": "opls",
    "improper_style": "harmonic; no impropers in bounded graph",
    "torsion_convention": "sum_n K_n/2*(1+(-1)^(n+1)*cos(n*phi)); n=1..4",
    "special_bonds": "lj/coul 0 0 0.5",
    "exclusions_12": "LJ=0; Coulomb=0",
    "exclusions_13": "LJ=0; Coulomb=0",
    "coulomb_14_scaling": "0.5",
    "lj_14_scaling": "0.5",
    "lj_mixing": "geometric_sigma_epsilon",
    "pair_cutoff_angstrom": "11",
    "lj_inner_angstrom": "9",
    "switching_function": "CHARMM r-squared energy switch for LJ 9..11; no Coulomb switch",
    "kspace_style": "pppm",
    "reciprocal_tolerance": "0.0001",
    "boundary_conditions": "p p p",
    "periodic_image_policy": None,
    "precision": None,
    "neighbor_list_policy": None,
    "long_range_solver": "pppm",
    "pressure_convention": None,
    "virial_convention": None,
    "pair_table_bits": None,
    "erfc_kernel": None,
    "lj_tail_correction": None,
}
REQUIRED_SETTINGS = tuple(sorted(SOURCE_VALUES))


def signed(value):
    return {**value, "identity": digest(value)}


def row(value, classification, evidence, **extra):
    return {
        "value": value,
        "classification": classification,
        "evidence": evidence,
        **extra,
    }


def _source(canonical_preparation):
    """Rehash the explicit ZIP, then cite exact lines and proven include scope."""
    from island.forcefields.pdms_opls import load_source

    request = canonical_preparation._backend.request
    options = request.validate()
    load_source(options)  # Checks raw archive and all pinned companions anew.
    raw_archive = Path(options.archive_path).read_bytes()
    require(
        sha256(raw_archive).hexdigest() == ARCHIVE_SHA256,
        "source_incompatible",
        "source_archive_mutation_during_audit",
    )
    with ZipFile(BytesIO(raw_archive)) as z:
        texts = {}
        inventory = {}
        for name in z.namelist():
            if name.endswith("/"):
                continue
            raw = z.read(name)
            relative = name.split("/", 1)[-1]
            inventory[relative] = sha256(raw).hexdigest()
            if name.endswith((".in", ".init", ".inp", "README.md")):
                texts[relative] = raw.decode()
    init_name = "param/OPLS/opls_system.in.init"
    init = texts[init_name]
    commands = {}
    for n, raw in enumerate(init.splitlines(), 1):
        words = raw.split("#", 1)[0].split()
        if words:
            commands[words[0]] = {
                "source_file": init_name,
                "source_sha256": PINS[init_name],
                "line": n,
                "raw_text": raw,
                "scope": "exact_pinned_OPLS_include_chain",
            }
    command_for = {
        "lammps_units": "units",
        "pair_style": "pair_style",
        "bond_style": "bond_style",
        "angle_style": "angle_style",
        "dihedral_style": "dihedral_style",
        "improper_style": "improper_style",
        "special_bonds": "special_bonds",
        "lj_mixing": "pair_modify",
        "pair_cutoff_angstrom": "pair_style",
        "lj_inner_angstrom": "pair_style",
        "kspace_style": "kspace_style",
        "reciprocal_tolerance": "kspace_style",
        "long_range_solver": "kspace_style",
    }
    derived = {
        **{
            k: ("units", "doc/src/units.rst")
            for k in ("distance_units", "energy_units", "force_units", "charge_units")
        },
        **{
            k: ("pair_style", "doc/src/pair_charmm.rst")
            for k in (
                "lj_functional_form",
                "lj_epsilon_convention",
                "lj_sigma_convention",
                "switching_function",
            )
        },
        "torsion_convention": ("dihedral_style", "doc/src/dihedral_opls.rst"),
        **{
            k: ("special_bonds", "doc/src/special_bonds.rst")
            for k in (
                "exclusions_12",
                "exclusions_13",
                "coulomb_14_scaling",
                "lj_14_scaling",
            )
        },
    }
    settings = {}
    for name, value in sorted(SOURCE_VALUES.items()):
        if name == "boundary_conditions":
            settings[name] = row(value, CLASSIFICATIONS[0], [PAPER_PBC])
        elif name in command_for:
            settings[name] = row(
                value, CLASSIFICATIONS[0], [commands[command_for[name]]]
            )
        elif name in derived:
            command, engine_doc = derived[name]
            settings[name] = row(
                value,
                CLASSIFICATIONS[3],
                [
                    commands[command],
                    {"file": engine_doc, **ENGINE_EVIDENCE[engine_doc]},
                ],
                formula="literal command interpreted by pinned LAMMPS definition",
                input_rows=[command],
                independent_reconstruction="J47/J49 source audit",
            )
        else:
            occurrences = []
            tokens = {
                "dielectric": "dielectric",
                "boundary_conditions": "boundary",
                "neighbor_list_policy": "neigh_modify",
                "pressure_convention": "thermo_style",
                "virial_convention": "compute",
            }
            for file, text in sorted(texts.items()):
                for n, raw in enumerate(text.splitlines(), 1):
                    if raw.strip().startswith(
                        tokens.get(name, "__not_a_command__") + " "
                    ):
                        occurrences.append(
                            {
                                "source_file": file,
                                "source_sha256": inventory[file],
                                "line": n,
                                "raw_text": raw,
                                "scope": "not_linked_to_pinned_OPLS_DFT_include_chain",
                            }
                        )
            settings[name] = row(
                None,
                CLASSIFICATIONS[5],
                occurrences,
                missing_evidence="OPLS DFT/NPA-scoped author declaration of " + name,
                searched_include_chain=[
                    "param/OPLS/opls_system.in",
                    init_name,
                    "param/OPLS/opls_system.in.settings",
                ],
                absence_evidence_sha256=digest(inventory),
            )
    return signed(
        {
            "kind": LAYERS[0],
            "settings": settings,
            "archive_sha256": ARCHIVE_SHA256,
            "paper_sha256": PAPER_SHA256,
            "companion_hashes": PINS,
            "archive_file_inventory": inventory,
            "cross_family_imports": [],
            "other_input_scope": "generic system.in.init/system.in.charges; no proven OPLS DFT include binding",
        }
    )


def _converter():
    # Explicit converter equations, not new parameter rows or source defaults.
    package = Path(__file__).parents[1]
    pins = {n: sha256((package / n).read_bytes()).hexdigest() for n in CONVERTER_PINS}
    require(
        pins == CONVERTER_PINS,
        "source_incompatible",
        "converter_implementation_mutation",
    )
    conversions = {
        "distance": "angstrom*0.1=nm",
        "energy": "kcal/mol*4.184=kJ/mol",
        "force": "kcal/(mol angstrom)*4.184=kJ/(mol angstrom)",
        "bond": "K_OpenMM=2*4.184*100*K_LAMMPS; r0_nm=0.1*r0_A",
        "angle": "K_OpenMM=2*4.184*K_LAMMPS; theta0_rad=theta0_deg*pi/180",
        "torsion": "K_OpenMM=4.184*K_LAMMPS; retain explicit 1/2",
        "lj": "epsilon_kJ=4.184*epsilon_kcal; sigma_nm=0.1*sigma_A; no r-star substitution",
        "charges": "verbatim serialized e; no normalization or redistribution",
        "pressure_proxy": "virial_kJ=P_config_atm*V_A3*4.184/68568.415; diagnostic only",
    }
    return signed(
        {
            "kind": LAYERS[1],
            "version": "equation_preserving_real_to_OpenMM_v1",
            "implementation_hashes": pins,
            "engine_definitions": ENGINE_EVIDENCE,
            "settings": {
                k: row(v, CLASSIFICATIONS[1], [{"implementation_hashes": pins}])
                for k, v in sorted(conversions.items())
            },
        }
    )


def _backend(source, policy):
    settings = json.loads(canonical(source["settings"]))
    choices = {
        "dielectric": "1",
        "coulomb_prefactor": "332.06371 kcal angstrom/(mol e^2)",
        "boundary_conditions": "p p p",
        "periodic_image_policy": "3D conducting periodic Ewald sum; minimum-image real space; no slab correction",
        "precision": "serial CPU double; output %.17g",
        "neighbor_list_policy": "neighbor 2 bin; every 1 delay 0 check no; atom sorting disabled",
        "pressure_convention": "configurational proxy: compute pressure NULL virial; atm; no kinetic term",
        "virial_convention": "xx yy zz xy xz yz; proxy_atm*V_A3*4.184/68568.415; not scientifically qualified",
        "pair_table_bits": "12",
        "erfc_kernel": "original pinned LAMMPS polynomial",
        "lj_tail_correction": "no",
    }
    if policy == POLICIES[1]:
        choices.update(
            kspace_style="ewald",
            long_range_solver="ewald",
            reciprocal_tolerance="0.0000001",
            pair_table_bits="0",
            erfc_kernel="J49 explicitly modified std::erfc kernel; force equations unchanged",
        )
    for name, value in choices.items():
        settings[name] = row(
            value,
            CLASSIFICATIONS[4],
            [{"policy": policy}],
            source_value=source["settings"][name]["value"],
            source_authority=False,
        )
    return signed(
        {
            "kind": LAYERS[2],
            "policy": policy,
            "settings": settings,
            "engine_sha256": AUTHOR_ENGINE if policy == POLICIES[0] else EXACT_ENGINE,
            "engine_definitions": ENGINE_EVIDENCE,
            "ewald_alpha_formula": None
            if policy == POLICIES[0]
            else "sqrt(-log(2*1e-7))/11 angstrom^-1",
            "ewald_kmax": None if policy == POLICIES[0] else [20, 20, 20],
            "newton": "on",
            "processes": 1,
            "NPT_authorized": False,
            "source_authority": False,
        }
    )


def _science():
    return signed(
        {
            "kind": LAYERS[3],
            "energy_bound_kj_mol": 3e-4,
            "force_bound_kj_mol_angstrom": 3e-4,
            "FD_bound_kj_mol_angstrom": 1e-4,
            "FD_steps_angstrom": [2e-4, 1e-4, 5e-5],
            "settings_classification": CLASSIFICATIONS[4],
            "reference": "J47-J49 unchanged numerical acceptance; no source authority",
            "prior_receipt_sha256": "8bcff997fc807c7fc45cadb9aca3a5357638d844bebfb79de3c8aa954fb5c4e9",
            "prior_verified_cases_only": [
                "HMDS",
                "capped_DP2",
                "capped_DP4",
                "capped_DP6",
            ],
            "current_graph_numerically_verified_by_preparation": False,
            "source_replay_parity": False,
            "exact_erfc_bounded_diagnostic_parity": True,
            "research_MD_authorized": False,
            "production_candidate_available": False,
            "unpassed_gates": [
                "OPLS_DFT_source_conventions",
                "author_numerical_reconciliation",
                "polymer_Ewald_PME",
                "backend_pressure_and_virial",
                "cutoff_recovery",
                "sampling",
                "long_NPT",
            ],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


def _document(system, assignment, mode, policy):
    assignment.validate_integrity(system)
    require(
        assignment.payload["family"] == "opls"
        and assignment.payload["selected_profile"] == PROFILE,
        "source_incompatible",
        "external_OPLS_PDMS_profile_required",
    )
    source = _source(assignment)
    records = {
        LAYERS[0]: source,
        LAYERS[1]: _converter(),
        LAYERS[2]: _backend(source, policy),
        LAYERS[3]: _science(),
    }
    unresolved = [
        {"layer": LAYERS[0], "setting": k, **v}
        for k, v in source["settings"].items()
        if v["classification"] == CLASSIFICATIONS[5]
    ]
    from island.electrostatic_audit import bind_inputs

    return signed(
        {
            "schema": SCHEMA,
            "family": "opls",
            "selected_profile": PROFILE,
            "source_kind": "external_literature",
            "charge_origin": assignment.payload["charge_origin"],
            "domain": assignment.payload["domain"],
            "assignment_identity": assignment.identity,
            "input_binding": bind_inputs(system, assignment).payload,
            "execution_mode": mode,
            "backend_policy": policy,
            **records,
            "unresolved_settings": unresolved,
            "required_settings": list(REQUIRED_SETTINGS),
            "default_MD_enabled": False,
            "research_MD_available": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "limitations": LIMITS,
        }
    )


@dataclass(frozen=True)
class ConventionRecord:
    json_text: str

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return self.payload["identity"]


@dataclass(frozen=True)
class ConventionPreparation(ConventionRecord):
    assignment: object = field(repr=False)
    mode: str = "diagnostic"
    policy: str = "source_replay_v1"

    def validate_integrity(self, system):
        require(
            self.mode in MODES and self.policy in POLICIES,
            "source_incompatible",
            "unknown_execution_mode_or_policy",
        )
        expected = _document(system, self.assignment, self.mode, self.policy)
        differences = []
        for layer in (LAYERS[0], LAYERS[2]):
            recorded = self.payload.get(layer, {}).get("settings", {})
            for setting in REQUIRED_SETTINGS:
                actual = recorded.get(setting)
                wanted = expected[layer]["settings"][setting]
                if actual != wanted:
                    differences.append(
                        {
                            "layer": layer,
                            "setting": setting,
                            "recorded": actual,
                            "expected": wanted,
                        }
                    )
        require(
            self.json_text == canonical(expected),
            "source_incompatible",
            "stale_or_mutated_convention_record",
            result_identity=self.identity,
            changed_settings=differences,
        )
        return self.assignment.validate_integrity(system)

    def records(self):
        return {
            kind: ConventionRecord(canonical(self.payload[kind])) for kind in LAYERS
        }

    def require_execution(self, system, action="static", *, research_opt_in=False):
        self.validate_integrity(system)
        require(
            action in ("static", "minimize", "nvt", "npt"),
            "source_incompatible",
            "unsupported_execution_action",
        )
        p = self.payload
        backend = p[LAYERS[2]]
        missing = [
            {
                "setting": k,
                "missing_evidence": "resolved selected-backend value",
                "evidence": backend["settings"].get(k, {}).get("evidence", []),
            }
            for k in REQUIRED_SETTINGS
            if backend["settings"].get(k, {}).get("value") is None
        ]
        require(
            not missing,
            "parameter_incomplete",
            "unresolved_backend_settings",
            missing=missing,
        )
        require(
            self.mode != "production_candidate",
            "scientifically_unqualified",
            "production_candidate_unavailable",
            missing=p[LAYERS[3]]["unpassed_gates"],
        )
        if self.mode == "diagnostic":
            require(
                action == "static",
                "scientifically_unqualified",
                "diagnostic_mode_does_not_authorize_MD",
                action=action,
            )
        else:
            require(
                research_opt_in is True,
                "scientifically_unqualified",
                "explicit_research_opt_in_required",
            )
            require(
                not p["unresolved_settings"],
                "parameter_incomplete",
                "source_conventions_unresolved",
                missing=p["unresolved_settings"],
            )
            require(
                p[LAYERS[3]]["research_MD_authorized"],
                "scientifically_unqualified",
                "author_convention_reconciliation_incomplete",
                missing=p[LAYERS[3]]["unpassed_gates"],
            )
        require(
            action != "npt",
            "scientifically_unqualified",
            "NPT_pressure_gate_unavailable",
        )
        return backend


def prepare_polymer_forcefield(
    system,
    forcefield="opls",
    *,
    execution_mode="diagnostic",
    backend_policy="source_replay_v1",
    **kwargs,
):
    """Automatic canonical chemistry, explicit four-layer policy; no expert profile needed."""
    require(forcefield == "opls", "source_incompatible", "OPLS_family_required")
    require(
        execution_mode in MODES and backend_policy in POLICIES,
        "source_incompatible",
        "unknown_execution_mode_or_policy",
    )
    assignment = resolve(system, forcefield=forcefield, **kwargs)
    p = _document(system, assignment, execution_mode, backend_policy)
    return ConventionPreparation(
        canonical(p), assignment, execution_mode, backend_policy
    )


__all__ = ["ConventionPreparation", "ConventionRecord", "prepare_polymer_forcefield"]
