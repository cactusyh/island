"""Reproducible target *requests*. No engine discovery, geometry generation or QM."""

from copy import deepcopy
from math import isfinite

from island.pcff_qm.targets import PROTOCOL as HISTORICAL_PROPOSED_PROTOCOL

from . import FAMILY, FIXTURES, FLAGS, PROFILE, digest, require, signed
from .inventory import chemical_identity

ROLES = ("training", "validation", "test")
PROTOCOL = {
    **deepcopy(HISTORICAL_PROPOSED_PROTOCOL),
    "authority": "proposed_extension_target_contract; no native/source authority",
    "software": "psi4",
    "software_version": None,
    "version_status": "required_before_execution_or_output_intake",
    "geometry_optimization": "equilibrium/constrained optimization requests need explicit future execution authorization",
    "charge": 0,
    "multiplicity": 1,
    "units": {
        "coordinates": "angstrom",
        "energy": "hartree",
        "gradient": "hartree/bohr",
        "hessian": "hartree/bohr^2",
    },
}


def plan(graphs, inventory):
    first = {}
    for i in inventory["instances"]:
        key = i["parameter_id"]
        if key not in first or len(graphs[i["fixture"]]["sites"]) < len(
            graphs[first[key]["fixture"]]["sites"]
        ):
            first[key] = i
    specs = []
    for t in inventory["terms"]:
        inst = first[t["id"]]
        f = t["family"]
        kind = (
            "constrained_bond_scan"
            if f == "quartic_bond"
            else "constrained_angle_scan"
            if f == "quartic_angle"
            else "torsion_or_methyl_rotation_scan"
            if f == "torsion_3"
            else "charge_LJ_reference_design_required"
            if f in ("partial_charge", "lj_9_6")
            else "coupled_two_dimensional_deformations"
        )
        specs.append(
            {
                "parameter_id": t["id"],
                "fixture": inst["fixture"],
                "sites": inst["sites"],
                "kind": kind,
                "equilibrium_dependencies": inst["equilibrium_dependencies"],
                "grid": {
                    "offsets": [-2, -1, 0, 1, 2],
                    "distance_step_angstrom": 0.02,
                    "angle_step_radian": 0.04,
                    "torsion_step_radian": 0.2,
                    "coupled_sign_combinations": "all four signs, independent axes; AAT at independent torsion blocks",
                },
                "required_properties": [
                    "energy",
                    "Cartesian_gradient",
                    "projected_curvature",
                ],
                "design_status": "not_identifiable_until_complete_orthogonal_design_is_supplied",
            }
        )
    for name in FIXTURES:
        for kind in (
            "equilibrium_geometry",
            "stationary_point_gradient",
            "Hessian_or_projected_curvature",
        ):
            specs.append(
                {
                    "parameter_id": None,
                    "fixture": name,
                    "sites": [],
                    "kind": kind,
                    "equilibrium_dependencies": [],
                    "required_properties": [kind],
                    "design_status": "reference_not_supplied",
                }
            )
    targets = []
    for spec in specs:
        name = spec["fixture"]
        payload = graphs[name]
        order = sorted(a["id"] for a in payload["sites"])
        coordinates = [payload["coordinates"][i] for i in order]
        for role in ROLES:
            spec_id = digest({"spec": spec, "role": role})
            request = {
                "id": "j78-" + spec_id[:24],
                "family": FAMILY,
                "profile": PROFILE,
                "dataset_role": role,
                "specification": spec,
                "graph_identity": chemical_identity(payload),
                "atom_order": order,
                "coordinate_order": order,
                "atom_mapping": [[i, i] for i in order],
                "seed_geometry": coordinates,
                "seed_geometry_sha256": digest(coordinates),
                "seed_origin": "historical_J54_unoptimized_input; not equilibrium or new_QM",
                "qm_protocol": deepcopy(PROTOCOL),
                "constraints": {
                    "sites": spec["sites"],
                    "grid": spec.get("grid"),
                    "actual_coordinates": None,
                    "status": "requires_independent_geometry_and_explicit_constraint_values",
                },
                "conformer_group": role + ":" + spec_id,
                "geometry_role_policy": "independent basin/geometry blocks; no identical/rigid/symmetry-equivalent geometries across roles",
            }
            targets.append(
                signed(
                    {
                        **request,
                        "input_hash": digest(request),
                        "output_hash": None,
                        "output_status": "not_executed",
                        "execution_authority": False,
                        "physical_observation": False,
                    }
                )
            )
    return signed(
        {
            "schema": "island.pdms-extension-QM-targets.v1",
            "targets": targets,
            "planned_target_requests": len(targets),
            "completed_QM_targets": 0,
            "native_engine_input_files": 0,
            "QM_executed": False,
            "execution_profile": None,
            "protocol": PROTOCOL,
            "scope": "four exact declared finite explicit-H neutral singlet graph fixtures",
            "seed_geometries_are_not_dataset_observations": True,
            **FLAGS,
        }
    )


def validate_target(target, expected):
    require(target == expected, "target_protocol_mapping_or_hash_mutation")
    require(
        target["input_hash"]
        == digest(
            {
                k: v
                for k, v in target.items()
                if k
                not in {
                    "identity",
                    "input_hash",
                    "output_hash",
                    "output_status",
                    "execution_authority",
                    "physical_observation",
                }
            }
        ),
        "input_hash_mutation",
    )
    require(
        target["atom_order"] == target["coordinate_order"], "coordinate_order_mutation"
    )
    return target


def validate_observation(observation, target, *, expected_protocol):
    """Intake only; original outputs and complete approved protocol are mandatory."""
    require(
        observation["qm_protocol"] == expected_protocol == target["qm_protocol"],
        "QM_protocol_mutation",
    )
    require(bool(expected_protocol["software_version"]), "QM_software_version_missing")
    for key in (
        "atom_order",
        "coordinate_order",
        "atom_mapping",
        "graph_identity",
        "dataset_role",
        "input_hash",
    ):
        require(
            observation[key] == target[key], "observation_binding_mutation", field=key
        )
    for key in ("output_hash", "input_hash"):
        value = observation.get(key)
        require(
            type(value) is str
            and len(value) == 64
            and all(c in "0123456789abcdef" for c in value),
            "native_input_output_hash_required",
        )
    require(
        observation.get("original_output_verified") is True,
        "original_output_not_verified",
    )
    require(observation.get("converged") is True, "unconverged_QM_output")
    x = observation["coordinates"]
    require(
        len(x) == len(target["atom_order"])
        and all(
            len(r) == 3 and all(type(v) in (int, float) and isfinite(v) for v in r)
            for r in x
        ),
        "complete_coordinates_required",
    )
    require(observation.get("geometry_sha256") == digest(x), "geometry_hash_mutation")


def split_manifest(targets):
    return signed(
        {
            "schema": "island.pdms-extension-split.v1",
            "roles": {
                role: [t["id"] for t in targets["targets"] if t["dataset_role"] == role]
                for role in ROLES
            },
            "physical_samples": {r: [] for r in ROLES},
            "split_before_fitting": True,
            "validation_reused_for_fit": False,
            "test_reused_for_fit": False,
            "independence_verified": False,
            "independence_status": "no_completed_geometries; role-specific_seed_copies_are_not_independent_data",
            "historical_J68_observations": "not_reused_as_independent_J78_heldout_data",
            "reserved_validation_activated": False,
        }
    )
