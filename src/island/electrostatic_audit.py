"""Independent coordinate-wise numerical audit; never grants MD authority.

Callbacks return raw kJ/mol components. Coordinates are angstrom, forces are
kJ/(mol angstrom). No source parser, charge repair or energy correction is used.
"""

import json
from dataclasses import dataclass

import numpy as np

from island.forcefields.opls_registry.api import require
from island.forcefields.selection.inventory import canonical, digest
from island.pdms import box_lengths, special_pairs

FD_BOUND = 1e-4  # Inherited, not fitted to this audit.
COMPONENTS = ("bonded", "lj", "coulomb", "total")


@dataclass(frozen=True)
class InputBinding:
    json_text: str

    @property
    def payload(self):
        return json.loads(self.json_text)

    @property
    def identity(self):
        return digest(self.payload)

    def validate(self, system, preparation):
        expected = bind_inputs(system, preparation)
        require(
            self.json_text == expected.json_text,
            "source_incompatible",
            "electrostatic_input_binding_mutation",
        )
        return expected.payload


def bind_inputs(system, preparation):
    """Recompute chemical/source identity before inspecting the numerical contract."""
    backend = preparation.validate_integrity(system).validate_integrity(system)
    p = backend.payload
    atoms = p["mapping"]
    value = {
        "schema": "island.electrostatic-input-binding.v1",
        "family": "opls",
        "selected_profile": preparation.payload["selected_profile"],
        "canonical_assignment_identity": preparation.identity,
        "source_identity": digest(p["source"]),
        "graph_identity": p["graph_identity"],
        "charge_origin": preparation.payload["charge_origin"],
        "domain": preparation.payload["domain"],
        "atoms": atoms,
        "atom_order": [a["site_id"] for a in atoms],
        "coordinates_angstrom": [
            list(map(float, system.coordinates.get(a["site_id"]))) for a in atoms
        ],
        "box_angstrom": box_lengths(system),
        "special_pairs": [
            [a, b, scale] for (a, b), scale in sorted(special_pairs(p).items())
        ],
        "source_model": p["model"],
        "numerical_contract": {
            "units": "angstrom, e, dalton; kcal/mol*4.184=kJ/mol",
            "precision": "double; %.17g engine output",
            "dielectric": 1,
            "boundary": "p p p; conducting Ewald boundary",
            "coulomb_cutoff_angstrom": 11,
            "coulomb_switching": "none",
            "lj_switching": "CHARMM r-squared 9 to 11 angstrom",
            "lj_mixing": "geometric sigma and epsilon",
            "coulomb_14": 0.5,
            "lj_14": 0.5,
            "exclusions": "graph distance 1 and 2 excluded",
            "OpenMM": {"solver": "Ewald", "tolerance": 1e-7, "platform": "Reference"},
            "LAMMPS": {"solver": "PPPM", "tolerance": 1e-4, "table": "default"},
            "matched_reciprocal_contract": "not supplied by source; diagnostic refinements grant no MD authority",
        },
    }
    return InputBinding(canonical(value))


def finite_difference(
    energy,
    coordinates,
    analytical_forces,
    *,
    probes=None,
    steps=(2e-4, 1e-4, 5e-5),
    bound=FD_BOUND,
):
    """Central differences with raw signed energies, across independent step sizes.

    `analytical_forces` maps each component to an (N,3) array. The caller must
    identify whether a component force is measured directly or reconstructed.
    Bound overrides are deliberately rejected. A failed sample never disappears
    into an average or an error-corrected result.
    """
    require(bound == FD_BOUND, "source_incompatible", "acceptance_bound_changed")
    xyz = np.asarray(coordinates, dtype=float)
    require(
        xyz.ndim == 2 and xyz.shape[1] == 3 and np.isfinite(xyz).all(),
        "source_incompatible",
        "invalid_fd_coordinates",
    )
    steps = tuple(steps)
    require(
        len(set(steps)) >= 2 and all(np.isfinite(h) and h > 0 for h in steps),
        "source_incompatible",
        "multiple_positive_fd_steps_required",
    )
    force = {k: np.asarray(v, dtype=float) for k, v in analytical_forces.items()}
    require(
        set(force) == set(COMPONENTS)
        and all(f.shape == xyz.shape and np.isfinite(f).all() for f in force.values()),
        "source_incompatible",
        "incomplete_fd_force_components",
    )
    probes = list(probes) if probes is not None else list(np.ndindex(xyz.shape))
    require(
        len(set(probes)) == len(probes)
        and all(0 <= n < len(xyz) and axis in (0, 1, 2) for n, axis in probes),
        "source_incompatible",
        "invalid_fd_probes",
    )
    records = []
    for n, axis in probes:
        for h in steps:
            values = []
            for sign in (1, -1):
                displaced = xyz.copy()
                displaced[n, axis] += sign * h
                row = energy(displaced)
                require(
                    set(row) == set(COMPONENTS)
                    and all(np.isfinite(v) for v in row.values()),
                    "source_incompatible",
                    "incomplete_fd_energy_components",
                )
                values.append({k: float(v) for k, v in row.items()})
            components = {}
            for k in COMPONENTS:
                fd = -(values[0][k] - values[1][k]) / (2 * h)
                observed = float(force[k][n, axis])
                components[k] = {
                    "positive_energy_kj_mol": values[0][k],
                    "negative_energy_kj_mol": values[1][k],
                    "central_difference_force": fd,
                    "analytical_force": observed,
                    "absolute_error": abs(fd - observed),
                    "sign_reversal_error": abs(fd + observed),
                    "nm_to_angstrom_error": abs(fd - observed * 0.1),
                    "kcal_to_kj_error": abs(fd - observed * 4.184),
                    "passed": abs(fd - observed) <= FD_BOUND,
                }
            records.append(
                {
                    "atom_index": n,
                    "axis": axis,
                    "h_angstrom": h,
                    "components": components,
                }
            )
    maxima = {
        k: max((r["components"][k]["absolute_error"] for r in records), default=0.0)
        for k in COMPONENTS
    }
    require(bool(records), "source_incompatible", "empty_fd_audit")
    return {
        "schema": "island.raw-central-difference.v1",
        "equation": "F_i=-(E(x+h*e_i)-E(x-h*e_i))/(2*h)",
        "coordinate_unit": "angstrom",
        "energy_unit": "kJ/mol",
        "force_unit": "kJ/(mol angstrom)",
        "steps_angstrom": list(steps),
        "bound": FD_BOUND,
        "comparisons_per_component": len(records),
        "max_error_by_component": maxima,
        "passed_by_component": {k: v <= FD_BOUND for k, v in maxima.items()},
        "passed": all(v <= FD_BOUND for v in maxima.values()),
        "raw_records": records,
    }
