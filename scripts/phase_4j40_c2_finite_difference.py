"""Separate ordinary and crosslink-associated bonded central differences.

Analytic forces are never used to construct numerical derivatives. Every atom
and Cartesian direction is tested. Native Ewald settings and original acceptance
bounds are unchanged. Different components may have different usable intervals.
"""

from copy import deepcopy
from itertools import pairwise

import numpy as np

from island._periodic_openmm import build
from island.analysis.equilibration import require
from island.analysis.npt import FLAGS, _signed

STEPS_NM = (1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8, 3e-9)
BOUND = 1e-5  # kJ / (mol Angstrom), unchanged J22/J40 bound


def audit(payload, xyz, atom_ids, *, partition):
    import openmm as mm
    from openmm import unit

    xyz = np.asarray(xyz, dtype=float)
    require(
        xyz.shape == (len(atom_ids), 3) and np.isfinite(xyz).all(),
        "Invalid coordinate mapping",
    )
    require(
        atom_ids == sorted(s["id"] for s in payload["numerical_data"]["sites"]),
        "Changed atom ordering",
    )
    require(partition in {"crosslink", "ordinary"}, "Unknown term partition")
    edges = {frozenset(r["sites"]) for r in payload["system"]["metadata"]["crosslinks"]}
    data = deepcopy(payload["numerical_data"])
    counts = {}
    for key in ("bonds", "angles", "torsions", "impropers"):
        rows = data.get(key, [])

        def touches(row):
            return any(frozenset(pair) in edges for pair in pairwise(row["sites"]))

        data[key] = [r for r in rows if touches(r) == (partition == "crosslink")]
        counts[key] = len(data[key])
    model = build(data, payload["config"], mm)
    # Nonbonded terms are tested in the complete-model audit, not attributed
    # to a nonexistent separate crosslink force-field functional form.
    for i in reversed(range(model.getNumForces())):
        if model.getForce(i).getForceGroup() >= 30:
            model.removeForce(i)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(xyz)
    groups = {}
    for force in model.getForces():
        groups.setdefault(force.getForceGroup(), []).append(force.getName())
    results = {}
    for group, names in groups.items():
        context.setPositions(xyz)
        analytical = (
            context.getState(getForces=True, groups=1 << group)
            .getForces(asNumpy=True)
            .value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
            / 10
        )
        require(np.isfinite(analytical).all(), "Nonfinite analytic force")
        rows = []
        for h in STEPS_NM:
            worst = (-1.0, None)
            for atom in range(len(xyz)):
                for axis in range(3):
                    energy = []
                    for sign in (-1, 1):
                        x = xyz.copy()
                        x[atom, axis] += sign * h
                        context.setPositions(x)
                        energy.append(
                            context.getState(getEnergy=True, groups=1 << group)
                            .getPotentialEnergy()
                            .value_in_unit(unit.kilojoule_per_mole)
                        )
                    require(np.isfinite(energy).all(), "Nonfinite displaced energy")
                    numerical = -(energy[1] - energy[0]) / (20 * h)
                    error = abs(numerical - analytical[atom, axis])
                    if error > worst[0]:
                        worst = (
                            float(error),
                            {
                                "atom_id": atom_ids[atom],
                                "axis": axis,
                                "numerical": float(numerical),
                                "analytical": float(analytical[atom, axis]),
                            },
                        )
            rows.append(
                {
                    "step_nm": h,
                    "maximum_absolute_error_kj_mol_angstrom": worst[0],
                    "worst": worst[1],
                    "passed": worst[0] <= BOUND,
                }
            )
        # At least two adjacent displacements must pass; a single lucky point is insufficient.
        intervals = [
            [a["step_nm"], b["step_nm"]]
            for a, b in pairwise(rows)
            if a["passed"] and b["passed"]
        ]
        results["+".join(names)] = {
            "rows": rows,
            "converged_intervals_nm": intervals,
            "passed": bool(intervals),
        }
    del context, integrator
    lengths = np.asarray(payload["config"]["box_lengths"]) / 10
    pairs = []
    cutoff = payload["config"]["nonbonded_cutoff"] / 10
    for i in range(len(xyz)):
        for j in range(i):
            displacement = xyz[i] - xyz[j]
            image = displacement - lengths * np.rint(displacement / lengths)
            pairs.append(
                (
                    float(np.min(np.abs(np.abs(image) - lengths / 2))),
                    abs(float(np.linalg.norm(image)) - cutoff),
                )
            )
    data = payload["numerical_data"]
    angles = []
    for row in data["angles"]:
        a, b, c = [atom_ids.index(s) for s in row["sites"]]
        u, v = xyz[a] - xyz[b], xyz[c] - xyz[b]
        u -= lengths * np.rint(u / lengths)
        v -= lengths * np.rint(v / lengths)
        angles.append(
            float(
                np.linalg.norm(np.cross(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))
            )
        )
    return _signed(
        dict(
            schema="island.j40-c2-partition-finite-difference.v1",
            partition=partition,
            term_counts=counts,
            crosslink_edges=sorted(sorted(e) for e in edges),
            backend_identity=payload["identity"],
            atom_ids=atom_ids,
            coordinates_nm=xyz.tolist(),
            units="kJ/(mol Angstrom)",
            sign="force = negative central energy gradient",
            bound=BOUND,
            steps_nm=STEPS_NM,
            engine=mm.version.full_version,
            platform="Reference",
            constraints=model.getNumConstraints(),
            exclusion_count=len(data["pairs"]),
            neighbor_list_policy="fresh state energy for each displaced coordinate; native cutoff unchanged",
            minimum_image_boundary_margin_nm=min(p[0] for p in pairs),
            minimum_cutoff_margin_nm=min(p[1] for p in pairs),
            minimum_angle_sine=min(angles),
            improper_status="tested"
            if data.get("impropers")
            else "not_present_in_source_assignment",
            components=results,
            passed=all(x["passed"] for x in results.values()),
            **FLAGS,
        )
    )
