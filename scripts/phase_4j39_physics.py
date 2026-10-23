"""Pinned Foyer.apply/ParmEd independent periodic OPLS energy and force check."""

from pathlib import Path

import numpy as np

from scripts.phase_4j32_real_fixture import XML
from scripts.phase_4j38_contract import require
from scripts.validate_oplsaa_energy import reference, reference_value


def check(backend, frames):
    import openmm as mm

    from island._periodic_openmm import build
    from island.periodic import periodic_chemical_view

    model, witness = reference(periodic_chemical_view(backend.system), XML)
    # Independent upstream model, explicit full RESP vector. Rebuild Coulomb
    # exceptions from independent bond distances; retain Foyer/ParmEd LJ terms.
    from collections import deque

    from openmm import unit

    ids = sorted(backend.system.topology.sites)
    rows = backend.payload["parameter_assignment"]["native_payload"]["record"][
        "resp_charges"
    ]["rows"]
    q = {r["atom_id"]: float.fromhex(r["ieee754_hex"]) for r in rows}
    require(
        all(
            q[r["id"]].hex() == r["charge"].hex()
            for r in backend.payload["numerical_data"]["sites"]
        ),
        "Numerical charges differ from independently audited native vector",
    )
    adj = {i: set() for i in ids}
    for a, b in backend.system.topology.bonds:
        adj[a].add(b)
        adj[b].add(a)
    nb = next(f for f in model.getForces() if isinstance(f, mm.NonbondedForce))
    for j, sid in enumerate(ids):
        _, sigma, eps = nb.getParticleParameters(j)
        nb.setParticleParameters(j, q[sid], sigma, eps)
    for j in range(nb.getNumExceptions()):
        a, b, _, sigma, eps = nb.getExceptionParameters(j)
        queue = deque([(ids[a], 0)])
        seen = {ids[a]}
        distance = None
        while queue:
            current, depth = queue.popleft()
            if current == ids[b]:
                distance = depth
                break
            for nxt in adj[current] - seen:
                seen.add(nxt)
                queue.append((nxt, depth + 1))
        require(distance in (1, 2, 3), "Unexpected upstream exception topology")
        nb.setExceptionParameters(
            j,
            a,
            b,
            q[ids[a]]
            * q[ids[b]]
            * (0.5 if distance == 3 else 0.0)
            * unit.elementary_charge**2,
            sigma,
            eps,
        )
    p = backend.payload
    cfg = p["config"]
    lengths = cfg["box_lengths"]
    model.setDefaultPeriodicBoxVectors(
        *[
            mm.Vec3(*(lengths[j] / 10 if i == j else 0 for j in range(3)))
            for i in range(3)
        ]
    )
    for f in model.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setNonbondedMethod(mm.NonbondedForce.Ewald)
            f.setCutoffDistance(cfg["nonbonded_cutoff"] / 10)
            f.setEwaldErrorTolerance(cfg["pme_tolerance"])
            f.setExceptionsUsePeriodicBoundaryConditions(True)
        elif isinstance(f, mm.CustomNonbondedForce):
            f.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
            f.setCutoffDistance(cfg["nonbonded_cutoff"] / 10)
        else:
            f.setUsesPeriodicBoundaryConditions(True)
    production = build(p["numerical_data"], cfg, mm)
    for f in production.getForces():
        f.setForceGroup({"bond": 0, "angle": 1, "rb_proper": 2}.get(f.getName(), 3))
    snapshots = []
    for i in (0, len(frames) - 1):
        f = frames[i]
        box = [mm.Vec3(*v) for v in f["box_vectors_nm"]]
        model.setDefaultPeriodicBoxVectors(*box)
        production.setDefaultPeriodicBoxVectors(*box)
        a = reference_value(model, np.array(f["positions_nm"]) * 10)
        b = reference_value(production, np.array(f["positions_nm"]) * 10)
        error = abs(a[0] - b[0])
        force = float(np.max(abs(a[1] - b[1])))
        components = {k: abs(a[2][k] - b[2][k]) for k in a[2]}
        require(
            error <= 1e-5 and force <= 1e-5 and max(components.values()) <= 1e-5,
            f"Independent periodic reference mismatch: {error}, {force}, {components}",
        )
        snapshots.append(
            {
                "frame_index": i,
                "energy_error_kj_mol": error,
                "force_error_kj_mol_angstrom": force,
                "component_errors_kj_mol": components,
            }
        )
    return {
        "passed": True,
        "reference": "pinned Foyer.apply -> ParmEd OPLS terms with explicit exact-graph RESP charges; independent Coulomb exceptions; periodic Ewald vs ISLAND; same OpenMM Reference engine",
        "independent_engine": False,
        "lammps_kspace_executed": False,
        "source": str(Path(XML)),
        "witness": witness,
        "charge_origin": "RESP/Psi4",
        "native_opls_charges": False,
        "snapshots": snapshots,
        "energy_tolerance_kj_mol": 1e-5,
        "force_tolerance_kj_mol_angstrom": 1e-5,
    }
