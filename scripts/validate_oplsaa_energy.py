"""Predeclared source-resolved OPLS parameter/energy acceptance, no MD."""

import argparse
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.forcefields.oplsaa import load_oplsaa_source, parameterize_oplsaa
from island.forcefields.oplsaa.adapter import check_foyer_installation
from island.workflows import storage

CASES = [
    ("ethane", "CC", 0),
    ("butane", "CCCC", 0),
    ("ethanol", "CCO", 0),
    ("dimethyl_ether", "COC", 0),
    ("benzene", "c1ccccc1", 0),
    ("pe3", "[*:1]CC[*:2]", 3),
    ("peo3", "[*:1]CCO[*:2]", 3),
    ("pe50", "[*:1]CC[*:2]", 50),
    ("ps3", "[*:1]CC(c1ccccc1)[*:2]", 3),
]
PLAN = {
    "schema": "island_opls_energy_acceptance_v1",
    "cases": CASES,
    "coordinate_seed": 20261003,
    "polymer_template_seed": 2026,
    "polymer_assembly_seed": 2026,
    "perturbation_seed": 71823,
    "perturbation_sigma_angstrom": [0.015, 0.03, 0.05],
    "pe50_frames": 1,
    "platform": "Reference",
    "energy_atol_kj_mol": 1e-5,
    "energy_rtol": 2e-10,
    "force_atol_kj_mol_angstrom": 1e-5,
    "force_rtol": 2e-10,
    "finite_difference_steps_angstrom": [1e-4, 5e-5],
    "fd_atol": 2e-3,
    "fd_rtol": 2e-5,
    "reference": "pinned Foyer.apply -> ParmEd geometric createSystem; no ISLAND resolution or pair helpers",
    "components": ["bond", "angle", "rb_proper", "nonbonded"],
    "expected_rejection": ["ps3"],
    "coordinate_nm_factor": 0.1,
    "force_angstrom_factor": 0.1,
    "outer_retries": 0,
    "typing_check": "exact stable-ID types against independent Foyer.apply",
}


def make_system(smiles, dp):
    if dp:
        from island.builders import build_linear_polymer

        return build_linear_polymer(
            smiles,
            dp=dp,
            coordinate_method="local_templates",
            template_seed=2026,
            assembly_seed=2026,
        )
    from island.chemistry import from_smiles

    return from_smiles(smiles, random_seed=20261003)


def reference(system, xml):
    """Independent upstream force selection and geometric mixing implementation."""
    import openmm as mm
    import parmed
    from foyer import Forcefield
    from openmm import app

    p = parmed.Structure()
    ids = sorted(system.topology.sites)
    mapping = {}
    for sid in ids:
        atom = system.topology.sites[sid]
        a = parmed.Atom(
            name=atom.element, atomic_number=atom.atomic_number, mass=atom.mass
        )
        p.add_atom(a, "M", 1)
        mapping[sid] = a
    for b in system.topology.bonds.values():
        p.bonds.append(parmed.Bond(mapping[b.site1], mapping[b.site2], order=b.order))
    p.coordinates = [system.coordinates.get(s) for s in ids]
    ff = Forcefield(forcefield_files=str(xml))
    p = ff.apply(
        p,
        use_residue_map=False,
        assert_bond_params=True,
        assert_angle_params=True,
        assert_dihedral_params=True,
        assert_improper_params=False,
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        removeCMMotion=False,
    )
    assert p.combining_rule == "geometric"
    model = p.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        removeCMMotion=False,
    )
    assert model.getNumConstraints() == 0 and not any(
        model.isVirtualSite(i) for i in range(len(ids))
    )
    for i, s in enumerate(ids):
        model.setParticleMass(i, system.topology.sites[s].mass)
    group = {
        mm.HarmonicBondForce: 0,
        mm.HarmonicAngleForce: 1,
        mm.RBTorsionForce: 2,
        mm.NonbondedForce: 3,
        mm.CustomNonbondedForce: 3,
    }
    nb = next(f for f in model.getForces() if type(f) is mm.NonbondedForce)
    custom = next(f for f in model.getForces() if type(f) is mm.CustomNonbondedForce)
    # Inspect actual upstream model. ParmEd stores sqrt(sigma), 2*sqrt(epsilon).
    assert "sigc=sigma1*sigma2" in custom.getEnergyFunction()
    assert custom.getNonbondedMethod() == mm.CustomNonbondedForce.NoCutoff
    for force in model.getForces():
        assert type(force) in group
        force.setForceGroup(group[type(force)])
        if type(force) is mm.NonbondedForce:
            assert force.getNonbondedMethod() == mm.NonbondedForce.NoCutoff
            force.setUseDispersionCorrection(False)
            force.setUseSwitchingFunction(False)
        if type(force) is mm.CustomNonbondedForce:
            force.setUseLongRangeCorrection(False)
            force.setUseSwitchingFunction(False)
    from openmm import unit

    for i in range(len(ids)):
        assert (
            nb.getParticleParameters(i)[2].value_in_unit(unit.kilojoule_per_mole) == 0
        )
        eps, sigma = custom.getParticleParameters(i)
        assert np.isclose(sigma**2, p.atoms[i].sigma * 0.1)
        assert np.isclose(eps**2 / 4, p.atoms[i].epsilon * 4.184)
    # Verify every LJ exception uses geometric sigma, not the raw Foyer LB value.
    for i in range(nb.getNumExceptions()):
        a, b, _q, sigma, eps = nb.getExceptionParameters(i)
        if eps.value_in_unit(unit.kilojoule_per_mole):
            assert np.isclose(
                sigma.value_in_unit(unit.nanometer),
                np.sqrt(p.atoms[a].sigma * p.atoms[b].sigma) * 0.1,
            )
    return model, {
        "forces": [type(f).__name__ for f in model.getForces()],
        "lj_expression": custom.getEnergyFunction(),
        "exceptions": nb.getNumExceptions(),
        "types": [atom.type for atom in p.atoms],
    }


def reference_value(model, xyz):
    import openmm as mm
    from openmm import unit

    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        context.setPositions(xyz * 0.1)
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces = (
            state.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer
            )
            * 0.1
        )
        components = {
            name: context.getState(getEnergy=True, groups=1 << i)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
            for i, name in enumerate(PLAN["components"])
        }
        return energy, forces, components
    finally:
        del context, integrator


def execute(root, xml, plan):
    source = load_oplsaa_source(xml)
    rng = np.random.default_rng(71823)
    rows = []
    for name, smiles, dp in CASES:
        row = {"case": name, "status": "failed"}
        start = time.perf_counter()
        try:
            system = make_system(smiles, dp)
            before = deepcopy(system.to_dict())
            try:
                parameters = parameterize_oplsaa(system, source)
            except Exception as error:
                if name == "ps3" and "Native component neutrality failed" in str(error):
                    row.update(status="expected_rejection", failure=str(error))
                    continue
                raise
            if name == "ps3":
                raise ValueError("Expected native neutrality rejection did not occur")
            parameters.validate_integrity(system, source)
            storage.publish(
                root / f"{name}-parameters.json", parameters.json_text.encode()
            )
            from island.workflows.bundle import system_data

            storage.publish(
                root / f"{name}-system.json",
                storage.json_bytes(storage.encode(system_data(system))),
            )
            snapshot = parameters.to_parameterized_system(system, source)
            evaluator = OPLSSinglePointEvaluator.from_parameterized_system(
                snapshot, source
            )
            ref, inspection = reference(system, xml)
            ids = sorted(system.topology.sites)
            assigned_sites = parameters.payload["resolved"]["sites"]
            assert inspection["types"] == [
                assigned_sites[sid]["type"] for sid in ids
            ], "Independent upstream type mapping differs"
            base = np.array([system.coordinates.get(s) for s in ids])
            frames = []
            fd = []
            for frame, amplitude in enumerate(
                PLAN["perturbation_sigma_angstrom"][: 1 if name == "pe50" else 3]
            ):
                xyz = base + rng.normal(0, amplitude, base.shape)
                coords = dict(zip(ids, xyz.tolist(), strict=True))
                result = evaluator.evaluate(coords)
                energy, forces, parts = reference_value(ref, xyz)
                actual_parts = {
                    k: result.energy_components[k]
                    for k in ("bond", "angle", "rb_proper")
                }
                actual_parts["nonbonded"] = (
                    result.energy_components["lj"] + result.energy_components["coulomb"]
                )
                diff = np.array([result.forces[s] for s in ids]) - forces
                assert np.isclose(
                    result.potential_energy, energy, atol=1e-5, rtol=2e-10
                )
                assert np.allclose(
                    [result.forces[s] for s in ids], forces, atol=1e-5, rtol=2e-10
                )
                assert all(
                    np.isclose(actual_parts[k], parts[k], atol=1e-5, rtol=2e-10)
                    for k in parts
                )
                frames.append(
                    {
                        "coordinate_fingerprint": result.coordinate_fingerprint,
                        "energy": result.potential_energy,
                        "energy_abs_error": abs(result.potential_energy - energy),
                        "force_max_abs_error": float(abs(diff).max()),
                        "component_abs_errors": {
                            k: abs(actual_parts[k] - parts[k]) for k in parts
                        },
                    }
                )
                if frame == 0 and name in ("butane", "ethanol", "benzene"):
                    for site_index, axis in ((0, 0), (len(ids) // 2, 1)):
                        errors = []
                        for h in PLAN["finite_difference_steps_angstrom"]:
                            energies = []
                            for sign in (-1, 1):
                                trial = xyz.copy()
                                trial[site_index, axis] += sign * h
                                energies.append(
                                    evaluator.evaluate(
                                        dict(zip(ids, trial.tolist(), strict=True))
                                    ).potential_energy
                                )
                            numeric = -(energies[1] - energies[0]) / (2 * h)
                            analytic = result.forces[ids[site_index]][axis]
                            assert np.isclose(numeric, analytic, atol=0.002, rtol=2e-5)
                            errors.append(abs(numeric - analytic))
                        fd.append(
                            {
                                "site": ids[site_index],
                                "axis": axis,
                                "absolute_errors": errors,
                            }
                        )
            assert system.to_dict() == before
            row.update(
                status="passed",
                sites=len(ids),
                coverage=parameters.payload["resolved"]["coverage"],
                parameter_identity=parameters.identity,
                reference_inspection=inspection,
                frames=frames,
                finite_differences=fd,
            )
        except Exception as error:  # noqa: BLE001 -- truthful bounded acceptance gate
            row["failure"] = f"{type(error).__name__}: {error}"
        finally:
            row["seconds"] = time.perf_counter() - start
            rows.append(row)
            print(name, row["status"], row.get("failure", ""), flush=True)
            storage.publish(
                root / "outcomes.json",
                storage.json_bytes(
                    {
                        "declaration": plan,
                        "cases": rows,
                        "complete": len(rows) == len(CASES)
                        and all(
                            r["status"]
                            == (
                                "expected_rejection" if r["case"] == "ps3" else "passed"
                            )
                            for r in rows
                        ),
                        "production_validated": False,
                        "simulation_readiness": "not_established",
                    }
                ),
                replace=True,
            )
    return all(
        r["status"] == ("expected_rejection" if r["case"] == "ps3" else "passed")
        for r in rows
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute-declared", action="store_true")
    args = parser.parse_args()
    plan = {
        **PLAN,
        "source": load_oplsaa_source(args.xml).identity,
        "dependencies": check_foyer_installation()[3],
    }
    plan = storage.decode(storage.encode(plan))
    if not args.execute_declared:
        args.output.mkdir(parents=True, exist_ok=False)
        storage.publish(
            args.output / "declaration.json", storage.json_bytes(storage.encode(plan))
        )
        return 0
    assert storage.decode(storage.read_json(args.output / "declaration.json")) == plan
    assert {p.name for p in args.output.iterdir()} == {"declaration.json"}
    return 0 if execute(args.output, args.xml, plan) else 1


if __name__ == "__main__":
    raise SystemExit(main())
