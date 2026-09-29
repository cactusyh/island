"""Independent whole-system acceptance using retained Phase 4D2.1 source files."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.evaluation import OpenMMSinglePointEvaluator
from island.forcefields import import_amber_prmtop


def reference(prmtop, mapping, coordinates):
    """Path B reads original Amber data through OpenMM's own reader."""
    import openmm as mm
    from openmm import app, unit

    system = app.AmberPrmtopFile(str(prmtop)).createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        rigidWater=False,
        implicitSolvent=None,
        removeCMMotion=False,
        hydrogenMass=None,
    )
    groups = {
        "HarmonicBondForce": 0,
        "HarmonicAngleForce": 1,
        "PeriodicTorsionForce": 2,
        "NonbondedForce": 3,
    }
    for force in system.getForces():
        force.setForceGroup(groups[type(force).__name__])
        if isinstance(force, mm.NonbondedForce):
            force.setUseSwitchingFunction(False)
            force.setUseDispersionCorrection(False)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(system, integrator, mm.Platform.getPlatformByName("Reference"))
    try:
        context.setPositions(
            np.array([coordinates[mapping[i]] for i in range(len(mapping))])
            * 0.1
            * unit.nanometer
        )
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        force_array = (
            state.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer
            )
            * 0.1
        )
        components = {
            name: context.getState(getEnergy=True, groups=1 << i)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
            for i, name in enumerate(("bond", "angle", "torsion", "nonbonded"))
        }
        return (
            energy,
            {mapping[i]: force_array[i] for i in range(len(mapping))},
            components,
        )
    finally:
        del context, integrator


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Refusing to overwrite existing acceptance output")
    manifest = json.loads(args.preparation_manifest.read_text())
    records = []
    for row in manifest["cases"]:
        directory = args.preparation_manifest.parent / row["artifact_subdirectory"]
        for name in ("input_system.json", "result.prmtop"):
            if (
                hashlib.sha256((directory / name).read_bytes()).hexdigest()
                != row["retained_file_sha256"][name]
            ):
                raise RuntimeError(f"{name} checksum mismatch")
        payload = json.loads((directory / "input_system.json").read_text())
        graph = Topology()
        for site in payload["sites"]:
            graph.add_site(
                AtomSite(**{key: value for key, value in site.items() if key != "kind"})
            )
        for bond in payload["bonds"]:
            graph.add_bond(**bond)
        xyz = {int(key): value for key, value in payload["coordinates"].items()}
        system = MolecularSystem(graph, Coordinates(xyz), metadata=payload["metadata"])
        mapping = {
            int(key): value
            for key, value in row["record"]["lineage"][
                "prmtop_index_to_site_id"
            ].items()
        }
        prmtop = directory / "result.prmtop"
        imported = import_amber_prmtop(
            system,
            prmtop,
            mapping,
            source=row["case"],
            charge_tolerance=row["record"]["charge_validation_tolerance_e"],
        )
        evaluator = OpenMMSinglePointEvaluator(system, imported)
        rng = np.random.default_rng(20260929)
        for frame in range(3):
            coordinates = {
                site: np.array(xyz[site]) + rng.normal(0, 0.025, 3)
                for site in sorted(xyz)
            }
            result = evaluator.evaluate(coordinates)
            energy, forces, components = reference(prmtop, mapping, coordinates)
            actual_components = dict(result.energy_components)
            actual_components["torsion"] = actual_components.pop(
                "proper_torsion"
            ) + actual_components.pop("periodic_improper")
            np.testing.assert_allclose(
                result.potential_energy, energy, rtol=2e-10, atol=2e-7
            )
            for name, value in components.items():
                np.testing.assert_allclose(
                    actual_components[name], value, rtol=2e-10, atol=2e-7
                )
            max_force_error = 0
            for site in sorted(xyz):
                np.testing.assert_allclose(
                    result.forces[site], forces[site], rtol=2e-9, atol=2e-6
                )
                max_force_error = max(
                    max_force_error,
                    float(np.max(np.abs(np.array(result.forces[site]) - forces[site]))),
                )
            records.append(
                {
                    "case": row["case"],
                    "frame": frame,
                    "source_prmtop_sha256": imported.source_sha256,
                    "coordinates_angstrom": {
                        site: coords.tolist() for site, coords in coordinates.items()
                    },
                    "energy_kj_mol": result.potential_energy,
                    "reference_energy_kj_mol": energy,
                    "components_kj_mol": dict(result.energy_components),
                    "reference_components_kj_mol": components,
                    "maximum_force_error_kj_mol_angstrom": max_force_error,
                    "coordinate_fingerprint": result.coordinate_fingerprint,
                    "parameter_fingerprint": result.parameter_fingerprint,
                    "model_fingerprint": result.model_fingerprint,
                    "evaluation_fingerprint": result.evaluation_fingerprint,
                    "backend_version": result.backend_version,
                    "platform": result.platform,
                    "settings": dict(result.settings),
                }
            )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(
            {
                "schema": "island_singlepoint_acceptance_v1",
                "preparation_manifest_sha256": hashlib.sha256(
                    args.preparation_manifest.read_bytes()
                ).hexdigest(),
                "generator_sha256": hashlib.sha256(
                    Path(__file__).read_bytes()
                ).hexdigest(),
                "cases": records,
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"{len(records)} independent whole-system comparisons passed; {args.output}")


if __name__ == "__main__":
    main()
