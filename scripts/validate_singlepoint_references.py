"""Independent whole-system acceptance using retained Phase 4D2.1 source files."""

import argparse
import hashlib
import json
from copy import deepcopy
from dataclasses import replace
from math import fsum
from pathlib import Path

import numpy as np

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.evaluation import OpenMMSinglePointEvaluator
from island.exceptions import EvaluationInputError
from island.forcefields import import_amber_prmtop
from island.forcefields.ambertools.models import (
    LEGACY_PREPARATION_SCHEMA,
    PREPARATION_SCHEMA,
    AmberToolsPreparationResult,
    digest,
)


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


def restore_preparation_record(raw):
    """Restore only v2/v3 schema-defined integer keys; never rewrite signed values."""
    record = deepcopy(raw)
    if record.get("schema") not in (PREPARATION_SCHEMA, LEGACY_PREPARATION_SCHEMA):
        raise EvaluationInputError("Unsupported archived preparation schema")

    def integer_keys(mapping):
        restored = {}
        for key, value in mapping.items():
            if not isinstance(key, str) or str(int(key)) != key:
                raise EvaluationInputError(
                    "Noncanonical integer mapping key in archive"
                )
            restored[int(key)] = value
        return restored

    for name in ("typed_mol2_index_to_site_id", "prmtop_index_to_site_id"):
        record["lineage"][name] = integer_keys(record["lineage"][name])
    record["expected_cip_by_site"] = integer_keys(record["expected_cip_by_site"])
    # Coordinate keys and generated atom-name keys are strings in the signed schema.
    return record


def _child(directory, name):
    if not isinstance(name, str) or Path(name).name != name or name in (".", ".."):
        raise EvaluationInputError("Archive paths must be single relative names")
    return directory / name


def load_archived_preparation(root, row):
    """Reconstruct and verify the existing signatures, never sign a replacement."""
    record = restore_preparation_record(row["record"])
    if digest(record) != row["record_signature"]:
        raise EvaluationInputError("Archived preparation record signature mismatch")
    directory = _child(root, row["artifact_subdirectory"])
    required = {
        **record["artifact_sha256"],
        "input.mol2": record["input_mol2_sha256"],
        "lineage.json": record["input_lineage_sha256"],
    }
    retained = row["retained_file_sha256"]
    if "input_system.json" not in retained:
        raise EvaluationInputError("Archive lacks authoritative input-system checksum")
    for name, checksum in required.items():
        if retained.get(name) != checksum:
            raise EvaluationInputError(
                f"{name} checksum contradicts signed preparation"
            )
    for name, checksum in retained.items():
        path = _child(directory, name)
        if hashlib.sha256(path.read_bytes()).hexdigest() != checksum:
            raise EvaluationInputError(f"{name} checksum mismatch")
    payload = json.loads((directory / "input_system.json").read_text())
    if payload.get("representation") != "atomistic" or payload.get("box") is not None:
        raise EvaluationInputError(
            "Archive requires the original nonperiodic atomistic input"
        )
    graph = Topology()
    for site in payload["sites"]:
        if site.get("kind") != "AtomSite":
            raise EvaluationInputError("Archived input contains a non-atomistic site")
        graph.add_site(
            AtomSite(**{key: value for key, value in site.items() if key != "kind"})
        )
    for bond in payload["bonds"]:
        graph.add_bond(**bond)
    system = MolecularSystem(
        graph,
        Coordinates({int(k): v for k, v in payload["coordinates"].items()}),
        metadata=payload["metadata"],
    )
    imported = import_amber_prmtop(
        system,
        directory / "result.prmtop",
        record["lineage"]["prmtop_index_to_site_id"],
        source=f"AmberTools {record['requested_force_field']} generated prmtop",
        force_field=record["requested_force_field"],
        charge_method="provided"
        if record["charge_method"] == "provided"
        else "AM1-BCC",
        charge_tolerance=record["charge_validation_tolerance_e"],
    )
    if (
        row["imported_signature"] != record["imported_result_signature"]
        or imported.source_sha256 != row["source_sha256"]
    ):
        raise EvaluationInputError(
            "Archived imported signature/source checksum mismatch"
        )
    shared = {
        key: value
        for key, value in record.items()
        if key != "imported_result_signature"
    }
    imported = replace(
        imported,
        provenance={**dict(imported.provenance), "ambertools_preparation": shared},
        result_signature=record["imported_result_signature"],
    )
    preparation = AmberToolsPreparationResult(imported, record, row["record_signature"])
    # This validates reconstructed parameter content against the historical hash,
    # then shared provenance, original coordinates, mapping and stage semantics.
    # A different parser/converter may fail; do not re-sign to hide incompatibility.
    preparation.validate_integrity(system)
    return system, preparation, directory


def charge_report(system, preparation, directory):
    imported = preparation.imported_result
    lines = (directory / "typed.mol2").read_text().splitlines()
    begin, end = lines.index("@<TRIPOS>ATOM") + 1, lines.index("@<TRIPOS>BOND")
    typed_total = fsum(
        float(line.split()[8]) for line in lines[begin:end] if line.strip()
    )
    formal = sum(site.formal_charge for site in system.topology.sites.values())
    return {
        "charge_unit": "elementary_charge",
        "charge_tolerance_e": preparation.record["charge_validation_tolerance_e"],
        "formal_charge_e": formal,
        "typed_mol2_total_charge_e": typed_total,
        "typed_mol2_charge_residual_e": typed_total - formal,
        "final_total_charge_e": imported.charge_result.total_assigned_charge,
        "final_charge_residual_e": imported.charge_result.total_charge_residual,
        "charges_modified": False,
        "residual_cause": "not established; reported totals are measured from checksum-verified files",
    }


def finite_difference_report(evaluator, coordinates, result):
    """Selected Cartesian derivatives in angstroms, with two independent steps."""
    site = min(coordinates)
    checks = []
    for axis in range(3):
        previous = None
        for step in (1e-4, 5e-5):
            left, right = deepcopy(coordinates), deepcopy(coordinates)
            left[site][axis] -= step
            right[site][axis] += step
            force = -(
                evaluator.evaluate(right).potential_energy
                - evaluator.evaluate(left).potential_energy
            ) / (2 * step)
            np.testing.assert_allclose(
                force, result.forces[site][axis], rtol=2e-6, atol=2e-5
            )
            if previous is not None:
                np.testing.assert_allclose(force, previous, rtol=2e-6, atol=2e-5)
            previous = force
            checks.append(
                {
                    "site_id": site,
                    "axis": axis,
                    "step_angstrom": step,
                    "finite_difference_force": force,
                    "backend_force": result.forces[site][axis],
                    "absolute_error": abs(force - result.forces[site][axis]),
                    "force_unit": "kJ/(mol*angstrom)",
                }
            )
    return checks


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
        system, preparation, directory = load_archived_preparation(
            args.preparation_manifest.parent, row
        )
        imported = preparation.imported_result
        xyz = {site: system.coordinates.get(site) for site in system.topology.sites}
        mapping = dict(imported.mapping)
        prmtop = directory / "result.prmtop"
        evaluator = OpenMMSinglePointEvaluator.from_preparation(system, preparation)
        charge_evidence = charge_report(system, preparation, directory)
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
                    "historical_preparation_signature": preparation.record_signature,
                    "historical_imported_signature": imported.result_signature,
                    "historical_signatures_verified": True,
                    "charge_report": charge_evidence,
                    "finite_difference_checks": finite_difference_report(
                        evaluator, coordinates, result
                    )
                    if frame == 0
                    else [],
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
