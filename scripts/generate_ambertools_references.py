"""Reproduce real GAFF/GAFF2 fixtures; never substitutes mocked outputs."""

import argparse
import json
import platform
import sys
from hashlib import sha256
from importlib.metadata import version
from math import atan2, cos, isclose, pi
from pathlib import Path

import numpy as np

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields import (
    AmberToolsOptions,
    AmberToolsParameterizationEngine,
)


def phenol_charges(system):
    """Explicit software-test charges, not a claimed physical charge model."""
    charges = {site_id: 0.0 for site_id in system.topology.sites}
    oxygen = next(i for i, site in system.topology.sites.items()
                  if site.element == "O")
    hydroxyl_h = next(i for i in system.topology.neighbors(oxygen)
                      if system.topology.sites[i].element == "H")
    charges[oxygen], charges[hydroxyl_h] = -0.2, 0.2
    return charges


def _ordered_phi() -> float:
    """Nondegenerate signed torsion from four ordered Cartesian points."""
    p0, p1, p2, p3 = (
        np.array((0.2, 0.1, -0.1)), np.array((1.1, 0.0, 0.3)),
        np.array((1.5, 1.0, 0.2)), np.array((2.4, 1.1, 1.2)),
    )
    axis = p2 - p1
    axis /= np.linalg.norm(axis)
    first = p0 - p1
    last = p3 - p2
    first -= np.dot(first, axis) * axis
    last -= np.dot(last, axis) * axis
    return atan2(np.dot(np.cross(axis, first), last), np.dot(first, last))


def independent_conversion_checks(result, prmtop: Path) -> dict:
    """Compare source-side formulas with converted records; not an MD evaluator."""
    import parmed as pmd

    source = pmd.load_file(str(prmtop))
    imported = result.imported_result
    mapping = imported.mapping
    checks = {}

    def compare(family: str, source_value: float, converted_value: float) -> None:
        if not isclose(source_value, converted_value, rel_tol=2e-5, abs_tol=1e-6):
            raise RuntimeError(
                f"{family} source/converted energy mismatch: "
                f"{source_value} vs {converted_value}"
            )
        checks[family] = {"source_kj_mol": source_value,
                          "converted_kj_mol": converted_value}

    if source.bonds:
        bond = source.bonds[0]
        key = tuple(sorted((mapping[bond.atom1.idx], mapping[bond.atom2.idx])))
        parameter = imported.bond_assignments[key].parameter
        radius = bond.type.req + 0.11  # Å, deliberately away from equilibrium
        compare("bond", bond.type.k * 0.11**2 * 4.184,
                0.5 * parameter.force_constant *
                (radius / 10 - parameter.equilibrium_length) ** 2)
    if source.angles:
        angle = source.angles[0]
        key = tuple(mapping[a.idx] for a in (angle.atom1, angle.atom2, angle.atom3))
        key = min(key, key[::-1])
        parameter = imported.angle_assignments[key].parameter
        displacement = 11 * pi / 180
        compare("angle", angle.type.k * displacement**2 * 4.184,
                0.5 * parameter.force_constant * displacement**2)
    phi = _ordered_phi()
    def ordered_key(dihedral, improper: bool) -> tuple[int, ...]:
        sites = tuple(mapping[a.idx] for a in (
            dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4,
        ))
        return sites if improper else min(sites, sites[::-1])

    for improper, family in ((False, "proper"), (True, "improper")):
        source_term = next((d for d in source.dihedrals if d.improper == improper), None)
        if source_term is None:
            checks[family] = {"status": "not_present_in_source"}
            continue
        ordered = tuple(mapping[a.idx] for a in (
            source_term.atom1, source_term.atom2, source_term.atom3,
            source_term.atom4,
        ))
        key = ordered if improper else min(ordered, ordered[::-1])
        source_rows = [d for d in source.dihedrals
                       if d.improper == improper and ordered_key(d, improper) == key]
        parameter = (imported.improper_assignments if improper else
                     imported.proper_torsion_assignments)[key].parameter
        source_energy = sum(
            row.type.phi_k * (1 + cos(row.type.per * phi
                                      - row.type.phase * pi / 180)) * 4.184
            for row in source_rows
        )
        converted_energy = sum(
            term.force_constant * (1 + cos(term.periodicity * phi
                                           - term.phase * pi / 180))
            for term in parameter.terms
        )
        compare(family, source_energy, converted_energy)
        checks[family]["ordered_phi_radians"] = phi
        checks[family]["stable_site_order"] = list(ordered)
    positive = next((a for a in source.atoms if a.epsilon > 0), None)
    if positive is not None:
        index = positive.nb_idx
        ntypes = source.ptr("NTYPES")
        pair_index = source.parm_data["NONBONDED_PARM_INDEX"][
            (index - 1) * ntypes + index - 1
        ] - 1
        coefficient_a = source.parm_data["LENNARD_JONES_ACOEF"][pair_index]
        coefficient_b = source.parm_data["LENNARD_JONES_BCOEF"][pair_index]
        radius = 2 * positive.rmin * 1.08  # Å, non-minimum separation
        parameter = imported.site_assignments[mapping[positive.idx]].parameter
        source_lj = (coefficient_a / radius**12
                     - coefficient_b / radius**6) * 4.184
        target_lj = 4 * parameter.epsilon * (
            (parameter.sigma / (radius / 10)) ** 12
            - (parameter.sigma / (radius / 10)) ** 6
        )
        compare("lj", source_lj, target_lj)
    excluded = set(imported.source_exclusions)
    electrostatic_pair = next((
        (left, right) for left in source.atoms for right in source.atoms
        if left.idx < right.idx
        and tuple(sorted((mapping[left.idx], mapping[right.idx]))) not in excluded
        and abs(left.charge * right.charge) > 1e-8
    ), None)
    if electrostatic_pair is not None:
        left, right = electrostatic_pair
        source_coul = 332.06371 * left.charge * right.charge / 5.0 * 4.184
        target_coul = 138.935455864 * (
            imported.charge_result.assignments[mapping[left.idx]].charge
            * imported.charge_result.assignments[mapping[right.idx]].charge
        ) / 0.5
        compare("coulomb", source_coul, target_coul)
    else:
        checks["coulomb"] = {"status": "no_nonexcluded_nonzero_charge_pair"}
    checks["zero_lj_site_count"] = sum(
        atom.epsilon == 0 and atom.rmin == 0 for atom in source.atoms
    )
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / "references.json"
    if manifest.exists():
        parser.error(f"Refusing to overwrite existing {manifest}")
    phenol = from_smiles(
        "Oc1ccccc1", site_ids=(101 + 7 * i for i in range(13)),
        random_seed=2026,
    )
    pe = build_linear_polymer(
        "[*]CC[*]", dp=3, coordinate_method="local_templates",
        template_seed=2026, assembly_seed=2026,
    )
    chiral = from_smiles(
        "F[C@H](Cl)Br", site_ids=(301 + 11 * i for i in range(5)),
        random_seed=2026,
    )
    cases = (
        ("phenol_gaff_am1bcc", phenol,
         AmberToolsOptions("gaff", "am1bcc", work_root=output,
                           charge_tolerance=0.002,
                           retain_success_artifacts=True)),
        ("phenol_gaff2_am1bcc", phenol,
         AmberToolsOptions("gaff2", "am1bcc", work_root=output,
                           charge_tolerance=0.002,
                           retain_success_artifacts=True)),
        ("phenol_gaff2_provided", phenol,
         AmberToolsOptions("gaff2", "provided", phenol_charges(phenol),
                           work_root=output, retain_success_artifacts=True)),
        ("pe_dp3_gaff2_provided", pe,
         AmberToolsOptions("gaff2", "provided",
                           {site_id: 0.0 for site_id in pe.topology.sites},
                           work_root=output, retain_success_artifacts=True)),
        ("halomethane_gaff2_provided", chiral,
         AmberToolsOptions("gaff2", "provided",
                           {site_id: 0.0 for site_id in chiral.topology.sites},
                           work_root=output, retain_success_artifacts=True)),
    )
    records = []
    for name, system, options in cases:
        original = system.to_dict()
        result = AmberToolsParameterizationEngine().parameterize(system, options)
        result.validate_integrity(system)
        snapshot = result.to_parameterized_system(system)
        assert system.to_dict() == original
        assert snapshot.metadata["aggregate"]["production_validated"] is False
        assert snapshot.metadata["aggregate"]["simulation_readiness"] == "not_established"
        directory = Path(result.record["artifact_dir"])
        (directory / "input_system.json").write_text(
            json.dumps(original, indent=2, sort_keys=True) + "\n"
        )
        if name.startswith("phenol") and len(result.imported_result.improper_assignments) == 0:
            raise RuntimeError(f"{name} lacks the required periodic improper")
        checks = independent_conversion_checks(
            result, Path(result.record["artifact_dir"]) / "result.prmtop"
        )
        records.append({
            "case": name,
            "artifact_subdirectory": directory.name,
            "retained_file_sha256": {
                path.name: sha256(path.read_bytes()).hexdigest()
                for path in sorted(directory.iterdir()) if path.is_file()
            },
            "input_graph_site_count": len(system.topology.sites),
            "input_formal_charge": sum(
                site.formal_charge for site in system.topology.sites.values()
            ),
            "record_signature": result.record_signature,
            "imported_signature": result.imported_result.result_signature,
            "source_sha256": result.imported_result.source_sha256,
            "record": dict(result.record),
            "improper_count": len(result.imported_result.improper_assignments),
            "source_14_pair_count": len(result.imported_result.source_14_pairs),
            "source_exclusion_count": len(result.imported_result.source_exclusions),
            "expected_cip_by_site": result.record["expected_cip_by_site"],
            "independent_conversion_checks": checks,
        })
        print(name, result.record_signature)
    manifest.write_text(json.dumps({
        "schema": "island_real_ambertools_references_v1",
        "status": "generated_with_actual_ambertools",
        "generator_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "python": {"executable": sys.executable, "version": platform.python_version()},
        "python_packages": {name: version(name) for name in ("numpy", "rdkit", "parmed")},
        "source_sha256": {
            str(path.relative_to(Path(__file__).parents[1])):
                sha256(path.read_bytes()).hexdigest()
            for path in sorted((Path(__file__).parents[1] / "src/island").rglob("*.py"))
            if ".ipynb_checkpoints" not in path.parts
        },
        "cases": records,
        "scientific_validation": False,
        "simulation_readiness": "not_established",
    }, indent=2, sort_keys=True) + "\n")
    print("manifest:", manifest)


if __name__ == "__main__":
    main()
