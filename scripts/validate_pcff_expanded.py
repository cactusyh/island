"""Declared partial J1 milestone against real pinned source and LAMMPS.

No case is silently removed. --require-full-source keeps the full-source gate
nonzero while any source row/rule/family remains unresolved.
"""

import argparse
import time
from collections import Counter
from decimal import Decimal
from pathlib import Path

import numpy as np
from validate_pcff_singlepoint import (
    GROUPS,
    TOL,
    aa_reference_coefficients,
    check_inventory,
    data_sections,
    digest,
    lammps_input,
    reference_inputs,
    run,
)

from island.charge_references.records import unpack
from island.chemistry import from_smiles
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    assign_pcff_source_types,
    define_pcff_model,
    load_pcff_source,
    save_pcff_automatic_record,
    save_pcff_model,
    save_pcff_parameters,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.coverage import (
    attach_fixture_evidence,
    pcff_coverage_ledger,
)
from island.forcefields.pcff.expanded import PROFILE_NAME
from island.forcefields.pcff.source import PIN
from island.workflows.storage import json_bytes, publish

# Independently authored expected labels follow the explicit SMILES heavy-atom
# traversal in these fixed input fixtures. H labels below use the documented
# parent classes, not the production environment/typing functions.
CASES = [
    ("cyclohexane", "C1CCCCC1", ["c2"] * 6, True),
    ("cyclic_ether", "C1COCCO1", ["c2", "c2", "oc", "c2", "c2", "oc"], True),
    ("epoxide", "C1CO1", ["c3h", "c3h", "o3e"], True),
    ("benzene", "c1ccccc1", ["cp"] * 6, True),
    ("ethene", "C=C", ["c="] * 2, True),
    ("methanethiol", "CS", ["c3", "sh"], True),
    ("acetal", "CCOC(C)OCC", ["c3", "c2", "oc", "coh", "c3", "oc", "c2", "c3"], False),
    ("acetone", "CC(=O)C", ["c3", "c_0", "o_1", "c3"], False),
    ("methylamine", "CN", ["c3", "na"], False),
    ("chloromethane", "CCl", ["c3", "cl"], False),
    ("amino_alcohol", "CCNCCO", ["c3", "c2", "na", "c2", "c2", "oh"], False),
    ("ester", "CC(=O)OC", ["c3", "c_1", "o_1", "o_2", "c3"], False),
    ("acetate", "CC(=O)[O-]", ["c3", "c-", "o-", "o-"], False),
    ("tetramethylammonium", "C[N+](C)(C)C", ["c3", "n4", "c3", "c3", "c3"], False),
]
H_TYPES = {"oh": "ho", "na": "hn2", "sh": "hs", "n4": "hn"}


def independent_inventory(raw):
    counts = Counter()
    atom_labels = set()
    sections = []
    current = None
    for n, line in enumerate(raw.decode().splitlines(), 1):
        words = line.split("!", 1)[0].split()
        if not words:
            continue
        if words[0].startswith("#"):
            current = (words[0][1:], words[1] if len(words) > 1 else "")
            sections.append(
                {"family": current[0], "namespace": current[1], "line": n, "rows": 0}
            )
        elif current and words[0][0].isdigit():
            counts[current] += 1
            if current == ("atom_types", "cff91"):
                atom_labels.add(words[2])
            sections[-1]["rows"] += 1
    expected = {
        "atom_types": 133,
        "equivalence": 134,
        "auto_equivalence": 108,
        "bond_increments": 564,
    }
    actual = {
        k: sum(v for (name, _), v in counts.items() if name == k) for k in expected
    }
    if actual != expected or len(atom_labels) != 133:
        raise ValueError("Independent inventory mismatch")
    return {
        "counts": actual,
        "distinct_atom_types": len(atom_labels),
        "sections": sections,
    }


def labels_for(system, heavy):
    ids = sorted(system.topology.sites)
    heavy_ids = [i for i in ids if system.topology.sites[i].element != "H"]
    if len(heavy_ids) != len(heavy):
        raise ValueError("Declared fixture inventory mismatch")
    labels = dict(zip(heavy_ids, heavy))
    for i in ids:
        if i in labels:
            continue
        neighbors = [
            b.site2 if b.site1 == i else b.site1
            for b in system.topology.bonds.values()
            if i in b.key
        ]
        if len(neighbors) != 1:
            raise ValueError("Fixture H count mismatch")
        labels[i] = H_TYPES.get(labels[neighbors[0]], "hc")
    return labels


def independent_charges(raw, system, labels):
    """Decimal endpoint sum from independently scanned direct/ordinary FRC rows.

    msi2lmp does not assign charges; its supplied charges are checked here first.
    Missing or conflicting source rows are failures, not synthetic zero charges.
    """
    equiv = {}
    increments = {}
    section = None
    for line in raw.decode().splitlines():
        w = line.split("!", 1)[0].split()
        if not w:
            continue
        if w[0].startswith("#"):
            section = w[0][1:]
            continue
        if not w[0][0].isdigit():
            continue
        if section == "equivalence" and (
            w[2] not in equiv or Decimal(w[0]) > equiv[w[2]][0]
        ):
            equiv[w[2]] = (Decimal(w[0]), w[4])
        if section == "bond_increments":
            key = tuple(w[2:4])
            value = (Decimal(w[0]), tuple(map(Decimal, w[4:6])))
            if key not in increments or value[0] > increments[key][0]:
                increments[key] = value
            elif value[0] == increments[key][0] and value[1] != increments[key][1]:
                raise ValueError("Ambiguous independent source rows")
    total = {i: Decimal(0) for i in labels}
    selected = []
    for a, b in sorted(system.topology.bonds):
        direct = (labels[a], labels[b])
        ordinary = tuple(equiv[t][1] for t in direct)
        values = None
        for query in (direct, ordinary):
            if query in increments:
                values = increments[query][1]
                break
            if query[::-1] in increments:
                values = increments[query[::-1]][1][::-1]
                break
        if values is None:
            raise ValueError(f"Independent missing charge row: {direct}")
        total[a] += values[0]
        total[b] += values[1]
        selected.append(
            {"sites": [a, b], "query": query, "increments": list(map(str, values))}
        )
    return {i: float(v) for i, v in total.items()}, selected


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--converter", required=True)
    p.add_argument("--lammps", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args(argv)
    out = Path(a.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    source_path = Path(a.source).resolve()
    converter = Path(a.converter).resolve()
    lammps = Path(a.lammps).resolve()
    declaration = {
        "schema": "island_pcff_j1_acceptance_v1",
        "source": PIN,
        "profile": PROFILE_NAME,
        "cases": CASES,
        "seed": 2026,
        "coordinates": "RDKit fixed seed; no minimization",
        "perturbation": "0.015*sin(arange(3*N).reshape(N,3)+0.3) angstrom",
        "policies": {"lj": [0, 0, 1], "coulomb": [0, 0, 1]},
        "charge_atol_e": 1e-12,
        "tolerances": TOL,
        "timeout_seconds": 180,
        "converter_sha256": digest(converter),
        "lammps_sha256": digest(lammps),
        "reference_AA_correction": "H5 retained explicit ABC,ABD,CBD roles from independent converter angles",
        "gate": "declared required numerical cases; full-source is separate and remains unmet",
    }
    publish(out / "declaration.json", json_bytes(declaration))
    source = load_pcff_source(source_path)
    publish(
        out / "independent-inventory.json",
        json_bytes(independent_inventory(source.raw)),
    )
    outcomes = []
    for name, smiles, expected, required in CASES:
        work = out / name
        work.mkdir()
        start = time.perf_counter()
        row = {
            "name": name,
            "required_numerical": required,
            "status": "failed",
            "stages": {},
        }
        try:
            system = from_smiles(smiles, random_seed=2026)
            labels = labels_for(system, expected)
            typing = type_pcff_atoms(system, source, profile=PROFILE_NAME)
            save_pcff_automatic_record(typing, work / "typing.json")
            if typing.assignments != labels:
                raise ValueError("Independent expected label mismatch")
            row["types"] = sorted(set(labels.values()))
            row["stages"]["typing"] = True
            explicit = assign_pcff_source_types(
                system,
                source,
                labels,
                provenance="Declared manual fixture labels from source descriptions, separately checked against graph",
            )
            charges = assign_automatic_pcff_charges(system, typing)
            save_pcff_automatic_record(charges, work / "charges.json")
            row["stages"]["charges"] = charges.complete
            row["charge_diagnostics"] = charges.payload["native_charge_record"][
                "diagnostics"
            ]
            if not charges.complete:
                row["status"] = "charge_coverage_missing"
                continue
            q, contrib = independent_charges(source.raw, system, labels)
            publish(
                work / "independent-charges.json",
                json_bytes({"charges": q, "contributions": contrib}),
            )
            if any(abs(q[i] - charges.charges[i]) > 1e-12 for i in q):
                raise ValueError("Independent charge mismatch")
            if (
                charges.charges
                != assign_automatic_pcff_charges(system, explicit).charges
            ):
                raise ValueError("Explicit/automatic divergence")
            assignment = assign_pcff_parameters(system, typing, charges)
            save_pcff_parameters(assignment, work / "parameters.json")
            spec = define_pcff_model(
                assignment,
                special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
            )
            save_pcff_model(spec, work / "model.json")
            payload = unpack(spec.json_text)
            row["stages"]["model"] = payload["model_definition_complete"]
            row["model_diagnostics"] = payload["diagnostics"]
            row["coverage"] = unpack(assignment.json_text)["coverage"]
            if not payload["model_definition_complete"]:
                row["status"] = "parameter_coverage_missing"
                continue
            ids = sorted(labels)
            xyz = np.array([system.coordinates.get(i) for i in ids])
            ids, neighbors = reference_inputs(system, labels, q, xyz, work)
            run(
                [
                    converter,
                    "reference",
                    "-class",
                    "II",
                    "-frc",
                    source_path,
                    "-p",
                    "3",
                    "-nocenter",
                ],
                work,
                "converter",
            )
            sections = data_sections(work / "reference.data")
            row["inventory"] = check_inventory(sections, ids, neighbors, q)
            aa = aa_reference_coefficients(sections)
            publish(work / "aa-overrides.json", json_bytes(aa))
            evaluator = PCFFSinglePointEvaluator(system, spec)
            row["comparisons"] = []
            for k, frame in enumerate(
                (
                    xyz,
                    xyz + 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3),
                )
            ):
                directory = work / f"frame-{k}"
                directory.mkdir()
                (directory / "reference.data").write_bytes(
                    (work / "reference.data").read_bytes()
                )
                (directory / "input.lmp").write_text(
                    lammps_input(frame, [0, 0, 1], [0, 0, 1], aa)
                )
                run([lammps, "-in", "input.lmp"], directory, "lammps")
                text = (directory / "lammps.log").read_text()
                line = [v for v in text.splitlines() if v.startswith("ENERGIES ")][-1]
                ref = np.array(list(map(float, line.split()[1:]))) * 4.184
                dump = (
                    (directory / "forces.dump")
                    .read_text()
                    .split("ITEM: ATOMS id x y z fx fy fz\n")[1]
                )
                forces = (
                    np.array(
                        [
                            list(map(float, l.split()[4:7]))
                            for l in dump.splitlines()
                            if l.strip()
                        ]
                    )
                    * 4.184
                )
                coords = dict(zip(ids, frame))
                result = evaluator.evaluate(coords)
                groups = {**GROUPS, "improper": ["angle-angle", "wilson_out_of_plane"]}
                actual = np.array(
                    [result.potential_energy]
                    + [
                        sum(result.energy_components[x] for x in terms)
                        for terms in groups.values()
                    ]
                )
                actual_forces = np.array([result.forces[i] for i in ids])
                delta = {
                    "energy_max_abs": float(np.max(np.abs(actual - ref))),
                    "forces_max_abs": float(np.max(np.abs(actual_forces - forces))),
                }
                row["comparisons"].append(delta)
                np.testing.assert_allclose(
                    actual, ref, atol=TOL["energy_atol_kj_mol"], rtol=TOL["rtol"]
                )
                np.testing.assert_allclose(
                    actual_forces,
                    forces,
                    atol=TOL["force_atol_kj_mol_angstrom"],
                    rtol=TOL["rtol"],
                )
                with evaluator.open_session() as session:
                    repeated = session.evaluate(coords)
                    np.testing.assert_allclose(
                        [repeated.forces[i] for i in ids],
                        actual_forces,
                        atol=1e-10,
                        rtol=1e-12,
                    )
            row["stages"]["independent_numerical"] = True
            row["status"] = "numerically_verified"
        except Exception as error:  # noqa: BLE001 - durable per-case failure evidence
            row["error"] = str(error)
        finally:
            row["seconds"] = time.perf_counter() - start
            row["hashes"] = {
                str(p.relative_to(work)): digest(p)
                for p in work.rglob("*")
                if p.is_file()
            }
            outcomes.append(row)
            publish(work / "outcome.json", json_bytes(row))
            print(name, row["status"], flush=True)
    ledger = attach_fixture_evidence(pcff_coverage_ledger(source), outcomes)
    publish(out / "coverage.json", json_bytes(ledger))
    passed = all(
        r["status"] == "numerically_verified"
        for r in outcomes
        if r["required_numerical"]
    )
    report = {
        "cases": outcomes,
        "declared_milestone_passed": passed,
        "full_source_complete": False,
        "source_unchanged": digest(source_path) == PIN["sha256"],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(out / "report.json", json_bytes(report))
    return (
        0 if passed and report["source_unchanged"] and not a.require_full_source else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
