"""Independent source and LAMMPS checks for the declared final-polymer matrix.

Uses retained systems/assignments, without construction, typing, or preparation.
The reference selects raw FRC rows independently and converts its own topology.
"""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from pcff_j5_reference import raw_inventory, read_source, resolve_raw
from validate_pcff_expanded import independent_charges
from validate_pcff_profile import numerical
from validate_pcff_singlepoint import (
    aa_reference_coefficients,
    check_inventory,
    data_sections,
    reference_inputs,
    run,
)

from island.charge_references.records import pack, unpack
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields import adopt_forcefield
from island.forcefields.pcff import (
    PCFFClass2Result,
    PCFFFinalGraphPreparation,
    define_pcff_model,
    load_pcff_source,
    select_pcff_profile,
    special_pair_policy,
)
from island.forcefields.pcff.registry import LINKED_BENZENOID
from island.workflows import storage
from island.workflows.bundle import system_from


def independent_selections(system, model, source_path):
    # Authored from the declared ring/methyl graph domain and source rows cp/c3/hc.
    # No production typer or repeat labels participate in this reference.
    labels = {
        i: "hc" if a.element == "H" else "cp" if a.metadata.get("aromatic") else "c3"
        for i, a in system.topology.sites.items()
    }
    raw = read_source(source_path.read_bytes())
    assignment = unpack(model.assignment.json_text)

    def key(f, sites):
        s = tuple(sites)
        if f == "angle-angle":
            s = min(s, (s[3], s[1], s[2], s[0]))
        elif f == "wilson_out_of_plane":
            arms = sorted((s[0], s[2], s[3]))
            s = (arms[0], s[1], arms[1], arms[2])
        else:
            s = min(s, s[::-1])
        return f, s

    expected = Counter(key(f, s) for f, s in raw_inventory(system))
    actual = Counter(key(a["family"], a["sites"]) for a in assignment["assignments"])
    if expected != actual:
        raise ValueError("Independent graph interaction mismatch")
    adjacency = {i: set() for i in labels}
    for bond in system.topology.bonds.values():
        adjacency[bond.site1].add(bond.site2)
        adjacency[bond.site2].add(bond.site1)
    checks = []
    for entry in assignment["assignments"]:
        if entry["status"] == "not_applicable":
            center = entry["sites"][1]
            if (
                entry["family"] != "wilson_out_of_plane"
                or labels[center] != "c3"
                or len(adjacency[center]) != 4
            ):
                raise ValueError("Unjustified independent applicability exclusion")
            checks.append(
                {
                    "family": entry["family"],
                    "sites": entry["sites"],
                    "independent": "tetrahedral methyl center: Wilson not applicable; all angle-angle requests retained",
                }
            )
            continue
        selected = resolve_raw(
            raw, entry["family"], [labels[i] for i in entry["sites"]]
        )
        if selected["status"] != entry["status"]:
            raise ValueError("Independent source resolution mismatch")
        if entry["status"] == "assigned":
            if selected["selected_ids"] != sorted(
                {r["record_id"] for r in entry["selected"]}
            ):
                raise ValueError("Independent source row mismatch")
            np.testing.assert_allclose(
                selected["values"], entry["normalized_values"], atol=1e-12, rtol=0
            )
        checks.append(
            {
                "family": entry["family"],
                "sites": entry["sites"],
                "independent": selected,
            }
        )
    charges, increments = independent_charges(source_path.read_bytes(), system, labels)
    np.testing.assert_allclose(
        [charges[i] for i in sorted(labels)],
        [
            {r["site"]: r["charge"] for r in unpack(model.json_text)["nonbonded"]}[i]
            for i in sorted(labels)
        ],
        atol=1e-12,
        rtol=0,
    )
    return (
        labels,
        charges,
        {
            "source_selection": checks,
            "increments": increments,
            "native_charge_total": sum(charges.values()),
        },
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument(
        "--reference-declaration",
        default="docs/evidence/phase_4j9_vertical_declaration.json",
    )
    a = p.parse_args()
    matrix = storage.read_json(a.matrix / "report.json")
    d = storage.read_json(a.reference_declaration)
    if storage.checksum(a.source.read_bytes()) != d["source"]["sha256"]:
        raise ValueError("Wrong declared source hash")
    for key in ("converter_path", "lammps_path"):
        if (
            storage.checksum(Path(d["reference"][key]).read_bytes())
            != d["reference"][key + "_sha256"]
        ):
            raise ValueError("Wrong declared reference executable")
    a.output.mkdir(parents=True, exist_ok=False)
    inputs = {}
    for case in matrix["cases"]:
        work = a.matrix / case["case"]
        inputs[case["case"]] = {
            n: storage.checksum((work / n).read_bytes())
            for n in ("system.json", "diagnostic.json", "final-graph.json")
        }
    declaration = {
        "schema": "island_j9_final_graph_numerical_declaration_v1",
        "input_hashes": inputs,
        "reference": d["reference"],
        "tolerances": d["tolerances"],
        "geometries": [
            "retained initial",
            "initial+0.015*sin(arange(3*N)+0.3) angstrom",
        ],
        "budgets": {
            "fresh_evaluations_per_case": 2,
            "external_seconds_per_command": 180,
            "new_minimization": 0,
            "new_dynamics": 0,
        },
        "source": d["source"],
        "full_source_complete": False,
    }
    storage.publish(a.output / "declaration.json", storage.json_bytes(declaration))
    source = load_pcff_source(a.source)
    selection = select_pcff_profile(LINKED_BENZENOID, sha256=source.identity["sha256"])
    outcomes = []
    for case in matrix["cases"]:
        work = a.matrix / case["case"]
        ref = a.output / case["case"]
        ref.mkdir()
        outcome = {"case": case["case"], "passed": False}
        try:
            for n, sha in inputs[case["case"]].items():
                if storage.checksum((work / n).read_bytes()) != sha:
                    raise ValueError("Changed retained input")
            system = system_from(
                storage.decode(storage.read_json(work / "system.json"))
            )
            diagnostic = storage.decode(storage.read_json(work / "diagnostic.json"))
            assignment = PCFFClass2Result(pack(diagnostic["assignment"]), source)
            assignment.validate_integrity(system)
            model = define_pcff_model(
                assignment,
                special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
            )
            prepared = adopt_forcefield(
                system, "pcff", model, pcff_profile=selection.profile()
            )
            evidence = PCFFFinalGraphPreparation(
                (work / "final-graph.json").read_text(), prepared
            )
            evidence.validate_integrity(system)
            if (
                model.identity != case["native"]
                or prepared.identity != case["prepared"]
            ):
                raise ValueError("Retained native/facade identity mismatch")
            labels, charges, checks = independent_selections(system, model, a.source)
            storage.publish(ref / "source-checks.json", storage.json_bytes(checks))
            ids = sorted(system.topology.sites)
            xyz = np.array([system.coordinates.get(i) for i in ids])
            _, neighbors = reference_inputs(system, labels, charges, xyz, ref)
            run(
                [
                    str(Path(d["reference"]["converter_path"]).resolve()),
                    "reference",
                    "-class",
                    "II",
                    "-frc",
                    str(a.source.resolve()),
                    "-p",
                    "3",
                    "-nocenter",
                ],
                ref,
                "converter",
            )
            sections = data_sections(ref / "reference.data")
            check_inventory(sections, ids, neighbors, charges)
            aa = aa_reference_coefficients(sections)
            storage.publish(ref / "aa-overrides.json", storage.json_bytes(aa))
            evaluator = PCFFSinglePointEvaluator(system, model)
            perturbed = xyz + 0.015 * np.sin(
                np.arange(xyz.size).reshape(xyz.shape) + 0.3
            )
            outcome["comparisons"] = [
                numerical(
                    evaluator,
                    system,
                    x,
                    ref,
                    name,
                    Path(d["reference"]["lammps_path"]).resolve(),
                    aa,
                    d["tolerances"],
                )
                for name, x in (("initial", xyz), ("perturbed", perturbed))
            ]
            outcome.update(
                passed=True, model=model.identity, prepared=prepared.identity
            )
        except Exception as exc:  # noqa: BLE001 -- fail closed, retaining reference failures
            outcome["error"] = type(exc).__name__ + ": " + str(exc)
        outcomes.append(outcome)
        storage.publish(ref / "outcome.json", storage.json_bytes(outcome))
    storage.publish(
        a.output / "report.json",
        storage.json_bytes({"cases": outcomes, "full_source_complete": False}),
    )
    print([(r["case"], r["passed"], r.get("error")) for r in outcomes])
    return 0 if all(r["passed"] for r in outcomes) else 1


if __name__ == "__main__":
    raise SystemExit(main())
