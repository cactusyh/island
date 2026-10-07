"""J17 declared external polyethylene: raw source, compiled converter and LAMMPS.

The independent labels and Decimal charge vector are established before calling
any production typing/charge/assignment API. Original converter output is kept.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from pcff_j5_reference import charge_raw, read_source, resolve_raw
from pcff_j13_reference import converter_zero_check, selection_check
from pcff_j14_reference import aa_raw_reference_coefficients
from validate_pcff_profile import numerical, reference
from validate_pcff_singlepoint import (
    check_inventory,
    data_sections,
    reference_inputs,
    run,
)
from validate_prepared_workflow import hashes

from island.builders import build_linear_polymer
from island.exceptions import PCFFError
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldSources,
    create_evaluator,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff import (
    assign_pcff_source_types,
    assign_typed_pcff_charges,
    bind_pcff_types,
    load_pcff_source,
    provide_pcff_charges,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY
from island.workflows import storage
from island.workflows.bundle import system_data


def independent_inputs(system, source):
    """Bounded source descriptions, checked without the production recognizer."""
    atoms = system.topology.sites
    adj = {i: [] for i in atoms}
    for b in system.topology.bonds.values():
        assert b.order == 1 and not b.aromatic
        adj[b.site1].append(b.site2)
        adj[b.site2].append(b.site1)
    labels = {}
    for i, a in atoms.items():
        assert a.formal_charge == 0 and not a.metadata.get("aromatic")
        assert not a.metadata.get("isotope") and not a.metadata.get("radical_electrons")
        if a.element == "C":
            assert len(adj[i]) == 4 and all(
                atoms[j].element in {"C", "H"} for j in adj[i]
            )
            labels[i] = "c"
        else:
            assert (
                a.element == "H"
                and len(adj[i]) == 1
                and atoms[adj[i][0]].element == "C"
            )
            labels[i] = "h"
    descriptors = []
    section = None
    for n, line in enumerate(source.read_text().splitlines(), 1):
        w = line.split()
        if not w:
            continue
        if w[0].startswith("#"):
            section = w[0]
        elif section == "#atom_types" and len(w) >= 6 and w[2] in {"c", "h"}:
            descriptors.append({"line": n, "raw": line, "label": w[2]})
    assert {r["label"] for r in descriptors} == {"c", "h"}
    assert "generic SP3 carbon" in next(
        r["raw"] for r in descriptors if r["label"] == "c"
    )
    assert "hydrogen bound to C" in next(
        r["raw"] for r in descriptors if r["label"] == "h"
    )
    charge = charge_raw(read_source(source.read_bytes()), system, labels)
    assert charge["complete"] and all(c["formal"] == 0 for c in charge["components"])
    return (
        labels,
        {i: float(q) for i, q in charge["partial_charges"].items()},
        descriptors,
        charge,
    )


def dependencies(model, source, labels):
    raw = read_source(source.read_bytes())
    checked = []
    for entry in model.assignment.payload["assignments"]:
        for dep in entry["dependencies"]:
            family = dep["assignment_id"].split(":")[0]
            resolved = resolve_raw(
                raw, family, [labels[i] for i in dep["sites"]], guarded_msi=True
            )
            assert resolved["status"] == dep["status"] == "assigned"
            np.testing.assert_allclose(
                resolved["values"][0], dep["equilibrium_value"], atol=1e-12, rtol=0
            )
            assert resolved["selected_ids"] == sorted(set(dep["source_rows"]))
            checked.append(
                {
                    "interaction": entry["id"],
                    "dependency": dep["assignment_id"],
                    "source_rows": resolved["selected_ids"],
                    "value": resolved["values"][0],
                }
            )
    return checked


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    root = a.output
    root.mkdir(parents=True, exist_ok=False)
    declaration = storage.read_json("docs/evidence/phase_4j17_declaration.json")
    historical = storage.read_json("docs/evidence/phase_4j13_declaration.json")
    assert storage.checksum(a.source.read_bytes()) == declaration["source_sha256"]
    refs = historical["reference"]
    for key in ("converter_path", "lammps_path"):
        assert storage.checksum(Path(refs[key]).read_bytes()) == refs[key + "_sha256"]
    tol = historical["tolerances"]
    storage.publish(
        root / "declaration.json",
        storage.json_bytes(
            {
                "parent_sha256": storage.checksum(
                    Path("docs/evidence/phase_4j17_declaration.json").read_bytes()
                ),
                "reference": refs,
                "tolerances": tol,
                "reference_correction": "Retained J14 center/shared-arm raw-source AA correction; original converter result also evaluated",
            }
        ),
    )
    out = {
        "passed": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        system = build_linear_polymer(
            **{k: v for k, v in declaration["polymer"].items() if k != "types"}
        )
        labels, charges, descriptors, independent_charge = independent_inputs(
            system, a.source
        )
        storage.publish(
            root / "system.json",
            storage.json_bytes(storage.encode(system_data(system))),
        )
        storage.publish(
            root / "independent-inputs.json",
            storage.json_bytes(
                storage.encode(
                    {
                        "labels": labels,
                        "descriptors": descriptors,
                        "charges": independent_charge,
                        "chemical_scope": "source-described generic sp3 C/H family; suitability not established",
                    }
                )
            ),
        )
        source = load_pcff_source(a.source)
        types = bind_pcff_types(
            system,
            source,
            labels,
            provenance="Independent final-graph saturated C/H checks and pinned FRC atom descriptions",
            evidence_references=[f"atom_types:cff91:{r['line']}" for r in descriptors],
        )
        try:
            assign_pcff_source_types(
                system,
                source,
                labels,
                provenance="Independent source-family assignment",
            )
        except PCFFError as error:
            assert "contradict" in str(error)
            out["legacy_rejection"] = str(error)
        else:
            raise AssertionError(
                "Legacy agreement unexpectedly accepted alternative labels"
            )
        native = assign_typed_pcff_charges(
            system, types, resolution_policy=COMPATIBILITY_POLICY
        )
        assert native.complete
        np.testing.assert_allclose(
            [native.charges[i] for i in sorted(labels)],
            [charges[i] for i in sorted(labels)],
            atol=1e-12,
            rtol=0,
        )
        provided = provide_pcff_charges(
            system,
            types,
            charges,
            unit="elementary_charge",
            provenance="Independent Decimal raw-FRC bond-increment reconstruction (not native ISLAND output)",
            evidence_references=[
                "independent-inputs.json",
                "scripts/pcff_j5_reference.py:charge_raw",
            ],
            component_totals={
                c["sites"][0]: c["formal"] for c in independent_charge["components"]
            },
            total_charge=0,
            resolution_policy=COMPATIBILITY_POLICY,
        )
        assert provided.charges == charges
        prepared_modes = {}
        for mode, q in (("native", native), ("provided", provided)):
            prepared_modes[mode] = prepare_forcefield(
                system,
                ForceFieldRequest(
                    "pcff",
                    PCFFOptions(
                        a.source,
                        (0, 0, 1),
                        (0, 0, 1),
                        resolution_policy=COMPATIBILITY_POLICY,
                        typed_graph=types,
                        graph_charges=q,
                    ),
                ),
            )
        prepared = prepared_modes["provided"]
        model = prepared.native_result
        _, _, selections = selection_check(
            system, model, a.source, reference_labels=labels
        )
        selections["dependencies"] = dependencies(model, a.source, labels)
        storage.publish(root / "source-checks.json", storage.json_bytes(selections))
        ids = sorted(labels)
        xyz = np.array([system.coordinates.get(i) for i in ids])
        _, neighbors = reference_inputs(system, labels, charges, xyz, root)
        command = [
            str(Path(refs["converter_path"]).resolve()),
            "reference",
            "-class",
            "II",
            "-frc",
            str(a.source.resolve()),
            "-p",
            "3",
            "-nocenter",
        ]
        run(command, root, "converter")
        sections = data_sections(root / "reference.data")
        out["converter_inventory"] = check_inventory(sections, ids, neighbors, charges)
        out["converter_bb13"] = converter_zero_check(sections)
        aa, aa_evidence = aa_raw_reference_coefficients(sections, labels, a.source)
        storage.publish(
            root / "aa-reference-correction.json",
            storage.json_bytes({"commands": aa, "evidence": aa_evidence}),
        )
        evaluator = create_evaluator(system, prepared)
        out["comparisons"] = []
        out["raw_converter_comparisons"] = []
        for name, x in (
            ("initial", xyz),
            (
                "perturbed",
                xyz + 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3),
            ),
        ):
            # Keep any converter mismatch; a corrected reference never overwrites it.
            try:
                raw_check = numerical(
                    evaluator,
                    system,
                    x,
                    root,
                    name + "-raw",
                    Path(refs["lammps_path"]).resolve(),
                    [],
                    tol,
                )
                raw_check["passed"] = True
            except AssertionError as error:
                raw_check = {"passed": False, "error": str(error)}
            out["raw_converter_comparisons"].append(raw_check)
            out["comparisons"].append(
                numerical(
                    evaluator,
                    system,
                    x,
                    root,
                    name + "-corrected",
                    Path(refs["lammps_path"]).resolve(),
                    aa,
                    tol,
                )
            )
        baseline = evaluator.evaluate_fresh()
        native_result = create_evaluator(
            system, prepared_modes["native"]
        ).evaluate_fresh()
        np.testing.assert_allclose(
            baseline.potential_energy,
            native_result.potential_energy,
            atol=1e-12,
            rtol=0,
        )
        # Differentiate the bound energy; compare directly to independent LAMMPS forces.
        _, independent_force = reference(
            Path(refs["lammps_path"]).resolve(), root / "fd-reference", xyz, ids, aa
        )
        errors = []
        with evaluator.open_session() as session:
            for h in tol["finite_difference_angstrom"]:
                fd = np.zeros_like(xyz)
                for i in range(len(ids)):
                    for axis in range(3):
                        plus, minus = xyz.copy(), xyz.copy()
                        plus[i, axis] += h
                        minus[i, axis] -= h
                        ep = session.evaluate(
                            dict(zip(ids, plus, strict=True))
                        ).potential_energy
                        em = session.evaluate(
                            dict(zip(ids, minus, strict=True))
                        ).potential_energy
                        fd[i, axis] = -(ep - em) / (2 * h)
                errors.append(float(np.max(np.abs(fd - independent_force))))
        assert errors[-1] < tol["force_atol_kj_mol_angstrom"]
        out["finite_difference_errors"] = errors
        save_prepared_forcefield(
            system,
            prepared,
            root / "bundle",
            sources=PreparedForceFieldSources(pcff_frc=a.source),
        )
        out.update(
            passed=True,
            bundle_hashes=hashes(root / "bundle"),
            prepared=prepared.identity,
            model=model.identity,
            typed_graph=types.identity,
            native_charges=native.identity,
            provided_charges=provided.identity,
            native_prepared=prepared_modes["native"].identity,
            converter_command=command,
            source_rows_checked=len(selections["rows"]),
            dependencies_checked=len(selections["dependencies"]),
        )
    except Exception as error:  # noqa: BLE001 -- preserve failures
        out.update(error_type=type(error).__name__, error=str(error))
    storage.publish(root / "outcome.json", storage.json_bytes(out))
    print(json.dumps(out, indent=2))
    return int(not out["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
