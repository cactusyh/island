"""J2 declared real-source assignment audit and bounded independent term checks.

No converter auto-equivalence or typing fidelity claim. Full-source remains an
explicit failing gate. Raw-source term oracles do not use production selection.
"""

import argparse
import json
import math
import subprocess
import time
from collections import Counter
from pathlib import Path

import numpy as np
from validate_pcff_expanded import CASES, independent_inventory
from validate_pcff_terms import XYZ, fmt, write_data

from island.charge_references.records import unpack
from island.chemistry import from_smiles
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldSources,
    create_evaluator,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.coverage import pcff_coverage_ledger
from island.forcefields.pcff.fallbacks import POLICY
from island.workflows.storage import checksum, json_bytes, publish

EXTRA = [
    ("water", "O"),
    ("silane", "[SiH4]"),
    ("dichlorine", "ClCl"),
    ("formaldehyde", "C=O"),
]


def raw_row(raw, family, types):
    """Independent exact raw-row oracle for declared term fixtures; no equivalences."""
    found = []
    section = None
    for line, text in enumerate(raw.decode().splitlines(), 1):
        words = text.split("!", 1)[0].split()
        if not words:
            continue
        if words[0].startswith("#"):
            section = words[0][1:]
        elif (
            section == family
            and words[0][0].isdigit()
            and words[2 : 2 + len(types)] == types
        ):
            found.append(
                {
                    "line": line,
                    "raw": text,
                    "values": list(map(float, words[2 + len(types) :])),
                }
            )
    if len(found) != 1:
        raise ValueError(f"Independent source row is not unique: {family} {types}")
    return found[0]


def lammps(exe, directory, xyz, kind, style, coefficients):
    directory.mkdir()
    write_data(directory / "data.lmp", xyz, kind, [list(range(len(xyz)))])
    # Deliberately handwritten LAMMPS style inputs; no production model exporter.
    script = f"""units real
atom_style full
boundary f f f
{kind}_style {style}
pair_style zero 30
special_bonds lj/coul 1 1 1
read_data data.lmp
pair_coeff * *
{kind}_coeff 1 {coefficients}
"""
    if kind == "improper":
        script += "improper_coeff 1 aa 0 0 0 109.5 109.5 109.5\n"
    script += 'thermo_style custom step pe\nthermo_modify format float %.17g\nrun 0\nprint "$(pe:%.17g)" file energy.txt\nwrite_dump all custom forces.dump id fx fy fz modify sort id format float %.17g\n'
    (directory / "in.lmp").write_text(script)
    cmd = [str(exe), "-in", "in.lmp", "-log", "log.lammps"]
    start = time.monotonic()
    p = subprocess.run(
        cmd, cwd=directory, capture_output=True, text=True, timeout=180, check=False
    )
    (directory / "stdout.txt").write_text(p.stdout)
    (directory / "stderr.txt").write_text(p.stderr)
    if p.returncode:
        raise RuntimeError(p.stdout[-1800:] + p.stderr)
    energy = float((directory / "energy.txt").read_text()) * 4.184
    lines = (
        (directory / "forces.dump")
        .read_text()
        .split("ITEM: ATOMS id fx fy fz\n")[1]
        .splitlines()
    )
    force = np.array([[float(x) for x in line.split()[1:]] for line in lines]) * 4.184
    return energy, force, {"command": cmd, "wall_seconds": time.monotonic() - start}


def compiled(family, values, xyz):
    import openmm as mm
    from openmm import unit

    from island.evaluation.pcff import build_model

    ids = list(range(len(xyz)))
    data = {
        "schema": "island_pcff_source_model_v2",
        "terms": [
            {"family": family, "sites": ids, "coefficients": values, "equilibria": []}
        ],
        "nonbonded": [{"site": i, "rmin": 1, "epsilon": 0, "charge": 0} for i in ids],
        "special_pair_inventory": [],
    }
    model = build_model(data, dict.fromkeys(ids, 12.0), mm)
    integrator = mm.VerletIntegrator(0.001)
    context = mm.Context(model, integrator, mm.Platform.getPlatformByName("Reference"))
    context.setPositions(np.array(xyz) * 0.1)
    state = context.getState(getEnergy=True, getForces=True)
    return state.getPotentialEnergy().value_in_unit(
        unit.kilojoule_per_mole
    ), state.getForces(asNumpy=True).value_in_unit(
        unit.kilojoule_per_mole / unit.nanometer
    ) * 0.1


def compare(e, f, ref, rf):
    np.testing.assert_allclose(e, ref, atol=1e-5, rtol=2e-10)
    np.testing.assert_allclose(f, rf, atol=1e-5, rtol=2e-10)
    return {
        "energy_error": float(abs(e - ref)),
        "force_max_error": float(np.max(np.abs(f - rf))),
    }


def terms(source, exe, root):
    from island.forcefields.pcff.terms import (
        FallbackClass2Term,
        finite_difference_forces,
    )

    # Source keys independently declared, not generated by the production resolver.
    fixtures = [
        ("quadratic_bond", ["cl_", "cl_"], "bond", "harmonic", 2),
        ("quadratic_angle", ["*2", "c_", "h_"], "angle", "harmonic", 3),
        ("torsion_1", ["*", "c_", "c_", "*"], "dihedral", "fourier", 4),
        ("wilson_out_of_plane", ["*", "c'_", "*", "*"], "improper", "class2", 4),
    ]
    results = []
    for family, key, kind, style, n in fixtures:
        r = raw_row(source.raw, family, key)
        p = r["values"]
        xyz = np.array(XYZ[:n])
        if family == "quadratic_bond":
            coeff = [p[0], p[1] * 4.184]
            command = fmt([p[1], p[0]])
        elif family == "quadratic_angle":
            coeff = [math.radians(p[0]), p[1] * 4.184]
            command = fmt([p[1], p[0]])
        elif family == "torsion_1":
            coeff = [p[0] * 4.184, p[1], math.radians(p[2])]
            command = "1 " + fmt(p)
        else:
            coeff = [p[0] * 4.184, math.radians(p[1])]
            command = fmt(p)
        from island.forcefields.pcff import resolve_pcff_source_record

        original_types = {
            "quadratic_bond": ["cl", "cl"],
            "quadratic_angle": ["hc", "c3", "hc"],
            "torsion_1": ["c3", "c2", "c2", "c3"],
            "wilson_out_of_plane": ["c3", "c_0", "o_1", "c3"],
        }[family]
        selection = resolve_pcff_source_record(
            source,
            family,
            original_types,
            namespace="cff91_auto",
            resolution_policy=POLICY,
        )
        if selection["status"] != "assigned" or selection["normalized_values"] != coeff:
            raise ValueError("Independent source selection mismatch: " + family)
        if not all(
            c["record_id"].endswith(":" + str(r["line"])) for c in selection["selected"]
        ):
            raise ValueError("Independent source record identity mismatch")
        ref, rf, meta = lammps(exe, root / family, xyz, kind, style, command)
        e, f = compiled(family, coeff, xyz)
        kernel = FallbackClass2Term(family, tuple(range(n)), tuple(coeff))
        coords = dict(enumerate(xyz))
        fd = []
        for h in (1e-4, 1e-5, 1e-6):
            forces = finite_difference_forces(kernel, coords, displacement=h)
            fd.append(float(np.max(np.abs(np.array(list(forces.values())) - f))))
        if fd[-1] > 1e-5:
            raise ValueError("Finite difference gate failed")
        results.append(
            {
                "family": family,
                "kind": "real_source_coefficients_on_declared_asymmetric_term_geometry",
                "source_row": r,
                "selection": selection,
                **compare(e, f, ref, rf),
                "finite_difference_errors": fd,
                **meta,
            }
        )
    # Nonzero phase exposes sign/phase conventions; explicitly synthetic coefficients.
    p = [2.3, 3, 37.0]
    xyz = np.array(XYZ)
    ref, rf, meta = lammps(
        exe, root / "nonzero_phase", xyz, "dihedral", "fourier", "1 " + fmt(p)
    )
    e, f = compiled("torsion_1", [p[0] * 4.184, p[1], math.radians(p[2])], xyz)
    results.append(
        {
            "family": "torsion_1",
            "kind": "synthetic nonzero-phase software fixture",
            **compare(e, f, ref, rf),
            **meta,
        }
    )
    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--lammps", type=Path)
    ap.add_argument("--require-full-source", action="store_true")
    ap.add_argument(
        "--cases", nargs="+", choices=[x[0] for x in CASES] + [x[0] for x in EXTRA]
    )
    args = ap.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = load_pcff_source(args.source)
    before = checksum(args.source.read_bytes())
    declaration = json.loads(
        Path("docs/evidence/phase_4j2_declaration.json").read_text()
    )
    declaration.update(command=__import__("sys").argv, source_identity=source.identity)
    publish(root / "declaration.json", json_bytes(declaration))
    report = {
        "inventory": independent_inventory(source.raw),
        "cases": [],
        "j2_acceptance": False,
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    for name, smiles in [(x[0], x[1]) for x in CASES] + EXTRA:
        if args.cases and name not in args.cases:
            continue
        started = time.monotonic()
        out = {"case": name, "smiles": smiles, "status": "failed"}
        work = root / name
        work.mkdir()
        try:
            system = from_smiles(smiles, random_seed=2026)
            typing = type_pcff_atoms(
                system, source, profile="island_pcff_source_graph_v2"
            )
            out["typed"] = typing.complete
            publish(work / "typing.json", typing.json_text.encode())
            charges = assign_automatic_pcff_charges(
                system, typing, resolution_policy=POLICY
            )
            out["charges_complete"] = charges.complete
            publish(work / "charges.json", charges.json_text.encode())
            if not charges.complete:
                out.update(
                    status="charge_missing",
                    diagnostics=charges.payload["native_charge_record"]["diagnostics"],
                )
                continue
            a = assign_pcff_parameters(
                system, typing, charges, resolution_policy=POLICY
            )
            publish(work / "parameters.json", a.json_text.encode())
            m = define_pcff_model(
                a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
            )
            publish(work / "model.json", m.json_text.encode())
            data = unpack(a.json_text)
            model = unpack(m.json_text)
            out.update(
                model_complete=model["model_definition_complete"],
                diagnostics=model["diagnostics"],
                source_coverage=data["coverage"],
                model_identity=m.identity,
                fallback_assignments=[
                    x for x in data["assignments"] if x.get("namespace") == "cff91_auto"
                ],
            )
            out["status"] = "complete" if out["model_complete"] else "parameter_missing"
            if name == "dichlorine" and out["model_complete"]:
                # Public fresh preparation is intentionally exercised, with unchanged graph.
                p = prepare_forcefield(
                    system,
                    ForceFieldRequest(
                        "pcff",
                        PCFFOptions(
                            args.source,
                            (0, 0, 1),
                            (0, 0, 1),
                            typing_profile="island_pcff_source_graph_v2",
                            resolution_policy=POLICY,
                        ),
                    ),
                )
                evaluator = create_evaluator(system, p)
                coords = np.array(
                    [system.coordinates.get(i) for i in sorted(system.topology.sites)]
                )
                out["numerical"] = []
                raw = raw_row(source.raw, "quadratic_bond", ["cl_", "cl_"])
                # Validate the actual selected identity/value against independently read row.
                assigned = next(
                    x for x in data["assignments"] if x["family"] == "quadratic_bond"
                )
                assert assigned["selected"][0]["record_id"].endswith(
                    ":" + str(raw["line"])
                )
                assert assigned["normalized_values"] == [
                    raw["values"][0],
                    raw["values"][1] * 4.184,
                ]
                with evaluator.open_session() as session:
                    for frame, xyz in enumerate(
                        [
                            coords,
                            coords
                            + 0.015
                            * np.sin(
                                np.arange(coords.size).reshape(coords.shape) + 0.3
                            ),
                        ]
                    ):
                        values = dict(zip(sorted(system.topology.sites), xyz))
                        result = session.evaluate(values)
                        fresh = evaluator.evaluate_fresh(values)
                        assert (
                            result.evaluation_fingerprint
                            == fresh.evaluation_fingerprint
                        )
                        assert result == fresh
                        if args.lammps:
                            ref, rf, meta = lammps(
                                args.lammps.resolve(),
                                work / f"reference{frame}",
                                xyz,
                                "bond",
                                "harmonic",
                                fmt(raw["values"][::-1]),
                            )
                            comparison = compare(
                                result.potential_energy,
                                np.array(list(result.forces.values())),
                                ref,
                                rf,
                            )
                            assert (
                                abs(result.energy_components["quadratic_bond"] - ref)
                                < 1e-5
                            )
                            out["numerical"].append({**comparison, **meta})
                save_prepared_forcefield(
                    system,
                    p,
                    work / "bundle",
                    sources=PreparedForceFieldSources(pcff_frc=args.source),
                )
                out["bundle_identity"] = p.identity
                out["bundle_hashes"] = {
                    str(f.relative_to(work / "bundle")): checksum(f.read_bytes())
                    for f in sorted((work / "bundle").rglob("*"))
                    if f.is_file()
                }
        except Exception as e:  # noqa: BLE001 -- durable failed experiments
            out.update(error=type(e).__name__ + ": " + str(e))
        finally:
            out["wall_seconds"] = time.monotonic() - started
            publish(work / "outcome.json", json_bytes(out))
            report["cases"].append(out)
            print(name, out["status"], flush=True)
    if args.lammps:
        try:
            report["term_checks"] = terms(source, args.lammps.resolve(), root)
        except Exception as e:  # noqa: BLE001 -- durable failed experiments
            report["term_failure"] = type(e).__name__ + ": " + str(e)
        report["executable_sha256"] = checksum(args.lammps.read_bytes())
    else:
        report["term_failure"] = "LAMMPS unavailable"
    report["source_unchanged"] = before == checksum(args.source.read_bytes())
    report["counts"] = {
        "denominator": len(report["cases"]),
        **dict(Counter(x["status"] for x in report["cases"])),
    }
    report["assignment_and_term_gate"] = (
        all(
            len(x.get("numerical", [])) == 2
            for x in report["cases"]
            if x["case"] == "dichlorine"
        )
        and len(report.get("term_checks", [])) == 5
    )
    report["workflow_gate"] = "must run validate_pcff_fallback_workflow.py separately"
    ledger = pcff_coverage_ledger(source)
    ledger["j2_fixture_evidence"] = report["counts"]
    ledger["j2_resolution_policy"] = POLICY
    ledger["j2_completion_note"] = (
        "New executable fallback does not complete unimplemented chemical rules or missing cross terms"
    )
    publish(root / "coverage.json", json_bytes(ledger))
    publish(root / "report.json", json_bytes(report))
    return (
        0 if report["assignment_and_term_gate"] and not args.require_full_source else 1
    )


if __name__ == "__main__":
    raise SystemExit(main())
