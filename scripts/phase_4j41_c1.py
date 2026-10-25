"""Execute the C1 mixed-Si/P source audit; blocked models never run dynamics."""

import argparse
import os
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from island.chemistry import from_smiles
from island.forcefields.pcff.mixed_sip import (
    assess_pcff_mixed_sip,
    load_pcff_sip_distribution,
    prepare_pcff_mixed_sip,
)
from island.workflows import storage
from island.workflows.bundle import system_data
from scripts.pcff_j41_c1_reference import compiled_source_audit, mixed_inventory
from scripts.pcff_j41_reference import audit_queries, require
from scripts.phase_4j41 import digest, reference_inputs, write
from scripts.validate_phase_4j41_receipt import validate_study as validate_baseline

CASES = [
    ("silane_control", "[SiH4]", "mixed_silane"),
    ("methylsilane", "C[SiH3]", "mixed_silane"),
    ("disiloxane", "[SiH3]O[SiH3]", "siloxane"),
    ("hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C", "siloxane"),
    ("sila_oligomer", "C[SiH2]C[SiH2]C", "organosilicon_polymer"),
    ("si_n", "N[SiH3]", "mixed_silane"),
    ("si_s", "S[SiH3]", "mixed_silane"),
    ("si_f", "F[SiH3]", "mixed_silane"),
    ("si_cl", "Cl[SiH3]", "mixed_silane"),
    ("si_br", "Br[SiH3]", "mixed_silane"),
    ("si_i", "I[SiH3]", "mixed_silane"),
    ("phosphate_ester", "COP(=O)(OC)OC", "phosphorus_audit"),
    ("phosphonate", "CP(=O)(OC)OC", "phosphorus_audit"),
    ("phosphine", "P", "phosphorus_audit"),
    ("phosphine_oxide", "CP(C)(C)=O", "phosphorus_audit"),
    ("pf5", "FP(F)(F)(F)F", "phosphorus_audit"),
    ("mixed_sip", "CO[Si](C)(C)OP(=O)(OC)OC", "mixed_sip_audit"),
]


def execute(args):
    root = args.output.resolve()
    root.mkdir(exist_ok=False)
    distribution = load_pcff_sip_distribution(args.distribution)
    for name in distribution.hashes:
        shutil.copyfile(args.distribution / name, root / name)
    shutil.copyfile(args.matcher, root / "matcher-leaf.so")
    command = [
        sys.executable,
        str(Path(__file__).with_name("phase_4j41.py").resolve()),
        "--distribution",
        str(args.distribution.resolve()),
        "--converter",
        str(args.converter.resolve()),
        "--lammps",
        str(args.lammps.resolve()),
        "--output",
        str(root / "baseline"),
    ]
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    (root / "baseline-execution.log").write_text(result.stdout + result.stderr)
    write(
        root / "baseline-command.json",
        {"argv": command, "exit_code": result.returncode},
    )
    require(result.returncode == 1, "Expected retained J41 source-incomplete gate")
    baseline = validate_baseline(root / "baseline")
    write(
        root / "mixed-source-inventory.json", mixed_inventory(distribution.source.raw)
    )
    write(
        root / "compiled-source-audit.json",
        compiled_source_audit(distribution.source, args.matcher.resolve()),
    )
    report = {
        "schema": "island.j41-c1-execution.v1",
        "argv": sys.argv,
        "distribution_hashes": distribution.hashes,
        "baseline_validation": baseline,
        "converter": {
            "path": str(args.converter.resolve()),
            "sha256": digest(args.converter),
        },
        "lammps": {"path": str(args.lammps.resolve()), "sha256": digest(args.lammps)},
        "matcher_sha256": digest(args.matcher),
        "cases": [],
        "acceptance_gate": 1,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "ewald_status": "ewald_unvalidated",
        "pressure_status": "pressure_unavailable",
        "long_npt_status": "blocked",
    }
    for name, smiles, suffix in CASES:
        folder = root / name
        folder.mkdir()
        system = from_smiles(smiles, random_seed=2026)
        ids = sorted(system.topology.sites)
        xyz = np.array([system.coordinates.get(i) for i in ids])
        xyz += 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3)
        for i, pos in zip(ids, xyz):
            system.coordinates.set(i, pos)
        write(folder / "system.json", system_data(system))
        a = assess_pcff_mixed_sip(
            system, distribution, profile="island_pcff_" + suffix + "_c1_v1"
        )
        (folder / "assessment.json").write_text(a.json_text)
        p = a.payload
        row = {
            "name": name,
            "smiles": smiles,
            "seed": 2026,
            "atom_count": len(ids),
            "assessment_identity": a.identity,
            "graph_identity": p["graph_identity"],
            "source_complete": p["complete"],
            "charge_origin": p["charge_origin"],
            "phosphorus_environments": p["phosphorus_environments"],
            "typing_status": "complete"
            if p["typing"] and p["typing"]["coverage"]["complete"]
            else "blocked",
            "missing_interactions": dict(
                Counter(
                    d["interaction"]["family"]
                    for d in p["diagnostics"]
                    if d["reason"] == "required_source_row_or_dependency_unresolved"
                )
            ),
            "diagnostic_reasons": dict(Counter(d["reason"] for d in p["diagnostics"])),
            "model_published": False,
            "new_dynamics_executed": False,
        }
        if p["inspection"]:
            audit = audit_queries(system, a, distribution.source.raw)
            write(folder / "independent-audit.json", audit)
            labels = p["typing"]["assignments"]
            charges = {
                i: float(q) for i, q in audit["charges"]["partial_charges"].items()
            }
            reference_inputs(system, labels, charges, xyz, folder)
            cmd = [
                str(args.converter.resolve()),
                "reference",
                "-class",
                "II",
                "-frc",
                str(root / "pcff.frc"),
                "-p",
                "3",
                "-nocenter",
            ]
            converted = subprocess.run(
                cmd,
                cwd=folder,
                capture_output=True,
                text=True,
                check=False,
                timeout=180,
            )
            (folder / "converter.log").write_text(converted.stdout + converted.stderr)
            row["converter"] = {
                "argv": cmd,
                "exit_code": converted.returncode,
                "assembly_is_not_chemical_typing": True,
            }
        else:
            row["converter"] = {
                "status": "not_executed",
                "reason": "No independently justified complete chemical type mapping; no guessed input labels",
            }
        if p["complete"]:
            require(
                name == "silane_control",
                "New fixture needs independent numerical qualification",
            )
            model = prepare_pcff_mixed_sip(system, a)
            baseline_report = storage.decode(
                storage.read_json(root / "baseline" / "report.json")
            )
            require(
                model.identity == baseline_report["cases"][0]["model_identity"],
                "Changed retained finite silane model",
            )
            row["model_identity"] = model.identity
            row["workflow_status"] = "fresh_J41_control_replay_passed"
        else:
            row["workflow_status"] = "blocked_before_model_publication"
        report["cases"].append(row)
        write(root / "report.json", report)
        print(name, row["workflow_status"], flush=True)
    return 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--distribution", type=Path)
    p.add_argument("--converter", type=Path)
    p.add_argument("--lammps", type=Path)
    p.add_argument("--matcher", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--reconstruct", type=Path)
    a = p.parse_args()
    if a.reconstruct:
        from validate_pcff_dynamics import reconstruct

        from scripts.validate_phase_4j41_c1_receipt import validate_study

        print(validate_study(a.reconstruct))
        for kind in ("nve", "baoab"):
            reconstruct(a.reconstruct / "baseline" / ("silane-" + kind))
        print("fresh-process reconstruction passed", os.getpid())
        return 0
    require(
        all((a.distribution, a.converter, a.lammps, a.matcher, a.output)),
        "All execution paths required",
    )
    return execute(a)


if __name__ == "__main__":
    raise SystemExit(main())
