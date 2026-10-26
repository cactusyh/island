"""Execute C4 source-bound converter resolution audits from final graphs.

CAR/MDF are generated only for the independent oracle. This script publishes
no runtime profile unless a separate complete numerical qualification exists.
The full-source acceptance command exits 1 for classified blocked results.
"""

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.exceptions import PCFFError
from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from island.forcefields.pcff.sip_completion import assess_pcff_sip_completion
from island.workflows.bundle import system_data, system_from
from scripts.pcff_j41_c3_reference import audit_inventory, compass_labels
from scripts.pcff_j41_c3_silane import silane_numerics
from scripts.pcff_j41_c4_reference import reconstruct, trace_audit
from scripts.pcff_j41_c4_tools import build, command, source_contract
from scripts.pcff_j41_reference import require
from scripts.phase_4j41_c2 import digest, read, write
from scripts.phase_4j41_c3 import sources
from scripts.validate_phase_4j41_c3_receipt import extract

BASE = "7e1150a443b72d518474dc050a02a1192bc486a0"
BRANCH = "codex/phase-4j41-c4-executable-msi2lmp-resolution"
CASES = (
    ("silane_control", "[SiH4]", "mixed_silane"),
    ("methylsilane", "C[SiH3]", "mixed_silane"),
    ("disiloxane", "[SiH3]O[SiH3]", "siloxane"),
    ("hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C", "siloxane"),
    ("phosphazene_probe", "[PH2]1=N[PH2]=N[PH2]=N1", "phosphorus_audit"),
    ("phosphate_ester", "COP(=O)(OC)OC", "phosphorus_audit"),
    ("phosphonate", "CP(=O)(OC)OC", "phosphorus_audit"),
    ("phosphine_oxide", "CP(C)(C)=O", "phosphorus_audit"),
    ("protonated_p", "[PH4+]", "phosphorus_audit"),
    ("pf5", "FP(F)(F)(F)F", "phosphorus_audit"),
    ("mixed_sip", "CO[Si](C)(C)OP(=O)(OC)OC", "mixed_sip_audit"),
    ("unsupported_si", "[SiH3]", "mixed_silane"),
)
BLOCKS = {
    "ewald": "ewald_unvalidated",
    "pressure": "pressure_unavailable",
    "cutoff_recovery": "blocked",
    "long_npt": "blocked",
}


def mapping(system, variant, domain, folder):
    if variant.filename == "pcff.frc":
        decision = assess_pcff_sip_completion(
            system, load_pcff_sip_distribution(folder), domain=domain
        )
        p = decision.payload["c1_assessment"]
        labels = p["typing"]["assignments"] if p["inspection"] else None
        return labels, {
            "c2_decision_identity": decision.identity,
            "source_evidence": decision.payload["source_evidence"],
            "phosphorus_decisions": decision.payload["phosphorus_decisions"],
            "assessment": p,
            "authority": "unchanged source-bound C1/C2 final-graph assessment; converter supplies no typing",
        }
    trial = compass_labels(system, variant)
    return (trial["labels"] if not trial["blockers"] else None), trial


def charge_vector(audit):
    charge = audit["charges"]
    if not charge or not charge["complete"]:
        return None
    return {
        i: float(q)
        for i, q in charge.get("partial_charges", charge.get("vector_e")).items()
    }


def convert(case, system, labels, charges, variant, tools):
    from scripts.phase_4j41 import reference_inputs

    ids = sorted(system.topology.sites)
    xyz = np.array([system.coordinates.get(i) for i in ids])
    runs = {}
    for kind in ("original", "trace"):
        folder = case / kind
        folder.mkdir()
        reference_inputs(system, labels, charges, xyz, folder)
        argv = [
            str(tools / kind / "msi2lmp"),
            "reference",
            "-class",
            "II",
            "-frc",
            str(case.parent / variant.filename),
            "-p",
            "3",
            "-nocenter",
        ]
        run = command(argv, folder)
        write(folder / "command.json", run)
        (folder / "converter.log").write_text(run["stdout"] + run["stderr"])
        runs[kind] = run
    require(
        runs["original"]["exit_code"] == runs["trace"]["exit_code"],
        "Trace changed converter result",
    )
    # The added trace is diagnostic only. Original stdout and stderr agree
    # exactly after removing trace lines (no error is swallowed or ignored).
    require(
        runs["original"]["stdout"] == runs["trace"]["stdout"],
        "Trace changed converter stdout",
    )
    stripped = (
        "\n".join(
            s
            for s in runs["trace"]["stderr"].splitlines()
            if not s.startswith("C4_MATCH ")
        )
        + "\n"
    )
    require(stripped == runs["original"]["stderr"], "Trace changed converter stderr")
    if runs["original"]["exit_code"] == 0:
        require(
            digest(case / "original/reference.data")
            == digest(case / "trace/reference.data"),
            "Trace changed coefficients",
        )
    checked = trace_audit(runs["trace"]["stderr"], variant)
    write(case / "full-converter-searches.json", checked)
    return {
        "exit_code": runs["original"]["exit_code"],
        "status": "conversion_completed"
        if runs["original"]["exit_code"] == 0
        else "converter_failure",
        "first_failure": next(
            (
                s.strip()
                for s in runs["original"]["stdout"].splitlines()
                if "Unable" in s
            ),
            None,
        ),
        "full_converter_calls_checked": len(checked),
        "trace_preserves_output": True,
        "chemical_typing_authority": False,
    }


def execute(a):
    root = a.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island.j41-c4-execution.v1",
        "base_commit": BASE,
        "cases": [],
        "acceptance_gate": 1,
        "new_execution_profiles": [],
        "production_validated": False,
        "simulation_readiness": "not_established",
        "periodic": BLOCKS,
        "public_input": "PSMILES/final graph; internally generated CAR/MDF oracle only",
    }
    # Missing tools/source inputs are classified, rather than a traceback being
    # reported as successful blocked evidence.
    if not a.source.is_dir() or not a.lammps.is_file():
        report.update(
            status="converter_unavailable",
            blocking_stage="authorized_converter_source_or_lammps_missing",
        )
        write(root / "report.json", report)
        return 1
    try:
        variants = sources(a.distribution)
    except (OSError, PCFFError, ValueError) as error:
        report.update(
            status="source_incomplete",
            blocking_stage="authorized_source_integrity",
            diagnostic=str(error),
        )
        write(root / "report.json", report)
        return 1
    try:
        tools = build(a.source.resolve(), root / "tools", a.lammps.resolve())
    except (OSError, subprocess.TimeoutExpired, ValueError) as error:
        tools = {"status": "converter_unavailable", "diagnostic": str(error)}
    if tools["status"] != "built":
        report.update(
            status="converter_unavailable",
            blocking_stage="converter_build",
            diagnostic=tools,
        )
        write(root / "report.json", report)
        return 1
    write(root / "converter-source-contract.json", source_contract(a.source))
    extract(a.c3_evidence / "execution.tar.xz", root / "retained-c3")
    old = root / "retained-c3/study"
    for v in variants:
        dest = root / v.filename.removesuffix(".frc")
        dest.mkdir()
        (dest / v.filename).write_bytes(v.raw)
        for n, raw, _ in v.companions:
            (dest / n).write_bytes(raw)
        (dest / "inventory.json").write_text(pack(v.payload))
        write(dest / "independent-inventory.json", audit_inventory(v))
        for name, smiles, domain in CASES:
            case = dest / name
            case.mkdir()
            prior = old / "pcff" / name / "system.json"
            system = (
                system_from(read(prior))
                if prior.exists()
                else from_smiles(smiles, random_seed=2026)
            )
            write(case / "system.json", system_data(system))
            labels, typing = mapping(system, v, domain, dest)
            audit = reconstruct(system, v, labels, typing)
            write(case / "raw-reconstruction.json", audit)
            row = {
                "source": v.filename,
                "name": name,
                "status": audit["status"],
                "system_sha256": digest(case / "system.json"),
                "audit_sha256": digest(case / "raw-reconstruction.json"),
                "source_sha256": v.expected_sha256,
                "companion_hashes": v.payload["companion_hashes"],
                "converter_sha256": tools["builds"]["original"]["executable_sha256"],
                "lammps_sha256": tools["lammps"]["sha256"],
                "execution_profile": None,
                "dynamics_executed": False,
                "missing_charge_increment": audit.get("missing_charge_increment"),
                "blocking_contexts": audit.get("missing_couplings", []),
                "blocking_stage": "chemical_typing"
                if labels is None
                else "required_source_interactions",
            }
            charges = charge_vector(audit)
            if labels is not None and charges is not None:
                row["converter"] = convert(
                    case, system, labels, charges, v, root / "tools"
                )
                if row["converter"]["exit_code"] == 0:
                    if name == "silane_control":
                        row["retained_numerical_control"] = silane_numerics(
                            system,
                            v,
                            case / "original",
                            a.lammps.resolve(),
                            execute=True,
                        )
                        row["blocking_stage"] = "retained_control_only_no_new_profile"
                    else:
                        row.update(
                            status="parameter_incomplete"
                            if audit["status"] == "source_incomplete"
                            else "scientifically_unqualified",
                            blocking_stage="full_independent_qualification_required",
                        )
                else:
                    row["blocking_stage"] = "converter_required_source_key"
                    if audit["status"] != "source_incomplete":
                        row["status"] = "parameter_incomplete"
            else:
                row["converter"] = {
                    "status": "not_executed",
                    "reason": "typing_blocked"
                    if labels is None
                    else "missing_charge_increment",
                }
                if labels is not None:
                    row["blocking_stage"] = "missing_charge_increment"
            report["cases"].append(row)
            print(
                v.filename, name, row["status"], row["converter"]["status"], flush=True
            )
    report.update(
        status="source_incomplete", blocking_stage="no_complete_new_mixed_domain"
    )
    write(root / "report.json", report)
    return 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--distribution", type=Path)
    p.add_argument("--source", type=Path)
    p.add_argument("--lammps", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument(
        "--c3-evidence", type=Path, default=Path("docs/evidence/phase_4j41_c3")
    )
    p.add_argument("--reconstruct", type=Path)
    a = p.parse_args()
    if a.reconstruct:
        from scripts.validate_phase_4j41_c4_receipt import validate_study

        print(json.dumps(validate_study(a.reconstruct, replay=True), indent=2))
        return 0
    if not all((a.distribution, a.source, a.lammps, a.output)):
        p.error(
            "explicit authorized source, converter source, LAMMPS and output paths required"
        )
    return execute(a)


if __name__ == "__main__":
    raise SystemExit(main())
