"""Execute authorized PCFF/COMPASS source-variant audits; incomplete models stop."""

import argparse
import json
import shutil
import subprocess
import tarfile
from pathlib import Path

import numpy as np

from island.charge_references.records import pack, unpack
from island.chemistry import from_smiles
from island.forcefields.pcff.charges import identity
from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from island.forcefields.pcff.sip_completion import assess_pcff_sip_completion
from island.forcefields.pcff.sip_variants import load_sip_source_variant
from island.workflows.bundle import system_data, system_from
from scripts.pcff_j41_c2_reference import compiled_searches, source_audit
from scripts.pcff_j41_c3_reference import audit_inventory, compass_trial, trial_searches
from scripts.pcff_j41_reference import require
from scripts.phase_4j41_c2 import digest, read, write

CASES = (
    ("silane_control", "[SiH4]", "mixed_silane"),
    ("methylsilane", "C[SiH3]", "mixed_silane"),
    ("disiloxane", "[SiH3]O[SiH3]", "siloxane"),
    ("hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C", "siloxane"),
    ("phosphate_ester", "COP(=O)(OC)OC", "phosphorus_audit"),
    ("phosphazene_probe", "[PH2]1=N[PH2]=N[PH2]=N1", "phosphorus_audit"),
    ("mixed_sip", "CO[Si](C)(C)OP(=O)(OC)OC", "mixed_sip_audit"),
)


def sources(distribution):
    """Explicit local request scope; never search the network for force fields."""
    from island.forcefields.pcff.mixed_sip import DISTRIBUTION_PINS

    folder = Path(distribution)
    specs = (
        ("pcff.frc", "PCFF", DISTRIBUTION_PINS["pcff.frc"], "4.0 / 8-October-13"),
        (
            "compass_published.frc",
            "COMPASS",
            "28156c9c95346a9efde17d0ae0aff8053410e1de86b05f8d284ef781c0bbadfe",
            "1.1 / 30-Jun-09; printed references 1-10",
        ),
    )
    result = []
    for filename, family, checksum, version in specs:
        companions = (
            [
                (folder / n, DISTRIBUTION_PINS[n])
                for n in ("pcff.rlb", "pcff_templates.dat")
            ]
            if family == "PCFF"
            else []
        )
        result.append(
            load_sip_source_variant(
                folder / filename,
                expected_sha256=checksum,
                provenance={
                    "family": family,
                    "version": version,
                    "repository": "local LAMMPS msi2lmp distribution",
                    "revision": "complete original file identified by SHA-256",
                    "source_path": "tools/msi2lmp/frc_files/" + filename,
                    "authorization": "authorized_local",
                    "authorization_basis": "Phase 4J41-C3 user-authorized local PCFF/COMPASS audit; distribution README retained; no proprietary download",
                },
                companions=companions,
            )
        )
    return result


def source_comparison(variants):
    """Separate source inventories, not a synthesized force-field profile."""
    return {
        "schema": "island.j41-c3-source-comparison.v1",
        "sources": [
            {
                "identity": v.identity,
                "filename": v.filename,
                "sha256": v.expected_sha256,
                "sip_labels": v.payload["sip_labels"],
                "family_counts": v.payload["family_counts"],
            }
            for v in variants
        ],
        "row_mixing_authorized": False,
        "selection_policy": "complete single-file source plus its declared companions; new numerical and chemical qualification required for every execution profile",
        "phosphorus_boundary": "PCFF p/p= connection/template conflict retained; COMPASS p4= reference 4 is phosphazene evidence, not phosphate/phosphonate authority",
    }


def execute(a):
    from scripts.pcff_j41_c3_silane import silane_numerics
    from scripts.phase_4j41 import reference_inputs

    root = a.output.resolve()
    root.mkdir(exist_ok=False)
    # Preserve original archived bytes; all C2 decisions are independently reloaded.
    with tarfile.open(a.c2_evidence / "execution.tar.xz") as archive:
        archive.extractall(root / "retained-c2", filter="data")
    retained = root / "retained-c2/study"
    variants = sources(a.distribution)
    report = {
        "schema": "island.j41-c3-execution.v1",
        "cases": [],
        "acceptance_gate": 1,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "ewald_status": "ewald_unvalidated",
        "pressure_status": "pressure_unavailable",
        "long_npt_status": "blocked",
        "new_execution_profiles": [],
        "converter_sha256": digest(a.converter),
        "matcher_sha256": digest(a.matcher),
        "lammps_sha256": digest(a.lammps),
    }
    write(root / "source-comparison.json", source_comparison(variants))
    shutil.copyfile(a.distribution / "README", root / "distribution-README")
    write(
        root / "source-scope.json",
        {
            "directory_entries": sorted(
                p.name for p in a.distribution.iterdir() if p.is_file()
            ),
            "selected_sources": [v.filename for v in variants],
            "excluded_families": "cff91/cvff/clayff/oplsaa are not PCFF/COMPASS source variants",
            "remote_searches": False,
        },
    )
    for v in variants:
        dest = root / v.filename.removesuffix(".frc")
        dest.mkdir()
        (dest / v.filename).write_bytes(v.raw)
        for n, raw, _ in v.companions:
            (dest / n).write_bytes(raw)
        (dest / "inventory.json").write_text(pack(v.payload))
        write(dest / "raw-source-audit.json", audit_inventory(v))
        queries = {}
        for name, smiles, domain in CASES:
            case = dest / name
            case.mkdir()
            prior_path = retained / "c1-replay" / name / "system.json"
            system = (
                system_from(read(prior_path))
                if prior_path.exists()
                else from_smiles(smiles, random_seed=2026)
            )
            write(case / "system.json", system_data(system))
            ids = sorted(system.topology.sites)
            xyz = np.array([system.coordinates.get(i) for i in ids])
            row = {
                "source": v.filename,
                "variant_identity": v.identity,
                "name": name,
                "atom_order": ids,
                "system_sha256": digest(case / "system.json"),
                "source_complete": False,
                "model_published": False,
                "dynamics_executed": False,
                "execution_profile": None,
            }
            labels, charges = None, None
            if v.filename == "pcff.frc":
                d = load_pcff_sip_distribution(dest)
                decision = assess_pcff_sip_completion(system, d, domain=domain)
                (case / "decision.json").write_text(decision.json_text)
                trial = source_audit(system, decision)
                row.update(
                    decision_identity=decision.identity,
                    source_complete=decision.payload["source_complete"],
                    preserved_c2_identity=None,
                )
                if (retained / name / "decision.json").exists():
                    require(
                        decision.identity
                        == identity(
                            unpack((retained / name / "decision.json").read_text())
                        ),
                        "Changed C2 identity",
                    )
                    row["preserved_c2_identity"] = decision.identity
                qs = trial.get("unique_searches", [])
                p = decision.payload["c1_assessment"]
                if p["inspection"]:
                    labels = p["typing"]["assignments"]
                    charges = {
                        i: float(q)
                        for i, q in trial["charges"]["partial_charges"].items()
                    }
                row["status"] = trial["status"]
            else:
                trial = compass_trial(system, v)
                qs = trial_searches(trial)
                row["source_complete"] = trial["source_complete"]
                row["status"] = (
                    "source_complete_unqualified"
                    if trial["source_complete"]
                    else (
                        "typing_blocked"
                        if trial["mapping"]["blockers"]
                        else "source_incomplete"
                    )
                )
                if trial["charges"] and trial["charges"]["complete"]:
                    labels = trial["mapping"]["labels"]
                    charges = {
                        i: float(q) for i, q in trial["charges"]["vector_e"].items()
                    }
            write(case / "raw-reconstruction.json", trial)
            for q in qs:
                queries[q["family"], q["namespace"], tuple(q["types"])] = q
            if labels is not None and name != "mixed_sip":
                reference_inputs(system, labels, charges, xyz, case)
                cmd = [
                    str(a.converter.resolve()),
                    "reference",
                    "-class",
                    "II",
                    "-frc",
                    str(dest / v.filename),
                    "-p",
                    "3",
                    "-nocenter",
                ]
                r = subprocess.run(
                    cmd,
                    cwd=case,
                    capture_output=True,
                    text=True,
                    check=False,
                    timeout=180,
                )
                (case / "converter.log").write_text(r.stdout + r.stderr)
                row["converter"] = {
                    "argv": cmd,
                    "exit_code": r.returncode,
                    "chemical_typing_authority": False,
                }
                write(case / "converter-command.json", row["converter"])
                if r.returncode == 0:
                    require(
                        name == "silane_control",
                        "New complete fixture requires qualification",
                    )
                    row["numerical_control"] = silane_numerics(
                        system, v, case, a.lammps.resolve(), execute=True
                    )
            else:
                row["converter"] = {
                    "status": "not_executed",
                    "blocking_stage": "unqualified_typing_or_incomplete_charges",
                }
            report["cases"].append(row)
            print(v.filename, name, row["status"], flush=True)
        compiled = compiled_searches(
            v.raw, a.matcher.resolve(), [queries[k] for k in sorted(queries)]
        )
        write(dest / "compiled-searches.json", compiled)
    write(root / "report.json", report)
    return 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in (
        "distribution",
        "converter",
        "matcher",
        "lammps",
        "output",
        "reconstruct",
    ):
        p.add_argument("--" + name, type=Path)
    p.add_argument(
        "--c2-evidence", type=Path, default=Path("docs/evidence/phase_4j41_c2")
    )
    a = p.parse_args()
    if a.reconstruct:
        from validate_pcff_dynamics import reconstruct

        from scripts.validate_phase_4j41_c3_receipt import validate_study

        print(json.dumps(validate_study(a.reconstruct), indent=2))
        for kind in ("nve", "baoab"):
            reconstruct(
                a.reconstruct
                / "retained-c2/study/c1-replay/baseline"
                / ("silane-" + kind)
            )
        print("fresh-process retained silane checkpoint/RNG reconstruction passed")
        return 0
    require(
        all((a.distribution, a.converter, a.matcher, a.lammps, a.output)),
        "All explicit execution paths required",
    )
    return execute(a)


if __name__ == "__main__":
    raise SystemExit(main())
