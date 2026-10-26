"""Validate C4 raw converter reconstruction, receipts and historical hashes."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import numpy as np

from island.charge_references.records import unpack
from island.chemistry import from_smiles
from island.workflows import storage
from island.workflows.bundle import system_data, system_from
from scripts.pcff_j41_c3_reference import audit_inventory
from scripts.pcff_j41_c3_silane import silane_numerics
from scripts.pcff_j41_c4_reference import reconstruct, trace_audit
from scripts.pcff_j41_c4_tools import source_contract
from scripts.phase_4j41_c2 import read
from scripts.phase_4j41_c3 import sources
from scripts.phase_4j41_c4 import BASE, BLOCKS, CASES, charge_vector, mapping
from scripts.validate_phase_4j41_c3_receipt import extract, reload_variant
from scripts.validate_phase_4j41_c3_receipt import validate_study as validate_c3
from scripts.validate_phase_4j41_receipt import require, sha, signed


def validate_case(case, variant, domain):
    system = system_from(read(case / "system.json"))
    labels, typing = mapping(system, variant, domain, case.parent)
    expected = reconstruct(system, variant, labels, typing)
    require(
        read(case / "raw-reconstruction.json") == expected,
        "Contradictory source, graph, atom ordering, labels, charges, dependencies or qualification",
    )
    return system, expected


def replay_converter(root, case, variant, expected_code):
    """Rerun the original executable after relocation; never modify the study."""
    with tempfile.TemporaryDirectory(prefix="island-c4-oracle-") as temp:
        import shutil

        folder = Path(temp)
        for n in ("reference.car", "reference.mdf"):
            shutil.copyfile(case / "original" / n, folder / n)
        cmd = [
            str(root / "tools/original/msi2lmp"),
            "reference",
            "-class",
            "II",
            "-frc",
            str(case.parent / variant.filename),
            "-p",
            "3",
            "-nocenter",
        ]
        result = subprocess.run(
            cmd, cwd=folder, capture_output=True, text=True, timeout=180, check=False
        )
        require(
            result.returncode == expected_code, "Relocated converter result changed"
        )
        original = read(case / "original/command.json")
        before = [s for s in original["stdout"].splitlines() if "Unable" in s]
        after = [s for s in result.stdout.splitlines() if "Unable" in s]
        require(before == after, "Relocated converter missing keys changed")
        if expected_code == 0:
            require(
                sha(folder / "reference.data") == sha(case / "original/reference.data"),
                "Relocated converter coefficients changed",
            )


def validate_study(root, *, replay=False):
    from scripts.phase_4j41 import reference_inputs

    root = Path(root).resolve()
    report = read(root / "report.json")
    require(
        report["schema"] == "island.j41-c4-execution.v1"
        and report["base_commit"] == BASE,
        "Wrong C4 base or schema",
    )
    require(
        report["acceptance_gate"] == 1
        and report["new_execution_profiles"] == []
        and report["production_validated"] is False
        and report["simulation_readiness"] == "not_established"
        and report["periodic"] == BLOCKS
        and report["status"] == "source_incomplete",
        "Promoted C4 acceptance",
    )
    require(
        report["public_input"]
        == "PSMILES/final graph; internally generated CAR/MDF oracle only",
        "CAR/MDF runtime requirement introduced",
    )
    retained = validate_c3(root / "retained-c3/study")
    tools = read(root / "tools/build.json")
    contract = source_contract(root / "tools/original")
    require(
        read(root / "converter-source-contract.json") == contract,
        "Changed algorithm provenance",
    )
    for kind in ("original", "trace"):
        record = read(root / "tools" / (kind + "-build.json"))
        require(
            record == tools["builds"][kind] and record["exit_code"] == 0,
            "Converter build failed",
        )
        require(
            sha(root / "tools" / kind / "msi2lmp") == record["executable_sha256"],
            "Changed converter executable",
        )
        for name, h in record["source_hashes"].items():
            require(sha(root / "tools" / kind / name) == h, "Changed converter source")
        if kind == "original":
            require(
                record["source_hashes"] == tools["source_hashes"],
                "Original converter was modified",
            )
    # Retained C3 copies of converter sources independently bind the revision.
    for n, h in tools["source_hashes"].items():
        require(
            sha(root / "retained-c3/converter" / n) == h,
            "Converter differs from C3 audited sources",
        )
    variants = [
        reload_variant(root / n.removesuffix(".frc"), n)
        for n in ("pcff.frc", "compass_published.frc")
    ]
    with tempfile.TemporaryDirectory(prefix="island-c4-source-") as temp:
        folder = Path(temp)
        for v in variants:
            (folder / v.filename).write_bytes(v.raw)
            for n, raw, _ in v.companions:
                (folder / n).write_bytes(raw)
        require(
            [v.identity for v in sources(folder)] == [v.identity for v in variants],
            "Unauthorized or mixed source",
        )
    require(
        [(r["source"], r["name"]) for r in report["cases"]]
        == [(v.filename, n) for v in variants for n, _, _ in CASES],
        "Incomplete or duplicate target matrix",
    )
    counts = {
        "converter_runs": 0,
        "source_row_calls_checked": 0,
        "typing_blocked": 0,
        "source_incomplete": 0,
        "scientifically_unqualified": 0,
    }
    findings = []
    for variant in variants:
        dest = root / variant.filename.removesuffix(".frc")
        require(
            unpack((dest / "inventory.json").read_text()) == variant.payload,
            "Changed source inventory",
        )
        require(
            read(dest / "independent-inventory.json") == audit_inventory(variant),
            "Changed raw source audit",
        )
        for row, (name, smiles, domain) in zip(
            [r for r in report["cases"] if r["source"] == variant.filename],
            CASES,
            strict=True,
        ):
            case = dest / name
            system, audit = validate_case(case, variant, domain)
            require(
                row["audit_sha256"] == sha(case / "raw-reconstruction.json")
                and row["system_sha256"] == sha(case / "system.json")
                and row["source_sha256"] == variant.expected_sha256
                and row["companion_hashes"] == variant.payload["companion_hashes"]
                and row["converter_sha256"]
                == tools["builds"]["original"]["executable_sha256"]
                and row["lammps_sha256"] == tools["lammps"]["sha256"],
                "Changed graph, source or executable binding",
            )
            prior = root / "retained-c3/study/pcff" / name / "system.json"
            if prior.exists():
                require(sha(prior) == row["system_sha256"], "Retained graph changed")
            else:
                encoded = storage.json_bytes(
                    storage.encode(system_data(from_smiles(smiles, random_seed=2026)))
                )
                from hashlib import sha256

                require(
                    sha256(encoded).hexdigest() == row["system_sha256"],
                    "Changed final graph or coordinates",
                )
            require(
                row["status"] == audit["status"]
                and row["execution_profile"] is None
                and row["dynamics_executed"] is False,
                "Promoted blocked target",
            )
            require(
                row["blocking_contexts"] == audit.get("missing_couplings", [])
                and row["missing_charge_increment"]
                == audit.get("missing_charge_increment"),
                "Changed exact missing interactions",
            )
            charges = charge_vector(audit)
            if audit["labels"] is None or charges is None:
                require(
                    row["converter"]
                    == {
                        "status": "not_executed",
                        "reason": "typing_blocked"
                        if audit["labels"] is None
                        else "missing_charge_increment",
                    },
                    "Guessed converter inputs",
                )
                require(not (case / "original").exists(), "Blocked target converted")
                require(
                    row["blocking_stage"]
                    == (
                        "chemical_typing"
                        if audit["labels"] is None
                        else "missing_charge_increment"
                    ),
                    "Wrong blocked stage",
                )
            else:
                with tempfile.TemporaryDirectory(prefix="island-c4-inputs-") as temp:
                    folder = Path(temp)
                    xyz = np.array(
                        [system.coordinates.get(i) for i in audit["atom_order"]]
                    )
                    reference_inputs(system, audit["labels"], charges, xyz, folder)
                    for kind in ("original", "trace"):
                        for n in ("reference.car", "reference.mdf"):
                            require(
                                sha(folder / n) == sha(case / kind / n),
                                "Changed atom/type/charge mapping",
                            )
                runs = {
                    kind: read(case / kind / "command.json")
                    for kind in ("original", "trace")
                }
                for kind, run in runs.items():
                    require(
                        "-ignore" not in run["argv"]
                        and "II" in run["argv"]
                        and "-nocenter" in run["argv"]
                        and Path(run["argv"][0]).parent.name == kind
                        and Path(run["argv"][run["argv"].index("-frc") + 1]).name
                        == variant.filename,
                        "Ignored or mixed-source converter",
                    )
                    require(
                        (case / kind / "converter.log").read_text()
                        == run["stdout"] + run["stderr"],
                        "Changed raw converter log",
                    )
                require(
                    runs["original"]["exit_code"]
                    == runs["trace"]["exit_code"]
                    == row["converter"]["exit_code"],
                    "Trace changed failure",
                )
                require(
                    runs["original"]["stdout"] == runs["trace"]["stdout"],
                    "Trace changed stdout",
                )
                stripped = (
                    "\n".join(
                        s
                        for s in runs["trace"]["stderr"].splitlines()
                        if not s.startswith("C4_MATCH ")
                    )
                    + "\n"
                )
                require(stripped == runs["original"]["stderr"], "Trace changed stderr")
                checked = trace_audit(runs["trace"]["stderr"], variant)
                require(
                    read(case / "full-converter-searches.json") == checked
                    and row["converter"]["full_converter_calls_checked"]
                    == len(checked),
                    "Changed full-converter reconstruction",
                )
                failure = next(
                    (
                        s.strip()
                        for s in runs["original"]["stdout"].splitlines()
                        if "Unable" in s
                    ),
                    None,
                )
                require(
                    row["converter"]["first_failure"] == failure
                    and row["converter"]["trace_preserves_output"] is True
                    and row["converter"]["chemical_typing_authority"] is False,
                    "Changed converter classification",
                )
                if runs["original"]["exit_code"] == 0:
                    require(
                        name == "silane_control" and failure is None,
                        "Unqualified new execution",
                    )
                    require(
                        sha(case / "original/reference.data")
                        == sha(case / "trace/reference.data"),
                        "Changed trace coefficients",
                    )
                    require(
                        row["retained_numerical_control"]
                        == silane_numerics(
                            system, variant, case / "original", None, execute=False
                        ),
                        "Changed raw energy/force control",
                    )
                    require(
                        row["converter"]["status"] == "conversion_completed"
                        and row["blocking_stage"]
                        == "retained_control_only_no_new_profile",
                        "Wrong control stage",
                    )
                else:
                    require(
                        failure
                        and audit["status"] == "source_incomplete"
                        and row["converter"]["status"] == "converter_failure"
                        and row["blocking_stage"] == "converter_required_source_key",
                        "Unclassified converter failure",
                    )
                if replay:
                    replay_converter(root, case, variant, runs["original"]["exit_code"])
                counts["converter_runs"] += 2
                counts["source_row_calls_checked"] += len(checked)
            counts[row["status"]] += 1
            findings.append(
                {
                    "source": variant.filename,
                    "name": name,
                    "status": row["status"],
                    "blocking_stage": row["blocking_stage"],
                    "converter": row["converter"],
                }
            )
    if replay:
        from scripts.validate_pcff_dynamics import reconstruct as continuation

        for kind in ("nve", "baoab"):
            continuation(
                root
                / "retained-c3/study/retained-c2/study/c1-replay/baseline"
                / ("silane-" + kind)
            )
    return {
        "independent_reconstruction": True,
        "acceptance_gate": 1,
        "counts": counts,
        "cases": findings,
        "silane_model_identity": retained["silane_model_identity"],
        "c2_decisions_preserved": retained["c2_decisions_preserved"],
        "c3_compiled_searches_preserved": retained["compiled_searches"],
    }


def validate(directory, *, repo=None):
    directory = Path(directory)
    repo = Path(repo) if repo else Path(__file__).resolve().parents[1]
    receipt = json.loads((directory / "receipt.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    require(
        receipt == signed(receipt) and receipt["schema"] == "island.j41-c4-receipt.v1",
        "Invalid receipt",
    )
    require(
        manifest == signed(manifest)
        and manifest["schema"] == "island.j41-c4-manifest.v1"
        and receipt["manifest_sha256"] == sha(directory / "manifest.json"),
        "Invalid manifest",
    )
    require(
        receipt["base_commit"] == manifest["base_commit"] == BASE
        and receipt["tested_implementation_commit"]
        == manifest["tested_implementation_commit"],
        "Wrong implementation binding",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["acceptance_gate"] == 1
        and receipt["new_execution_profiles"] == []
        and receipt["periodic"] == BLOCKS,
        "Promoted readiness",
    )
    require(
        set(manifest["files"]) == {"execution.tar.xz", "validation.tar.xz"},
        "Incomplete archive set",
    )
    for n, h in manifest["files"].items():
        require(sha(directory / n) == h, "Changed archive " + n)
    for key in ("historical_files", "implementation_files"):
        require(bool(manifest[key]), "Missing immutable inventory")
        for n, h in manifest[key].items():
            require(sha(repo / n) == h, "Changed " + n)
    expected = {
        "focused",
        "regressions",
        "full",
        "ruff",
        "compileall",
        "diff_check",
        "pip_supported",
        "pip_foyer",
        "pip_analysis",
        "pip_amber",
        "pip_quantum",
    }
    require(set(receipt["verification"]) == expected, "Missing verification")
    with tempfile.TemporaryDirectory(prefix="island-c4-receipt-") as temp:
        folder = Path(temp)
        extract(directory / "execution.tar.xz", folder / "raw")
        raw = folder / "raw"
        bundle = json.loads((raw / "bundle-manifest.json").read_text())
        require(
            bundle == signed(bundle)
            and set(bundle["files"])
            == {
                str(p.relative_to(raw))
                for p in raw.rglob("*")
                if p.is_file() and p != raw / "bundle-manifest.json"
            },
            "Changed portable inventory",
        )
        for n, h in bundle["files"].items():
            require(sha(raw / n) == h, "Changed raw artifact " + n)
        require(
            validate_study(raw / "study") == receipt["independent_validation"],
            "Changed study reconstruction",
        )
        extract(directory / "validation.tar.xz", folder / "checks")
        checks = folder / "checks/validation"
        for name, result in receipt["verification"].items():
            require(
                result["exit_code"] == 0
                and result["tested_implementation_commit"]
                == receipt["tested_implementation_commit"],
                "Failed or stale check " + name,
            )
            command = json.loads((checks / (name + "-command.json")).read_text())
            require(
                all(result[k] == v for k, v in command.items()),
                "Changed check command " + name,
            )
            if name in ("focused", "regressions", "full"):
                import xml.etree.ElementTree as ET

                counts = {k: 0 for k in ("tests", "failures", "errors", "skipped")}
                for suite in (
                    ET.parse(checks / (name + ".xml")).getroot().iter("testsuite")
                ):
                    for k in counts:
                        counts[k] += int(suite.attrib.get(k, 0))
                counts["passed"] = counts["tests"] - sum(
                    counts[k] for k in ("failures", "errors", "skipped")
                )
                require(
                    counts == result["junit"]
                    and counts["passed"]
                    and counts["errors"] == counts["failures"] == 0,
                    "Failed JUnit " + name,
                )
        mutations = json.loads((checks / "mutation-audit.json").read_text())
        require(
            mutations == receipt["audits"]["negative_controls"]
            and mutations["all_rejected"],
            "Missing mutation rejection",
        )
        require(
            receipt["audits"]["relocation"][
                "all_study_files_identical_after_relocation"
            ]
            is True,
            "Missing relocation",
        )
    return {
        "receipt_valid": True,
        "historical_files": len(manifest["historical_files"]),
        "acceptance_gate": 1,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--study", action="store_true")
    a = p.parse_args()
    print(
        json.dumps(
            validate_study(a.directory) if a.study else validate(a.directory), indent=2
        )
    )
