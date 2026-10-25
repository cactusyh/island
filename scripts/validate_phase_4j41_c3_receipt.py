"""Reconstruct C3 single-source decisions, raw numerics and immutable evidence."""

import argparse
import json
import tarfile
import tempfile
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import unpack
from island.forcefields.pcff.charges import identity
from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from island.forcefields.pcff.sip_completion import assess_pcff_sip_completion
from island.forcefields.pcff.sip_variants import (
    load_sip_source_variant,
    validate_sip_variant_inventory,
)
from island.workflows.bundle import system_from
from scripts.pcff_j5_reference import match, read_source
from scripts.pcff_j41_c2_reference import source_audit
from scripts.pcff_j41_c3_reference import audit_inventory, compass_trial, trial_searches
from scripts.pcff_j41_c3_silane import silane_numerics
from scripts.phase_4j41_c2 import read
from scripts.phase_4j41_c3 import CASES, source_comparison, sources
from scripts.validate_phase_4j41_c2_receipt import validate_study as validate_c2
from scripts.validate_phase_4j41_receipt import require, sha, signed


def reload_variant(folder, filename):
    p = unpack((folder / "inventory.json").read_text())
    return load_sip_source_variant(
        folder / filename,
        expected_sha256=p["sha256"],
        provenance=p["provenance"],
        companions=[(folder / n, h) for n, h in p["companion_hashes"].items()],
    )


def validate_study(root):
    import numpy as np

    from island.chemistry import from_smiles
    from island.workflows import storage
    from island.workflows.bundle import system_data
    from scripts.phase_4j41 import reference_inputs

    root = Path(root)
    report = read(root / "report.json")
    require(report["schema"] == "island.j41-c3-execution.v1", "Wrong study schema")
    require(
        report["acceptance_gate"] == 1
        and report["production_validated"] is False
        and report["simulation_readiness"] == "not_established"
        and report["new_execution_profiles"] == [],
        "Promoted C3 readiness",
    )
    require(
        report["ewald_status"] == "ewald_unvalidated"
        and report["pressure_status"] == "pressure_unavailable"
        and report["long_npt_status"] == "blocked",
        "Promoted periodic gates",
    )
    baseline = validate_c2(root / "retained-c2/study")
    require(len(baseline["cases"]) == 17, "Lost historical decisions")
    variants = [
        reload_variant(root / n.removesuffix(".frc"), n)
        for n in ("pcff.frc", "compass_published.frc")
    ]
    # Re-hashing a replacement does not authorize it for this execution matrix.
    with tempfile.TemporaryDirectory(prefix="island-c3-sources-") as temp:
        folder = Path(temp)
        for v in variants:
            (folder / v.filename).write_bytes(v.raw)
            for n, raw, _ in v.companions:
                (folder / n).write_bytes(raw)
        require(
            [v.identity for v in variants] == [v.identity for v in sources(folder)],
            "Changed authorized source hashes, provenance or companion files",
        )
    require(
        read(root / "source-comparison.json") == source_comparison(variants),
        "Changed source comparison",
    )
    require(
        [(r["source"], r["name"]) for r in report["cases"]]
        == [(v.filename, name) for v in variants for name, _, _ in CASES],
        "Incomplete or duplicate source/fixture matrix",
    )
    findings, searched = [], 0
    for v in variants:
        dest = root / v.filename.removesuffix(".frc")
        validate_sip_variant_inventory(dest / "inventory.json", v)
        require(
            read(dest / "raw-source-audit.json") == audit_inventory(v),
            "Changed independent source audit",
        )
        queries = {}
        for row, (name, smiles, domain) in zip(
            [r for r in report["cases"] if r["source"] == v.filename],
            CASES,
            strict=True,
        ):
            case = dest / name
            system = system_from(read(case / "system.json"))
            require(
                row["system_sha256"] == sha(case / "system.json")
                and row["atom_order"] == sorted(system.topology.sites)
                and row["variant_identity"] == v.identity,
                "Changed source, graph, coordinates or atom order",
            )
            prior = root / "retained-c2/study/c1-replay" / name / "system.json"
            if prior.exists():
                require(
                    sha(prior) == row["system_sha256"], "Changed retained final graph"
                )
            else:
                expected = storage.json_bytes(
                    storage.encode(system_data(from_smiles(smiles, random_seed=2026)))
                )
                require(
                    sha(case / "system.json") == sha256(expected).hexdigest(),
                    "Changed C3 phosphazene graph or coordinates",
                )
            labels, charges = None, None
            if v.filename == "pcff.frc":
                d = load_pcff_sip_distribution(dest)
                decision = assess_pcff_sip_completion(system, d, domain=domain)
                require(
                    (case / "decision.json").read_text() == decision.json_text
                    and row["decision_identity"] == decision.identity,
                    "Changed pinned decision",
                )
                old = root / "retained-c2/study" / name / "decision.json"
                if old.exists():
                    require(
                        identity(unpack(old.read_text()))
                        == decision.identity
                        == row["preserved_c2_identity"],
                        "Changed C2 identity",
                    )
                trial = source_audit(system, decision)
                qs = trial.get("unique_searches", [])
                source_complete = decision.payload["source_complete"]
                status = trial["status"]
                p = decision.payload["c1_assessment"]
                if p["inspection"]:
                    labels = p["typing"]["assignments"]
                    charges = {
                        i: float(q)
                        for i, q in trial["charges"]["partial_charges"].items()
                    }
            else:
                trial = compass_trial(system, v)
                qs = trial_searches(trial)
                source_complete = trial["source_complete"]
                status = (
                    "source_complete_unqualified"
                    if source_complete
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
            require(
                read(case / "raw-reconstruction.json") == trial
                and row["source_complete"] == source_complete
                and row["status"] == status,
                "Changed raw reconstruction or blocking stage",
            )
            require(
                row["model_published"] is False
                and row["dynamics_executed"] is False
                and row["execution_profile"] is None,
                "Unqualified model published",
            )
            for q in qs:
                queries[q["family"], q["namespace"], tuple(q["types"])] = q
            if labels is not None and name != "mixed_sip":
                xyz = np.array([system.coordinates.get(i) for i in sorted(labels)])
                with tempfile.TemporaryDirectory(prefix="island-c3-input-") as temp:
                    reference_inputs(system, labels, charges, xyz, Path(temp))
                    for n in ("reference.car", "reference.mdf"):
                        require(
                            sha(Path(temp) / n) == sha(case / n),
                            "Changed converter input labels, charges or graph mapping",
                        )
                command = read(case / "converter-command.json")
                require(
                    command == row["converter"]
                    and command["chemical_typing_authority"] is False,
                    "Changed converter command",
                )
                require(
                    Path(command["argv"][command["argv"].index("-frc") + 1]).name
                    == v.filename
                    and "-ignore" not in command["argv"],
                    "Mixed or ignored source in converter",
                )
                log = (case / "converter.log").read_text()
                if command["exit_code"] == 0:
                    require(
                        name == "silane_control" and "Error" not in log,
                        "Incomplete conversion promoted",
                    )
                    require(
                        row["numerical_control"]
                        == silane_numerics(system, v, case, None, execute=False),
                        "Changed independent energy/force evidence",
                    )
                else:
                    require(
                        not source_complete and ("Unable" in log or "ERROR" in log),
                        "Converter failure missing raw diagnostic",
                    )
            else:
                require(
                    row["converter"]
                    == {
                        "status": "not_executed",
                        "blocking_stage": "unqualified_typing_or_incomplete_charges",
                    },
                    "Guessed converter input",
                )
                require(
                    not (case / "reference.data").exists(), "Blocked fixture converted"
                )
            findings.append(
                {
                    "source": v.filename,
                    "name": name,
                    "status": status,
                    "source_complete": source_complete,
                }
            )
        raw = read_source(v.raw)
        indexed = {r["id"]: r for rs in raw.values() for r in rs}
        compiled = read(dest / "compiled-searches.json")
        require(
            [(q["family"], q["namespace"], tuple(q["types"])) for q in compiled]
            == sorted(queries),
            "Incomplete compiled source searches",
        )
        for q in compiled:
            oracle = match(
                raw, q["family"], q["namespace"], q["types"], guarded_msi=True
            )
            require(q["independent"] == oracle, "Changed independent match")
            selected = q["compiled_row"]
            require(
                (selected is None and oracle["status"] == "missing")
                or selected in oracle.get("selected_ids", oracle["candidate_ids"]),
                "Wrong compiled selection",
            )
            require(
                q["compiled_values"]
                == (list(map(str, indexed[selected]["values"])) if selected else None),
                "Changed source values",
            )
        searched += len(compiled)
    return {
        "independent_reconstruction": True,
        "cases": findings,
        "compiled_searches": searched,
        "c2_decisions_preserved": 17,
        "silane_model_identity": baseline["silane_model_identity"],
        "acceptance_gate": 1,
    }


def extract(archive, target):
    with tarfile.open(archive) as tar:
        names = set()
        for member in tar.getmembers():
            require(
                member.name not in names
                and (member.isfile() or member.isdir())
                and not Path(member.name).is_absolute()
                and ".." not in Path(member.name).parts,
                "Unsafe or duplicate archive entry",
            )
            names.add(member.name)
        tar.extractall(target, filter="data")


def validate(directory, *, repo=None):
    directory = Path(directory)
    repo = Path(repo) if repo else Path(__file__).resolve().parents[1]
    receipt = json.loads((directory / "receipt.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    require(
        receipt == signed(receipt) and receipt["schema"] == "island.j41-c3-receipt.v1",
        "Invalid receipt",
    )
    require(
        manifest == signed(manifest)
        and manifest["schema"] == "island.j41-c3-manifest.v1",
        "Invalid manifest",
    )
    require(
        receipt["manifest_sha256"] == sha(directory / "manifest.json"),
        "Changed manifest",
    )
    require(
        receipt["acceptance_gate"] == 1
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    require(
        set(manifest["files"]) == {"execution.tar.xz", "validation.tar.xz"},
        "Incomplete evidence archives",
    )
    for n, h in manifest["files"].items():
        require(sha(directory / n) == h, "Changed evidence " + n)
    for key in ("historical_files", "implementation_files"):
        require(bool(manifest[key]), "Empty immutable inventory")
        for n, h in manifest[key].items():
            require(sha(repo / n) == h, "Changed " + key + ": " + n)
    require(
        set(receipt["verification"])
        == {
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
        and all(v["exit_code"] == 0 for v in receipt["verification"].values()),
        "Missing or failed verification",
    )
    with tempfile.TemporaryDirectory(prefix="island-j41-c3-") as temp:
        root = Path(temp)
        extract(directory / "execution.tar.xz", root / "execution")
        raw = root / "execution"
        bundle = json.loads((raw / "bundle-manifest.json").read_text())
        require(
            bundle == signed(bundle)
            and set(bundle["files"])
            == {
                str(p.relative_to(raw))
                for p in raw.rglob("*")
                if p.is_file() and p != raw / "bundle-manifest.json"
            },
            "Changed raw bundle inventory",
        )
        for n, h in bundle["files"].items():
            require(sha(raw / n) == h, "Changed raw artifact " + n)
        require(
            validate_study(raw / "study") == receipt["independent_validation"],
            "Changed independent reconstruction",
        )
        extract(directory / "validation.tar.xz", root / "checks")
        checks = root / "checks/validation"
        for n in ("focused", "regressions", "full"):
            import xml.etree.ElementTree as ET

            counts = {k: 0 for k in ("tests", "failures", "errors", "skipped")}
            for suite in ET.parse(checks / (n + ".xml")).getroot().iter("testsuite"):
                for k in counts:
                    counts[k] += int(suite.attrib.get(k, 0))
            counts["passed"] = counts["tests"] - sum(
                counts[k] for k in ("failures", "errors", "skipped")
            )
            require(
                counts == receipt["verification"][n]["junit"]
                and counts["passed"] > 0
                and counts["failures"] == counts["errors"] == 0,
                "Invalid JUnit results",
            )
    return {
        "receipt_valid": True,
        "historical_files": len(manifest["historical_files"]),
        "acceptance_gate": 1,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--study", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            validate_study(args.directory) if args.study else validate(args.directory),
            indent=2,
        )
    )
