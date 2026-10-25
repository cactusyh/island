"""Validate C2 raw decisions, compiled searches, blocked gates and immutable bytes."""

import argparse
import json
import tarfile
import tempfile
from pathlib import Path

from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from island.forcefields.pcff.sip_completion import load_pcff_sip_completion
from island.workflows.bundle import system_from
from scripts.pcff_j5_reference import match, read_source
from scripts.pcff_j41_c2_reference import source_audit
from scripts.phase_4j41_c2 import read
from scripts.validate_phase_4j41_c1_receipt import validate_study as validate_c1
from scripts.validate_phase_4j41_receipt import require, sha, signed


def validate_study(root):
    root = Path(root)
    d = load_pcff_sip_distribution(root)
    report = read(root / "report.json")
    require(report["schema"] == "island.j41-c2-execution.v1", "Wrong study schema")
    require(report["distribution_hashes"] == d.hashes, "Changed source hashes")
    require(
        report["acceptance_gate"] == 1
        and report["production_validated"] is False
        and report["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    require(
        report["ewald_status"] == "ewald_unvalidated"
        and report["pressure_status"] == "pressure_unavailable"
        and report["long_npt_status"] == "blocked",
        "Promoted periodic gates",
    )
    baseline = validate_c1(root / "c1-replay")
    previous = read(root / "c1-replay/report.json")
    expected = {r["name"]: r for r in previous["cases"]}
    require(
        set(expected) == {r["name"] for r in report["cases"]}
        and len(expected) == len(report["cases"]),
        "Changed retained matrix",
    )
    queries, results = {}, []
    for row in report["cases"]:
        name = row["name"]
        system = system_from(read(root / "c1-replay" / name / "system.json"))
        decision = load_pcff_sip_completion(
            root / name / "decision.json", d, system=system
        )
        require(
            row["decision_identity"] == decision.identity
            and row["c1_identity"]
            == decision.payload["c1_identity"]
            == expected[name]["assessment_identity"],
            "Changed identity",
        )
        require(
            row["graph_identity"] == decision.payload["c1_assessment"]["graph_identity"]
            and row["source_complete"]
            == decision.payload["source_complete"]
            == expected[name]["source_complete"],
            "Changed source completeness",
        )
        require(
            row["execution_profile"] == decision.payload["execution_profile"]
            and row["blocking_interaction_count"]
            == len(decision.payload["blocking_interactions"]),
            "Changed source decision",
        )
        audit = source_audit(system, decision)
        require(
            read(root / name / "raw-source-audit.json") == audit
            and row["typing_status"] == audit["status"],
            "Changed raw-source reconstruction",
        )
        require(
            row["converter"] == expected[name]["converter"],
            "Changed converter evidence",
        )
        require(
            row["model_published"] is False
            and row["new_domain_dynamics_executed"] is False
            and not (root / name / "model.json").exists(),
            "Incomplete domain published",
        )
        require(
            row["workflow_status"]
            == (
                "unchanged_silane_fresh_reference_replay_passed"
                if name == "silane_control"
                else "blocked_before_model_publication"
            ),
            "Changed workflow state",
        )
        if name == "silane_control":
            require(
                row["model_identity"] == baseline["silane_model_identity"],
                "Changed silane identity",
            )
        for q in audit.get("unique_searches", []):
            queries[q["family"], q["namespace"], tuple(q["types"])] = q
        results.append(
            {
                "name": name,
                "identity": decision.identity,
                "source_complete": row["source_complete"],
            }
        )
    compiled = read(root / "compiled-completion-searches.json")
    require(
        len(compiled) == len(queries) == report["compiled_search_count"],
        "Incomplete compiled search matrix",
    )
    require(
        sha(root / "matcher-leaf.so") == report["matcher_sha256"],
        "Changed compiled matcher",
    )
    raw = read_source(d.source.raw)
    indexed = {r["id"]: r for rs in raw.values() for r in rs}
    require(
        [(q["family"], q["namespace"], tuple(q["types"])) for q in compiled]
        == sorted(queries),
        "Changed compiled query ordering",
    )
    for q in compiled:
        oracle = match(raw, q["family"], q["namespace"], q["types"], guarded_msi=True)
        require(oracle == q["independent"], "Changed independent match")
        selected = q["compiled_row"]
        require(
            (selected is None and oracle["status"] == "missing")
            or selected
            in oracle[
                "selected_ids" if oracle["status"] == "assigned" else "candidate_ids"
            ],
            "Wrong compiled selection",
        )
        require(
            q["compiled_values"]
            == (list(map(str, indexed[selected]["values"])) if selected else None),
            "Changed compiled coefficient values",
        )
    return {
        "independent_reconstruction": True,
        "silane_model_identity": baseline["silane_model_identity"],
        "cases": results,
        "compiled_searches": len(compiled),
        "acceptance_gate": 1,
    }


def validate(directory, *, repo=None):
    directory = Path(directory)
    repo = Path(repo) if repo else Path(__file__).resolve().parents[1]
    receipt = json.loads((directory / "receipt.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    require(
        receipt == signed(receipt) and receipt["schema"] == "island.j41-c2-receipt.v1",
        "Invalid receipt",
    )
    require(
        manifest == signed(manifest)
        and manifest["schema"] == "island.j41-c2-manifest.v1",
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
        "Missing evidence",
    )
    for name, checksum in manifest["files"].items():
        require(sha(directory / name) == checksum, "Changed evidence " + name)
    for group in ("historical_files", "implementation_files"):
        require(bool(manifest[group]), "Empty preservation inventory")
        for name, checksum in manifest[group].items():
            require(sha(repo / name) == checksum, "Changed " + group + ": " + name)
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
    require(
        set(receipt["verification"]) == expected
        and all(v["exit_code"] == 0 for v in receipt["verification"].values()),
        "Required verification failed",
    )
    require(
        all(
            receipt["verification"][name]["junit"]["passed"] > 0
            and not receipt["verification"][name]["junit"]["failures"]
            and not receipt["verification"][name]["junit"]["errors"]
            for name in ("focused", "regressions", "full")
        ),
        "Invalid JUnit results",
    )
    with tempfile.TemporaryDirectory(prefix="island-j41-c2-") as temp:
        with tarfile.open(directory / "execution.tar.xz") as archive:
            for member in archive.getmembers():
                require(
                    (member.isfile() or member.isdir())
                    and not member.issym()
                    and not member.islnk()
                    and not Path(member.name).is_absolute()
                    and ".." not in Path(member.name).parts,
                    "Unsafe archive",
                )
            archive.extractall(temp, filter="data")
        root = Path(temp)
        bundle = json.loads((root / "bundle-manifest.json").read_text())
        require(bundle == signed(bundle), "Changed portable manifest")
        require(
            set(bundle["files"])
            == {
                str(p.relative_to(root))
                for p in root.rglob("*")
                if p.is_file() and p.name != "bundle-manifest.json"
            },
            "Incomplete raw inventory",
        )
        for name, checksum in bundle["files"].items():
            require(sha(root / name) == checksum, "Changed raw artifact")
        require(
            validate_study(root / "study") == receipt["independent_validation"],
            "Changed independent audit",
        )
    return {
        "receipt_valid": True,
        "historical_files": len(manifest["historical_files"]),
        "acceptance_gate": 1,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--study", action="store_true")
    a = p.parse_args()
    print(
        json.dumps(
            validate_study(a.directory) if a.study else validate(a.directory), indent=2
        )
    )


if __name__ == "__main__":
    main()
