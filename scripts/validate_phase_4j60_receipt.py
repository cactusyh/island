"""Seal/replay J60 numerical evidence with the existing history/JUnit helpers."""

import argparse
import json
import lzma
from collections import Counter
from hashlib import sha256
from pathlib import Path

from island.pcff_completion import digest, pin, require, signed, write
from scripts.phase_4j60 import replay
from scripts.seal_phase_4j51_receipt import junit
from scripts.validate_phase_4j47_receipt import history
from scripts.validate_phase_4j47_receipt import pin as historical_pin

BASE = "b4aecb4c601ea4b7d24cec9eba17d92747280297"
BRANCH = "codex/phase-4j60-pcff-hamiltonian-identifiable-pilot-fit"


def read(path):
    r = json.loads(Path(path).read_text())
    require(
        r["identity"] == digest({k: v for k, v in r.items() if k != "identity"}),
        "record_identity_mutation",
        file=str(path),
    )
    return r


def source_checks(root):
    definition = read(root / "hamiltonian-definition.json")
    for source in definition["sources"].values():
        path = root / "sources" / (source["file"].replace("/", "__") + ".xz")
        require(
            sha256(lzma.decompress(path.read_bytes())).hexdigest() == source["sha256"],
            "source_replacement",
        )
    for convention in definition["conventions"]:
        require(
            convention["identity"]
            == digest({k: v for k, v in convention.items() if k != "identity"}),
            "convention_mutation",
        )
        for evidence in convention["evidence"]:
            raw = lzma.decompress(
                (
                    root / "sources" / (evidence["file"].replace("/", "__") + ".xz")
                ).read_bytes()
            )
            lo, hi = evidence["line_range"]
            require(
                sha256(raw).hexdigest() == evidence["sha256"]
                and "\n".join(raw.decode().splitlines()[lo - 1 : hi])
                == evidence["raw_text"],
                "source_location_mutation",
            )
    require(
        definition["author_confirmations"]["original_author_artifact_status"]
        == "not_yet_supplied",
        "author_provenance_promoted",
    )
    require(
        not definition["registered_profile"] and not definition["native_PCFF_claim"],
        "profile_promoted",
    )
    return definition


def implementation(repo):
    return {
        str(p.relative_to(repo)): pin(p)
        for p in sorted(
            list((repo / "src/island/pcff_qm_pilot").glob("hamiltonian*.py"))
            + list((repo / "scripts").glob("*phase_4j60*.py"))
            + [repo / "tests/test_phase_4j60.py", repo / "docs/phase_4j60.md"]
        )
    }


def artifacts(root):
    return {
        str(p.relative_to(root)): pin(p)
        for p in sorted(root.rglob("*"))
        if p.is_file() and p.name not in ("receipt.json", "manifest.json")
    }


def audit(root, repo, do_replay=True):
    definition = source_checks(root)
    old = history(repo, BASE)
    require(
        all(historical_pin(repo / name) == h for name, h in old.items()),
        "historical_bytes_changed",
    )
    criteria = read(root / "plan/reference-plan.json")
    require(
        criteria["budget"]["maximum_gradient_jobs"] == 35, "compute_budget_mutation"
    )
    design = read(root / "joint-target-design.json")
    require(
        design["independent_coefficient_count"] == 20
        and design["selected_complete_design"]["rank"] == 20,
        "partial_design",
    )
    require(
        not {s["geometry_identity"] for s in design["selected_training"]}
        & {s["geometry_identity"] for s in design["held_out"]},
        "held_out_leakage",
    )
    decision = read(root / "fitting-decision.json")
    reference = read(root / "reference-convergence.json")
    expected_gates = {
        c["effect"]: c["predeclared_gate_passed"]
        for c in reference["comparisons"]
        if c["predeclared_gate_passed"] is not None
    }
    require(
        decision["reference_gates"] == expected_gates
        and decision["basis_gate_numbers"]
        == next(c for c in reference["comparisons"] if c["effect"] == "cardinal_QZ_5Z"),
        "reference_decision_mutation",
    )
    fit_path = root / "experimental-joint-fit.json"
    require(decision["fit_completed"] == fit_path.exists(), "fit_completion_mutation")
    if fit_path.exists():
        fit = read(fit_path)
        before_fit = read(root / "plan/reference-gate-before-fit.json")
        require(
            before_fit["reference_identity"] == reference["identity"]
            and before_fit["reference_gates_passed"] is True
            and before_fit["fitting_values_used_for_gate"] is False
            and before_fit["fit_native_results_present"] == 0
            and all(
                pin(root / k) == v for k, v in before_fit["evidence_sha256"].items()
            ),
            "pre_fit_reference_authority_mutation",
        )
        require(
            set(fit["coefficients"]) == set(design["columns"])
            and decision["fitted_coefficients"] == list(fit["coefficients"])
            and all(expected_gates.values()),
            "fit_parameter_or_gate_mutation",
        )
    else:
        require(decision["fitted_coefficients"] == [], "absent_fit_claim")
    require(
        decision["parameter_completeness"] == "incomplete"
        and decision["production_validated"] is False
        and decision["simulation_readiness"] == "not_established",
        "readiness_promoted",
    )
    checks = json.loads((root / "validation/checks.json").read_text())
    require(all(r["returncode"] == 0 for r in checks), "required_software_check_failed")
    suites = json.loads((root / "validation/junit-files.json").read_text())
    counts, nodes = junit(root, suites)
    require(len(set(nodes)) == len(nodes), "repeated_suite_nodes")
    collection = json.loads((root / "collection.json").read_text())["nodes"]
    require(Counter(nodes) == Counter(collection), "suite_coverage_mismatch")
    if do_replay:
        replay(repo, root)
    return {
        "sources": definition["identity"],
        "historical_files": old,
        "software_counts": counts,
        "focused": sum(n.startswith("tests/test_phase_4j60.py::") for n in nodes),
        "J54_J59_regressions": sum(
            any(n.startswith(f"tests/test_phase_4j{j}") for j in range(54, 60))
            for n in nodes
        ),
        "decision_identity": decision["identity"],
        "QM_jobs": {
            "primary": len(
                json.loads((root / "plan/campaign.json").read_text())["jobs"]
            ),
            **{
                p.parent.parent.name: len(json.loads(p.read_text())["jobs"])
                for p in sorted((root / "stages").glob("*/plan/campaign.json"))
            },
        },
    }


def seal(root, repo):
    a = audit(root, repo)
    manifest = signed(
        {
            "schema": "island.j60-manifest.v1",
            "base_commit": BASE,
            "historical_files": a["historical_files"],
            "implementation_files": implementation(repo),
            "artifacts": artifacts(root),
        }
    )
    write(root / "manifest.json", manifest)
    receipt = signed(
        {
            "schema": "island.j60-receipt.v1",
            "base_commit": BASE,
            "J59_prerequisite": "9681e4926d92f5a38ea9dbc329c15c1207d819cd",
            "branch": BRANCH,
            "manifest_sha256": pin(root / "manifest.json"),
            "hamiltonian_identity": a["sources"],
            "fitting_decision_identity": a["decision_identity"],
            "software_counts": a["software_counts"],
            "focused_tests": a["focused"],
            "J54_J59_regressions": a["J54_J59_regressions"],
            "fitted_coefficients": read(root / "fitting-decision.json")[
                "fitted_coefficients"
            ],
            "accepted_coefficients": [],
            "QM_jobs": a["QM_jobs"],
            "new_energy_gradient_calculations": sum(a["QM_jobs"].values()),
            "fit_completed": True,
            "fit_validated": False,
            "fit_decision": "experimental_fit_rejected_by_untouched_validation",
            "user_attested_confirmations": 28,
            "original_author_artifact_status": "not_yet_supplied",
            "new_registered_profiles": [],
            "PCFF_MD_executed": False,
            "parameter_completeness": "incomplete",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    write(root / "receipt.json", receipt)
    return {
        "receipt_sha256": pin(root / "receipt.json"),
        "software_counts": a["software_counts"],
        "historical_files": len(a["historical_files"]),
    }


def validate(root, repo, expected=None):
    if expected:
        require(
            pin(root / "receipt.json") == expected, "independent_receipt_hash_mismatch"
        )
    r = read(root / "receipt.json")
    m = read(root / "manifest.json")
    require(
        r["base_commit"] == m["base_commit"] == BASE and r["branch"] == BRANCH,
        "base_branch_mutation",
    )
    require(pin(root / "manifest.json") == r["manifest_sha256"], "manifest_mutation")
    require(
        m["artifacts"] == artifacts(root)
        and m["implementation_files"] == implementation(repo),
        "artifact_implementation_mutation",
    )
    a = audit(root, repo)
    require(
        m["historical_files"] == a["historical_files"]
        and r["software_counts"] == a["software_counts"],
        "history_software_mutation",
    )
    require(
        r["fitting_decision_identity"] == a["decision_identity"]
        and r["hamiltonian_identity"] == a["sources"],
        "result_identity_mutation",
    )
    decision = read(root / "fitting-decision.json")
    require(
        r["fitted_coefficients"] == decision["fitted_coefficients"]
        and r["accepted_coefficients"] == decision["accepted_coefficients"] == []
        and bool(r["fitted_coefficients"]) == decision["fit_completed"]
        and r["parameter_completeness"] == "incomplete"
        and r["cross_engine_equivalence"] == "not_claimed"
        and r["scientific_validation"] == "not_established"
        and r["production_validated"] is False
        and r["simulation_readiness"] == "not_established"
        and r["new_registered_profiles"] == []
        and r["PCFF_MD_executed"] is False
        and r["focused_tests"] == a["focused"]
        and r["J54_J59_regressions"] == a["J54_J59_regressions"],
        "receipt_claim_mutation",
    )
    require(
        r["QM_jobs"] == a["QM_jobs"]
        and r["new_energy_gradient_calculations"] == sum(a["QM_jobs"].values()) == 87
        and r["fit_completed"] == decision["fit_completed"] is True
        and r["fit_validated"] == decision["fit_validated"] is False
        and r["fit_decision"] == decision["decision"]
        and r["user_attested_confirmations"] == 28
        and r["original_author_artifact_status"] == "not_yet_supplied",
        "numerical_or_author_claim_mutation",
    )
    return {
        "status": "passed",
        "receipt_sha256": pin(root / "receipt.json"),
        "historical_files": len(a["historical_files"]),
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--seal", action="store_true")
    p.add_argument("--expected-sha256")
    a = p.parse_args()
    print(
        json.dumps(
            seal(a.directory.resolve(), a.repo.resolve())
            if a.seal
            else validate(a.directory.resolve(), a.repo.resolve(), a.expected_sha256),
            sort_keys=True,
        )
    )
