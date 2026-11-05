"""Fenced offline PCFF-compatible PDMS extension plans and blocked qualification."""

import argparse
import json
import subprocess
from pathlib import Path

from island.charge_references.records import unpack
from island.pcff_pdms.extension import FAMILY, FIXTURES, FLAGS, PROFILE
from island.pcff_pdms.extension.fitting import contract
from island.pcff_pdms.extension.inventory import inventory
from island.pcff_pdms.extension.qualification import candidate, reports
from island.pcff_pdms.extension.targets import plan, split_manifest
from scripts.phase_4j73_source_intake import (
    digest,
    exact,
    pin,
    read,
    require,
    signed,
    verify_signed,
    write,
)
from scripts.phase_4j76_source_resolution import scope

BASE = "b4354693d6f888fda442f6944d6af7700860e639"
BRANCH = "codex/phase-4j78-pcff-compatible-pdms-extension"
ROOT = Path("docs/evidence/phase_4j78")
FREEZE_SHA = "43a97d4f2bb3ea1c296c6e011d15c5d0794f714c4ba27661c9886a37d6f92b19"
IMPLEMENTATION = (
    "src/island/pcff_pdms/extension/__init__.py",
    "src/island/pcff_pdms/extension/inventory.py",
    "src/island/pcff_pdms/extension/targets.py",
    "src/island/pcff_pdms/extension/fitting.py",
    "src/island/pcff_pdms/extension/qualification.py",
    "scripts/phase_4j78_extension.py",
    "scripts/validate_phase_4j78_receipt.py",
    "tests/test_phase_4j78.py",
    "docs/phase_4j78.md",
)
INPUTS = (
    *(
        f"docs/evidence/phase_4j54/audit/{n}/{f}.json"
        for n in FIXTURES
        for f in ("system", "report", "manifest")
    ),
    "docs/evidence/phase_4j76/receipt.json",
    "docs/evidence/phase_4j76/native-pcff-coverage.json",
    "docs/evidence/phase_4j75/source-file-inventory.json",
    "docs/evidence/phase_4j72/receipt.json",
    "docs/evidence/phase_4j72_c1/quarantine-manifest.json",
)
SCOPE = {
    "J78_new_QM_attempts": 0,
    "J78_MD_steps": 0,
    "LAMMPS_executions": 0,
    "OpenMM_executions": 0,
    "converter_attempts": 0,
    "msi2lmp_executions": 0,
    "scientific_fitting_executions": 0,
    "scientific_parameter_generation_executions": 0,
    "scientific_evaluator_executions": 0,
    "reserved_validation_executions": 0,
    "fitted_coefficients": 0,
    "chemical_transfer": "not_established",
    "scientific_validation": "not_established",
    "native_PCFF_unresolved_instances": 378,
    "native_PCFF_exact_rows_resolved": 0,
    "J72_unintended_external_QM_attempts": 4,
    "J72_unintended_MD_steps_verified_minimum": 46,
    "J72_unintended_MD_steps_possible_maximum": 48,
    "J72_zero_execution_scope_met": False,
    "synthetic_matrix_controls_are_scientific_evidence": False,
    "author_confirmed_reported": 28,
    "author_confirmed_audited": 0,
    "original_author_artifact_status": "not_yet_supplied",
    **FLAGS,
}


def capture(repo, root):
    scope()
    require(not root.exists(), "do_not_overwrite_capture")
    original = {n: pin(repo / n) for n in INPUTS}
    graphs = {}
    for n in FIXTURES:
        directory = repo / f"docs/evidence/phase_4j54/audit/{n}"
        text = (directory / "system.json").read_text()
        unpack(text)  # validates the original data-only encoding/hash
        graphs[n] = {
            "original_system_text": text,
            "original_system_sha256": pin(directory / "system.json"),
            "original_report_identity": read(directory / "report.json")["identity"],
        }
    rows = subprocess.check_output(
        ["git", "ls-tree", "-r", "--full-tree", BASE], cwd=repo, text=True
    )
    history = {}
    inspected = {}
    replay_sources = []
    for row in rows.splitlines():
        meta, name = row.split("\t", 1)
        p = repo / name
        history[name] = (
            digest(str(p.readlink())) if meta.startswith("120000") else pin(p)
        )
        if name.startswith(
            (
                "src/island/pcff_pdms/",
                "src/island/pcff_qm/",
                "src/island/pcff_qm_pilot/",
                "src/island/forcefields/pcff/",
            )
        ) and name.endswith((".py", ".json")):
            inspected[name] = pin(p)
        if (name.startswith("src/island/") and name.endswith((".py", ".json"))) or (
            name.startswith("scripts/") and name.endswith(".py")
        ):
            replay_sources.append(name)
    sources = {
        "native_PCFF": {
            "status": "context_only_not_imported_as_extension_parameters",
            "inventory_sha256": original[
                "docs/evidence/phase_4j75/source-file-inventory.json"
            ],
            "J76_receipt_sha256": original["docs/evidence/phase_4j76/receipt.json"],
        },
        "external_literature": "not_imported",
        "analog_values": "not_imported",
        "COMPASS_OPLS_GAFF_CFF91_overlay": "forbidden",
        "target_origin": "J78_proposed_QM_fitted_extension; no_fitted_values",
    }
    write(
        root / "graph-inputs.json",
        signed({"graphs": graphs, "source_identities": sources}),
    )
    write(root / "historical-files.json", history)
    write(
        root / "implementation-reuse.json",
        signed(
            {
                "inspected_modules": inspected,
                "direct_reuse": {
                    "scope": "island.pcff_pdms.domain_diagnostics + core.Topology.rebuild_derived_interactions",
                    "records": "island.pcff_qm.records.Record and record; pcff_completion hash/require helpers",
                    "QM_protocol": "island.pcff_qm.targets.PROTOCOL as explicitly proposed context only",
                    "solver": "island.pcff_qm_pilot.hmds_fit.lsq_linear",
                    "conditioning": "island.pcff_qm_pilot.hmds_design.metrics",
                    "relocation": "island.pcff_qm_pilot.bundle data-only pack/unpack",
                },
                "adapter_rationale": "Historical fit entry points bind 58 columns/27 atoms or 9 atoms and source-bound evaluators. New staged matrix adapter reuses solver/diagnostics; no duplicate molecular fitter or evaluator. Historical module bytes and identity globs remain unchanged.",
                "native_registry_changed": False,
                "replay_source_files": replay_sources,
            }
        ),
    )
    write(
        root / "frozen-inputs.json",
        signed(
            {
                "schema": "island.j78-input-freeze.v1",
                "base": BASE,
                "fetched_origin_main": BASE,
                "original_inputs": original,
                "capture_files": {
                    n: pin(root / n)
                    for n in (
                        "graph-inputs.json",
                        "historical-files.json",
                        "implementation-reuse.json",
                    )
                },
                "scope": SCOPE,
            }
        ),
    )
    print(
        json.dumps(
            {
                "freeze_sha256": pin(root / "frozen-inputs.json"),
                "historical_files": len(history),
            }
        )
    )


def inputs(root):
    scope()
    require(
        pin(root / "frozen-inputs.json") == FREEZE_SHA,
        "fabricated_or_resigned_source_freeze",
    )
    f = read(root / "frozen-inputs.json")
    verify_signed(f)
    exact(f["scope"], SCOPE, "scientific_scope_mutation")
    require(f["base"] == f["fetched_origin_main"] == BASE, "base_mutation")
    for n, sha in f["capture_files"].items():
        require(pin(root / n) == sha, "captured_source_hash_mutation")
    captured = read(root / "graph-inputs.json")
    graphs = {
        n: unpack(captured["graphs"][n]["original_system_text"]) for n in FIXTURES
    }
    return graphs, captured["source_identities"]


def reconstruct(root):
    graphs, sources = inputs(root)
    inv = inventory(graphs)
    targets = plan(graphs, inv)
    split = split_manifest(targets)
    fitting = contract(inv)
    bundle = candidate(inv, targets, split, fitting, sources)
    return {
        "unique-term-inventory.json": inv,
        "qm-target-manifest.json": targets,
        "training-validation-split.json": split,
        "fitting-contract.json": fitting,
        "candidate-parameter-bundle.json": bundle,
        **reports(inv, targets, bundle, native_immutable=True),
        "capability-report.json": signed(
            {
                "forcefield_selector": "pcff",
                "profile": PROFILE,
                "family": FAMILY,
                "native_PCFF_registry_modified": False,
                "registered_derived_profiles": [],
                "selection_status": "fail_closed_not_qualified",
                "runtime_adapter": None,
                "training_scope_not_transfer_authority": True,
                **FLAGS,
            }
        ),
    }


def replay(root, *, build=False):
    documents = reconstruct(root)
    for name, value in documents.items():
        if build:
            require(not (root / name).exists(), "do_not_overwrite_attempt")
            write(root / name, value)
        exact(read(root / name), value, "resigned_extension_semantic_mutation:" + name)
    q = documents["qualification-report.json"]
    return {
        "decision": q["decision"],
        "unresolved_terms": q["unresolved_unique_terms"],
        "unresolved_scalar_parameters": q["unresolved_scalar_parameters"],
        "planned_QM_targets": q["planned_QM_targets"],
        "scope": SCOPE,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("capture", "build", "replay"))
    p.add_argument("--directory", type=Path, default=ROOT)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    args = p.parse_args()
    try:
        scope()
        result = (
            capture(args.repo, args.directory)
            if args.mode == "capture"
            else replay(args.directory, build=args.mode == "build")
        )
        print(json.dumps(result, sort_keys=True))
        return 0
    except (ValueError, RuntimeError, OSError, KeyError, TypeError) as e:
        print(json.dumps({"status": "blocked", "reason": str(e), "scope": SCOPE}))
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
