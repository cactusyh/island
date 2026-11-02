"""Additive quarantine reconstruction; no literature decisions or science run.

The independent checks use retained bytes, explicit scientific dependency
inventories and original producer mappings. Audit ledgers may describe the
quarantine; scientific inputs may not depend on it.
"""

import hashlib
import json
import lzma
import re
import tarfile
from collections import Counter
from pathlib import Path

from island.pcff_completion import pin, require, signed
from island.pcff_qm_pilot.source_adjudication import verify_signed

BASE = REQUIRED = "30737605b8252e2bf0f2265352f0937ee931fd6d"
BRANCH = "codex/phase-4j72-c1-execution-scope-correction"
ROOT = Path("docs/evidence/phase_4j72_c1")
PARENT = Path("docs/evidence/phase_4j72")
DEVIATION_SHA = "b1de0425346d6fbfe67f5eb66f2453750269da3cbbfbed6db705972627fb46dd"
ARCHIVE_SHA = "43149b96791fb7963b89405616636b1e8dc161a97897b5e28bb431a98ca20236"
# Filled once when the additive input freeze is written, before verification.
CONTRACT_SHA = "038da2cc8a212c21de3548360a10dce203c81e4f628bc969661ae63560e1ced3"

DECLARATIONS = {
    "authorized_new_QM": 0,
    "authorized_MD_steps": 0,
    "HMDS_new_QM": 0,
    "reserved_validation_QM": 0,
    "unintended_external_QM_attempts": 4,
    "unintended_MD_steps_verified_minimum": 46,
    "unintended_MD_steps_possible_maximum": 48,
    "zero_execution_scope_met": False,
    "scientific_evidence_contaminated": False,
    "quarantine_complete": True,
    "source_adjudication_excluded_quarantine": True,
    "C1_new_native_QM_attempts": 0,
    "C1_MD_steps": 0,
    "C1_scientific_fits": 0,
    "C1_converter_attempts": 0,
    "C1_new_papers_read": 0,
    "decision": "source_conflict_unresolved",
    "accepted_coefficients": 0,
    "automatic_assignment": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "parameter_completeness": "incomplete",
    "chemical_transfer": "not_established",
    "scientific_validation": "not_established",
    "original_author_artifact_status": "not_yet_supplied",
    "author_confirmations": "28 user-attested",
    "reserved_validation_activated": False,
    "converter_complete": False,
    "static_complete_closure": False,
}

SCIENTIFIC_FILES = (
    "docs/evidence/phase_4j70/completeness-matrix.json",
    "docs/evidence/phase_4j71/parameter-records.json",
    "docs/evidence/phase_4j71/extraction-ledger.json",
    "docs/evidence/phase_4j71/reading-status.json",
    "docs/evidence/phase_4j71/source-evidence-table.json",
    "docs/evidence/phase_4j71/partial-literature-bundle.json",
    "docs/evidence/phase_4j71/mapping-report.json",
    "docs/evidence/phase_4j72/comparison-normalized-records.json",
    "docs/evidence/phase_4j72/authority-hierarchy.json",
    "docs/evidence/phase_4j72/source-conflict-matrix.json",
    "docs/evidence/phase_4j72/row-adjudication.json",
    "docs/evidence/phase_4j72/mathematical-resolutions.json",
    "docs/evidence/phase_4j72/convention-adjudications.json",
    "docs/evidence/phase_4j72/unresolved-evidence-requirements.json",
    "docs/evidence/phase_4j72/candidate-source-bundles.json",
    "docs/evidence/phase_4j72/static-closure-report.json",
    "docs/evidence/phase_4j72/strict-converter-report.json",
    "docs/evidence/phase_4j72/canonical-capability.json",
    "docs/evidence/phase_4j67/fitting-dataset.json",
    "docs/evidence/phase_4j67/independent-native-checks.json",
    "docs/evidence/phase_4j67/compute-accounting.json",
    "docs/evidence/phase_4j67/validation-activation.json",
)


def read(path):
    return json.loads(Path(path).read_text())


def typed_equal(actual, expected, reason):
    require(type(actual) is type(expected), reason)
    if isinstance(expected, dict):
        require(set(actual) == set(expected), reason)
        for key in expected:
            typed_equal(actual[key], expected[key], reason + ":" + key)
    elif isinstance(expected, list):
        require(len(actual) == len(expected), reason)
        for a, b in zip(actual, expected, strict=True):
            typed_equal(a, b, reason)
    else:
        require(actual == expected, reason)


def quarantine(repo):
    directory = repo / PARENT / "scope-deviation"
    require(
        pin(directory / "execution-deviation.json") == DEVIATION_SHA,
        "C1_original_scope_violation_modified",
    )
    require(
        pin(directory / "quarantined-j31-regressions.tar.xz") == ARCHIVE_SHA,
        "C1_original_quarantine_modified",
    )
    original = read(directory / "execution-deviation.json")
    raw, members = {}, {}
    with tarfile.open(directory / original["archive"], "r:xz") as archive:
        for member in archive:
            require(
                member.isfile()
                and member.name not in raw
                and not member.name.startswith("/")
                and ".." not in Path(member.name).parts,
                "C1_unsafe_or_duplicate_quarantine_member",
            )
            data = archive.extractfile(member).read()
            raw[member.name] = data
            members[member.name] = {
                "sha256": hashlib.sha256(data).hexdigest(),
                "bytes": len(data),
            }
    typed_equal(members, original["files"], "C1_quarantine_members")
    require(len(raw) == 192, "C1_quarantine_file_count")
    attempts = []
    for name in sorted(raw):
        if not name.endswith("/sqm.out"):
            continue
        prefix = name.rsplit("/", 1)[0]
        stdout = raw[prefix + "/antechamber.stdout.log"].decode()
        commands = [
            line.removeprefix("Running: ")
            for line in stdout.splitlines()
            if line.startswith("Running: ") and "/sqm " in line
        ]
        require(
            len(commands) == 1 and b"Calculation Completed" in raw[name],
            "C1_unintended_attempt_reconstruction",
        )
        attempts.append(
            {
                "case": int(re.search(r"source_c(\d+)/", name)[1]),
                "command_verbatim_from_stdout": commands[0],
                "process_id": None,
                "process_id_status": "not_retained_in_original_evidence",
                "native_executable_sha256": None,
                "native_executable_hash_status": "not_retained_at_execution_time",
                "stdout_command_log": prefix + "/antechamber.stdout.log",
                "parent_stderr": prefix + "/antechamber.stderr.log",
                "native_output": name,
                "native_input": prefix + "/sqm.in",
                "native_separate_stderr_status": "not_independently_captured",
                "execution_namespace": "unintended_J31_CHO_AM1_BCC_test_execution",
                "authorized": False,
                "scientific_use": False,
            }
        )
    require(
        [a["case"] for a in attempts] == [2, 3, 6, 7], "C1_four_unintended_attempts"
    )
    md = []
    for case in (1, 2, 3, 5, 6, 7):
        prefix = f"test_real_final_graph_source_c{case}/acceptance/"
        maxima = {}
        record_paths = []
        for name, data in raw.items():
            if (
                name.startswith(prefix)
                and "/generation-" in name
                and name.endswith("/record.json")
            ):
                branch = name[len(prefix) :].split("/")[0]
                frames = json.loads(data)["payload"]["frames"]
                maxima[branch] = max(
                    maxima.get(branch, 0), *(f["step"] for f in frames)
                )
                record_paths.append(name)
        require(
            maxima == {"full": 4, "relocated": 2 if case == 7 else 4, "split": 2},
            "C1_MD_journal_reconstruction",
        )
        md.append(
            {
                "case": case,
                "recorded_steps": 4 + 2 + maxima["relocated"] - 2,
                "uncommitted_steps_possible": 2 if case == 7 else 0,
                "journal_paths": sorted(record_paths),
            }
        )
    require(sum(m["recorded_steps"] for m in md) == 46, "C1_MD_minimum")
    return signed(
        {
            "schema": "island.j72-c1-quarantine-manifest.v1",
            "archive_path": str(PARENT / "scope-deviation" / original["archive"]),
            "archive_sha256": ARCHIVE_SHA,
            "original_deviation_sha256": DEVIATION_SHA,
            "members": members,
            "member_count": 192,
            "unintended_attempts": attempts,
            "MD_cases": md,
            "verified_MD_steps": 46,
            "possible_MD_steps": 48,
            "interruption_uncertainty": original["MD_step_uncertainty"],
            "complete_scope": "all 192 originally retained quarantine members; no claim that unrecorded PIDs, native stderr or uncommitted frames can be recovered",
            "quarantine_complete": True,
            "allowed_role": "integrity_audit_only",
            "original_contract_sha256": pin(repo / PARENT / "contract.json"),
            "original_receipt_sha256": pin(repo / PARENT / "receipt.json"),
            "original_scope_report_sha256": pin(
                repo / PARENT / "execution-scope-report.json"
            ),
            "original_launch_records": {
                str(PARENT / "validation" / name): pin(
                    repo / PARENT / "validation" / name
                )
                for name in ("start_suite.py", "suite-results.json", "start-suite.log")
            },
        }
    )


def reject_quarantine(value, manifest):
    """Reject quarantined paths, byte hashes and derived output identities.

    This function applies to scientific inputs only. The corrective receipt,
    execution-scope ledger and integrity manifest intentionally name quarantine.
    """
    forbidden_hashes = {v["sha256"] for v in manifest["members"].values() if v["bytes"]}
    forbidden_hashes.add(manifest["archive_sha256"])

    def walk(item):
        if isinstance(item, dict):
            for key, data in item.items():
                walk(key)
                walk(data)
        elif isinstance(item, list):
            for data in item:
                walk(data)
        elif isinstance(item, str):
            require(
                not any(
                    s in item.replace("\\", "/").lower()
                    for s in (
                        "scope-deviation/",
                        "quarantined-j31",
                        "test_real_final_graph_source_c",
                        "island-ambertools-",
                        "island-j72-validation/pytest-terminal-0",
                    )
                ),
                "C1_quarantine_path_in_scientific_evidence",
            )
            hashes = set(re.findall(r"[0-9a-f]{64}", item))
            require(
                not hashes.intersection(forbidden_hashes),
                "C1_quarantine_hash_in_scientific_evidence",
            )

    walk(value)


def exclusion(repo, manifest, documents=None):
    docs = (
        documents
        if documents is not None
        else {n: read(repo / n) for n in SCIENTIFIC_FILES}
    )
    require(set(docs) == set(SCIENTIFIC_FILES), "C1_scientific_dependency_inventory")
    for value in docs.values():
        reject_quarantine(value, manifest)
    records = docs["docs/evidence/phase_4j71/parameter-records.json"]["records"]
    normalized = docs[str(PARENT / "comparison-normalized-records.json")]["records"]
    require(len(records) == len(normalized) == 467, "C1_literature_records")
    original = {r["record_id"]: r for r in records}
    for record in normalized:
        require(
            record["J71_identity"] == original[record["record_id"]]["identity"]
            and record["source_file_sha256"]
            == original[record["record_id"]]["source_file_sha256"],
            "C1_original_literature_producer",
        )
        parent = original[record["record_id"]]
        for field, inherited in (
            ("raw_value", "raw_value"),
            ("raw_unit", "raw_unit"),
            ("raw_ordered_sites", "ordered_sites"),
            ("raw_energy_equation", "scalar_energy_equation"),
            ("citation", "citation"),
            ("location", "location"),
        ):
            typed_equal(
                record[field], parent[inherited], "C1_original_source_cells:" + field
            )
    adjudication = docs[str(PARENT / "row-adjudication.json")]
    summary = adjudication["summary"]
    require(
        len(adjudication["instances"]) == 378
        and summary["instances_unresolved"] == 378
        and summary["real_duplicate_equivalences_replayed"] == 16
        and summary["documented_unequal_source_assignment_conflicts"] == 2
        and summary["exact_target_instances_resolved"] == 0,
        "C1_preserved_scientific_result",
    )
    targets = docs["docs/evidence/phase_4j70/completeness-matrix.json"]["instances"]
    rows = adjudication["instances"]
    require(
        {r["id"] for r in rows} == {t["id"] for t in targets},
        "C1_exact_target_inventory",
    )
    target_by_id = {t["id"]: t for t in targets}
    comparisons = {r["record_id"]: r for r in normalized}
    require(
        Counter(r["final_J72_status"] for r in rows)
        == {
            "analog_unqualified": 294,
            "convention_unresolved": 60,
            "source_missing": 24,
        },
        "C1_row_classifications",
    )
    for row in rows:
        require(
            row["stable_sites"] == target_by_id[row["id"]]["stable_sites"]
            and row["J70_identity"] == target_by_id[row["id"]]["identity"]
            and row["selected_source"] is None
            and not row["resolved"],
            "C1_exact_row_dependency",
        )
        for candidate in row["candidate_records"]:
            require(
                candidate["J71_identity"]
                == original[candidate["record_id"]]["identity"]
                and candidate["comparison_identity"]
                == comparisons[candidate["record_id"]]["identity"],
                "C1_source_authority_dependency",
            )
    for name in (
        "canonical-capability.json",
        "strict-converter-report.json",
        "static-closure-report.json",
    ):
        document = docs[str(PARENT / name)]
        require(
            not document.get("executed", False)
            and not document.get("complete", False)
            and not document.get("automatic_assignment", False),
            "C1_execution_gates",
        )
    bundles = docs[str(PARENT / "candidate-source-bundles.json")]
    require(
        bundles["execution_profile"] is None and bundles["resolved_instances"] == 0,
        "C1_bundle_excludes_execution",
    )
    receipt = read(repo / PARENT / "receipt.json")
    decision = read(repo / PARENT / "decision.json")
    for ledger in (receipt, decision):
        require(
            ledger["new_QM"] == 4
            and ledger["HMDS_new_QM"] == 0
            and ledger["MD_steps"] is None
            and ledger["MD_steps_verified_minimum"] == 46
            and ledger["MD_steps_possible_maximum"] == 48
            and ledger["zero_execution_scope_met"] is False
            and ledger["reserved_validation_activated"] is False,
            "C1_historical_counts_cannot_be_laundered",
        )
    # Complete dependency closure is byte-bound to the pre-C1 base, rather
    # than relying only on a scan for recognizable quarantine names.
    return signed(
        {
            "schema": "island.j72-c1-contamination-exclusion.v1",
            "scientific_inputs": {n: pin(repo / n) for n in SCIENTIFIC_FILES},
            "checks": {
                role: "quarantine_excluded"
                for role in (
                    "literature_records",
                    "conflict_matrices",
                    "row_adjudication",
                    "source_authority",
                    "parameter_bundles",
                    "static_closure",
                    "capability_registration",
                    "HMDS_observations_and_producers",
                    "scientific_counts",
                    "QM_MD_evidence_summaries",
                )
            },
            "audit_metadata_policy": "quarantine references permitted only in integrity/scope ledgers; never scientific dependencies",
            "J71_original_record_identities_checked": 467,
            "J70_unresolved_instances_preserved": 378,
            "equal_duplicate_equivalences_preserved": 16,
            "unequal_COMPASS_assignments_unresolved": 2,
            "scientific_evidence_contaminated": False,
            "source_adjudication_excluded_quarantine": True,
            "no_new_literature_adjudication": True,
        }
    )


def reconstruct(repo, root):
    require(
        pin(root / "frozen-inputs.json") == CONTRACT_SHA, "C1_frozen_input_mutation"
    )
    frozen = read(root / "frozen-inputs.json")
    verify_signed(frozen)
    require(
        frozen["base_commit"] == BASE and frozen["required_J72"] == REQUIRED,
        "C1_ancestry",
    )
    for path, digest in frozen["inputs"].items():
        require(pin(repo / path) == digest, "C1_parent_input_mutation", path=path)
    q = quarantine(repo)
    report = exclusion(repo, q)
    return {
        "quarantine-manifest.json": q,
        "contamination-exclusion.json": report,
        "scope-correction.json": signed(
            {
                "schema": "island.j72-c1-scope-correction.v1",
                **DECLARATIONS,
                "historical_scope": "J72 violation retained; C1 adds no actual scientific execution",
                "quarantine_manifest_identity": q["identity"],
                "contamination_exclusion_identity": report["identity"],
                "scientific_result": {
                    "unresolved_instances": 378,
                    "duplicate_equivalences": 16,
                    "unresolved_unequal_COMPASS_assignments": 2,
                },
            }
        ),
    }


def verify_events(directories):
    """Check every allowed native launch against its frozen hash allowlist."""
    totals = Counter()
    for directory in directories:
        policy = read(directory / "policy.json")
        require(
            policy["schema"] == "island.j72-c1-execution-policy.v1"
            and policy["authorized_new_QM"] == policy["authorized_MD_steps"] == 0
            and policy["native_default"] == policy["shell_default"] == "deny",
            "C1_execution_policy_semantics",
        )
        digest = pin(directory / "policy.json")
        ledger = directory / "events.jsonl"
        # Lossless compression retains every original ledger byte without
        # making execution-time absolute paths a replay dependency.
        if ledger.is_file():
            raw_events = ledger.read_text()
        else:
            raw_events = lzma.decompress(
                (directory / "events.jsonl.xz").read_bytes()
            ).decode()
        events = [json.loads(line) for line in raw_events.splitlines()]
        require(
            events
            and (
                any(e["action"] == "installed" for e in events)
                or all(e["action"] in ("rejected", "exec_missing") for e in events)
            ),
            "C1_fence_not_installed",
        )
        for event in events:
            require(
                event["schema"] == "island.j72-c1-scope-event.v1"
                and event["timestamp_utc"],
                "C1_event_schema",
            )
            totals[event["action"]] += 1
            if event["action"] == "exec_allowed":
                identity = event["executable_identity"]
                entry = policy["executables"].get(identity["resolved"])
                require(
                    entry
                    and identity["sha256"] == entry["sha256"]
                    and entry["kind"] in ("python", "git", "ruff"),
                    "C1_unapproved_native_launch",
                )
                executable_name = Path(identity["resolved"]).name
                require(
                    (
                        entry["kind"] == "python"
                        and re.fullmatch(r"python(?:\d+(?:\.\d+)?)?", executable_name)
                    )
                    or (
                        entry["kind"] in ("git", "ruff")
                        and executable_name == entry["kind"]
                    ),
                    "C1_mislabeled_scientific_executable",
                )
                if entry["kind"] == "python":
                    require(
                        event["environment_policy"]["ISLAND_C1_POLICY_SHA256"]
                        == digest,
                        "C1_unfenced_python_launch",
                    )
            if event["action"] == "rejected":
                require(
                    {
                        "command",
                        "executable_identity",
                        "parent_process",
                        "environment_policy",
                        "reason",
                    }
                    <= set(event),
                    "C1_rejection_record",
                )
    return {
        "event_counts": dict(totals),
        "C1_new_native_QM_attempts": 0,
        "C1_MD_steps": 0,
    }
