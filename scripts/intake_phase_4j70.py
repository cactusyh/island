"""Frozen-source missing-row intake; native converter only, no QM/MD/fitting."""

import argparse
import lzma
import shutil
import subprocess
import tempfile
import time
import zipfile
from collections import Counter
from hashlib import sha256
from pathlib import Path

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_hamiltonian import load_model
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.source_intake import (
    CATEGORIES,
    ROW_KEYS,
    completeness,
    missing_rows,
)
from scripts.audit_phase_4j69 import inputs as parent_inputs

ROOT = Path("docs/evidence/phase_4j70")
BASE = "efacd4b798a1f345ea05f928905be88ad49acf53"
REQUIRED = "9f379f85f58807d9ae6502da1db712d8fc5680c9"
FROZEN = "162cc2d03c8619c9b753152d198aba342dc6c9a15e8162d32ea52219d6a4cd00"


def inputs(repo, root):
    require(pin(root / "plan.json") == FROZEN, "J70_frozen_plan")
    p = read(root / "plan.json")
    require(
        pin(root / "source-inventory.json") == p["source_inventory_sha256"],
        "J70_frozen_sources",
    )
    require(
        pin(repo / "docs/evidence/phase_4j69/receipt.json")
        == p["parent_receipt_sha256"],
        "J70_parent_receipt",
    )
    return p


def parse_native(raw, model):
    from island.forcefields.pcff.class2 import inspect_pcff_class2, resolve
    from island.forcefields.pcff.source import PCFFSource, records

    source = PCFFSource(raw, sha256(raw).hexdigest())
    catalog = inspect_pcff_class2(source)["records"]
    equivalence = records(source.inventory, "equivalence")
    rows = missing_rows(model)
    decisions = {
        r["id"]: resolve(
            r["family"],
            r["ordered_types"],
            catalog,
            equivalence,
        )
        for r in rows
    }
    return signed(
        {
            "status_counts": dict(Counter(d["status"] for d in decisions.values())),
            "queries": decisions,
            "wildcards_used": False,
            "zero_fill": False,
        }
    )


def discover_parse(repo, root):
    inputs(repo, root)
    inv = read(root / "source-inventory.json")
    # Hash every original before any new parser sees its contents.
    require(
        all(pin(r["path"]) == r["sha256"] for r in inv["sources"]),
        "J70_discovered_source_mutation",
    )
    model = load_model(repo)
    review = []
    archives = []
    parsed = []
    old = read(repo / "docs/evidence/phase_4j44/source-audit.json")
    literature = read(repo / "docs/evidence/phase_4j54/audit/literature-context.json")
    for r in inv["sources"]:
        path = Path(r["path"])
        item = {
            **r,
            "intake_status": "not_a_missing_PCFF_row_source",
            "eligible_native_overlay": False,
        }
        if path.suffix == ".frc":
            raw = path.read_bytes()
            # Other force-field families are inventoried, never parsed as PCFF.
            if "pcff" in path.name.lower():
                dest = Path("source-bytes") / (r["sha256"] + ".frc.xz")
                (root / dest).parent.mkdir(exist_ok=True)
                (root / dest).write_bytes(lzma.compress(raw))
                item.update(archive=str(dest), archive_sha256=pin(root / dest))
                try:
                    result = parse_native(raw, model)
                    parsed.append(
                        {
                            "path": r["path"],
                            "sha256": r["sha256"],
                            "category": r["category"],
                            "result": result,
                        }
                    )
                    item["intake_status"] = "parsed_separately_no_cross_source_overlay"
                except (ValueError, RuntimeError) as e:
                    item["intake_status"] = "parser_rejected"
                    item["parser_error"] = str(e)
            else:
                item["intake_status"] = "excluded_other_force_field_family"
        if path.suffix == ".zip":
            with zipfile.ZipFile(path) as z:
                names = [n for n in z.namelist() if not n.endswith("/")]
                require(
                    len(names) == len(set(names))
                    and all(
                        not Path(n).is_absolute() and ".." not in Path(n).parts
                        for n in names
                    ),
                    "unsafe_archive",
                )
                members = [
                    {
                        "name": n,
                        "sha256": sha256(z.read(n)).hexdigest(),
                        "bytes": len(z.read(n)),
                    }
                    for n in names
                ]
                archives.append(
                    {
                        "path": r["path"],
                        "sha256": r["sha256"],
                        "members": members,
                        "intake": "COMPASS/OPLS/other models remain separate; no exact missing PCFF rows supplied"
                        if "PDMS_benchmark" in path.name
                        else "separately inventoried external archive",
                    }
                )
                if path.name == "LAMMPS-pcff-main.zip":
                    raw = z.read("LAMMPS-pcff-main/pcff.frc")
                    require(
                        sha256(raw).hexdigest()
                        == "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c",
                        "engineering_FRC_replacement",
                    )
                    item["intake_status"] = "byte_identical_native_FRC_adds_no_rows"
        if path.is_relative_to(repo / "island_ffref"):
            relative = str(path.relative_to(repo / "island_ffref"))
            item["retained_review"] = [
                d for d in old["dependency_findings"] if d["source"] == relative
            ]
            item["retained_citation"] = [
                d for d in literature["sources"] if d["file"] == relative
            ]
            item["retained_bibliography_and_license"] = [
                d
                for d in old["files"]
                if d["filename"] == relative and d["sha256"] == r["sha256"]
            ]
            item["review_scope"] = (
                "Exact retained J44/J54 PDF row/location reviews reused with original producer identity; not a new PDF transcription"
            )
        review.append(item)
    write(
        root / "source-review.json",
        signed(
            {
                "schema": "island.j70-source-review.v1",
                "originals_checked": len(review),
                "sources": review,
                "archive_members": archives,
                "parsed_PCFF_files": parsed,
                "historical_literature_review": {
                    "file": "docs/evidence/phase_4j54/audit/literature-context.json",
                    "sha256": pin(
                        repo / "docs/evidence/phase_4j54/audit/literature-context.json"
                    ),
                },
                "author_artifact_status": "author_artifact_not_supplied",
                "accepted_rows": 0,
                "network_downloads": 0,
            }
        ),
    )


def reports(repo, root):
    inputs(repo, root)
    model = load_model(repo)
    rows = missing_rows(model)
    require(not completeness(model, rows), "unexpected_complete_source")
    matrix = signed(
        {
            "schema": "island.j70-completeness-matrix.v1",
            "graph_identity": model["graph_identity"],
            "atom_order": model["atom_order"],
            "instances": rows,
            "family_counts": dict(Counter(r["family"] for r in rows)),
            "enumerated_instances": len(rows),
            "symmetry_reduced_physical_couplings": 56,
            "native_rows_ingested": 0,
            "author_primary_rows_ingested": 0,
            "literature_external_rows_ingested": 0,
            "research_values_ingested": 0,
            "source_complete": False,
        }
    )
    targets = []
    for r in rows:
        if r["status"] == "qm_fit_required":
            targets.append(
                signed(
                    {
                        "schema": "island.j70-frozen-fit-target.v1",
                        "id": r["id"],
                        "graph_identity": r["graph_identity"],
                        "stable_sites": r["stable_sites"],
                        "types": r["ordered_types"],
                        "equation": r["contract"],
                        "dependencies": r["dependencies"],
                        "coefficient_initial_values": None,
                        "input_template": {
                            "coordinates": "exact retained HMDS geometry required; no new geometry generated",
                            "energy_scope": "complete conservative energy/forces; coupling isolation not assumed",
                            "producer": None,
                            "QM_protocol": None,
                            "license_and_provenance": None,
                            "selection_rule": None,
                            "identifiability_and_validation_contract": None,
                        },
                        "status": "not_executable_missing_frozen_fit_protocol",
                        "native_PCFF_authority_on_fit": False,
                        "QM_authorized": False,
                        "fitting_authorized": False,
                    }
                )
            )
    decision = signed(
        {
            "schema": "island.j70-decision.v1",
            "decision": "parameter_incomplete",
            "source_complete": False,
            "converter_complete": False,
            "parameter_assignment_capability": False,
            "execution_profile_created": False,
            "registered_profiles": [],
            "accepted_coefficients": 0,
            "automatic_assignment": False,
            "new_QM": 0,
            "MD_steps": 0,
            "production_fits": 0,
            "reserved_validation_activated": False,
            "reserved_validation_consumed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
            "scientific_validation": "not_established",
            "chemical_transfer": "not_established",
            "native_PCFF_equivalence": "not_established",
            "original_author_artifact_status": "not_yet_supplied",
            "author_confirmations": "28 user-attested",
            "J68_bound_free_incompatibility": "preserved; native values for already free columns alone cannot remove incompatibility",
            "next_dependency": "Exact hash-bound authoritative numerical rows and native Phi/AAT/pair conventions; no generic new QM or automatic fitting allocation",
        }
    )
    return {
        "parameter-bundle-schema.json": signed(
            {
                "schema": "island.j70-bundle-schema.v1",
                "row_schema": "island.pcff-source-row.v1",
                "required_row_fields": sorted(ROW_KEYS),
                "source_categories": list(CATEGORIES),
                "instance_count": 378,
                "missing_coefficients": None,
                "missing_is_not_zero": True,
                "immutability": "frozen serialized SourceRow; payload copies; identity plus reparse plus graph/dependency semantic checks",
                "completeness_gates": read(root / "plan.json")["success_gates"],
                "registration_authority": False,
            }
        ),
        "completeness-matrix.json": matrix,
        "missing-row-and-convention-report.json": signed(
            {
                "schema": "island.j70-blockers.v1",
                "missing_instances": 378,
                "statuses": dict(Counter(r["status"] for r in rows)),
                "convention_unresolved_instances": [
                    r["id"]
                    for r in rows
                    if r["contract"]["equation_status"] == "convention_unresolved"
                ],
                "source_unavailable": [r["id"] for r in rows],
                "author_primary": "author_artifact_not_supplied",
                "native_versus_external": "No external/research values substituted; complete-source authority unavailable",
                "full_inventory_closure": "blocked: missing rows, native Phi/AAT and pair-convention authority; empty torsion-torsion section is not numeric zero",
            }
        ),
        "qm-fit-targets.json": signed({"targets": targets, "executed": False}),
        "parameter-bundle.json": signed(
            {
                "schema": "island.pcff-missing-row-bundle.v1",
                "matrix_identity": matrix["identity"],
                "status": "parameter_incomplete",
                "rows": [r["identity"] for r in rows],
                "source_complete": False,
                "execution_authorized": False,
                "source_categories_separate": True,
                "coefficients_supplied": 0,
                "profile": None,
            }
        ),
        "decision.json": decision,
    }


def converters(repo, root):
    inputs(repo, root)
    out = root / "converter"
    require(not out.exists(), "J70_no_overwrite_converter")
    out.mkdir()
    parent = repo / "docs/evidence/phase_4j69"
    _, raw = parent_inputs(repo, parent)
    results = []
    with tempfile.TemporaryDirectory(prefix="J70-strict-") as temp:
        base = Path(temp)
        exe = base / "msi2lmp"
        shutil.copyfile(parent / "converter/msi2lmp", exe)
        exe.chmod(0o700)
        frc = base / "pcff.frc"
        frc.write_bytes(raw["tools/msi2lmp/frc_files/pcff.frc"])
        for i in (1, 2):
            attempt = out / f"attempt-{i}"
            attempt.mkdir()
            work = base / f"relocated-{i}"
            shutil.copytree(parent / "converter/input", work)
            command = [
                str(exe),
                "system",
                "-class",
                "II",
                "-frc",
                str(frc),
                "-print",
                "3",
                "-nocenter",
            ]
            start = time.monotonic()
            try:
                r = subprocess.run(
                    command, cwd=work, capture_output=True, timeout=60, check=False
                )
                code = r.returncode
                stdout = r.stdout
                stderr = r.stderr
                failure = None
            except subprocess.TimeoutExpired as e:
                code = None
                stdout = e.stdout or b""
                stderr = e.stderr or b""
                failure = "interrupted_timeout"
            for name, blob in (("stdout", stdout), ("stderr", stderr)):
                (attempt / (name + ".log.xz")).write_bytes(lzma.compress(blob))
            result = signed(
                {
                    "schema": "island.j70-converter-attempt.v1",
                    "attempt": i,
                    "argv": command,
                    "returncode": code,
                    "failure": failure,
                    "runtime_seconds": time.monotonic() - start,
                    "executable_sha256": pin(exe),
                    "source_sha256": pin(frc),
                    "input_files": {
                        n: pin(work / n)
                        for n in (
                            "system.car",
                            "system.mdf",
                            "system.json",
                            "manifest.json",
                        )
                    },
                    "stdout_sha256": sha256(stdout).hexdigest(),
                    "stderr_sha256": sha256(stderr).hexdigest(),
                    "complete": False,
                    "ignore_used": False,
                    "new_QM": 0,
                    "MD_steps": 0,
                }
            )
            write(attempt / "attempt.json", result)
            results.append(result)
    write(
        out / "report.json",
        signed(
            {
                "schema": "island.j70-converter-report.v1",
                "attempts": results,
                "complete": False,
                "successive_fresh_directory_returncodes": [
                    r["returncode"] for r in results
                ],
                "deterministic_failure": all(r["returncode"] == 17 for r in results)
                and len({r["stdout_sha256"] for r in results}) == 1,
                "deterministic_complete_output": "not_established",
                "missing": "endbonddihedral hc c3 sio c3",
                "full_emission": "blocked",
                "new_QM": 0,
                "MD_steps": 0,
            }
        ),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--discover-parse", action="store_true")
    p.add_argument("--converter", action="store_true")
    a = p.parse_args()
    repo = a.repo.resolve()
    root = a.root.resolve()
    if a.discover_parse:
        discover_parse(repo, root)
    if a.converter:
        converters(repo, root)
    for n, v in reports(repo, root).items():
        write(root / n, v)
    print("J70 parameter_incomplete;378 exact sites;zero QM/MD/fitting")


if __name__ == "__main__":
    main()
