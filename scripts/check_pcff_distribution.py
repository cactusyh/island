"""Recheck J16 artifacts and native semantics; never accept success booleans.

This performs offline consistency calculations, not force evaluation or external
commands. The current gate requires registered independently checked vertical
artifacts for any claimed executable polymer; converter output is insufficient.
"""

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from pcff_distribution_contract import frozen
from pcff_distribution_reference import labels_for, missing_request
from pcff_j5_reference import charge_raw, read_source
from validate_pcff_distribution import artifacts, assess_set, native_case
from validate_pcff_singlepoint import check_inventory, data_sections
from validate_pcff_urethanes import contract

from island.forcefields.pcff import (
    inspect_pcff_distribution,
    load_pcff_distribution_audit,
    load_pcff_source,
    pcff_domain_coverage,
)
from island.workflows import storage
from island.workflows.bundle import system_from

EVIDENCE = Path(__file__).parents[1] / "docs/evidence"


def checked_set(results, declaration):
    names = [c["name"] for c in declaration["cases"]]
    if (
        not names
        or len(names) != len(set(names))
        or Counter(r.get("case") for r in results) != Counter(names)
    ):
        raise ValueError("Missing, duplicate, empty or unexpected mandatory cases")


def artifact_check(root, expected):
    if not expected or any(Path(n).name != n for n in expected):
        raise ValueError("Invalid case artifact references")
    if artifacts(root) != expected:
        raise ValueError("Missing, stale, altered or unexpected case artifacts")


def check_coordinates(atoms, original):
    """Exact WriteDataFile.c:391 decimal serialization, not a force tolerance."""
    actual = np.array(
        [list(map(float, row[4:7])) for row in sorted(atoms, key=lambda x: int(x[0]))]
    )
    expected = np.array([[float(f"{v:.9f}") for v in xyz] for xyz in original])
    np.testing.assert_array_equal(actual, expected)
    return {
        "format": "%15.9f",
        "reference": "WriteDataFile.c:391",
        "exact_serialized_frame": True,
        "maximum_difference_from_original_angstrom": float(
            np.max(np.abs(actual - np.asarray(original)))
        ),
    }


def derive_gates(
    cases,
    verified_vertical_names,
    *,
    audit_validated,
    total_source_labels,
    verified_source_labels,
):
    """Inputs here are computed by the checker, never read as receipt booleans."""
    audit = audit_validated and bool(cases) and all(not r["errors"] for r in cases)
    polymers = [r for r in cases if r["polymer_target"]]
    operational = [
        r["case"]
        for r in polymers
        if r["native_model_complete"] and r["case"] in verified_vertical_names
    ]
    return {
        "audit_gate": audit,
        "converter_boundary_gate": audit,
        "executable_polymer_gate": audit and bool(operational),
        "verified_polymer_cases": operational,
        "full_source_complete": audit
        and total_source_labels > 0
        and verified_source_labels == total_source_labels
        and all(
            r["native_model_complete"] and r["case"] in verified_vertical_names
            for r in cases
        ),
    }


def check(root, *, frc, rlb, templates, converter):
    d = frozen("declaration")
    report = storage.decode(storage.read_json(root / "report.json"))
    expected_build = frozen("converter")
    if (
        storage.read_json(root / "converter-build.json") != expected_build
        or storage.checksum(converter.read_bytes()) != expected_build["sha256"]
    ):
        raise ValueError("Converter evidence identity mismatch")
    if report.get("declaration_sha256") != storage.checksum(
        (EVIDENCE / "phase_4j16_declaration.json").read_bytes()
    ):
        raise ValueError("Stale or mismatched declaration evidence")
    hashes = {
        k: d["distribution"][n]["sha256"]
        for k, n in [
            ("frc", "pcff.frc"),
            ("rlb", "pcff.rlb"),
            ("templates", "pcff_templates.dat"),
        ]
    }
    audit = inspect_pcff_distribution(
        frc_path=frc, rlb_path=rlb, templates_path=templates, expected_hashes=hashes
    )
    audit = load_pcff_distribution_audit(
        root / "distribution.json", frc_path=frc, rlb_path=rlb, templates_path=templates
    )
    if audit.identity != report.get("distribution_identity"):
        raise ValueError("Stale distribution evidence")
    source = load_pcff_source(frc)
    raw = read_source(source.raw)
    if report.get("source") != source.identity:
        raise ValueError("Wrong source evidence")
    checked_set(report["cases"], d)
    recorded = {r["case"]: r for r in report["cases"]}
    expected = contract()["cases"]
    links = frozen("reference_links")
    checks = []
    verified_vertical = set()
    actual_results = []
    for c in d["cases"]:
        r = recorded[c["name"]]
        if r.get("execution_error"):
            raise ValueError("Recorded execution failure: " + c["name"])
        work = root / c["name"]
        artifact_check(work, r["artifacts"])
        original = Path(c["system_path"])
        if storage.checksum(original.read_bytes()) != c["system_sha256"]:
            raise ValueError("Frozen graph mismatch: " + c["name"])
        s = system_from(storage.decode(storage.read_json(original)))
        current, typing, inspections = native_case(s, c, source, frc)
        if any(r.get(k) != v for k, v in current.items()):
            raise ValueError("Contradictory native stage or diagnostic: " + c["name"])
        if storage.decode(storage.read_json(work / "native.json")) != inspections:
            raise ValueError("Contradictory native inspection: " + c["name"])
        detail = {
            "case": c["name"],
            "native_model_complete": bool(current.get("model_complete")),
            "polymer_target": "psmiles" in c,
            "errors": [],
        }
        if current.get("charges_complete"):
            labels = labels_for(s, c)
            charges = charge_raw(raw, s, labels)
            if (
                labels != typing.assignments
                or not charges["complete"]
                or storage.decode(storage.read_json(work / "independent-charge.json"))
                != charges
            ):
                raise ValueError("Independent label/charge contradiction")
            current["independent_types_and_charges"] = True
            command = storage.read_json(work / "command.json")
            if len(command) != 9 or command[1:] != [
                "reference",
                "-class",
                "II",
                "-frc",
                str(frc.resolve()),
                "-p",
                "3",
                "-nocenter",
            ]:
                raise ValueError("Wrong converter invocation")
            if (
                storage.checksum(Path(command[0]).read_bytes())
                != expected_build["sha256"]
            ):
                raise ValueError("Wrong recorded converter executable")
            text = (work / "stdout.txt").read_text() + (work / "stderr.txt").read_text()
            code = r["converter"].get("returncode")
            if type(code) is not int:
                raise ValueError("Missing converter exit status")
            if code:
                detail["converter"] = missing_request(
                    raw, code, text, system=s, reference_labels=labels
                )
                current["converter"] = {
                    "returncode": code,
                    "independent_missing": detail["converter"],
                }
            else:
                if "Unable to find" in text:
                    raise ValueError(
                        "Success contradicts converter missing-parameter output"
                    )
                sections = data_sections(work / "reference.data")
                ids = sorted(labels)
                adj = {i: set() for i in ids}
                for b in s.topology.bonds.values():
                    adj[b.site1].add(b.site2)
                    adj[b.site2].add(b.site1)
                vector = {i: float(charges["partial_charges"][i]) for i in ids}
                detail["inventory"] = check_inventory(sections, ids, adj, vector)
                detail["coordinate_serialization"] = check_coordinates(
                    sections["Atoms"], [s.coordinates.get(i) for i in ids]
                )
                current["converter"] = {
                    "returncode": 0,
                    "charge_vector_preserved": True,
                    "output_sections": {k: len(v) for k, v in sections.items()},
                }
                if c["name"] not in links:
                    detail["errors"].append(
                        "Independent parameter/numerical reference unavailable"
                    )
                else:
                    link = links[c["name"]]
                    old = Path(link["root"])
                    for name, h in link["files"].items():
                        if storage.checksum((old / name).read_bytes()) != h:
                            raise ValueError("Historical independent evidence mismatch")
                    # Byte-identical coefficients/inventories/coordinates to the
                    # independently authored, previously numerically verified raw
                    # reference. Do not trust a success flag in a new report.
                    if (work / "reference.data").read_bytes() != (
                        old / "reference/reference.data"
                    ).read_bytes():
                        detail["errors"].append(
                            "Converter output differs from retained independently verified reference"
                        )
                    else:
                        detail["independent_reference_data_sha256"] = link["files"][
                            "reference/reference.data"
                        ]
                        detail["reference_status"] = (
                            "byte-identical retained independent reference; no new force calculation"
                        )
                        verified_vertical.add(c["name"])
            if current.get("model_complete") and c["name"] not in verified_vertical:
                detail["errors"].append("Executable model lacks independent evidence")
        else:
            current["converter"] = {
                "state": "blocked_by_native_typing_or_charge",
                "not_executed": True,
            }
            if r.get("converter") != current["converter"]:
                raise ValueError("Invalid prerequisite rejection")
        actual_results.append(current)
        checks.append(detail)
    failures = assess_set(actual_results, d, expected)
    if failures:
        raise ValueError(str(failures))
    ledger = pcff_domain_coverage(source, profile=d["typing_profile"])
    gates = derive_gates(
        checks,
        verified_vertical,
        audit_validated=True,
        total_source_labels=ledger["denominator"],
        verified_source_labels=ledger["global_verified_labels"],
    )
    return {
        "schema": "island_j16_checked_boundaries_v1",
        "cases": checks,
        "distribution_identity": audit.identity,
        "input_report_sha256": storage.checksum((root / "report.json").read_bytes()),
        "source": source.identity,
        "production_validated": False,
        "simulation_readiness": "not_established",
        **gates,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("input", "frc", "rlb", "templates", "converter", "output"):
        p.add_argument("--" + key, type=Path, required=True)
    p.add_argument(
        "--gate",
        choices=("audit", "converter", "polymer", "full-source"),
        default="polymer",
    )
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    try:
        receipt = check(
            a.input, frc=a.frc, rlb=a.rlb, templates=a.templates, converter=a.converter
        )
    except Exception as e:  # noqa: BLE001 -- retain failure, never a chemical success
        receipt = {
            "schema": "island_j16_checked_boundaries_v1",
            "error": type(e).__name__ + ": " + str(e),
        }
    storage.publish(a.output, storage.json_bytes(receipt))
    gate = "full-source" if a.require_full_source else a.gate
    key = {
        "audit": "audit_gate",
        "converter": "converter_boundary_gate",
        "polymer": "executable_polymer_gate",
        "full-source": "full_source_complete",
    }[gate]
    print({k: v for k, v in receipt.items() if k not in ("cases", "source")})
    return 0 if receipt.get(key) is True else 1


if __name__ == "__main__":
    raise SystemExit(main())
