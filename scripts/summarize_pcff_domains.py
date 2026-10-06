"""Publish J3 evidence and all-label gaps; full-source acceptance fails closed."""

import argparse
from pathlib import Path

from pcff_domain_gates import numerical_gate

from island.forcefields.pcff import load_pcff_source
from island.forcefields.pcff.domain_coverage import pcff_domain_coverage
from island.workflows.storage import checksum, json_bytes, publish, read_json


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--matrix", type=Path, required=True)
    p.add_argument("--numerical", type=Path, required=True)
    p.add_argument("--workflows", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--offline", type=Path, required=True)
    p.add_argument("--charge-readback", type=Path, required=True)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    source = load_pcff_source(a.source)
    ledger = pcff_domain_coverage(source)
    matrix = read_json(a.matrix / "report.json")
    numerical = read_json(a.numerical / "report.json")
    assert matrix["source"] == source.identity
    cases = matrix["cases"]
    assert len({c["name"] for c in cases}) == len(cases) == 37
    from island.charge_references.records import unpack

    gaps = []
    assignments = {}
    for c in cases:
        path = a.matrix / c["name"] / "assignment.json"
        assignments[c["name"]] = (
            unpack(path.read_text())["assignments"] if path.is_file() else []
        )
        by_id = {r["id"]: r for r in assignments[c["name"]]}
        for stage, key in [
            ("typing", "typing_diagnostics"),
            ("charges", "charge_diagnostics"),
            ("model", "diagnostics"),
        ]:
            for diagnostic in c.get(key, []):
                detail = by_id.get(diagnostic.get("id"), {})
                reason = diagnostic.get("reason", "")
                classification = (
                    "source_parameter_missing"
                    if "no source row" in reason or reason == "source_parameter_missing"
                    else "interpretation_unresolved"
                )
                gaps.append(
                    {
                        "case": c["name"],
                        "stage": stage,
                        "classification": classification,
                        **detail,
                        **diagnostic,
                    }
                )
        if c.get("error") or c.get("expected_label_discrepancy"):
            gaps.append(
                {
                    "case": c["name"],
                    "stage": "harness_or_expected_labels",
                    "error": c.get("error"),
                    "discrepancy": c.get("expected_label_discrepancy"),
                }
            )
    for row in ledger["type_rows"]:
        row["fixture_observations"] = [
            {
                "case": c["name"],
                "typing": c["typed"],
                "native_charges": c["charges_complete"],
                "model": c["model_complete"],
                "expected_labels_match": c.get("expected_labels_match"),
            }
            for c in cases
            if row["type"] in c.get("types", [])
        ]
    for row in ledger["type_rows"]:
        for family in row["families"]:
            family["observations"] = [
                {
                    "case": case,
                    "interaction": r["id"],
                    "status": r["status"],
                    "sites": r["sites"],
                }
                for case, rows in assignments.items()
                for r in rows
                if row["type"] in r.get("supplied_types", [])
                and family["family"] == r.get("family")
            ]
    required = ("hydrogen_sulfide", "heavy_water")
    checks = []
    for name in required:
        selected = [c for c in numerical["cases"] if c["case"] == name]
        outcome = read_json(a.workflows / name / "outcome.json")
        checks.append({"case": name, "numerical": selected, "workflow": outcome})
    domains = {
        c["domain"]
        for c in cases
        if c["domain"]
        in {
            "molecular_hydrogen",
            "sulfide_hydride",
            "isotope_water",
            "amide",
            "substituted_aromatic_n",
            "imine",
            "aromatic_amine",
            "phosphorus_oxo",
            "peroxide",
            "acid",
        }
        and c["typed"]
        and c["charges_complete"]
        and c.get("expected_labels_match") is True
        and not c.get("error")
    }
    gate = (
        len(domains) >= 4
        and numerical_gate(numerical["cases"], numerical.get("source_unchanged"))
        and all(
            c["workflow"].get("passed") is True and not c["workflow"].get("error")
            for c in checks
        )
    )
    for check in checks:
        name = check["case"]
        n = check["numerical"][0] if len(check["numerical"]) == 1 else {}
        declaration = read_json(a.workflows / name / "declaration.json")
        from validate_prepared_workflow import hashes

        gate = gate and (
            hashes(a.numerical / name / "bundle")
            == n.get("bundle_hashes")
            == declaration["input_hashes"]
        )
        gate = (
            gate
            and n.get("prepared_identity") == declaration["prepared_identity"]
            and n.get("model_identity") == declaration["native_identity"]
        )
        gate = gate and declaration["source_sha256"] == source.identity["sha256"]
    offline = read_json(a.offline)
    charge_readback = read_json(a.charge_readback)
    gate = (
        gate
        and offline.get("passed") is True
        and charge_readback["source_sha256"] == source.identity["sha256"]
    )
    for check in checks:
        name = check["case"]
        native = check["numerical"][0]
        q = [r for r in charge_readback["cases"] if r["case"] == name]
        b = [
            r
            for r in offline["bundles"]
            if r["identity"] == native["prepared_identity"]
            and r["native"] == native["model_identity"]
        ]
        w = [
            r
            for r in offline["workflows"]
            if r["identity"] == native["prepared_identity"]
            and r["status"] == "completed"
            and r["frames"] == 5
        ]
        gate = (
            gate
            and len(q) == len(b) == len(w) == 1
            and q[0]["model_identity"] == native["model_identity"]
            and q[0]["max_error"] <= 1e-12
        )
    counts = lambda rows: {
        "denominator": len(rows),
        "typed": sum(c["typed"] for c in rows),
        "native_charges": sum(c["charges_complete"] for c in rows),
        "complete_models": sum(c["model_complete"] for c in rows),
    }
    report = {
        "source": source.identity,
        "source_sha256": checksum(a.source.read_bytes()),
        "matrix_counts": counts(cases),
        "original18": counts([c for c in cases if c["domain"] == "original_j2"]),
        "new_domains_with_labels_and_native_charges": sorted(domains),
        "numerical_workflow": checks,
        "bounded_operational_gate": gate,
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "offline_receipt": {
            "sha256": checksum(a.offline.read_bytes()),
            "passed": offline.get("passed"),
        },
        "charge_readback_sha256": checksum(a.charge_readback.read_bytes()),
        "matrix_sha256": checksum((a.matrix / "report.json").read_bytes()),
        "numerical_sha256": checksum((a.numerical / "report.json").read_bytes()),
    }
    publish(a.output / "coverage.json", json_bytes(ledger))
    publish(
        a.output / "gaps.json",
        json_bytes(
            {
                "interactions": gaps,
                "unresolved_labels": [
                    r for r in ledger["type_rows"] if r["declared_predicate"] is None
                ],
            }
        ),
    )
    publish(a.output / "report.json", json_bytes(report))
    return 0 if gate and not a.require_full_source else 1


if __name__ == "__main__":
    raise SystemExit(main())
