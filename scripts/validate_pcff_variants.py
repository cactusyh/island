"""Declared J7 source comparison, independent raw-FRC checks and child reload."""

import argparse
import builtins
import hashlib
import shutil
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

from pcff_j5_reference import charge_raw, read_source, resolve_raw

from island.charge_references.records import pack, unpack
from island.forcefields.pcff import (
    audit_pcff_source_variants,
    compare_pcff_sources,
    load_pcff_resolution_assessment,
    load_pcff_source_comparison,
    load_pcff_source_variant_audit,
    save_pcff_source_comparison,
    save_pcff_source_variant_audit,
)
from island.forcefields.pcff.charges import identity
from island.forcefields.pcff.variants_cli import declared_sources
from island.workflows.bundle import system_from
from island.workflows.storage import child as contained
from island.workflows.storage import decode, publish, read_json


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    publish(path, pack(value).encode())


def scan_independent(raw):
    # Deliberately authored raw-line reader; no production parser/catalog usage.
    arity = {
        "atom_types": 1,
        "equivalence": 1,
        "auto_equivalence": 1,
        "bond_increments": 2,
        "quartic_bond": 2,
        "quadratic_bond": 2,
        "quartic_angle": 3,
        "quadratic_angle": 3,
        "torsion_3": 4,
        "torsion_1": 4,
        "wilson_out_of_plane": 4,
        "nonbond(9-6)": 1,
        "bond-bond": 3,
        "bond-angle": 3,
        "bond-bond_1_3": 4,
        "angle-angle": 4,
        "end_bond-torsion_3": 4,
        "middle_bond-torsion_3": 4,
        "angle-torsion_3": 4,
        "angle-angle-torsion_1": 4,
    }
    section = namespace = None
    records = {}
    sections = []
    for line, text in enumerate(raw.decode().splitlines(), 1):
        fields = text.split("!", 1)[0].split()
        if not fields:
            continue
        if fields[0].startswith("#"):
            section = fields[0][1:]
            namespace = " ".join(fields[1:])
            sections.append((section, namespace))
            continue
        if section in ("version", "reference", "end") or not fields[0][0].isdigit():
            continue
        if len(fields) < 2 or not fields[1].isdigit():
            raise ValueError(f"Independent row format {line}")
        k = arity.get(section, 0)
        tokens = fields[2 : 2 + k]
        for token in fields[2 + k :]:
            try:
                value = Decimal(token)
                if not value.is_finite():
                    raise ValueError("nonfinite independent source")
                tokens.append(str(value.normalize()) if value else "0")
            except ArithmeticError:
                tokens.append(token)
        v = Decimal(fields[0])
        version = str(v.normalize()) if v else "0"
        records[f"{section}:{namespace}:{line}"] = {
            "section": section,
            "namespace": namespace,
            "version": version,
            "reference": fields[1],
            "tokens": tokens,
        }
    return records, sections


def independent_compare(variant):
    expected, sections = scan_independent(variant.raw)
    p = variant.payload
    actual = {r["id"]: r["content"] for r in p["rows"]}
    if expected != actual:
        raise ValueError("Independent source-row inventory disagreement")
    if sections != [(s["name"], s["namespace"]) for s in p["sections"]]:
        raise ValueError("Independent section disagreement")
    return {
        "variant_identity": variant.identity,
        "source": variant.selection,
        "row_count": len(expected),
        "sections": len(sections),
        "independent_raw_inventory_verified": True,
    }


def check_audit(audit, raw_rows, system, labels, tolerance):
    p = unpack(audit.json_text)
    comparisons = []
    for variant in p["variant_results"]:
        rows = raw_rows[variant["source"]["sha256"]]
        count = 0
        for e in variant["entries"]:
            raw = resolve_raw(rows, e["family"], e["types"])
            q = e["query"]["forward"]
            if raw["status"] != q["status"]:
                raise ValueError(f"Independent status disagreement: {e['request_id']}")
            selected = sorted({r["record_id"] for r in q.get("selected", [])})
            if raw["status"] == "assigned":
                if selected != raw["selected_ids"]:
                    raise ValueError(
                        f"Independent selected-row disagreement: {e['request_id']}"
                    )
                if len(raw["values"]) != len(q["normalized_values"]) or any(
                    abs(a - b) > tolerance
                    for a, b in zip(raw["values"], q["normalized_values"], strict=True)
                ):
                    raise ValueError(
                        f"Independent coefficient disagreement: {e['request_id']}"
                    )
            elif (
                sorted({r["record_id"] for r in q.get("candidates", [])})
                != raw["candidate_ids"]
            ):
                raise ValueError(
                    f"Independent candidate disagreement: {e['request_id']}"
                )
            count += 1
        charge = charge_raw(rows, system, labels)
        production = variant["charge_audit"]
        if charge["complete"] != production["source_charge_calculation_complete"]:
            raise ValueError("Independent charge completeness disagreement")
        if any(
            Decimal(v) != Decimal(production["charges_e"][i])
            for i, v in charge["partial_charges"].items()
        ):
            raise ValueError("Independent exact charge disagreement")
        for c, pc in zip(charge["components"], production["components"], strict=True):
            if c["coverage_complete"] != pc["coverage_complete"] or Decimal(
                c["known_contribution_total"]
            ) != Decimal(pc["known_contribution_total_e"]):
                raise ValueError("Independent component disagreement")
        comparisons.append(
            {
                "source": variant["source"],
                "requests": count,
                "independent_selected_rows_and_coefficients": True,
                "independent_charge_decimal": charge,
                "model_complete": variant["native_model_complete"],
            }
        )
    return comparisons


def _load_comparison_identity(path, variants):
    # The loader fully validates this immutable JSON-owned result at the boundary.
    result = load_pcff_source_comparison(path, variants=variants)
    return identity(unpack(result.json_text))


def _load_audit_identity(path, assessment, variants):
    result = load_pcff_source_variant_audit(
        path, assessment=assessment, variants=variants
    )
    return identity(unpack(result.json_text))


def child(args):
    forbidden = {"openmm", "scipy", "rdkit", "parmed", "foyer"}
    if forbidden.intersection(sys.modules):
        raise RuntimeError("Scientific import before offline reconstruction")
    importer = builtins.__import__

    def blocked(name, *a, **kw):
        if name.split(".")[0] in forbidden:
            raise RuntimeError("Forbidden science import: " + name)
        return importer(name, *a, **kw)

    builtins.__import__ = blocked
    variants = declared_sources(args.declaration)
    manifest = unpack((args.child / "manifest.json").read_text())
    for name, h in manifest["hashes"].items():
        if sha(contained(args.child, name)) != h:
            raise ValueError("Relocated artifact changed")
    comparison_id = _load_comparison_identity(args.child / "comparison.json", variants)
    if comparison_id != manifest["comparison_identity"]:
        raise ValueError("Relocated comparison identity changed")
    pinned = next(v for v in variants if v.payload["native_assignment_authorized"])
    source = load_pcff_source_from_variant(pinned)
    results = []
    for row in manifest["cases"]:
        system = system_from(decode(read_json(contained(args.child, row["system"]))))
        assessment = load_pcff_resolution_assessment(
            contained(args.child, row["assessment"]), source, system=system
        )
        audit_id = _load_audit_identity(
            contained(args.child, row["audit"]), assessment, variants
        )
        if audit_id != row["identity"]:
            raise ValueError("Relocated cross-source audit identity changed")
        results.append({"case": row["case"], "identity": audit_id, "unchanged": True})
    write(
        args.child.parent / "child-reconstruction.json",
        {
            "comparison_identity": comparison_id,
            "cases": results,
            "scientific_imports_blocked": True,
            "separate_process": True,
            "sources_explicit": True,
        },
    )
    return 0


def load_pcff_source_from_variant(variant):
    from island.forcefields.pcff.source import PCFFSource

    source = PCFFSource(variant.raw, variant.expected_sha256)
    source.require_assignment()
    return source


def run(args):
    declaration = read_json(args.declaration)
    args.output.mkdir(parents=True, exist_ok=False)
    write(args.output / "declaration.json", declaration)
    variants = declared_sources(args.declaration)
    report = {
        "schema": "island_j7_multisource_acceptance_v1",
        "sources": [],
        "cases": [],
        "comparison_gate": False,
        "operational_gate": False,
        "full_source_gate": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    staging = args.output / "staging"
    staging.mkdir()
    c = compare_pcff_sources(variants)
    save_pcff_source_comparison(c, staging / "comparison.json")
    for v in variants:
        report["sources"].append(independent_compare(v))
    raw_rows = {v.expected_sha256: read_source(v.raw) for v in variants}
    # Predeclared conversion-only tolerance; physical reversal uses exact values.
    tolerance = read_json(args.tolerances)["coefficient_atol"]
    source = load_pcff_source_from_variant(
        next(v for v in variants if v.payload["native_assignment_authorized"])
    )
    manifest = {"comparison_identity": c.identity, "cases": [], "hashes": {}}
    for case in declaration["cases"]:
        name = case["name"]
        path = Path(case["retained_input"]["path"])
        try:
            if sha(path) != case["retained_input"]["sha256"]:
                raise ValueError("Historical system hash mismatch")
            system = system_from(decode(read_json(path)))
            original = args.assessments / (name + ".assessment.json")
            assessment = load_pcff_resolution_assessment(
                original, source, system=system
            )
            audit = audit_pcff_source_variants(assessment, variants)
            p = unpack(audit.json_text)
            checked = check_audit(
                audit,
                raw_rows,
                system,
                p["original_assessment"]["typing"]["assignments"],
                tolerance,
            )
            save_pcff_source_variant_audit(audit, staging / (name + ".audit.json"))
            shutil.copyfile(original, staging / (name + ".assessment.json"))
            shutil.copyfile(path, staging / (name + ".system.json"))
            row = {
                "case": name,
                "identity": identity(p),
                "independent_comparisons": checked,
                "source_variant_only": p["source_variant_only"],
                "per_variant": [
                    {
                        "source": v["source"],
                        "classification_counts": v["classification_counts"],
                        "required_interaction_blockers": v[
                            "required_interaction_blockers"
                        ],
                        "source_charge_calculation_complete": v["charge_audit"][
                            "source_charge_calculation_complete"
                        ],
                        "components": v["charge_audit"]["components"],
                        "native_model_complete": v["native_model_complete"],
                    }
                    for v in p["variant_results"]
                ],
            }
            manifest["cases"].append(
                {
                    "case": name,
                    "identity": row["identity"],
                    "audit": name + ".audit.json",
                    "assessment": name + ".assessment.json",
                    "system": name + ".system.json",
                }
            )
        except Exception as e:  # noqa: BLE001 -- durable truthful failed-case evidence
            row = {"case": name, "error": f"{type(e).__name__}: {e}"}
        report["cases"].append(row)
        print(name, row.get("per_variant", row.get("error")), flush=True)
    manifest["hashes"] = {
        str(p.relative_to(staging)): sha(p) for p in staging.iterdir() if p.is_file()
    }
    write(staging / "manifest.json", manifest)
    relocated = args.output / "relocated"
    staging.rename(relocated)
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--child",
        str(relocated.resolve()),
        "--declaration",
        str(args.declaration.resolve()),
    ]
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=False,
        timeout=declaration["budgets"]["external_child_seconds"],
    )
    (args.output / "child.stdout").write_text(result.stdout)
    (args.output / "child.stderr").write_text(result.stderr)
    report["child"] = {
        "command": cmd,
        "returncode": result.returncode,
        "relocated": True,
        "source_libraries_external": True,
    }
    report["comparison_gate"] = (
        result.returncode == 0
        and len(manifest["cases"]) == 6
        and {r["case"] for r in manifest["cases"]}
        == {
            "thioformaldehyde",
            "thioacetone",
            "pyridinium",
            "guanidinium",
            "alpha_lactam",
            "beta_lactam",
        }
        and len(report["sources"]) == 4
    )
    write(args.output / "report.json", report)
    return 0 if report["comparison_gate"] and not args.require_full_source else 1


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--declaration", type=Path, required=True)
    p.add_argument("--output", type=Path)
    p.add_argument("--assessments", type=Path)
    p.add_argument("--tolerances", type=Path)
    p.add_argument("--child", type=Path)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    if a.child:
        return child(a)
    if not all((a.output, a.assessments, a.tolerances)):
        p.error("--output, --assessments and --tolerances required")
    existed_before = a.output.exists()
    try:
        return run(a)
    except Exception as e:
        if (
            not existed_before
            and a.output.exists()
            and not (a.output / "failure.json").exists()
        ):
            write(
                a.output / "failure.json",
                {"error": f"{type(e).__name__}: {e}", "comparison_gate": False},
            )
        raise


if __name__ == "__main__":
    raise SystemExit(main())
