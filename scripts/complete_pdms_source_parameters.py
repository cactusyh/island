"""Reconstruct available PDMS rows from one pinned FRC; do not fill source gaps."""

import argparse
import json
from collections import Counter
from decimal import Decimal
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff.class2 import inspect_pcff_class2
from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from island.forcefields.pdms_source_completion import (
    DATA_ID,
    PDMSCompletionError,
    PDMSCompletionRecord,
    assess_pdms_source_completion,
    file_hash,
    source_row,
)
from island.forcefields.siloxane_literature import load_siloxane_parameters
from island.workflows.bundle import system_data, system_from
from scripts.assess_sun_siloxane_parameters import fixtures

PAPER_PIN = "d23244f11d6be8fb683ab01b1f586ee466f60c3b65d10b8da9379eda111961ca"
ANGLE_LINES = {
    "Si-Si-H": 2410,
    "H-Si-H": 2409,
    "Si-Si-Si": 2411,
    "H-Si-C": 2406,
    "Si-C-H": 2196,
    "Si-Si-C": 2407,
    "Cp-Cp-Si": 2282,
    "Si-C-C": 2158,
}


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def corroborate(supplement, distribution):
    """Corroborate whole FRC rows; never overwrite paper blanks or round FRC."""
    catalog = inspect_pcff_class2(distribution.source)["records"]
    rows = []
    for paper in supplement.payload["rows"]:
        if paper["missing_fields"] != ["k4"]:
            continue
        row = source_row(
            catalog[f"quartic_angle:cff91:{ANGLE_LINES[paper['labels']]}"],
            distribution.hashes["pcff.frc"],
        )
        for printed, value in zip(
            paper["original_values"][:3], row["original_values"][:3], strict=True
        ):
            if Decimal(value).quantize(Decimal(printed)) != Decimal(printed):
                raise PDMSCompletionError(
                    "source_incompatible", "paper_frc_rounding_mismatch"
                )
        if Decimal(row["original_values"][3]) != 0:
            raise PDMSCompletionError("source_incompatible", "unexpected_source_k4")
        rows.append(
            {
                "paper_row": paper["id"],
                "paper_evidence_checksum": paper["evidence_checksum"],
                "paper_original_values_unchanged": paper["original_values"],
                "complete_frc_row": row,
                "agreement": "theta/K2/K3 agree at printed paper precision; FRC explicitly prints k4=0.0000",
                "assignment_scope": "silane si, not siloxane sio; whole FRC tuple only",
            }
        )
    if len(rows) != 8:
        raise PDMSCompletionError(
            "source_incompatible", "eight_angle_evidence_not_reconstructed"
        )
    return {
        "paper_catalog_pin": PAPER_PIN,
        "paper_catalog_identity": supplement.identity,
        "status": "source_corroboration_only",
        "rows": rows,
        "execution_authorized": False,
    }


def run(distribution_path, catalog, references, output, *, reconstruct=False):
    distribution = load_pcff_sip_distribution(distribution_path)
    supplement = load_siloxane_parameters(catalog, PAPER_PIN, references=references)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    comparisons = corroborate(supplement, distribution)
    if reconstruct:
        if (
            json.loads((output / "eight-angle-source-corroboration.json").read_text())
            != comparisons
        ):
            raise PDMSCompletionError("source_incompatible", "angle_evidence_changed")
    else:
        write(output / "eight-angle-source-corroboration.json", comparisons)
    summaries, all_rows, gaps, rejected = [], {}, [], []
    for name, system in fixtures():
        if name == "si_o_c_negative":
            # This graph is outside a proven PDMS template environment. Existing
            # Si/O typing decisions remain conservative and are not widened.
            try:
                assess_pdms_source_completion(system, distribution)
            except PDMSCompletionError as error:
                rejected.append(
                    {
                        "fixture": name,
                        "status": error.status,
                        "diagnostics": error.diagnostics,
                    }
                )
            else:
                raise PDMSCompletionError(
                    "source_incompatible", "unsupported_si_o_c_was_not_rejected"
                )
            continue
        result = assess_pdms_source_completion(system, distribution)
        if reconstruct:
            rebuilt = system_from(unpack((output / f"{name}-graph.json").read_text()))
            stored = PDMSCompletionRecord(
                (output / f"{name}-parameters.json").read_text()
            )
            stored.validate_integrity(rebuilt, distribution)
            if stored != result:
                raise PDMSCompletionError(
                    "source_incompatible", "fixture_or_end_groups_changed"
                )
        else:
            (output / f"{name}-graph.json").write_text(pack(system_data(system)))
            (output / f"{name}-parameters.json").write_text(result.json_text)
        p = result.payload
        all_rows.update(p["source_rows"])
        gaps.extend({"fixture": name, **g} for g in p["missing_dependencies"])
        summaries.append(
            {
                "fixture": name,
                "status": p["status"],
                "record_identity": result.identity,
                "graph_identity": p["requirements"]["graph_identity"],
                "atoms": len(system.topology.sites),
                "bonds": len(system.topology.bonds),
                "typing_complete": p["typing"]["coverage"]["complete"],
                "charges_complete": p["native_charge_record"]["complete"],
                "selected_unique_source_rows": len(p["source_rows"]),
                "coverage": p["coverage"],
                "missing_interaction_requests": len(p["missing_dependencies"]),
            }
        )
    for row in comparisons["rows"]:
        r = row["complete_frc_row"]
        all_rows[r["id"]] = r
    # Exact site-specific gaps plus grouped source keys allow the next source
    # acquisition/fitting task to target coefficients rather than atom counts.
    grouped = Counter((g["family"], tuple(g.get("ordinary_types", []))) for g in gaps)
    report = {
        "data_identity": DATA_ID,
        "source_hashes": distribution.hashes,
        "paper_catalog_pin": PAPER_PIN,
        "eight_angle_rows_corroborated": len(comparisons["rows"]),
        "fixtures": summaries,
        "negative_controls": rejected,
        "unique_source_rows_including_silane_comparisons": len(all_rows),
        "missing_requests": len(gaps),
        "missing_source_keys": len(grouped),
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "native_pcff_acceptance_gate": 1,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    grouped_rows = [
        {"family": f, "ordinary_types": list(t), "affected_requests": n}
        for (f, t), n in sorted(grouped.items())
    ]
    objects = {
        "report.json": report,
        "source-rows.json": all_rows,
        "missing-interactions.json": gaps,
        "missing-source-keys.json": grouped_rows,
    }
    for name, value in objects.items():
        if reconstruct:
            if json.loads((output / name).read_text()) != value:
                raise PDMSCompletionError(
                    "source_incompatible", "report_or_source_rows_changed:" + name
                )
        else:
            write(output / name, value)
    manifest = {
        p.name: file_hash(p)
        for p in sorted(output.iterdir())
        if p.is_file() and p.suffix == ".json" and p.name != "manifest.json"
    }
    if reconstruct:
        if json.loads((output / "manifest.json").read_text()) != manifest:
            raise PDMSCompletionError("source_incompatible", "manifest_changed")
    else:
        write(output / "manifest.json", manifest)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distribution", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--reconstruct", action="store_true")
    parser.add_argument("--require-executable", action="store_true")
    args = parser.parse_args()
    try:
        report = run(
            args.distribution,
            args.catalog,
            args.references,
            args.output,
            reconstruct=args.reconstruct,
        )
    except (PDMSCompletionError, PCFFError) as error:
        print(
            json.dumps(
                {
                    "status": getattr(error, "status", "source_incompatible"),
                    "diagnostics": getattr(error, "diagnostics", [str(error)]),
                }
            )
        )
        return 1
    print(json.dumps(report, indent=2))
    return 1 if args.require_executable else 0


if __name__ == "__main__":
    raise SystemExit(main())
