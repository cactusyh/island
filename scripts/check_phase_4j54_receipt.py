"""Re-sign J54 audit mutations and prove raw-source replay still rejects them."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path

from island.forcefields.selection.inventory import canonical
from scripts.phase_4j54 import signed, write
from scripts.validate_phase_4j47_receipt import pin
from scripts.validate_phase_4j54_receipt import validate_audit


def controls(directory, *, distribution, archive):
    records = []
    for mutation in (
        "unit",
        "charge",
        "functional_form",
        "missing_coupling",
        "source_hash",
        "profile",
        "topology",
        "readiness",
        "author_compass",
    ):
        with tempfile.TemporaryDirectory(prefix="island-j54-resigned-") as tmp:
            root = Path(tmp) / "audit"
            shutil.copytree(directory, root)
            summary = json.loads((root / "summary.json").read_text())
            case = next(c for c in summary["fixtures"] if c["fixture"] == "hmds")
            report = json.loads((root / "hmds/report.json").read_text())
            row = next(
                r
                for r in report["source_rows"].values()
                if r["family"] == "quartic_bond"
            )
            if mutation == "unit":
                row["original_units"][0] = "nm"
            elif mutation == "charge":
                report["charge_vector"][next(iter(report["charge_vector"]))] += 0.01
            elif mutation == "functional_form":
                row["functional_form"] = "harmonic"
            elif mutation == "missing_coupling":
                report["interaction_inventory"] = [
                    r
                    for r in report["interaction_inventory"]
                    if r["family"] != "bond-bond"
                ]
            elif mutation == "source_hash":
                row["source_sha256"] = "0" * 64
            elif mutation == "profile":
                report["selected_profile"] = report["candidate_profile"]
            elif mutation == "topology":
                from island.charge_references.records import pack, unpack

                system = unpack((root / "hmds/system.json").read_text())
                system["sites"].reverse()
                (root / "hmds/system.json").write_text(pack(system))
            elif mutation == "readiness":
                summary["production_validated"] = True
            else:
                author = json.loads((root / "author-compass.json").read_text())
                author["active_coefficient_rows"][0]["original_units"][0] = "kJ/mol"
                write(root / "author-compass.json", author)
            report = signed({k: v for k, v in report.items() if k != "identity"})
            (root / "hmds/report.json").write_text(canonical(report))
            manifest = json.loads((root / "hmds/manifest.json").read_text())
            manifest["files"] = {n: pin(root / "hmds" / n) for n in manifest["files"]}
            manifest["report_identity"] = report["identity"]
            (root / "hmds/manifest.json").write_text(
                canonical(
                    signed({k: v for k, v in manifest.items() if k != "identity"})
                )
            )
            case["manifest_sha256"] = pin(root / "hmds/manifest.json")
            case["report_identity"] = report["identity"]
            write(
                root / "summary.json",
                signed({k: v for k, v in summary.items() if k != "identity"}),
            )
            try:
                validate_audit(root, distribution=distribution, archive=archive)
            except (ValueError, KeyError, TypeError) as error:
                records.append(
                    {"mutation": mutation, "rejected": True, "reason": str(error)}
                )
            else:
                raise AssertionError("Accepted re-signed mutation: " + mutation)
    return {"controls": records, "source_or_charge_edits_applied_to_live_assets": False}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--distribution", type=Path, required=True)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    write(
        a.output, controls(a.directory, distribution=a.distribution, archive=a.archive)
    )


if __name__ == "__main__":
    main()
