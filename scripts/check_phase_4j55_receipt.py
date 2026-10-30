"""Re-sign tampered J55 records and prove independent validation still rejects."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path

from island.opls_conventions import signed
from scripts.phase_4j47 import write
from scripts.validate_phase_4j47_receipt import pin, require
from scripts.validate_phase_4j55_receipt import validate


def sign(path, value):
    write(path, signed({k: v for k, v in value.items() if k != "identity"}))


def controls(directory, repo):
    original = Path(directory).resolve()
    result = []
    for mutation in (
        "readiness",
        "equivalence",
        "physical_convention",
        "charge",
        "profile",
        "source",
        "missing_row",
        "runtime",
        "graph",
    ):
        with tempfile.TemporaryDirectory(prefix="island-j55-receipt-") as tmp:
            root = Path(tmp) / "evidence"
            shutil.copytree(original, root)
            r = json.loads((root / "receipt.json").read_text())
            summary = json.loads((root / "audit/summary.json").read_text())
            case = root / "audit/dp0"
            p = json.loads((case / "report.json").read_text())
            if mutation == "readiness":
                r["production_validated"] = True
            elif mutation == "equivalence":
                r["cross_engine_equivalence"] = "passed"
            elif mutation == "physical_convention":
                p["conventions"]["dielectric"]["value"] = "1"
                p["missing_conventions"].remove("dielectric")
            elif mutation == "charge":
                p["charge_record"]["vector"][0][1] = "0"
            elif mutation == "profile":
                p["selected_profile"] = "pcff"
            elif mutation == "source":
                p["source_assignment"]["source_identity"] = "0" * 64
            elif mutation == "missing_row":
                p["interaction_inventory"].pop(0)
            elif mutation == "runtime":
                p["lammps_runtime_status"] = "failed"
            elif mutation == "graph":
                from island.charge_references.records import pack, unpack

                raw = unpack((case / "bundle/system.json").read_text())
                raw["coordinates"][next(iter(raw["coordinates"]))] = (1.0, 2.0, 3.0)
                (case / "bundle/system.json").write_text(pack(raw))
            sign(case / "report.json", p)
            shutil.copyfile(case / "report.json", case / "bundle/report.json")
            bundle = json.loads((case / "bundle/manifest.json").read_text())
            bundle["files"] = {n: pin(case / "bundle" / n) for n in bundle["files"]}
            sign(case / "bundle/manifest.json", bundle)
            c = summary["cases"][1]
            c["report_identity"] = json.loads((case / "report.json").read_text())[
                "identity"
            ]
            c["bundle_manifest_sha256"] = pin(case / "bundle/manifest.json")
            if mutation == "runtime":
                c["lammps_runtime_status"] = "failed"
            sign(root / "audit/summary.json", summary)
            r["summary_identity"] = json.loads(
                (root / "audit/summary.json").read_text()
            )["identity"]
            m = json.loads((root / "manifest.json").read_text())
            m["artifacts"] = {n: pin(root / n) for n in m["artifacts"]}
            sign(root / "manifest.json", m)
            r["manifest_sha256"] = pin(root / "manifest.json")
            sign(root / "receipt.json", r)
            try:
                validate(root, repo=repo)
            except (ValueError, KeyError, TypeError) as error:
                result.append(
                    {"mutation": mutation, "rejected": True, "reason": str(error)}
                )
            else:
                raise AssertionError("Accepted re-signed mutation: " + mutation)
    require(all(r["rejected"] for r in result), "Accepted tampering")
    return {
        "preliminary_receipt_sha256": pin(original / "receipt.json"),
        "controls": result,
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(argv)
    write(a.output, controls(a.directory, a.repo))


if __name__ == "__main__":
    main()
