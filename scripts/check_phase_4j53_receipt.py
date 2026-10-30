"""Adversarially re-sign modified receipts; semantic validation must still fail."""

import argparse
import json
import shutil
import tempfile
from pathlib import Path

from island.opls_conventions import signed
from scripts.phase_4j47 import write
from scripts.validate_phase_4j47_receipt import pin
from scripts.validate_phase_4j53_receipt import validate


def controls(directory, repo):
    rows = []
    cases = (
        "readiness",
        "MD",
        "new_profile",
        "status",
        "gate",
        "lost_blocker",
        "permission",
        "convention",
        "charge",
        "source_hash",
        "reopen",
        "schema",
    )
    for case in cases:
        with tempfile.TemporaryDirectory(prefix="island-j53-receipt-") as temp:
            root = Path(temp) / "evidence"
            shutil.copytree(directory, root)
            receipt = json.loads((root / "receipt.json").read_text())
            manifest = json.loads((root / "manifest.json").read_text())
            path, value = None, None
            if case == "readiness":
                receipt["production_validated"] = True
            elif case == "MD":
                receipt["MD_executed"] = True
            elif case == "new_profile":
                receipt["new_profiles"] = ["unauthorized"]
            elif case == "status":
                receipt["source_status"] = "source_complete_independently_verified"
            elif case == "gate":
                receipt["research_opt_in_acceptance_gate"] = 0
            elif case == "lost_blocker":
                receipt["permanent_blockers"].pop()
            elif case == "permission":
                path = root / "audit/source-intake.json"
                value = json.loads(path.read_text())
                value["automatically_selectable"] = True
            elif case == "convention":
                path = root / "audit/dp2/preparation.json"
                value = json.loads(path.read_text())
                value["source_authoritative_settings"]["settings"]["dielectric"][
                    "value"
                ] = "1"
            elif case == "charge":
                path = root / "audit/dp2/system.json"
                # Invalid packed bytes remain rejected even if every outer pin is updated.
                value = {"payload": "modified charge", "sha256": "0" * 64}
            elif case == "source_hash":
                manifest["source_hashes"][next(iter(manifest["source_hashes"]))] = (
                    "0" * 64
                )
            elif case == "reopen":
                path = root / "audit/blocker-reopen.json"
                value = json.loads(path.read_text())
                value["rows"][0]["reopened"] = True
            else:
                path = root / "audit/source-intake-schema.json"
                value = json.loads(path.read_text())
                value["additionalProperties"] = True
            if path is not None:
                if "identity" in value:
                    value = signed({k: v for k, v in value.items() if k != "identity"})
                write(path, value)
                manifest["artifacts"][str(path.relative_to(root))] = pin(path)
            manifest = signed({k: v for k, v in manifest.items() if k != "identity"})
            write(root / "manifest.json", manifest)
            receipt["manifest_sha256"] = pin(root / "manifest.json")
            receipt = signed({k: v for k, v in receipt.items() if k != "identity"})
            write(root / "receipt.json", receipt)
            try:
                validate(root, repo=repo)
            except (ValueError, KeyError, TypeError) as error:
                rows.append({"case": case, "rejected": True, "reason": str(error)})
            else:
                raise ValueError("Accepted re-signed mutation: " + case)
    return {"all_rejected": True, "cases": rows}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args(argv)
    write(args.output, controls(args.directory, args.repo))


if __name__ == "__main__":
    main()
