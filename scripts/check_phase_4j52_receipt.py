"""Re-sign adversarial J52 records; prove semantic rejection and receipt relocation."""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from island.opls_conventions import signed
from island.opls_final_resolution import NINE_CONVENTIONS
from scripts.phase_4j47 import write
from scripts.validate_phase_4j47_receipt import pin, require
from scripts.validate_phase_4j52_receipt import validate


def controls(directory, repo):
    root, repo = directory.resolve(), repo.resolve()
    validate(root, repo=repo)
    outcomes = []
    with tempfile.TemporaryDirectory(prefix="island-j52-receipt-") as tmp:
        moved = Path(tmp) / "moved"
        shutil.copytree(root, moved)
        receipt_bytes = (moved / "receipt.json").read_bytes()
        manifest_bytes = (moved / "manifest.json").read_bytes()
        for name in (
            *NINE_CONVENTIONS,
            "missing_blockers",
            "raw_FD_energy",
            "source_identity",
            "profile",
            "readiness",
        ):
            r = json.loads(receipt_bytes)
            m = json.loads(manifest_bytes)
            changed = None
            if name in NINE_CONVENTIONS:
                changed = moved / "audit/final-nine-conventions.json"
                value = json.loads(changed.read_text())
                row = next(v for v in value["rows"] if v["setting"] == name)
                row["resolution"], row["value"] = "source_resolved", "invented"
                replacement = signed({k: v for k, v in row.items() if k != "identity"})
                row.clear()
                row.update(replacement)
            elif name == "missing_blockers":
                changed = moved / "audit/permanent-blockers.json"
                value = json.loads(changed.read_text())
                value["rows"] = []
            elif name == "raw_FD_energy":
                changed = moved / "audit/dp0/exact_erfc_diagnostic_v1/components.json"
                value = json.loads(changed.read_text())
                value["finite_difference"]["raw_records"][0]["components"]["coulomb"][
                    "positive_energy_kj_mol"
                ] += 1
            elif name == "source_identity":
                r["source_audit_identity"] = "0" * 64
            elif name == "profile":
                r["profile_identity"] = "native_pcff"
            elif name == "readiness":
                r["production_validated"] = True
            original = changed.read_bytes() if changed else None
            if changed:
                if "identity" in value:
                    value = signed({k: v for k, v in value.items() if k != "identity"})
                write(changed, value)
                m["artifacts"][str(changed.relative_to(moved))] = pin(changed)
            write(
                moved / "manifest.json",
                signed({k: v for k, v in m.items() if k != "identity"}),
            )
            r["manifest_sha256"] = pin(moved / "manifest.json")
            write(
                moved / "receipt.json",
                signed({k: v for k, v in r.items() if k != "identity"}),
            )
            try:
                validate(moved, repo=repo)
            except (ValueError, KeyError) as error:
                outcomes.append(
                    {
                        "control": name,
                        "rejected": True,
                        "reason": str(error),
                        "records_and_manifest_resigned": True,
                    }
                )
            else:
                raise ValueError("Re-signed tampering accepted: " + name)
            finally:
                if changed:
                    changed.write_bytes(original)
                (moved / "receipt.json").write_bytes(receipt_bytes)
                (moved / "manifest.json").write_bytes(manifest_bytes)
        proc = subprocess.run(
            [
                sys.executable,
                "-m",
                "scripts.validate_phase_4j52_receipt",
                str(moved),
                "--repo",
                str(repo),
                "--sha256",
                pin(root / "receipt.json"),
            ],
            capture_output=True,
            text=True,
            check=False,
            env={**os.environ, "PYTHONHASHSEED": "99397"},
        )
        require(
            proc.returncode == 0,
            "Separate-process relocated receipt failed: " + proc.stderr,
        )
        relocation = {
            "verified": True,
            "separate_process": True,
            "returncode": proc.returncode,
            "validator_result": json.loads(proc.stdout),
        }
    return {"controls": outcomes, "relocated_receipt": relocation}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(argv)
    result = controls(a.directory, a.repo)
    write(a.output, result)
    print(
        json.dumps(
            {
                "rejected_controls": len(result["controls"]),
                "relocated_receipt_verified": True,
            }
        )
    )


if __name__ == "__main__":
    main()
