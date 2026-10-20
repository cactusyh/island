"""Read-only revalidation of retained J31 scientific/software assets.

Copies analysis bundles only to a temporary relocation directory. No dynamics,
parameter generation or publication of historical evidence is performed.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from hashlib import sha256
from pathlib import Path

from island.analysis import load_network_analysis_bundle
from island.analysis.network import FLAGS
from scripts.phase_4j31_contract import CONTRACT
from scripts.validate_phase_4j31_receipt import validate_receipt

REVIEWED = "6606895d643806c731b52cf0fafa00e338ce3adb"
RECEIPT = Path("docs/evidence/phase_4j31.json")
MANIFEST = Path("docs/paper/reproducibility_manifest.json")


def historical_hashes():
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", REVIEWED, "docs/evidence"], text=True
    ).splitlines()
    result = {}
    for path in paths:
        original = subprocess.check_output(["git", "show", f"{REVIEWED}:{path}"])
        expected = sha256(original).hexdigest()
        if sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError("Historical evidence modified: " + path)
        result[path] = expected
    for path in ("docs/paper/reproducibility_manifest.json", ".gitattributes"):
        original = subprocess.check_output(["git", "show", f"{REVIEWED}:{path}"])
        if Path(path).read_bytes() != original:
            raise ValueError(
                "Historical manifest/publication attributes modified: " + path
            )
    return result


def revalidate():
    before = historical_hashes()
    receipt = validate_receipt(
        RECEIPT,
        MANIFEST,
        expected_identity="32a5688d61cd367658c62dbab094794e969fac4046dc256c540f110ec9f91907",
    )
    cases = {}
    with tempfile.TemporaryDirectory(prefix="island-j31-c1-relocation-") as temporary:
        for name, row in receipt["force_field_source_completeness"]["matrix"].items():
            if row["status"] != "complete":
                continue
            source = (
                Path(receipt["acceptance_directory"])
                / name.replace("/", "-")
                / "analysis-relocated"
            )
            target = Path(temporary) / name.replace("/", "-")
            shutil.copytree(source, target)
            _, selection, record = load_network_analysis_bundle(
                target, expected_identity=row["execution"]["analysis_record_identity"]
            )
            code = 'from island.analysis import load_network_analysis_bundle; import sys; print(load_network_analysis_bundle(sys.argv[1],expected_identity=sys.argv[2])[2].to_json(),end="")'
            child = subprocess.run(
                [sys.executable, "-c", code, str(target), record.identity],
                capture_output=True,
                text=True,
                check=True,
                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
            )
            if child.stdout != record.to_json():
                raise ValueError("Fresh-process reconstruction mismatch: " + name)
            cases[name] = {
                "analysis_record_identity": record.identity,
                "selection_identity": selection.identity,
                "observables_match": record.payload["observables"]
                == row["execution"]["network_observables"],
                "same_process_relocation": "PASS",
                "fresh_process_reconstruction": "PASS",
                "archived_workflow_comparison": row["execution"]["comparison"],
            }
    if before != historical_hashes():
        raise ValueError("Historical evidence changed during revalidation")
    return {
        "contract": CONTRACT,
        "reviewed_commit": REVIEWED,
        "historical_receipt_identity": receipt["identity"],
        "historical_file_hashes": before,
        "revalidated_cases": cases,
        "new_dynamics_executed": False,
        **FLAGS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = revalidate()
    args.output.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(
        "Retained J31: six bundles, workflows, relocation and fresh reconstruction PASS"
    )


if __name__ == "__main__":
    main()
