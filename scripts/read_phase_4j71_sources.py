"""Optional verification of external originals; never claims automatic paper reading.

All literature reading was performed in the J71 session. This command verifies
source bytes and exact archive commands against the reviewed ledger. It neither
OCRs tables nor imports a PDF library. Full papers remain outside git.
"""

import argparse
import zipfile
from pathlib import Path

from island.pcff_completion import pin, require
from island.pcff_qm_pilot.hmds_precision import read
from scripts.extract_phase_4j71 import ROOT, inputs


def verify_originals(repo, root):
    _, cells, reading, _ = inputs(repo, root)
    for source in reading["files"]:
        require(
            pin(repo / source["relative_path"]) == source["sha256"],
            "J71_original_source_hash_mutation",
        )
    archives = read(root / "archive-members.json")["archives"]
    checked = 0
    for archive in archives:
        with zipfile.ZipFile(repo / archive["path"]) as z:
            require(
                sorted(n for n in z.namelist() if not n.endswith("/"))
                == sorted(m["path"] for m in archive["members"]),
                "J71_archive_members_mutation",
            )
            from hashlib import sha256

            for member in archive["members"]:
                raw = z.read(member["path"])
                require(
                    sha256(raw).hexdigest() == member["sha256"]
                    and len(raw) == member["bytes"],
                    "J71_archive_source_hash_mutation",
                )
            for cell in cells:
                loc = cell["location"]
                if (
                    "archive_member" not in loc
                    or cell["source"] != archive["sha256"][:12]
                ):
                    continue
                line = (
                    z.read(loc["archive_member"])
                    .decode()
                    .splitlines()[loc["line"] - 1]
                    .strip()
                )
                require(line == loc["raw_snippet"], "J71_archive_source_row_mutation")
                checked += 1
    return {
        "original_files": len(reading["files"]),
        "archive_numeric_rows_checked": checked,
        "full_text_reading_claim_added": False,
        "QM": 0,
        "MD": 0,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--root", type=Path, default=ROOT)
    a = p.parse_args()
    print(verify_originals(a.repo.resolve(), a.root.resolve()))
