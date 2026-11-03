"""Optional independent PDF table recheck, without OCR or full-text retention.

Requires the mounted J71 originals and a separately installed PyMuPDF. It does
not read new papers, alter extraction identities, or establish PCFF authority.
"""

import argparse
import hashlib
import json
import re
from decimal import Decimal
from pathlib import Path

NUMBER = r"[-+−]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


def numeric_row(text, label, count, command=""):
    pattern = r"(?<!\S)" + re.escape(command + label) + r"\s+"
    pattern += r"\s+".join("(" + NUMBER + ")" for _ in range(count))
    hits = re.findall(pattern, text)
    if len(hits) != 1:
        raise ValueError("ambiguous_or_missing_PDF_row: " + label)
    return list(hits[0])


def check(repo):
    import fitz  # Optional, isolated reader. Never imported by the ordinary API.

    root = repo / "docs/evidence/phase_4j71"
    plan = json.loads((repo / "docs/evidence/phase_4j72/contract.json").read_text())
    for name in ("parameter-records.json", "reading-status.json"):
        expected = plan["parent_pins"]["docs/evidence/phase_4j71/" + name]
        if hashlib.sha256((root / name).read_bytes()).hexdigest() != expected:
            raise ValueError("J71_pinned_input_mutation")
    records = json.loads((root / "parameter-records.json").read_text())["records"]
    sources = json.loads((root / "reading-status.json").read_text())["files"]
    checked, originals = [], []
    for prefix, page in (
        ("5372747c08f6", 8),
        ("1e95027a1202", 65),
        ("6e86aa640854", 4),
    ):
        source = next(
            s
            for s in sources
            if s["sha256"].startswith(prefix) and s["scope"] == "primary_directory"
        )
        path = repo / source["relative_path"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError("original_PDF_hash_mutation")
        with fitz.open(path) as document:
            text = document[page - 1].get_text(sort=True)
        originals.append(
            {
                "sha256": source["sha256"],
                "relative_path": source["relative_path"],
                "PDF_page": page,
            }
        )
        if prefix != "6e86aa640854" and "in kJ/mol" not in text:
            raise ValueError("printed_table_unit_unresolved")
        for r in records:
            if r["source_file_sha256"] != source["sha256"]:
                continue
            if prefix == "6e86aa640854":
                if r["interaction_family"] != "quartic_angle" or r["location"][
                    "row_label"
                ] not in ("c_si_o", "c_si_c"):
                    continue
                # PDF whitespace may wrap neighboring commands. Each exact
                # command/label followed by four numbers must be unique.
                flattened = re.sub(r"\s+", " ", text)
                tokens = numeric_row(
                    flattened, r["location"]["row_label"], 4, "angle_coeff "
                )
            else:
                if r["interaction_family"] != "RB_torsion":
                    continue
                tokens = numeric_row(text, "".join(r["ordered_sites"]), 6)
            if [Decimal(t.replace("−", "-")) for t in tokens] != [
                Decimal(v) for v in r["raw_value"]
            ]:
                raise ValueError(
                    "original_PDF_numeric_row_disagreement: " + r["record_id"]
                )
            checked.append(
                {
                    "record_id": r["record_id"],
                    "J71_identity": r["identity"],
                    "source_sha256": source["sha256"],
                    "PDF_page": page,
                    "table": r["location"]["table"],
                    "source_numeric_tokens": tokens,
                    "J71_raw_values": r["raw_value"],
                    "raw_units": r["raw_unit"],
                    "numeric_equality": True,
                }
            )
    if len(checked) != 34:
        raise ValueError("PDF_table_coverage")
    return {
        "schema": "island.j72-independent-PDF-recheck.v1",
        "reader": "PyMuPDF " + fitz.VersionBind,
        "originals": originals,
        "checked_rows": checked,
        "numeric_cells": sum(len(r["source_numeric_tokens"]) for r in checked),
        "OCR_used": False,
        "full_document_read_claim_added": False,
        "parameter_extraction_records_added": 0,
        "source_conflicts_repaired": 0,
        "native_PCFF_authority_established": False,
        "execution_count_scope": "this read-only checker invocation only; not phase totals",
        "new_QM": 0,
        "QM_fits": 0,
        "MD_steps": 0,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    a = p.parse_args()
    print(json.dumps(check(a.repo.resolve()), indent=2, ensure_ascii=False))
