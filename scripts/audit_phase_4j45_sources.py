"""Independent offline reconstruction of the selected J45 table and citations.

Uses the explicit local files only. No production coefficient parser, network,
parameter inference or backend is involved. Optional pypdf is offline tooling.
"""

import argparse
import json
from decimal import Decimal
from hashlib import sha256
from pathlib import Path


def reconstruct(references, selection, catalog):
    from pypdf import PdfReader, __version__

    assert (
        sha256(Path(selection).read_bytes()).hexdigest()
        == "64ecdd14316fcc740210853b96917f4f674dfb8ec63d0ebc244ec8030643ede9"
    ), "selection_evidence_mutated"
    evidence = json.loads(Path(selection).read_text())
    prior = json.loads(Path(catalog).read_text())
    # The independent reader verifies bytes directly, outside pinned_bytes/J43.
    assert (
        sha256(Path(catalog).read_bytes()).hexdigest() == evidence["j44_catalog_sha256"]
    )
    files = {f["filename"]: f for f in prior["files"]}
    for filename, entry in files.items():
        raw = (Path(references) / filename).read_bytes()
        assert (
            len(raw) == entry["size"] and sha256(raw).hexdigest() == entry["sha256"]
        ), filename
    pages = {}
    for record in evidence["candidate_rows"] + evidence["citations"]:
        filename = record["source_file"]
        assert record["source_sha256"] == files[filename]["sha256"]
        assert record["doi"] == files[filename]["bibliography"]["doi"]
        if filename not in pages:
            pages[filename] = [
                (page.extract_text() or "").splitlines()
                for page in PdfReader(Path(references) / filename).pages
            ]
        raw = pages[filename][record["page"] - 1][record["line"] - 1]
        assert raw == record["raw_text"], (filename, record["page"], record["line"])
        if "numeric_columns" in record:
            values = [raw.split()[i] for i in record["numeric_columns"]]
            assert values == record["original_values"]
            assert all(Decimal(v).is_finite() for v in values)
    # Independently reconstruct both charge columns, without assigning any atom.
    charges = {}
    for record in evidence["citations"]:
        if record["page"] == 6 and record["source_file"].endswith("si_001.pdf"):
            tokens = record["raw_text"].split()
            if len(tokens) == 3 and tokens[0] in ("Si", "O", "C", "H"):
                charges[tokens[0]] = tokens[1:]
                assert all(Decimal(v).is_finite() for v in tokens[1:])
    assert set(charges) == {"Si", "O", "C", "H"}
    assert len(evidence["candidate_rows"]) == 17
    assert evidence["accepted_rows"] == evidence["execution_profiles"] == []
    return {
        "raw_source_files_verified": len(files),
        "selected_coefficient_rows_reconstructed": 17,
        "selected_charge_rows_reconstructed": len(charges),
        "charge_columns_unassigned": charges,
        "citation_lines_reconstructed": len(evidence["citations"]),
        "pypdf_version": __version__,
        "acceptance_gate": 1,
        "execution_profiles": [],
        "backend_executions": [],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--selection", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(reconstruct(args.references, args.selection, args.catalog)))


if __name__ == "__main__":
    main()
