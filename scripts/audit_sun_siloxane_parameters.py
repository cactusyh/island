"""Independently reconstruct the supplement from the two original local PDFs.

Does not import the production loader or the transcription/generation script.
Only explicit PDF-OCR corrections recorded after visual inspection are applied.
"""

import argparse
import json
import re
from collections import Counter
from decimal import Decimal
from hashlib import sha256
from pathlib import Path


def audit(catalog, expected_sha256, references):
    from pypdf import PdfReader, __version__

    raw = Path(catalog).read_bytes()
    if sha256(raw).hexdigest() != expected_sha256:
        raise ValueError("catalog_hash_mismatch")
    data = json.loads(raw)
    readers = {}
    for key, source in data["sources"].items():
        path = Path(references) / source["file"]
        if sha256(path.read_bytes()).hexdigest() != source["sha256"]:
            raise ValueError("source_hash_mismatch:" + key)
        readers[key] = PdfReader(path)
    evidence = data["numeric_evidence"]
    text = (Path(catalog).parent / evidence["file"]).read_bytes()
    assert sha256(text).hexdigest() == evidence["sha256"]
    lines = text.decode().splitlines()
    correction_count = 0
    for row in data["rows"]:
        e = row["source_evidence"]
        page = readers[row["source"]].pages[e["pdf_page"] - 1]
        x0, y0, x1, y1 = e["bbox_pdf_points"]
        cells = []

        def visit(t, cm, tm, font, size, bounds=(x0, y0, x1, y1), cells=cells):
            x0, y0, x1, y1 = bounds
            x = tm[4] * cm[0] + tm[5] * cm[2] + cm[4]
            y = tm[4] * cm[1] + tm[5] * cm[3] + cm[5]
            if t.strip() and x0 <= round(x, 2) <= x1 and y0 <= round(y, 2) <= y1:
                cells.append((round(x, 2), t.strip()))

        page.extract_text(visitor_text=visit)
        text = " | ".join(t for _, t in sorted(cells, key=lambda x: x[0]))
        assert text == e["raw_pdf_text"], row["id"]
        for correction in e["ocr_corrections"]:
            assert correction["raw"] in text
            text = text.replace(
                correction["raw"], correction["page_image_transcription"]
            )
            correction_count += 1
        text = re.sub(r"-\s*-\s*", "-", text)
        text = re.sub(r"-\s+(?=\d)", "-", text)
        numbers = re.findall(r"-?\d+\.\d+", text)
        if row["source"] == "sun1995" and row["family"] == "nonbond(9-6)":
            numbers = [numbers[-1], numbers[-2]]
        assert numbers == e["numeric_decimals_reconstructed"], row["id"]
        assert numbers == [v for v in row["original_values"] if v is not None]
        assert all(Decimal(v).is_finite() for v in numbers)
        assert row["conversion"]["converted_values"] == row["original_values"]
        assert row["conversion"]["converted_units"] == row["source_units"]
        citation = row["numeric_record"]["citations"][0]
        assert lines[citation["line"] - 1] == citation["raw"]
        assert [citation["raw"].split()[i] for i in citation["columns"]] == numbers
        assert row["missing_fields"] == [
            f for f, v in zip(row["fields"], row["original_values"]) if v is None
        ]
    for quote in data["inheritance_evidence"]:
        text = readers[quote["source"]].pages[quote["page"] - 1].extract_text()
        assert quote["raw_text"] in text
    return {
        "source_pdf_hashes_verified": len(readers),
        "parameter_rows_reconstructed": len(data["rows"]),
        "rows_by_source": dict(Counter(r["source"] for r in data["rows"])),
        "partial_rows_not_zero_filled": sum(
            bool(r["missing_fields"]) for r in data["rows"]
        ),
        "explicit_page_image_ocr_corrections": correction_count,
        "pypdf_version": __version__,
        "execution_authorized": False,
        "scientific_qualification": "not_established",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    parser.add_argument("--references", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.catalog, args.sha256, args.references), indent=2))


if __name__ == "__main__":
    main()
