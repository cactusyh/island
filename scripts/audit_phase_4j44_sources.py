"""Independently reconstruct reviewed citations from explicit local assets.

Optional pypdf is required only by this offline extraction command. It never
fetches files, executes archives, installs packages, or grants MD authority.
"""

import argparse
import json
from decimal import Decimal
from hashlib import sha256
from pathlib import Path
from zipfile import ZipFile

from island.forcefields.literature_audit import (
    LiteratureSourceReview,
    inventory_literature,
)


def reconstruct(references, catalog):
    import pypdf

    review = LiteratureSourceReview(Path(catalog).read_text())
    result = review.validate_inventory(inventory_literature(references))
    pages, members, excerpts, headings = {}, 0, 0, 0
    for entry in review.payload["files"]:
        path = Path(references) / entry["filename"]
        if entry["file_type"] == "pdf":
            pdf = pypdf.PdfReader(path)
            assert {str(k): str(v) for k, v in (pdf.metadata or {}).items()} == entry[
                "pdf_metadata"
            ]
            texts = [page.extract_text() or "" for page in pdf.pages]
            assert len(texts) == entry["page_count"], entry["filename"]
            pages[entry["filename"]] = texts
            for evidence in entry["evidence"]:
                text = texts[evidence["page"] - 1]
                assert evidence["excerpt"] in text, entry["filename"]
                assert sha256(text.encode()).hexdigest() == evidence["text_sha256"]
                excerpts += 1
            for heading in entry["table_locations"]:
                assert heading["raw_heading"] in texts[heading["page"] - 1].splitlines()
                headings += 1
        else:
            with ZipFile(path) as archive:
                observed = [
                    {
                        "name": info.filename,
                        "size": info.file_size,
                        "sha256": sha256(archive.read(info)).hexdigest(),
                    }
                    for info in archive.infolist()
                    if not info.is_dir()
                ]
            assert observed == entry["archive_members"], entry["filename"]
            members += len(observed)
    for finding in review.payload["dependency_findings"]:
        if finding["source"] in pages:
            for raw in finding.get("source_rows", []):
                assert raw in "\n".join(pages[finding["source"]])
    zero_tokens = 0
    for row in review.payload["extracted_row_candidates"]:
        text = pages[row["source_file"]][row["page"] - 1].splitlines()[row["line"] - 1]
        assert text == row["raw_text"]
        # Independent direct Decimal reconstruction, outside the production/J43 reader.
        values = [text.split()[i] for i in row["numeric_columns"]]
        assert values == row["original_values"]
        zero_tokens += sum(Decimal(v) == 0 for v in values)
    return {
        **result,
        "raw_source_files_verified": len(review.payload["files"]),
        "pdf_files_reconstructed": len(pages),
        "archive_members_verified": members,
        "citation_excerpts_reconstructed": excerpts,
        "table_headings_reconstructed": headings,
        "explicit_zero_tokens_preserved_as_unqualified_evidence": zero_tokens,
        "pypdf_version": pypdf.__version__,
        "backend_executions": [],
        "acceptance_gate": 1,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--catalog", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(reconstruct(args.references, args.catalog)))


if __name__ == "__main__":
    main()
