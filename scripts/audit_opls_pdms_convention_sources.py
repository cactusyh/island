"""Independent PDF/decimal reconstruction; no ISLAND parameter parser imports."""

import argparse
import json
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

from pypdf import PdfReader, __version__


def audit(references, evidence):
    document = json.loads(Path(evidence).read_text())
    checksum = document.pop("evidence_checksum")
    assert (
        sha256(
            (
                json.dumps(
                    document, sort_keys=True, separators=(",", ":"), allow_nan=False
                )
                + "\n"
            ).encode()
        ).hexdigest()
        == checksum
    )
    readers = {}
    for name, pin in document["source_hashes"].items():
        path = Path(references) / name
        assert sha256(path.read_bytes()).hexdigest() == pin
        readers[name] = PdfReader(path)
    for citation in document["citations"]:
        raw = (
            readers[citation["source_file"]].pages[citation["page"] - 1].extract_text()
        )
        assert raw.splitlines()[citation["line"] - 1] == citation["raw_text"]
    fig = document["figure7"]
    images = readers[fig["source_file"]].pages[fig["page"] - 1].images
    image = next(image for image in images if image.name == fig["image_resource"])
    assert sha256(image.data).hexdigest() == fig["image_sha256"]
    energy = document["unit_evidence"]["figure_s1"]
    images = readers[energy["source_file"]].pages[energy["page"] - 1].images
    image = next(image for image in images if image.name == energy["image_resource"])
    assert sha256(image.data).hexdigest() == energy["image_sha256"]
    si = next(reader for name, reader in readers.items() if name.endswith("si_001.pdf"))
    table = si.pages[5].extract_text().splitlines()
    charges = {
        tokens[0]: tokens[1:]
        for line in table
        if len(tokens := line.split()) == 3 and tokens[0] in ("O", "Si", "C", "H")
    }
    assert charges == fig["Qeq_Charge"]
    q = {e: list(map(Decimal, values)) for e, values in charges.items()}
    monomer = q["Si"][0] + q["O"][0] + 2 * q["C"][0] + 6 * q["H"][0]
    end_pair = 2 * (q["Si"][1] + 3 * q["C"][1] + 9 * q["H"][1]) + q["O"][1]
    assert monomer == end_pair == 0
    return {
        "source_files_verified": len(readers),
        "text_citations_verified": len(document["citations"]),
        "figure7_image_verified": True,
        "charge_rows_reconstructed": len(charges),
        "M_unit_charge": str(monomer),
        "E_pair_charge": str(end_pair),
        "visual_fragment_transcription": "image pinned; independent human visual review still applicable",
        "pypdf_version": __version__,
        "acceptance_gate": 1,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.references, args.evidence)))
