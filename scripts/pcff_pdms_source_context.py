"""Reuse pinned local literature reviews; keep every numerical version separate."""

import json
import subprocess
from hashlib import sha256
from pathlib import Path

from scripts.phase_4j54 import signed

CATALOG = "docs/evidence/sun_siloxane_parameters/parameters.json"
CATALOG_PIN = "d23244f11d6be8fb683ab01b1f586ee466f60c3b65d10b8da9379eda111961ca"
SCOPE = {
    "pcff/ma00107a006.pdf": "Sun1995 silanes; not a complete same-version polysiloxane model; eight printed angle blanks preserved",
    "pcff/1-s2.0-S1386142597000139-main.pdf": "Sun1997 siloxane rows; referenced C/H base and complete couplings unavailable; selective inheritance is not a free numerical overlay",
    "pcff/main.pdf": "Patil/Mark2000 aromatic polysiloxane application; no new transferable coefficient rows in the reviewed tables",
    "pcff/jp047434r.pdf": "Smith2004 Buckingham exp-6 and sine-squared bend-torsion; incompatible with requested LJ9-6 Class-II domain",
    "pcff/COMPASS, An ab Initio Force-Field Optimized forCondensed-Phase Applications.pdf": "COMPASS organic development source; no complete source-bound methyl-PDMS domain from this asset",
    "pcff/compass_ff.pdf": "COMPASS functional form alone does not supply the exact missing graph couplings",
    "pcff/jp4c08471.pdf": "Xiang benchmark; models and source families remain separate",
    "pcff/jp4c08471_si_001.pdf": "Tables S2/S3 Class-II candidates have recorded missing/cross-dependency conflicts; S4 OPLS rows are excluded from PCFF assignment",
}


def review(repo, references, *, pdf_python):
    repo, references = Path(repo), Path(references)
    catalog = repo / CATALOG
    if sha256(catalog.read_bytes()).hexdigest() != CATALOG_PIN:
        raise ValueError("historical_literature_catalog_changed")
    data = json.loads(catalog.read_text())
    old = json.loads((repo / "docs/evidence/phase_4j44/source-audit.json").read_text())
    known = {r["filename"]: r for r in old["files"]}
    for s in data["sources"].values():
        known[s["file"]] = {
            **s,
            "bibliography": {k: s[k] for k in ("title", "citation", "doi")},
        }
    main = json.loads(
        (repo / "docs/evidence/pcff_main_pdms_review/review.json").read_text()
    )["source"]
    known["pcff/main.pdf"] = main
    sources = []
    for filename, scope in SCOPE.items():
        path = references / filename
        expected = known[filename]["sha256"]
        if sha256(path.read_bytes()).hexdigest() != expected:
            raise ValueError("literature_source_replacement:" + filename)
        sources.append(
            {
                "file": filename,
                "sha256": expected,
                "citation": known[filename].get(
                    "bibliography", known[filename].get("title")
                ),
                "scope": scope,
                "parameter_overlay_authorized": False,
            }
        )
    command = [
        str(pdf_python),
        str(repo / "scripts/audit_sun_siloxane_parameters.py"),
        "--catalog",
        str(catalog),
        "--sha256",
        CATALOG_PIN,
        "--references",
        str(references),
    ]
    run = subprocess.run(command, capture_output=True, text=True, check=True)
    return signed(
        {
            "sources": sources,
            "catalog": {"file": CATALOG, "sha256": CATALOG_PIN},
            "independent_original_PDF_reconstruction": json.loads(run.stdout),
            "source_row_mapping": [
                {
                    "id": r["id"],
                    "source": r["source"],
                    "family": r["family"],
                    "labels": r["labels"],
                    "location": r["source_evidence"],
                    "evidence_checksum": r["evidence_checksum"],
                }
                for r in data["rows"]
            ],
            "ClassII_diagnostics_from_J44": [
                r
                for r in old["dependency_findings"]
                if r.get("source") == "pcff/jp4c08471_si_001.pdf"
                and "OPLS" not in r.get("interaction", "")
            ],
            "review_policy": "known source hashes and prior evidence reused; no blind recursive literature inventory; no cross-version borrowing",
        }
    )
