"""Offline transcription/reconstruction of Sun 1995 and Sun--Rigby 1997 tables.

The seeds below were checked against the original page images. PDF coordinates,
verbatim OCR and explicit OCR corrections make every numerical row auditable.
This command produces data, never execution authority or a native FRC overlay.
"""

import argparse
import json
import re
from collections import Counter
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

SOURCES = {
    "sun1995": {
        "file": "pcff/ma00107a006.pdf",
        "sha256": "9c6d6315b91e0c15a3cba5f90f011a862da79ae193bfe66900bfbdc9c5c8fe48",
        "title": "Ab Initio Calculations and Force Field Development for Computer Simulation of Polysilanes",
        "citation": "H. Sun, Macromolecules 28 (1995) 701–712",
        "doi": "10.1021/ma00107a006",
    },
    "sun1997": {
        "file": "pcff/1-s2.0-S1386142597000139-main.pdf",
        "sha256": "65bda3a5b14cfd05bc778f2306aa6b7221dda8124d66b4d7c7fd817853d01fe7",
        "title": "Polysiloxanes: ab initio force field and structural, conformational and thermophysical properties",
        "citation": "H. Sun and D. Rigby, Spectrochimica Acta Part A 53 (1997) 1301–1323",
        "doi": None,
    },
}

# source, page, family, x limits, y limits, followed by literal label/value seeds.
# '?' is a printed blank, NEVER a zero. Coefficient units follow the equations;
# equilibrium angles are degrees and angular displacements are radians.
TABLES = [
    (
        "sun1995",
        4,
        "bond_increment",
        (330, 540),
        (250, 300),
        """
Si-H .126
Si-C .135
Si-Cp .117
C-H -.053
Cp-H -.140
""",
    ),
    (
        "sun1995",
        4,
        "nonbond(9-6)",
        (325, 565),
        (190, 210),
        """
Si 4.450 .190
H 2.940 .023
""",
    ),
    (
        "sun1995",
        11,
        "quartic_bond",
        (60, 300),
        (662, 696),
        """
Si-Si 2.338 114.22 -140.42 80.71
Si-H 1.478 202.78 -305.36 280.27
Si-C 1.899 189.65 -279.42 307.51
Si-Cp 1.863 233.24 -276.87 161.67
""",
    ),
    (
        "sun1995",
        11,
        "quartic_angle",
        (60, 300),
        (562, 645),
        """
Si-Si-H 112.09 22.51 -11.59 ?
H-Si-H 108.61 32.54 -8.32 ?
Si-Si-Si 114.27 24.95 -19.59 ?
H-Si-C 112.10 36.48 -12.81 ?
Si-C-H 112.04 28.77 -13.95 ?
Si-Si-C 113.00 19.47 -34.35 ?
C-Si-C 113.19 36.21 -20.39 20.02
Cp-Cp-Si 120.00 30.47 -23.54 ?
Cp-Si-H 109.59 41.95 -42.36 48.14
Si-C-C 112.67 39.52 -7.44 ?
""",
    ),
    (
        "sun1995",
        11,
        "torsion_3",
        (60, 300),
        (487, 545),
        """
*-Si-Si-* 0.00 0.00 -.07
*-Si-C-* 0.00 0.00 -.07
Si-Si-C-C -.35 0.00 -.07
Cp-Cp-Si-* 0.00 0.00 -.02
Si-Cp-Cp-Cp 0.00 4.33 0.00
Si-Cp-Cp-H 0.00 1.51 0.00
*-C-C-* 0.00 .05 -.14
""",
    ),
    (
        "sun1995",
        11,
        "wilson_out_of_plane",
        (90, 300),
        (458, 468),
        """
Cp-Cp-Si-Cp 5.37
""",
    ),
    (
        "sun1995",
        11,
        "bond-bond",
        (60, 180),
        (396, 428),
        """
Si-Si/Si-Si 6.07
Si-Si/Si-H 3.52
H-Si/Si-H 4.64
H-Si/Si-C 3.93
""",
    ),
    (
        "sun1995",
        11,
        "bond-bond",
        (182, 300),
        (396, 428),
        """
C-Si/Si-C 3.74
Cp-Cp/Cp-Si 21.39
Cp-Si/Si-H 3.93
Si-Si/Si-C 2.30
""",
    ),
    (
        "sun1995",
        11,
        "bond-angle",
        (60, 180),
        (334, 376),
        """
Si-Si-H 2.07 5.66
H-Si-H 9.35 9.35
Si-Si-Si 8.99 8.99
H-Si-C 7.41 13.40
Si-C-H 18.28 16.69
""",
    ),
    (
        "sun1995",
        11,
        "bond-angle",
        (182, 300),
        (340, 376),
        """
Si-Si-C 11.44 16.95
C-Si-C 18.58 18.58
Cp-Cp-Si 14.58 23.77
Cp-Si-H 22.59 8.78
""",
    ),
    (
        "sun1995",
        11,
        "angle-angle",
        (60, 180),
        (249, 315),
        """
H-Si-Si/Si-Si-H 1.61
H-Si-H/H-Si-H 2.07
Si-Si-Si/Si-Si-H 4.20
Si-Si-H/H-Si-Si 3.49
H-Si-H/H-Si-C 4.46
H-Si-C/C-Si-H 4.68
H-C-Si/Si-C-H 2.21
H-Si-Si/Si-Si-C -2.96
""",
    ),
    (
        "sun1995",
        11,
        "angle-angle",
        (182, 300),
        (249, 315),
        """
Si-Si-H/H-Si-C 3.48
C-Si-C/C-Si-H 3.38
C-Si-H/H-Si-C 2.80
Si-Si-Si/Si-Si-C 4.53
C-Si-Si/Si-Si-C 2.08
Si-Si-C/C-Si-Si -5.68
Si-Si-C/C-Si-C 1.35
""",
    ),
    (
        "sun1995",
        11,
        "middle_bond-torsion_3",
        (60, 300),
        (172, 239),
        """
H-Si-Si-H 0.00 0.00 -.63
H-Si-C-H 0.00 0.00 -.59
H-Si-Si-C 0.00 0.00 -.69
Si-Si-C-H 0.00 0.00 -.19
Si-Cp-Cp-Cp 0.00 11.16 0.00
Si-Cp-Cp-H 0.00 6.22 0.00
Cp-Cp-Si-H 0.00 0.00 -.31
""",
    ),
    (
        "sun1995",
        11,
        "angle-torsion_3",
        (60, 300),
        (130, 165),
        """
Si-Cp-Cp-Cp 0.00 4.33 0.00 0.00 -5.54 0.00
Si-Cp-Cp-H 0.00 1.111 0.00 0.00 4.59 0.00
Cp-Cp-Si-H 0.00 0.000 -.28 0.000 0.00 -.19
""",
    ),
    (
        "sun1995",
        11,
        "angle-angle-torsion_1",
        (60, 180),
        (91, 128),
        """
H-Si-Si-H -10.82
H-Si-Si-Si -12.29
H-Si-C-H -12.93
""",
    ),
    (
        "sun1995",
        11,
        "angle-angle-torsion_1",
        (182, 300),
        (91, 128),
        """
H-Si-Si-C -16.91
Si-Si-C-H -13.37
C-Si-C-H -17.58
""",
    ),
    (
        "sun1997",
        21,
        "bond_increment",
        (45, 220),
        (100, 150),
        """
o2-si4 -.2225
""",
    ),
    (
        "sun1997",
        21,
        "quartic_bond",
        (45, 450),
        (50, 110),
        """
o2-si4 1.640 350.1232 -517.3424 673.7067
""",
    ),
    (
        "sun1997",
        22,
        "quartic_angle",
        (45, 470),
        (560, 628),
        """
h1-o2-si4 122.8 23.7764 -19.8152 9.6331
si4-o2-si4 159.0 8.5000 -13.4188 -4.1785
c3a-si4-o2 114.9 23.0218 -31.3993 24.9814
c4-si4-o2 114.9 23.0218 -31.3993 24.9814
h1-si4-o2 107.4 57.6643 -10.6506 4.6274
o2-si4-o2 110.7 70.3069 -6.9375 0.0000
""",
    ),
    (
        "sun1997",
        22,
        "torsion_3",
        (45, 470),
        (445, 524),
        """
h1-o2-si4-c4 0.000 0.000 -.050
h1-o2-si4-h1 0.000 0.000 -.050
h1-o2-si4-o2 0.000 0.000 -.050
si4-o2-si4-c4 0.000 0.000 -.010
si4-o2-si4-h1 0.000 0.000 -.010
si4-o2-si4-o2 -.225 0.000 -.010
*-o2-si4-* 0.000 0.000 -.010
""",
    ),
    (
        "sun1997",
        22,
        "nonbond(9-6)",
        (45, 220),
        (378, 409),
        """
si4 4.405 .198
si4c 4.290 .131
o2 3.300 .080
""",
    ),
    (
        "sun1997",
        22,
        "bond-bond",
        (45, 470),
        (287, 342),
        """
h1-o2-si4 6.3820
si4-o2-si4 41.1143
c4-si4-o2 5.4896
h1-si4-o2 11.6183
o2-si4-o2 41.1143
""",
    ),
    (
        "sun1997",
        22,
        "bond-angle",
        (45, 470),
        (196, 250),
        """
h1-o2-si4 18.0902 31.0726
si4-o2-si4 28.6686 28.6686
c4-si4-o2 6.4278 20.5669
h1-si4-o2 6.4278 20.5669
o2-si4-o2 23.4380 23.4380
""",
    ),
]

FIELDS = {
    "quartic_bond": (
        ["r0", "k2", "k3", "k4"],
        [
            "angstrom",
            "kcal/(mol*angstrom^2)",
            "kcal/(mol*angstrom^3)",
            "kcal/(mol*angstrom^4)",
        ],
    ),
    "quartic_angle": (
        ["theta0", "k2", "k3", "k4"],
        ["degree", "kcal/(mol*radian^2)", "kcal/(mol*radian^3)", "kcal/(mol*radian^4)"],
    ),
    "torsion_3": (["k1", "k2", "k3"], ["kcal/mol"] * 3),
    "wilson_out_of_plane": (["k2"], ["kcal/(mol*radian^2)"]),
    "bond_increment": (["delta_ij"], ["elementary_charge"]),
    "nonbond(9-6)": (["r_min", "epsilon"], ["angstrom", "kcal/mol"]),
    "bond-bond": (["k"], ["kcal/(mol*angstrom^2)"]),
    "bond-angle": (["k_left", "k_right"], ["kcal/(mol*angstrom*radian)"] * 2),
    "angle-angle": (["k"], ["kcal/(mol*radian^2)"]),
    "middle_bond-torsion_3": (["k1", "k2", "k3"], ["kcal/(mol*angstrom)"] * 3),
    "angle-torsion_3": (
        ["k_left1", "k_left2", "k_left3", "k_right1", "k_right2", "k_right3"],
        ["kcal/(mol*radian)"] * 6,
    ),
    "angle-angle-torsion_1": (["k"], ["kcal/(mol*radian^2)"]),
}


def norm(text):
    return (
        text.strip()
        .replace("—", "-")
        .replace("–", "-")
        .replace(" ", "")
        .replace("--", "-")
    )


def chunks(page):
    result = []

    def visit(text, cm, tm, font, size):
        if text.strip():
            x = tm[4] * cm[0] + tm[5] * cm[2] + cm[4]
            y = tm[4] * cm[1] + tm[5] * cm[3] + cm[5]
            result.append((round(x, 2), round(y, 2), text.strip()))

    page.extract_text(visitor_text=visit)
    return result


def write(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n")


def extract(references, output):
    from pypdf import PdfReader

    from island.forcefields.external_pdms import LIMITS, digest

    output.mkdir(parents=True, exist_ok=True)
    readers = {}
    for sid, source in SOURCES.items():
        raw = (references / source["file"]).read_bytes()
        if sha256(raw).hexdigest() != source["sha256"]:
            raise ValueError("source_hash_mismatch:" + sid)
        readers[sid] = PdfReader(references / source["file"])
    rows, text_rows = [], []
    for sid, pn, family, xr, yr, seeds in TABLES:
        items = chunks(readers[sid].pages[pn - 1])
        seeds = [line.split() for line in seeds.strip().splitlines()]
        for index, (label, *values) in enumerate(seeds):
            # The 1997 OCR column order is damaged. Fixed row baselines preserve
            # the actual table roles; these coordinates are not chemical rules.
            if sid == "sun1997":
                baselines = {
                    (21, "bond_increment"): [129.0],
                    (21, "quartic_bond"): [86.25],
                    (22, "quartic_angle"): [
                        624.75,
                        612.75,
                        600.75,
                        588.75,
                        576.0,
                        564.0,
                    ],
                    (22, "torsion_3"): [
                        521.25,
                        508.5,
                        496.5,
                        484.5,
                        472.5,
                        460.5,
                        448.5,
                    ],
                    (22, "nonbond(9-6)"): [405.75, 393.75, 381.75],
                    (22, "bond-bond"): [338.25, 326.25, 314.25, 302.25, 290.25],
                    (22, "bond-angle"): [246.75, 234.75, 222.75, 210.75, 198.75],
                }
                y = baselines[pn, family][index]
            else:
                candidates = [
                    (x, y, t)
                    for x, y, t in items
                    if xr[0] <= x <= xr[1] and yr[0] <= y <= yr[1] and norm(t) == label
                ]
                if family == "wilson_out_of_plane":
                    candidates = [
                        (x, y, t)
                        for x, y, t in items
                        if t == "Cp-Cp-Si-" and yr[0] <= y <= yr[1]
                    ]
                if len(candidates) != 1:
                    raise ValueError(
                        f"label_not_unique:{sid}:{family}:{label}:{candidates}"
                    )
                y = candidates[0][1]
            selected = sorted(
                [
                    (x, y0, t)
                    for x, y0, t in items
                    if xr[0] <= x <= xr[1] and abs(y0 - y) <= 3.6
                ],
                key=lambda p: p[0],
            )
            raw_text = " | ".join(t for _, _, t in selected)
            corrected = raw_text
            corrections = []
            for old, new in [
                ("9.633 I", "9.6331"),
                ("3 1.3993", "31.3993"),
                ("i.198", "0.198"),
            ]:
                if old in corrected:
                    corrections.append({"raw": old, "page_image_transcription": new})
                    corrected = corrected.replace(old, new)
            corrected = re.sub(r"-\s*-\s*", "-", corrected)
            corrected = re.sub(r"-\s+(?=\d)", "-", corrected)
            found = re.findall(r"-?\d+\.\d+", corrected)
            printed = [v for v in values if v != "?"]
            if sid == "sun1995" and family == "nonbond(9-6)":
                # Table 5 prints initial epsilon,D then refined epsilon,D.
                found = [found[-1], found[-2]]
            if [Decimal(x) for x in found] != [Decimal(x) for x in printed]:
                raise ValueError(
                    f"numeric_mismatch:{sid}:{family}:{label}:{found}:{printed}:{raw_text}"
                )
            # Keep the paper's exact decimal spellings, including trailing zeros.
            # Seeds are visual checks, not an alternate numerical authority.
            printed = found
            exact = iter(printed)
            values = ["?" if v == "?" else next(exact) for v in values]
            fields, units = FIELDS[family]
            row_id = f"{sid}:{family}:{label}"
            numeric_line = " ".join([row_id, *printed])
            text_rows.append(numeric_line)
            values = [None if v == "?" else v for v in values]
            row = {
                "id": row_id,
                "source": sid,
                "family": family,
                "labels": label,
                "fields": fields,
                "original_values": values,
                "source_units": units,
                "unit_basis": "source equations: E=coefficient times bond/angle displacements; equilibrium angles in degrees, angular displacements in radians",
                "conversion": {
                    "policy": "identity_decimal_v1",
                    "converted_values": values,
                    "converted_units": units,
                    "code": "island.identity_decimal_v1",
                    "version": 1,
                },
                "missing_fields": [f for f, v in zip(fields, values) if v is None],
                "explicit_zero_fields": [
                    f
                    for f, v in zip(fields, values)
                    if v is not None and Decimal(v) == 0
                ],
                "source_evidence": {
                    "pdf_page": pn,
                    "location": "Table 5 (refined nonbond values only)"
                    if pn == 4
                    else "Appendix intrinsic/cross-coupling tables",
                    "bbox_pdf_points": [
                        xr[0],
                        round(y - 3.6, 2),
                        xr[1],
                        round(y + 3.6, 2),
                    ],
                    "raw_pdf_text": raw_text,
                    "ocr_corrections": corrections,
                    "numeric_decimals_reconstructed": found,
                },
                "numeric_record": {
                    "family": "literature_numeric_record",
                    "types": [row_id.replace("*", "wildcard")],
                    "values": printed,
                    "units": [u for u, v in zip(units, values) if v is not None],
                    "provenance": SOURCES[sid]["citation"] + f"; PDF p{pn}; {label}",
                    "citations": [
                        {
                            "file": "source-rows.tsv",
                            "line": len(text_rows),
                            "raw": numeric_line,
                            "columns": list(range(1, len(printed) + 1)),
                        }
                    ],
                },
                "execution_authorized": False,
            }
            row["evidence_checksum"] = digest(row)
            rows.append(row)
    assert len(rows) == 105
    text = ("\n".join(text_rows) + "\n").encode()
    (output / "source-rows.tsv").write_bytes(text)
    quotes = []
    for pn, needle in [
        (6, "parameters for the Si-C"),
        (7, "The carbon and hydrogen parameters"),
        (7, "two atom types"),
        (23, "(1995) 701"),
    ]:
        page = readers["sun1997"].pages[pn - 1].extract_text()
        start = page.index(needle)
        quotes.append(
            {
                "source": "sun1997",
                "page": pn,
                "raw_text": page[max(0, start - 70) : start + 430],
            }
        )
    payload = {
        "schema": "island_sun_siloxane_parameter_supplement_v1",
        "identity": "island_external_sun_siloxane_class2_data_v1",
        "classification": "external_literature_parameter_supplement; not an execution profile",
        "sources": SOURCES,
        "numeric_evidence": {
            "file": "source-rows.tsv",
            "sha256": sha256(text).hexdigest(),
        },
        "rows": rows,
        "inheritance_evidence": quotes,
        "inheritance_scope": "1997 Section 3.2 explicitly refers to 1995 [28c] for Si-C(methyl/phenyl) bonds and related angles/torsions; other transfers require separate justification. 1997 reoptimized Si/O nonbond values take precedence in its own model; 1995 nonbond records remain comparison-only.",
        "unresolved_policies": [
            "same-version hydrocarbon C/H base rows (1997 reference 38 is submitted work)",
            "blank 1995 angle k4 cells; no zero inference",
            "precise omitted cross-term applicability",
            "unlike-pair mixing and explicit 1-4 weights",
            "Si4/Si4c body versus broad Appendix type description",
            "chemical-state-specific inheritance of charges and remaining cross terms",
            "independent backend and numerical qualification",
        ],
        "source_rights": "user-authorized local extraction; ACS/Elsevier copyright notices retained; redistribution permission not established",
        "execution_profiles": [],
        "execution_authorized": False,
        "native_pcff_coverage": False,
        "msi2lmp_parity": False,
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    payload["checksum"] = digest(payload)
    write(output / "parameters.json", payload)
    print(
        json.dumps(
            {
                "rows": len(rows),
                "by_source": dict(Counter(r["source"] for r in rows)),
                "partial_rows": sum(bool(r["missing_fields"]) for r in rows),
                "catalog_sha256": sha256(
                    (output / "parameters.json").read_bytes()
                ).hexdigest(),
                "payload_identity": digest(payload),
            }
        )
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--references", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    extract(args.references, args.output)


if __name__ == "__main__":
    main()
