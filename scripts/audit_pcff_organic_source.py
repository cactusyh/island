"""Independent raw source audit; external FRC is never modified."""

import argparse
import hashlib
import json
from decimal import Decimal
from pathlib import Path

parser = argparse.ArgumentParser(
    description="Independent J4 raw-source counts and exact charge-row audit"
)
parser.add_argument("--source", type=Path, required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
raw = args.source.read_bytes()
assert (
    hashlib.sha256(raw).hexdigest()
    == "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c"
)

sections = {}
sec = None
for n, line in enumerate(raw.decode().splitlines(), 1):
    w = line.split("!", 1)[0].split()
    if not w:
        continue
    if w[0].startswith("#"):
        sec = w[0][1:]
        sections.setdefault(sec, [])
    elif sec and w[0][0].isdigit():
        sections[sec].append((n, w, line))
labels = {"h", "s'", "nh+", "c+", "nr", "n3n", "n4n"}
selected = {
    k: [{"line": n, "raw": line} for n, w, line in sections[k] if w[2] in labels]
    for k in ("atom_types", "equivalence", "auto_equivalence")
}
row = next(
    (n, w, line) for n, w, line in sections["bond_increments"] if w[2:4] == ["c+", "nr"]
)
value = 3 * (Decimal(row[1][4]) + Decimal(row[1][5]))
out = {
    "source_sha256": hashlib.sha256(raw).hexdigest(),
    "counts": {
        k: len(sections[k])
        for k in ("atom_types", "equivalence", "auto_equivalence", "bond_increments")
    },
    "distinct_atom_types": len({w[2] for n, w, line in sections["atom_types"]}),
    "selected_source_rows": selected,
    "guanidinium_independent_conservation": {
        "row": row[0],
        "raw": row[2],
        "bonds": 3,
        "other_bonds": "nr-H increments are equal-and-opposite and do not change molecular sum",
        "total_e": str(value),
        "formal_e": 1,
        "residual_e": str(value - 1),
        "status": "failed unchanged 1e-12 acceptance",
    },
    "source_domains_unresolved": {
        "sulfur_oxidation": "s/s1/s2/sh/sp/sio? inspect exact sulfur rows; all sulfur atom rows have degree 1 or 2; no blanket degree3/4 sulfur label",
        "sf": "source description sulfonate but connections1; cannot assign central tetracoordinate sulfur",
        "p=": "source phosphazene label connections5; common PCl5 degree5 alone is not phosphazene evidence",
        "imidazolium": "ci/ni/hi charge-equivalence and ring-H convention need complete independent audit",
        "aliases": "n+/n1/n2/nb etc remain ambiguous; no name substitution",
    },
    "production_validated": False,
    "simulation_readiness": "not_established",
}
out["source_domains_unresolved"]["sulfur_oxidation"] = (
    "All source S atom rows have connections1 or2; no declared support for sulfoxide/sulfone central sulfur coordination3/4."
)
with args.output.open("x") as file:
    file.write(json.dumps(out, indent=2) + "\n")
print(
    out["counts"],
    out["distinct_atom_types"],
    out["guanidinium_independent_conservation"],
)
