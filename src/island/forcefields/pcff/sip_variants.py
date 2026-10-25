"""Authorized, single-file Si/P source inventories. Audits confer no model authority.

COMPASS's published FRC omits the PCFF atom connectivity column and automatic
equivalences. Preserve that distinction instead of adapting it into pinned PCFF.
"""

import re
from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path

from island.charge_references.records import pack, unpack

from .charges import identity
from .class2 import FAMILIES
from .source import FLAGS, boundary, digest, require

SCHEMA = "island.pcff-sip-source-variant-c3.v1"
ARITIES = {f: v[0] for f, v in FAMILIES.items()}
ARITIES.update(
    {
        "quadratic_bond": 2,
        "quadratic_angle": 3,
        "torsion_1": 4,
        "bond_increments": 2,
        "atom_types": 1,
        "equivalence": 1,
        "auto_equivalence": 1,
    }
)
PROVENANCE = {
    "family",
    "version",
    "repository",
    "revision",
    "source_path",
    "authorization",
    "authorization_basis",
}


def scan(raw):
    """Lossless section/line inventory; no automatic chemical type assignment."""
    require(type(raw) is bytes, "Immutable source bytes required")
    lines = raw.decode("ascii").splitlines()
    require(
        lines and lines[0].split() == ["!BIOSYM", "forcefield", "1"],
        "Unsupported FRC header",
    )
    sections, rows, declarations = [], [], []
    current = None
    ended = False
    for line, text in enumerate(lines, 1):
        words = text.split("!", 1)[0].split()
        if text.lstrip().startswith("#"):
            require(not ended, "Content after FRC end")
            words = text.split()
            name = words[0][1:]
            current = {
                "name": name,
                "namespace": " ".join(words[1:]),
                "line": line,
                "raw": text,
                "row_ids": [],
            }
            sections.append(current)
            if name in ("version", "define"):
                declarations.append(current.copy())
            ended = name == "end"
        elif words and re.fullmatch(r"\d+(?:\.\d+)?", words[0]):
            require(current is not None and not ended, "Row outside FRC section")
            name, ns = current["name"], current["namespace"]
            arity = ARITIES.get(name, 0)
            require(len(words) >= 2 + arity, "Truncated FRC row")
            row = {
                "id": f"{name}:{ns}:{line}",
                "family": name,
                "namespace": ns,
                "line": line,
                "raw": text,
                "fields": words,
                "types": words[2 : 2 + arity],
            }
            if name == "atom_types":
                require(len(words) >= 5, "Truncated atom declaration")
                row.update(
                    element=words[4],
                    mass=words[3],
                    connections=int(words[5])
                    if words[5:6] and words[5].isdigit()
                    else None,
                )
            elif name in FAMILIES or name in (
                "quadratic_bond",
                "quadratic_angle",
                "torsion_1",
                "bond_increments",
            ):
                require(
                    all(Decimal(t).is_finite() for t in words[2 + arity :]),
                    "Nonfinite source coefficient",
                )
                row["values"] = words[2 + arity :]
            rows.append(row)
            current["row_ids"].append(row["id"])
    require(any(r["family"] == "atom_types" for r in rows), "No atom types")
    return {
        "sections": sections,
        "rows": rows,
        "declarations": declarations,
        "end_marker_present": ended,
    }


def template_blocks(raw):
    text = raw.decode("ascii")
    return [
        {"type": m[1], "line": text[: m.start()].count("\n") + 1, "raw": m[0]}
        for m in re.finditer(
            r"(?m)^type:\s*([^\s]+)[ \t]*\n[\s\S]*?^end_type[ \t]*$", text
        )
    ]


@dataclass(frozen=True)
class SiPSourceVariant:
    """One complete FRC and its own explicitly declared companion files."""

    filename: str
    raw: bytes
    expected_sha256: str
    provenance_json: str
    companions: tuple = ()

    @boundary
    def validate_integrity(self):
        require(Path(self.filename).name == self.filename, "Invalid source filename")
        require(
            digest(self.raw) == self.expected_sha256, "Variant source hash mismatch"
        )
        p = unpack(self.provenance_json)
        require(
            set(p) == PROVENANCE
            and all(type(v) is str and v.strip() for v in p.values()),
            "Complete source provenance and authorization required",
        )
        require(p["family"] in ("PCFF", "COMPASS"), "Unsupported source family")
        require(
            p["authorization"] in ("authorized_local", "explicitly_supplied"),
            "Source is not authorized",
        )
        inv = scan(self.raw)
        namespace = "cff91" if p["family"] == "PCFF" else "compass"
        require(
            any(
                s["name"] == "define" and s["namespace"].lower() == namespace
                for s in inv["sections"]
            ),
            "Source family contradicts FRC declaration",
        )
        names = {self.filename}
        for name, raw, checksum in self.companions:
            require(
                Path(name).name == name and name not in names,
                "Duplicate or invalid companion",
            )
            require(
                type(raw) is bytes and digest(raw) == checksum,
                "Companion source hash mismatch",
            )
            names.add(name)

    @property
    def payload(self):
        self.validate_integrity()
        inv = scan(self.raw)
        atoms = [r for r in inv["rows"] if r["family"] == "atom_types"]
        labels = {r["types"][0] for r in atoms if r["element"] in ("Si", "P")}
        aliases = set(labels)
        for r in inv["rows"]:
            if (
                r["family"] in ("equivalence", "auto_equivalence")
                and r["types"][0] in labels
            ):
                aliases.update(r["fields"][3:])
        selected = [
            r["id"]
            for r in inv["rows"]
            if set(r["types"]) & aliases or any(t.startswith("*") for t in r["types"])
        ]
        return {
            "schema": SCHEMA,
            "filename": self.filename,
            "sha256": self.expected_sha256,
            "provenance": unpack(self.provenance_json),
            **inv,
            "atom_types": atoms,
            "sip_labels": sorted(labels),
            "sip_aliases": sorted(aliases),
            "sip_row_ids": selected,
            "family_counts": dict(
                sorted(Counter(r["family"] for r in inv["rows"]).items())
            ),
            "sip_families": {
                f: [
                    r["id"]
                    for r in inv["rows"]
                    if r["family"] == f and r["id"] in selected
                ]
                for f in sorted(ARITIES)
            },
            "templates": {
                self.filename: template_blocks(self.raw),
                **{n: template_blocks(raw) for n, raw, _ in self.companions},
            },
            "companion_hashes": {n: h for n, _, h in self.companions},
            "execution_authorized": False,
            "policy": "single complete source; inventory is not a typing or execution profile; no cross-file row mixing",
            **FLAGS,
        }

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def load_sip_source_variant(path, *, expected_sha256, provenance, companions=()):
    v = SiPSourceVariant(
        Path(path).name,
        Path(path).read_bytes(),
        expected_sha256,
        pack(provenance),
        tuple((Path(p).name, Path(p).read_bytes(), h) for p, h in companions),
    )
    v.validate_integrity()
    return v


@boundary
def validate_sip_variant_inventory(path, variant):
    require(
        pack(unpack(Path(path).read_text())) == pack(variant.payload),
        "Contradictory variant rows, provenance, source or execution authority",
    )
    return variant.identity
