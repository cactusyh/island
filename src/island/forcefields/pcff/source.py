"""Local, section-aware FRC ingestion. No network or chemical backend imports."""

import json
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal
from functools import lru_cache, wraps
from hashlib import sha256
from math import isfinite
from pathlib import Path

from island.exceptions import PCFFError

PIN = json.loads(Path(__file__).with_name("pin.json").read_text())
FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
SEMANTIC = {"atom_types", "equivalence", "auto_equivalence", "bond_increments"}
EQUIVALENCE = ("nonbond", "bond", "angle", "torsion", "out_of_plane")
AUTO = (
    "nonbond",
    "bond_increment",
    "bond",
    "angle_end",
    "angle_apex",
    "torsion_end",
    "torsion_center",
    "out_of_plane_end",
    "out_of_plane_center",
)


def require(ok, message):
    if not ok:
        raise PCFFError(message)


def boundary(fn):
    @wraps(fn)
    def call(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except PCFFError:
            raise
        except Exception as error:
            raise PCFFError(f"{fn.__name__}: {error}") from error

    return call


def digest(raw):
    return sha256(raw).hexdigest()


def number(token):
    value = float(token)
    require(isfinite(value), "Nonfinite source number")
    return value


def _parse_uncached(raw):
    """Preserve every line, including unknown sections; decode four record families."""
    require(type(raw) is bytes, "Source must be immutable bytes")
    lines = raw.decode("utf-8").splitlines()
    require(
        lines and lines[0].strip() == "!BIOSYM forcefield 1", "Unsupported FRC header"
    )
    sections = []
    declarations = []
    section = None
    ended = False
    for index, line in enumerate(lines, 1):
        stripped = line.strip()
        if stripped.startswith("#"):
            require(not ended, f"Content after #end at line {index}")
            fields = stripped.split()
            name = fields[0][1:]
            if name == "version":
                require(len(fields) == 4, f"Malformed version declaration at {index}")
                declarations.append({"line": index, "raw": line, "fields": fields[1:]})
            section = {
                "name": name,
                "namespace": " ".join(fields[1:]),
                "line": index,
                "raw_header": line,
                "lines": [],
                "records": [],
                "interpretation": "semantic" if name in SEMANTIC else "raw_only",
            }
            sections.append(section)
            ended = name == "end"
            continue
        if section is not None:
            section["lines"].append({"line": index, "raw": line})
        if not stripped or stripped.startswith(("!", ">", "@")):
            continue
        require(not ended, f"Content after #end at line {index}")
        if section is None or section["name"] not in SEMANTIC:
            continue
        fields = stripped.split("!", 1)[0].split()
        name = section["name"]
        try:
            version = Decimal(fields[0])
            require(version.is_finite() and version >= 0, "Invalid record version")
            require(fields[1].isdigit(), "Invalid reference identifier")
            record = {
                "id": f"{name}:{section['namespace']}:{index}",
                "line": index,
                "raw": line,
                "section": name,
                "namespace": section["namespace"],
                "version": fields[0],
                "reference": fields[1],
            }
            if name == "atom_types":
                require(len(fields) >= 6, "Truncated atom type")
                mass = number(fields[3])
                require(mass > 0 and fields[5].isdigit(), "Invalid mass/connectivity")
                record["data"] = {
                    "type": fields[2],
                    "mass": mass,
                    "element": fields[4],
                    "connections": int(fields[5]),
                    "comment": " ".join(fields[6:]),
                    "mass_unit": "dalton",
                }
            elif name in ("equivalence", "auto_equivalence"):
                names = EQUIVALENCE if name == "equivalence" else AUTO
                require(
                    len(fields) == len(names) + 3, "Truncated/extra equivalence fields"
                )
                record["data"] = {
                    "type": fields[2],
                    "families": dict(zip(names, fields[3:])),
                }
            else:
                require(len(fields) == 6, "Truncated/extra bond increment fields")
                record["data"] = {
                    "types": fields[2:4],
                    "increments": [number(t) for t in fields[4:]],
                    "unit": "elementary_charge",
                }
            section["records"].append(record)
        except Exception as error:
            raise PCFFError(f"{name} line {index}: {error}") from error
    require(ended, "Truncated FRC: missing #end")
    for name in SEMANTIC:
        require(
            any(s["name"] == name and s["records"] for s in sections),
            f"Missing required section {name}",
        )
    return {"declarations": declarations, "sections": sections}


@lru_cache(maxsize=8)
def _parsed(raw, parser_context):
    return _parse_uncached(raw)


def _parser_context():
    return tuple(sorted(SEMANTIC)), EQUIVALENCE, AUTO, number, require, _parse_uncached


@boundary
def parse_frc(raw):
    require(type(raw) is bytes, "Source must be immutable bytes")
    return deepcopy(_parsed(raw, _parser_context()))


@dataclass(frozen=True)
class PCFFSource:
    """Immutable source bytes; unknown hashes may be inspected but not assigned."""

    raw: bytes
    expected_sha256: str

    @boundary
    def validate_integrity(self):
        require(
            type(self.expected_sha256) is str and len(self.expected_sha256) == 64,
            "Expected SHA256 required",
        )
        require(digest(self.raw) == self.expected_sha256, "PCFF source hash mismatch")
        _parsed(self.raw, _parser_context())

    @property
    def inventory(self):
        self.validate_integrity()
        return deepcopy(_parsed(self.raw, _parser_context()))

    @property
    def identity(self):
        self.validate_integrity()
        audited = self.expected_sha256 == PIN["sha256"]
        return {
            "family": "PCFF",
            "variant": "LAMMPS-distributed pcff.frc" if audited else "unreviewed FRC",
            "sha256": self.expected_sha256,
            "parser": PIN["parser"],
            "profile": PIN["profile"] if audited else None,
            "repository": PIN["repository"] if audited else None,
            "commit": PIN["commit"] if audited else None,
            "source_path": PIN["source_path"] if audited else None,
        }

    def require_assignment(self):
        self.validate_integrity()
        require(
            self.expected_sha256 == PIN["sha256"],
            "Source structurally parsed; charge semantics not audited for this hash",
        )


@boundary
def load_pcff_source(path, *, expected_sha256=None):
    source = PCFFSource(
        Path(path).read_bytes(),
        PIN["sha256"] if expected_sha256 is None else expected_sha256,
    )
    source.validate_integrity()
    return source


def records(inventory, name):
    namespace = (
        "cff91_auto" if name in ("auto_equivalence", "bond_increments") else "cff91"
    )
    return [
        r
        for s in inventory["sections"]
        if s["name"] == name and s["namespace"] == namespace
        for r in s["records"]
    ]


def select(candidates):
    """Highest decimal version; equal-version disagreements never use file order."""
    if not candidates:
        return None
    latest = max(Decimal(r["version"]) for r in candidates)
    top = [r for r in candidates if Decimal(r["version"]) == latest]
    require(
        all(r["data"] == top[0]["data"] for r in top),
        f"Conflicting source candidates: {[r['id'] for r in top]}",
    )
    return {
        "record": top[0],
        "candidate_ids": [r["id"] for r in candidates],
        "equivalent_selected_ids": [r["id"] for r in top],
        "policy": PIN["selection_policy"],
    }
