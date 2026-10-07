"""Data-only local distribution evidence; never an operational typing authority.

Only the envelope grammar of the supplied 1995 template dialect is interpreted.
Pattern and precedence trees are retained as text/structure, not executed. No
runtime preparation depends on this module or searches for companion files.
"""

import re
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows import storage

from .catalog import catalog_from_inventory
from .source import FLAGS, PIN, boundary, digest, parse_frc, records, require

SCHEMA = "island_pcff_local_distribution_audit_v1"
ROLES = ("frc", "rlb", "templates")


def _ascii(raw):
    require(type(raw) is bytes and b"\x00" not in raw, "Invalid auxiliary bytes")
    return raw.decode("ascii")


def _tree(text):
    """Parse parentheses without assigning chemical override semantics."""
    tokens = re.findall(r"\(|\)|[^\s()]+", text)
    root, stack = [], []
    for token in tokens:
        if token == "(":
            branch = []
            (stack[-1] if stack else root).append(branch)
            stack.append(branch)
        elif token == ")":
            require(bool(stack), "Unbalanced template precedence")
            stack.pop()
        else:
            require(bool(stack), "Bare precedence token")
            stack[-1].append(token)
    require(not stack, "Truncated template precedence")
    return root


@boundary
def inspect_pcff_templates(raw):
    """Preserve duplicate records, syntax variants and uninterpreted predicates."""
    text = _ascii(raw)
    rows, current, test, precedence = [], None, None, None
    diagnostics, header, body = [], [], []
    ended = False
    for line, original in enumerate(text.splitlines(), 1):
        s = original.split("!", 1)[0].strip()
        if not s:
            if current is None and precedence is None:
                header.append({"line": line, "text": original})
            continue
        if s == "precedence:":
            require(
                current is None and precedence is None and not ended,
                f"Unexpected precedence at {line}",
            )
            precedence = {"line": line}
        elif s == "end_precedence":
            require(precedence is not None and not ended, f"Unexpected end at {line}")
            precedence.update(end_line=line, tree=_tree("\n".join(body)))
            ended = True
        elif precedence is not None:
            require(not ended, f"Content after precedence at {line}")
            body.append(s)
        elif re.fullmatch(r"type(?::|\s)\s*\S+", s):
            require(current is None and test is None, f"Unclosed type at {line}")
            label = re.sub(r"^type(?::|\s)\s*", "", s)
            if not s.startswith("type:"):
                diagnostics.append(
                    {"line": line, "reason": "non_colon_type_syntax_retained"}
                )
            require(
                bool(label) and not any(x.isspace() for x in label),
                f"Malformed type at {line}",
            )
            current = {
                "id": f"template:{line}:{label}",
                "line": line,
                "label": label,
                "patterns": [],
                "tests": [],
                "unsupported": [],
                "interpretation": "structural_only_not_authorized_for_assignment",
            }
        elif s == "end_type":
            require(
                current is not None and test is None, f"Unexpected end_type at {line}"
            )
            current["end_line"] = line
            rows.append(current)
            current = None
        elif current is not None:
            m = re.match(r"template\s*(:?)\s*(.*)", s)
            if m:
                require(test is None and bool(m[2]), f"Malformed template at {line}")
                current["patterns"].append(
                    {"line": line, "text": m[2], "colon": bool(m[1])}
                )
            elif s.startswith("atom_test:"):
                require(test is None, f"Nested atom_test at {line}")
                index = s.split(":", 1)[1].strip()
                require(
                    index.isdigit() and int(index) > 0, f"Invalid test index at {line}"
                )
                test = {"line": line, "index": int(index), "conditions": []}
            elif s == "end_test":
                require(test is not None, f"Unexpected end_test at {line}")
                test["end_line"] = line
                current["tests"].append(test)
                test = None
            elif test is not None and ":" in s:
                key, value = s.split(":", 1)
                test["conditions"].append(
                    {"line": line, "key": key.strip(), "value": value.strip()}
                )
            else:
                current["unsupported"].append({"line": line, "text": original})
                diagnostics.append({"line": line, "reason": "uninterpreted_record"})
        else:
            diagnostics.append(
                {"line": line, "reason": "uninterpreted_top_level", "text": original}
            )
    require(
        current is None and test is None and ended, "Truncated template distribution"
    )
    require(bool(rows), "Empty template distribution")
    counts = Counter(r["label"] for r in rows)
    return {
        "schema": "island_pcff_template_envelope_v1",
        "sha256": digest(raw),
        "header": header,
        "records": rows,
        "precedence": precedence,
        "counts": {
            "records": len(rows),
            "labels": len(counts),
            "patterns": sum(len(r["patterns"]) for r in rows),
            "tests": sum(len(r["tests"]) for r in rows),
        },
        "duplicate_labels": {k: v for k, v in sorted(counts.items()) if v > 1},
        "diagnostics": diagnostics,
        "unsupported_semantics": [
            "pattern token and bracket semantics",
            "atom numbering inside patterns",
            "planar/nonplanar geometry predicates",
            "override/conflict and repeated precedence nodes",
            "formal-charge/isotope/radical completeness",
            "numerical parameter availability",
        ],
        "operational_authority": False,
    }


def _payload(raw, provenance):
    require(
        type(provenance) is dict and set(provenance) == set(ROLES),
        "Explicit file provenance required",
    )
    require(
        all(type(provenance[r]) is str and provenance[r] for r in ROLES),
        "Invalid file provenance",
    )
    inv = parse_frc(raw["frc"])
    catalog = catalog_from_inventory(inv, {"sha256": digest(raw["frc"])})
    templates = inspect_pcff_templates(raw["templates"])
    rlb = _ascii(raw["rlb"]).splitlines()
    atoms = records(inv, "atom_types")
    labels = {r["data"]["type"] for r in atoms}
    tlabels = {r["label"] for r in templates["records"]}
    return {
        "schema": SCHEMA,
        "files": {
            r: {
                "sha256": digest(raw[r]),
                "bytes": len(raw[r]),
                "encoding": "ASCII" if raw[r].isascii() else "UTF-8",
                "historical_resolved_path": provenance[r],
            }
            for r in ROLES
        },
        "frc": {
            "identical_to_historical_pin": digest(raw["frc"]) == PIN["sha256"],
            "declarations": inv["declarations"],
            "atom_labels": sorted(labels),
            "semantic_counts": {
                k: len(records(inv, k))
                for k in (
                    "atom_types",
                    "equivalence",
                    "auto_equivalence",
                    "bond_increments",
                )
            },
            "sections": [
                {
                    "name": s["name"],
                    "namespace": s["namespace"],
                    "line": s["line"],
                    "raw_sha256": digest(storage.json_bytes(s)),
                    "parsed_records": len(s["records"]),
                    "numerical_record_ids": [
                        r["id"]
                        for r in catalog["records"]
                        if r["section"] == s["name"]
                        and r["namespace"] == s["namespace"]
                    ],
                }
                for s in inv["sections"]
            ],
        },
        "rlb": {
            "lines": rlb,
            "interpretation": "VERSION/elib marker only; no chemical or numerical records"
            if rlb == ["VERSION", "elib"]
            else "unsupported_content",
            "operational_authority": False,
        },
        "templates": templates,
        "links": {
            "template_labels_without_frc_atom": sorted(tlabels - labels),
            "frc_labels_without_template": sorted(labels - tlabels),
            "by_source_label": {
                label: [r["id"] for r in templates["records"] if r["label"] == label]
                for label in sorted(labels)
            },
        },
        "runtime_authorization": False,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFDistributionAudit:
    """Owned inspection record, with no conversion to a model or typing profile."""

    json_text: str

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return digest(self.json_text.encode())

    @boundary
    def validate_integrity(self, *, frc_path, rlb_path, templates_path):
        p = self.payload
        require(p.get("schema") == SCHEMA, "Unsupported distribution audit schema")
        paths = dict(zip(ROLES, (frc_path, rlb_path, templates_path), strict=True))
        raw = {r: Path(paths[r]).read_bytes() for r in ROLES}
        require(
            all(digest(raw[r]) == p["files"][r]["sha256"] for r in ROLES),
            "Distribution dependency hash mismatch",
        )
        expected = _payload(
            raw, {r: p["files"][r]["historical_resolved_path"] for r in ROLES}
        )
        require(p == expected, "Contradictory distribution audit")


@boundary
def inspect_pcff_distribution(*, frc_path, rlb_path, templates_path, expected_hashes):
    """Read three explicit local files; no inferred paths or runtime promotion."""
    require(
        type(expected_hashes) is dict and set(expected_hashes) == set(ROLES),
        "Three explicit hashes required",
    )
    paths = dict(zip(ROLES, (frc_path, rlb_path, templates_path), strict=True))
    raw = {r: Path(paths[r]).read_bytes() for r in ROLES}
    require(
        all(digest(raw[r]) == expected_hashes[r] for r in ROLES),
        "Distribution dependency hash mismatch",
    )
    return PCFFDistributionAudit(
        pack(_payload(raw, {r: str(Path(paths[r]).resolve()) for r in ROLES}))
    )


@boundary
def save_pcff_distribution_audit(path, result, *, frc_path, rlb_path, templates_path):
    require(isinstance(result, PCFFDistributionAudit), "Distribution audit required")
    result.validate_integrity(
        frc_path=frc_path, rlb_path=rlb_path, templates_path=templates_path
    )
    storage.publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_distribution_audit(path, *, frc_path, rlb_path, templates_path):
    result = PCFFDistributionAudit(Path(path).read_text())
    result.validate_integrity(
        frc_path=frc_path, rlb_path=rlb_path, templates_path=templates_path
    )
    return result
