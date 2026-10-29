"""Explicit complete external XML sources; never an overlay on native OPLS."""

import json
import re
from dataclasses import dataclass
from decimal import Decimal
from functools import lru_cache
from hashlib import sha256
from math import isfinite
from pathlib import Path
from xml.etree import ElementTree as ET

from island.dynamics._checkpoint_data import strict_load
from island.exceptions import ForceFieldError
from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.pdms_literature import pinned_bytes
from island.workflows.storage import child

SCHEMA = "island_external_opls_source_v1"
UNITS = "openmm_xml_nm_kJ_mol_radian_elementary_charge_v1"
SOURCE_FAMILY = "external_opls_compatible"
FORMS = {
    "HarmonicBondForce": ("Bond", {"length", "k"}),
    "HarmonicAngleForce": ("Angle", {"angle", "k"}),
    "RBTorsionForce": ("Proper", {f"c{i}" for i in range(6)}),
    "OPLSTorsionForce": ("Proper", {f"k{i}" for i in range(1, 5)}),
}


class ExternalOPLSError(ForceFieldError):
    def __init__(self, status, diagnostics, *, result=None):
        self.status, self.diagnostics, self.result = status, diagnostics, result
        super().__init__(status + ": " + "; ".join(d["stage"] for d in diagnostics))


def decimal(value):
    require(type(value) is str, "source_decimal_string_required")
    v = Decimal(value)
    require(
        v.is_finite() and isfinite(float(v)) and (v == 0 or float(v) != 0),
        "nonfinite_or_unrepresentable_source_value",
    )
    return v


def inspect_xml(xml, document):
    # Keep source citations as exact file bytes/lines. No external XML entities.
    require(b"<!DOCTYPE" not in xml and b"<!ENTITY" not in xml, "xml_entities")
    root = ET.fromstring(xml)
    require(root.tag == "ForceField", "forcefield_root")
    require(
        set(root.attrib) == {"name", "version", "combining_rule"}
        and root.attrib["combining_rule"] == "geometric",
        "unsupported_mixing_or_metadata",
    )
    sections = [e.tag for e in root]
    require(len(sections) == len(set(sections)), "duplicate_force_section")
    require(
        set(sections) <= {"AtomTypes", "NonbondedForce", *FORMS},
        "unsupported_force_form_or_cross_terms",
    )
    require(
        {"AtomTypes", "NonbondedForce", "HarmonicBondForce", "HarmonicAngleForce"}
        <= set(sections),
        "missing_force_sections",
    )
    require(
        not ({"RBTorsionForce", "OPLSTorsionForce"} <= set(sections)),
        "mixed_torsion_conventions",
    )
    types = {}
    for row in root.find("AtomTypes"):
        a = row.attrib
        require(
            row.tag == "Type" and {"name", "class", "mass"} <= set(a),
            "atom_type_fields",
        )
        require(
            a["name"] not in types and ("def" not in a or a["def"].strip()),
            "duplicate_or_empty_type_rule",
        )
        require(decimal(a["mass"]) >= 0, "invalid_source_mass")
        if a.get("def"):
            require(
                a.get("element") and decimal(a["mass"]) > 0, "nonphysical_typable_site"
            )
        types[a["name"]] = dict(a)
    require(any(a.get("def") for a in types.values()), "missing_atom_typing_rules")
    nonbonded = root.find("NonbondedForce")
    require(
        set(nonbonded.attrib) == {"lj14scale", "coulomb14scale"},
        "missing_special_pair_policy",
    )
    for v in nonbonded.attrib.values():
        require(0 <= decimal(v) <= 1, "invalid_special_pair_scale")
    nb = {}
    for row in nonbonded:
        a = row.attrib
        require(
            row.tag == "Atom" and set(a) == {"type", "charge", "sigma", "epsilon"},
            "nonbonded_fields",
        )
        require(a["type"] in types and a["type"] not in nb, "nonbonded_type_conflict")
        decimal(a["charge"])
        require(decimal(a["sigma"]) >= 0 and decimal(a["epsilon"]) >= 0, "invalid_LJ")
        nb[a["type"]] = dict(a)
    require(set(nb) == set(types), "missing_nonbonded_or_charge_rows")
    for section, (tag, fields) in FORMS.items():
        for row in root.findall(section + "/" + tag):
            for name in fields:
                decimal(row.attrib[name])
        if root.find(section) is not None:
            require(
                all(e.tag == tag for e in root.find(section)),
                "unsupported_interaction_role",
            )
    evidence = _row_evidence(xml, document)

    return root, types, nb, evidence


def _row_evidence(xml, document):
    """Locate exact raw element bytes independently of ElementTree serialization."""
    from xml.parsers import expat

    parser = expat.ParserCreate()
    stack, counters, opened, evidence = [], {}, [], []
    row_tags = {"Type", "Atom", "Bond", "Angle", "Proper"}

    def start(tag, attrs):
        parent = stack[-1] if stack else None
        stack.append(tag)
        if tag not in row_tags:
            return
        key = f"{parent}/{tag}"
        index = counters.get(key, 0)
        counters[key] = index + 1
        opened.append(
            (
                tag,
                key,
                index,
                parser.CurrentByteIndex,
                parser.CurrentLineNumber,
                dict(attrs),
            )
        )

    def end(tag):
        stack.pop()
        if tag not in row_tags:
            return
        _, key, index, begin, line, attrs = opened.pop()
        finish = parser.CurrentByteIndex
        if not xml[begin:finish].rstrip().endswith(b"/>"):
            finish = xml.index(b">", finish) + 1
        item = {
            "source_record": f"{key}[{index}]",
            "source_file": document["xml_file"],
            "source_sha256": document["xml_sha256"],
            "line": line,
            "raw_text": xml[begin:finish].decode(),
            "source_attributes": attrs,
            "citation": document["provenance"]["citation"],
            "doi": document["provenance"]["doi"],
            "units": UNITS,
            "conversion_code": "external_opls_openmm_xml_v1",
        }
        evidence.append({**item, "evidence_checksum": digest(item)})

    parser.StartElementHandler, parser.EndElementHandler = start, end
    parser.Parse(xml, True)
    return evidence


@lru_cache(maxsize=8)
def _checked_xml(xml, document_text):
    # Cache only a successful assertion, never a mutable XML tree.
    inspect_xml(xml, strict_load(document_text))


@dataclass(frozen=True)
class ExternalOPLSSource:
    manifest_text: str
    xml: bytes
    root_path: Path
    manifest_file: str | None = None

    @property
    def document(self):
        return strict_load(self.manifest_text)

    @property
    def identity(self):
        return digest(
            {
                "manifest_sha256": sha256(self.manifest_text.encode()).hexdigest(),
                "xml_sha256": sha256(self.xml).hexdigest(),
            }
        )

    def validate_integrity(self):
        p = self.document
        if self.manifest_file is not None:
            require(
                pinned_bytes(
                    child(self.root_path, self.manifest_file),
                    sha256(self.manifest_text.encode()).hexdigest(),
                ).decode()
                == self.manifest_text,
                "stale_source_manifest",
            )
        require(
            p["schema"] == SCHEMA and p["source_family"] == SOURCE_FAMILY,
            "external_source_family_or_schema",
        )
        require(
            re.fullmatch(
                r"island_external_opls_[a-z0-9_]+_v[1-9][0-9]*", p["profile_id"] or ""
            ),
            "external_versioned_profile_required",
        )
        require(p["units"] == UNITS, "unresolved_or_wrong_units")
        require(p["limitations"] == LIMITS, "promoted_readiness")
        require(
            set(p)
            == {
                "schema",
                "source_family",
                "profile_id",
                "units",
                "xml_file",
                "xml_sha256",
                "provenance",
                "model",
                "limitations",
            },
            "unknown_source_fields_or_forged_authority",
        )
        require(
            pinned_bytes(child(self.root_path, p["xml_file"]), p["xml_sha256"])
            == self.xml,
            "stale_XML",
        )
        provenance = p["provenance"]
        require(
            set(provenance)
            == {
                "citation",
                "doi",
                "charge_method",
                "license_status",
                "license_evidence",
                "origin",
                "conversion",
            }
            and all(
                type(provenance[k]) is str and provenance[k].strip()
                for k in ("citation", "charge_method", "origin", "conversion")
            ),
            "source_provenance_required",
        )
        if provenance["license_status"] == "authorized_local_use":
            proof = provenance["license_evidence"]
            require(
                type(proof) is dict
                and set(proof) == {"file", "sha256", "line", "raw_text"},
                "license_row_evidence_required",
            )
            license_bytes = pinned_bytes(
                child(self.root_path, proof["file"]), proof["sha256"]
            )
            require(
                type(proof["line"]) is int
                and proof["line"] > 0
                and license_bytes.decode().splitlines()[proof["line"] - 1]
                == proof["raw_text"],
                "license_evidence_changed",
            )
        model = p["model"]
        require(
            set(model)
            == {
                "nonbonded_method",
                "periodic",
                "improper_policy",
                "chemical_state",
                "elements",
            },
            "model_policy_incomplete",
        )
        require(
            model["periodic"] is False
            and model["chemical_state"] == "neutral_closed_shell",
            "unsupported_model_scope",
        )
        require(
            type(model["elements"]) is list
            and model["elements"]
            and len(set(model["elements"])) == len(model["elements"]),
            "element_scope_required",
        )
        _checked_xml(self.xml, self.manifest_text)
        root = ET.fromstring(self.xml)
        require(
            set(model["elements"])
            <= {a.attrib.get("element") for a in root.find("AtomTypes")},
            "declared_element_without_source_types",
        )
        require(
            root.attrib["name"].upper() not in {"PCFF", "COMPASS", "GAFF", "GAFF2"},
            "cross_family_source",
        )

    def execution_diagnostics(self):
        p = self.document
        diagnostics = []
        if p["provenance"]["license_status"] != "authorized_local_use":
            diagnostics.append({"stage": "source_use_provenance_unestablished"})
        if p["model"]["nonbonded_method"] != "NoCutoff":
            diagnostics.append(
                {
                    "stage": "unsupported_nonbonded_method",
                    "requested": p["model"]["nonbonded_method"],
                    "supported": "finite_NoCutoff",
                }
            )
        if (
            p["model"]["improper_policy"]
            != "source_defines_no_impropers_for_declared_domain"
        ):
            diagnostics.append(
                {"stage": "improper_applicability_or_backend_incomplete"}
            )
        return diagnostics


def load_source(path, expected_sha256):
    """Only an explicit manifest path and an independent pin may select a source."""
    try:
        path = Path(path)
        raw = pinned_bytes(path, expected_sha256)
        p = strict_load(raw.decode())
        xml = pinned_bytes(child(path.parent, p["xml_file"]), p["xml_sha256"])
        result = ExternalOPLSSource(raw.decode(), xml, path.parent, path.name)
        result.validate_integrity()
        return result
    except FileNotFoundError as error:
        raise ExternalOPLSError(
            "asset_unavailable",
            [{"stage": "source_asset_missing", "reason": str(error)}],
        ) from error
    except ExternalOPLSError:
        raise
    except Exception as error:
        raise ExternalOPLSError(
            "source_incompatible",
            [{"stage": "source_validation_failed", "reason": str(error)}],
        ) from error


def write_source_manifest(xml, output, *, profile_id, provenance, model):
    """Publish an explicit local descriptor; this never confers scientific validity.

    The original XML must already be beside the requested manifest. No source is
    searched, downloaded, merged, edited or assigned a license by this helper.
    """
    xml, output = Path(xml), Path(output)
    require(
        not output.exists() and xml.parent.resolve() == output.parent.resolve(),
        "explicit_colocated_source_and_new_manifest_required",
    )
    p = {
        "schema": SCHEMA,
        "source_family": SOURCE_FAMILY,
        "profile_id": profile_id,
        "units": UNITS,
        "xml_file": xml.name,
        "xml_sha256": sha256(xml.read_bytes()).hexdigest(),
        "provenance": provenance,
        "model": model,
        "limitations": LIMITS,
    }
    raw = json.dumps(p, ensure_ascii=False, sort_keys=True, indent=2) + "\n"
    candidate = ExternalOPLSSource(raw, xml.read_bytes(), output.parent)
    candidate.validate_integrity()
    output.write_text(raw)
    return sha256(raw.encode()).hexdigest()
