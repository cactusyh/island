"""Pinned local Foyer OPLS-AA source; no network or scientific imports."""

import json
from dataclasses import dataclass
from functools import wraps
from hashlib import sha256
from math import isfinite
from pathlib import Path
from xml.etree import ElementTree as ET

from island.exceptions import OPLSAssignmentError

PIN = json.loads(Path(__file__).with_name("pin.json").read_text())
FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
CAPABILITIES = {
    "atom_typing": True,
    "native_charges": True,
    "complete_bonded_assignment": False,
    "energy_evaluation": False,
    "dynamics": False,
}


def require(ok, message):
    if not ok:
        raise OPLSAssignmentError(message)


def boundary(fn):
    @wraps(fn)
    def wrapped(*a, **kw):
        try:
            return fn(*a, **kw)
        except OPLSAssignmentError:
            raise
        except Exception as error:
            raise OPLSAssignmentError(f"{fn.__name__}: {error}") from error

    return wrapped


def digest(data):
    return sha256(data).hexdigest()


@dataclass(frozen=True)
class FoyerOPLSSource:
    """Owned exact XML bytes, verified on every access; source-specific identity."""

    xml: bytes

    @boundary
    def validate_integrity(self):
        require(
            type(self.xml) is bytes and digest(self.xml) == PIN["xml_sha256"],
            "OPLS source hash mismatch; install the pinned XML",
        )
        root = ET.fromstring(self.xml)
        require(
            root.attrib
            == {"name": "OPLS-AA", "version": "0.1.0", "combining_rule": "geometric"},
            "Source metadata mismatch",
        )
        require(
            not root.findall(".//Include"), "Includes are not supported by this pin"
        )

    @property
    def identity(self):
        self.validate_integrity()
        return {
            "family": "OPLS-AA",
            "variant": "Foyer-distributed OPLS-AA",
            "source_revision": PIN["revision"],
            "xml_sha256": PIN["xml_sha256"],
            "declared_version": "0.1.0",
            "typing_implementation": "pinned_foyer_v1",
            "charge_policy": "native_nonbonded_atom_charge_v1",
        }

    @property
    def entries(self):
        self.validate_integrity()
        root = ET.fromstring(self.xml)
        types = {r.attrib["name"]: dict(r.attrib) for r in root.find("AtomTypes")}
        charges = {}
        for r in root.find("NonbondedForce"):
            q = float(r.attrib["charge"])
            require(isfinite(q), "Nonfinite native charge")
            require(r.attrib["type"] not in charges, "Duplicate native charge")
            charges[r.attrib["type"]] = q
        return types, charges


@boundary
def load_oplsaa_source(path):
    source = FoyerOPLSSource(Path(path).read_bytes())
    source.validate_integrity()
    return source
