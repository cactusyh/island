"""Code-authorized PCFF profiles; catalog presence never authorizes execution.

The registry is data-only. JSON cannot promote a candidate or change a domain.
Runtime selection uses one hash and final-graph perception, never CAR/MDF.
"""

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .automatic import chemical_graph
from .charges import identity
from .source import FLAGS, PIN, boundary, require

LINKED_BENZENOID = "island_pcff_linked_benzenoid_polymer_v1"


def _linked_definition():
    from .operational_profile import DEFINITION

    result = deepcopy(DEFINITION)
    result.update(
        name=LINKED_BENZENOID,
        authorized_atom_labels=["cp", "hc", "c3"],
        chemical_rule="connected_six_carbon_aromatic_cycles_single_links_methyl_ends_v1",
        unsupported_domains=[
            "fused or non-six-member aromatic rings",
            "heteroatoms",
            "unsaturation outside rings",
            "charged/radical/isotopic states",
            "non-methyl aliphatic environments",
            "periodic/virtual-site systems",
            "PCFF-IFF comparison candidates",
        ],
        evidence="J9 predeclared final PSMILES linked-phenylene source audit; independent source and numerical receipts",
        scope="Final graph only: repeat, sequence, branch and crosslink names are provenance, never dispatch rules",
    )
    return result


def _entries():
    from .operational_profile import DEFINITION

    entries = []
    for definition in (DEFINITION, _linked_definition()):
        entries.append(
            {
                "schema": "island_pcff_profile_registry_entry_v1",
                "name": definition["name"],
                "state": "operational",
                "definition": definition,
                "definition_identity": identity(definition),
                "validation_status": "bounded_implementation_fidelity; not scientific validation",
                "capabilities": {
                    k: True
                    for k in (
                        "preparation",
                        "native_charges",
                        "model_assembly",
                        "evaluation",
                        "bundles",
                        "workflows",
                    )
                },
                **FLAGS,
            }
        )
    for name, sha, state, reason in (
        (
            "lammps_pcff_full_catalog_v1",
            PIN["sha256"],
            "audited",
            "Full source has unresolved typing, increments and cross terms",
        ),
        (
            "iff_pcff_interface_v1_5_candidate",
            "3ad5a1be7334c646ed6cb813b769d0e89a1aa03fe695013e940626b7df922693",
            "candidate",
            "Comparison-only; no authorized native interpretation",
        ),
        (
            "lunar_pcff_interface_v1_6_candidate",
            "66a7798e8f7acd676c29c380ca050ad299a9b7c49c3fecf37821e67219c4e80b",
            "candidate",
            "Comparison-only; no authorized native interpretation",
        ),
        (
            "pcff_pyridinium_unresolved_v1",
            PIN["sha256"],
            "unsupported",
            "Missing nh+ native increments; no nh substitution",
        ),
        (
            "pcff_guanidinium_unresolved_v1",
            PIN["sha256"],
            "unsupported",
            "Printed native component total 0.9999 e; no charge repair",
        ),
    ):
        entries.append(
            {
                "schema": "island_pcff_profile_registry_entry_v1",
                "name": name,
                "source_sha256": sha,
                "state": state,
                "blocking_reason": reason,
                "capabilities": {
                    k: False
                    for k in (
                        "preparation",
                        "native_charges",
                        "model_assembly",
                        "evaluation",
                        "bundles",
                        "workflows",
                    )
                },
                **FLAGS,
            }
        )
    for entry in entries:
        if entry["state"] != "operational":
            candidate = entry["state"] == "candidate"
            provenance = {k: PIN[k] for k in ("repository", "commit", "source_path")}
            if entry["name"].startswith("iff_"):
                provenance = {
                    "repository": "https://github.com/hendrikheinz/INTERFACE-force-field-and-surface-models",
                    "commit": "584179265906d93aa40f7d5b275143985871b654",
                    "source_path": "pcff_interface_v1_5.frc",
                }
            elif entry["name"].startswith("lunar_"):
                provenance = {
                    "repository": "https://github.com/CMMRLab/LUNAR",
                    "commit": "67dabeda9e6bd3cc8f968aa0c6a88730886138ef",
                    "source_path": "frc_files/pcff_interface_v1_6mBN.frc",
                }
            entry["definition"] = {
                "schema": "island_pcff_nonauthorizing_profile_v1",
                "name": entry["name"],
                "source_sha256": entry["source_sha256"],
                "source_profile": "island_pcff_frc_candidate_v1"
                if candidate
                else PIN["profile"],
                "provenance": provenance,
                "typing_profile": None,
                "authorized_atom_labels": [],
                "authorized_charge_policy": None,
                "charge_tolerance_e": None,
                "authorized_source_families": [],
                "unsupported_families": ["all runtime families in this registration"],
                "chemical_rule": "no operational graph authorization",
                "physical_model": None,
                "validation_status": entry["blocking_reason"],
                **FLAGS,
            }
            entry["definition_identity"] = identity(entry["definition"])
    return entries


def list_pcff_profiles():
    """Owned registry entries, available without scientific imports or files."""
    return deepcopy(_entries())


@boundary
def inspect_pcff_profile(name):
    matches = [p for p in _entries() if p["name"] == name]
    require(len(matches) == 1, "Unknown PCFF registry profile")
    return deepcopy(matches[0])


def operational_definition(name):
    entry = inspect_pcff_profile(name)
    require(entry["state"] == "operational", "PCFF profile is not operational")
    return deepcopy(entry["definition"])


@boundary
def select_pcff_profile(name, *, sha256):
    from .operational_profile import PCFFOperationalSelection

    entry = inspect_pcff_profile(name)
    require(
        entry["state"] == "operational",
        "Candidate/audited/unsupported profile cannot authorize preparation",
    )
    require(
        entry["definition"]["source_sha256"] == sha256, "Profile source hash mismatch"
    )
    return PCFFOperationalSelection(name, sha256)


@boundary
def validate_pcff_profile_system(system, selection):
    from .operational_profile import PCFFOperationalSelection

    require(
        type(selection) is PCFFOperationalSelection,
        "Exact registered selection required",
    )
    selection.profile().validate_system(system)


@boundary
def validate_linked_benzenoid_graph(system):
    graph = chemical_graph(system)
    require(
        graph["representation"] == "atomistic" and not graph["has_box"],
        "Finite nonperiodic atomistic graph required",
    )
    atoms = {a["id"]: a for a in graph["sites"]}
    require(bool(atoms), "Empty graph")
    neighbors = {i: set() for i in atoms}
    aromatic = {i for i, a in atoms.items() if a["metadata"]["aromatic"]}
    require(bool(aromatic), "No authorized benzenoid ring")
    for b in graph["bonds"]:
        i, j = b["sites"]
        require(b["order"] == (1.5 if b["aromatic"] else 1), "Unsupported bond order")
        require(
            not b["aromatic"] or {i, j}.issubset(aromatic), "Inconsistent aromatic bond"
        )
        neighbors[i].add(j)
        neighbors[j].add(i)
    for i, a in atoms.items():
        m = a["metadata"]
        require(
            a["element"] in ("C", "H")
            and a["atomic_number"] == (6 if a["element"] == "C" else 1)
            and a["formal_charge"] == 0
            and m["isotope"] == m["radical_electrons"] == 0
            and m["chiral_tag"] == "CHI_UNSPECIFIED"
            and m["cip_label"] is None,
            "Unsupported element/charge/isotope/radical/stereo state",
        )
        require(
            abs(a["mass"] - (12.01115 if a["element"] == "C" else 1.00797)) <= 0.02,
            "Unsupported mass",
        )
        require(a["element"] == "C" or i not in aromatic, "Aromatic H is invalid")
        n = Counter(atoms[j]["element"] for j in neighbors[i])
        if a["element"] == "H":
            require(n == {"C": 1}, "Explicit hydrogen parent mismatch")
        elif i in aromatic:
            require(
                len(neighbors[i]) == 3 and n in ({"C": 2, "H": 1}, {"C": 3}),
                "Aromatic valence/H coverage mismatch",
            )
        else:
            require(
                n == {"C": 1, "H": 3} and any(j in aromatic for j in neighbors[i]),
                "Only methyl end/substituent sp3 carbon authorized",
            )
    ring_neighbors = {i: set() for i in aromatic}
    for b in graph["bonds"]:
        if b["aromatic"]:
            i, j = b["sites"]
            ring_neighbors[i].add(j)
            ring_neighbors[j].add(i)
    remaining = set(aromatic)
    ring_component = {}
    while remaining:
        todo = [min(remaining)]
        seen = set()
        while todo:
            i = todo.pop()
            if i not in seen:
                seen.add(i)
                todo.extend(ring_neighbors[i] - seen)
        require(
            len(seen) == 6 and all(len(ring_neighbors[i]) == 2 for i in seen),
            "Only isolated six-carbon aromatic cycles authorized",
        )
        for site in seen:
            ring_component[site] = min(seen)
        remaining -= seen
    for bond in graph["bonds"]:
        i, j = bond["sites"]
        if not bond["aromatic"] and i in aromatic and j in aromatic:
            require(
                ring_component[i] != ring_component[j],
                "Single links inside one aromatic ring are not authorized",
            )
    todo, seen = [min(atoms)], set()
    while todo:
        i = todo.pop()
        if i not in seen:
            seen.add(i)
            todo.extend(neighbors[i] - seen)
    require(seen == set(atoms), "Disconnected operational graph")


@dataclass(frozen=True)
class PCFFProfileRegistration:
    """Persisted registry inspection; trust is in code, never a mutable state flag."""

    json_text: str

    @boundary
    def validate_integrity(self):
        data = unpack(self.json_text)
        require(
            pack(data) == pack(inspect_pcff_profile(data.get("name"))),
            "Unregistered or contradictory PCFF registration",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def save_pcff_profile_registration(name, path):
    result = PCFFProfileRegistration(pack(inspect_pcff_profile(name)))
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())
    return result


@boundary
def load_pcff_profile_registration(path):
    result = PCFFProfileRegistration(Path(path).read_text())
    result.validate_integrity()
    return result
