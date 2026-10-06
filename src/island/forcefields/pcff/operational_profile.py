"""Explicit operational authorization, separate from candidate source queries."""

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .automatic import chemical_graph
from .charges import identity
from .class2 import FAMILIES
from .source import FLAGS, PIN, boundary, require

SCHEMA = "island_pcff_operational_source_profile_v1"
NAME = "island_pcff_benzenoid_c6h6_v1"
DEFINITION = {
    "schema": SCHEMA,
    "name": NAME,
    "source_sha256": PIN["sha256"],
    "source_profile": PIN["profile"],
    "provenance": {k: PIN[k] for k in ("repository", "commit", "source_path")},
    "authorized_atom_labels": ["cp", "hc"],
    "authorized_charge_policy": "source_native_bond_increments_source_graph_v1",
    "charge_tolerance_e": 1e-12,
    "typing_profile": "island_pcff_source_graph_v1",
    "resolution_policy": None,
    "authorized_source_families": sorted(FAMILIES),
    "authorized_operations": {
        k: True for k in ("typing", "charges", "model_assembly", "evaluation")
    },
    "chemical_rule": "connected_neutral_explicit_standard_mass_C6H6_aromatic_cycle_v1",
    "unsupported_families": [
        "torsion-torsion_1",
        "nonzero_wilson_equilibrium",
        "quadratic_bond",
        "quadratic_angle",
        "torsion_1",
        "policy_derived_parameter_zero",
    ],
    "unsupported_domains": [
        "substituted/fused rings",
        "heteroaromatics",
        "charged/radical/isotopic states",
        "assigned stereochemistry",
        "acyclic compounds in this profile",
        "periodic/virtual-site systems",
        "PCFF-IFF comparison candidates",
    ],
    "physical_model": "island_lammps_pcff_source_graph_v1",
    "nonbonded": "9-6 minimum-distance, sixth-power mixing; explicit caller LJ/Coulomb shortest-path weights",
    "wilson": "source equilibrium zero only; signed LAMMPS Class II mean of three out-of-plane angles",
    "torsion_torsion": "excluded by named operational model scope, not an inferred source zero",
    "evidence": "Pinned FRC rows and pinned LAMMPS implementation; retained J1 nonplanar benzene whole-system evidence; J8 declared complete vertical slice",
    "scope": "Authorized bounded software model; all native integrity, graph, charge and coverage gates remain required",
    **FLAGS,
}


@dataclass(frozen=True)
class PCFFOperationalProfile:
    json_text: str

    @boundary
    def validate_integrity(self):
        p = unpack(self.json_text)
        require(
            pack(p) == pack(_definition(p.get("name"))),
            "Unregistered or contradictory PCFF operational profile",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @boundary
    def validate_system(self, system):
        self.validate_integrity()
        if self.payload["name"] != NAME:
            from .registry import FUSED_BENZENOID, validate_linked_benzenoid_graph

            validate_linked_benzenoid_graph(system, fused=self.payload["name"] == FUSED_BENZENOID)
            return
        graph = chemical_graph(system)
        require(
            graph["representation"] == "atomistic" and not graph["has_box"],
            "Operational profile requires finite nonperiodic atomistic system",
        )
        atoms = {s["id"]: s for s in graph["sites"]}
        require(
            Counter(a["element"] for a in atoms.values()) == {"C": 6, "H": 6},
            "Outside authorized C6H6 graph domain",
        )
        adjacency = {i: set() for i in atoms}
        for b in graph["bonds"]:
            i, j = b["sites"]
            aromatic = atoms[i]["element"] == atoms[j]["element"] == "C"
            require(
                b["aromatic"] == aromatic and b["order"] == (1.5 if aromatic else 1),
                "Unsupported aromatic/attachment bond environment",
            )
            adjacency[i].add(j)
            adjacency[j].add(i)
        for i, a in atoms.items():
            meta = a["metadata"]
            require(
                a["formal_charge"] == 0
                and meta["isotope"] == meta["radical_electrons"] == 0
                and meta["chiral_tag"] == "CHI_UNSPECIFIED"
                and meta["cip_label"] is None,
                "Unsupported charge/isotope/radical/stereochemistry",
            )
            require(
                meta["aromatic"] == (a["element"] == "C"),
                "Aromatic atom metadata mismatch",
            )
            require(
                abs(a["mass"] - (12.01115 if a["element"] == "C" else 1.00797)) <= 0.02,
                "Unsupported mass inventory",
            )
            require(
                Counter(atoms[j]["element"] for j in adjacency[i])
                == ({"C": 2, "H": 1} if a["element"] == "C" else {"C": 1}),
                "Unsupported explicit-H/ring connectivity",
            )
        todo = [next(iter(atoms))]
        seen = set()
        while todo:
            i = todo.pop()
            if i not in seen:
                seen.add(i)
                todo.extend(adjacency[i] - seen)
        require(seen == set(atoms), "Disconnected graph not authorized")

    @boundary
    def validate_native(self, system, native):
        from .model import PCFFModelSpecification

        require(
            type(native) is PCFFModelSpecification,
            "Existing native PCFF model required",
        )
        self.validate_system(system)
        native.validate_integrity(system)
        from .validation_cache import profile_validation_summary

        summary = profile_validation_summary(native)
        p = summary["model"]
        profile = self.payload
        require(
            p["source"]["sha256"] == profile["source_sha256"]
            and p["source"]["profile"] == profile["source_profile"],
            "Operational source hash/profile mismatch",
        )
        require(
            p["schema"] == "island_pcff_source_model_v1"
            and p["compatibility_profile"]["name"] == profile["physical_model"],
            "Unauthorized model interpretation profile",
        )
        require(
            p["model_definition_complete"] and not p["diagnostics"],
            "Incomplete operational model",
        )
        assignment = summary["assignment"]
        auto = assignment["charge_record"]["automatic_typing"]
        require(
            auto["profile"]["name"] == profile["typing_profile"]
            and (
                set(auto["assignments"].values())
                == set(profile["authorized_atom_labels"])
                if profile["name"] == NAME
                else set(auto["assignments"].values()).issubset(
                    profile["authorized_atom_labels"]
                )
            ),
            "Unauthorized typing/labels",
        )
        for term in p["terms"]:
            require(
                term["family"] in profile["authorized_source_families"]
                and term["origin"] == "source_row"
                and bool(term["source_rows"]),
                "Policy zero or unauthorized family cannot enter operational model",
            )
            if term["family"] == "wilson_out_of_plane":
                require(
                    term["coefficients"][1] == 0,
                    "Nonzero Wilson equilibrium not authorized",
                )
        if profile["name"] != NAME:
            require(
                assignment["parameter_coverage_complete"], "Missing source parameters"
            )
            return
        required = {
            "quartic_bond": 12,
            "quartic_angle": 18,
            "torsion_3": 24,
            "wilson_out_of_plane": 6,
            "angle-angle": 18,
        }
        counts = Counter(t["family"] for t in p["terms"])
        require(
            all(counts[k] == n for k, n in required.items()),
            "Incomplete aromatic interaction multiplicity",
        )
        for family in ("bond-bond", "bond-angle"):
            require(counts[family] == 18, "Required adjacent coupling missing")
        for family in (
            "end_bond-torsion_3",
            "middle_bond-torsion_3",
            "angle-torsion_3",
            "angle-angle-torsion_1",
            "bond-bond_1_3",
        ):
            require(counts[family] == 24, "Required torsional coupling missing")


@dataclass(frozen=True)
class PCFFOperationalSelection:
    name: str
    sha256: str

    def __post_init__(self):
        require(
            type(self.name) is str
            and self.name in _operational_names()
            and type(self.sha256) is str
            and self.sha256 == _definition(self.name)["source_sha256"],
            "Unregistered PCFF operational hash/profile selection",
        )

    def profile(self):
        self.__post_init__()
        result = PCFFOperationalProfile(pack(deepcopy(_definition(self.name))))
        result.validate_integrity()
        return result


@boundary
def load_pcff_operational_profile(path):
    profile = PCFFOperationalProfile(Path(path).read_text())
    profile.validate_integrity()
    return profile


@boundary
def save_pcff_operational_profile(profile, path):
    require(type(profile) is PCFFOperationalProfile, "Expected operational profile")
    profile.validate_integrity()
    publish(Path(path), profile.json_text.encode())


def _definition(name):
    # J8's exact signed payload is preserved. New authorizations are code-owned.
    if name == NAME:
        return DEFINITION
    from .registry import operational_definition

    return operational_definition(name)


def _operational_names():
    from .registry import list_pcff_profiles

    return {p["name"] for p in list_pcff_profiles() if p["state"] == "operational"}
