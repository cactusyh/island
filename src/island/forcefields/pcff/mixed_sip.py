"""C1 opt-in mixed Si/P assessments with pinned companion typing evidence.

An inventory or template match is not execution authority.  This tranche never
substitutes a wildcard or a converter-generated zero for a required Si/P row.
The historical J41 finite-silane model and all historical rules are unchanged.
"""

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack, unpack

from .automatic import chemical_graph, type_pcff_atoms
from .charges import identity
from .expanded import environments
from .fallbacks import MSI_POLICY
from .operational import inspect_pcff_operational_support
from .silicon_phosphorus import SPECIAL, assess_pcff_sip, prepare_pcff_sip
from .source import FLAGS, PIN, boundary, load_pcff_source, require

DISTRIBUTION_PINS = {
    "pcff.frc": PIN["sha256"],
    "pcff.rlb": "82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36",
    "pcff_templates.dat": "dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093",
}
PROFILES = {
    "island_pcff_mixed_silane_c1_v1": "mixed_silane",
    "island_pcff_siloxane_c1_v1": "siloxane",
    "island_pcff_organosilicon_polymer_c1_v1": "organosilicon_polymer",
    "island_pcff_phosphorus_audit_c1_v1": "phosphorus_audit",
    "island_pcff_mixed_sip_audit_c1_v1": "mixed_sip_audit",
}
POLICY = {
    "version": "island_pcff_mixed_sip_completeness_v1",
    "resolution_policy": MSI_POLICY,
    "si_p_rows": "exact_or_source_declared_ordinary_equivalence; no wildcard",
    "missing_cross_term": "reject; never substitute zero",
    "p_typing": "no execution while source connections/template conflict persists",
    "charge_origin": "PCFF/source_native_bond_increments",
    "periodic": "not_authorized",
}


@dataclass(frozen=True)
class PCFFSiPDistribution:
    source: object
    rlb: bytes
    templates: bytes

    @boundary
    def validate_integrity(self):
        self.source.require_assignment()
        require(self.hashes == DISTRIBUTION_PINS, "Si/P distribution hash mismatch")

    @property
    def hashes(self):
        return {
            "pcff.frc": sha256(self.source.raw).hexdigest(),
            "pcff.rlb": sha256(self.rlb).hexdigest(),
            "pcff_templates.dat": sha256(self.templates).hexdigest(),
        }


@boundary
def load_pcff_sip_distribution(directory):
    directory = Path(directory)
    result = PCFFSiPDistribution(
        load_pcff_source(directory / "pcff.frc"),
        (directory / "pcff.rlb").read_bytes(),
        (directory / "pcff_templates.dat").read_bytes(),
    )
    result.validate_integrity()
    return result


def phosphorus_environments(graph):
    """Describe chemical state without guessing a source type or altering a graph."""
    atoms, adj, bonds, env = environments(graph)
    result = []
    for i, atom in atoms.items():
        if atom["element"] != "P":
            continue
        e = env[i]
        orders = e["bond_orders"]
        doubles = [
            j
            for j in adj[i]
            if next(b["order"] for b in graph["bonds"] if set(b["sites"]) == {i, j})
            == 2
        ]
        if orders == [1, 1, 1]:
            kind = "phosphine_like_trivalent"
        elif orders == [1, 1, 1, 2] and len(doubles) == 1:
            if atoms[doubles[0]]["element"] == "O":
                kind = (
                    "phosphate_ester_pentavalent"
                    if e["neighbor_elements"] == {"O": 4}
                    else "phosphonate_pentavalent"
                    if e["neighbor_elements"] == {"C": 1, "O": 3}
                    else "other_oxo_pentavalent"
                )
            else:
                kind = "non_oxo_pentavalent"
        else:
            kind = "unsupported_connectivity_or_valence"
        # This conventional formal oxidation-number description is diagnostic,
        # never an atom-type rule. P-H assigns H +I; P-C/O/N/S/halogen assigns
        # the bonding electrons away from P. Unknown ligands stay unresolved.
        signs = {
            "H": -1,
            "Si": -1,
            "P": 0,
            "C": 1,
            "O": 1,
            "N": 1,
            "S": 1,
            "F": 1,
            "Cl": 1,
            "Br": 1,
            "I": 1,
        }
        oxidation = (
            atom["formal_charge"]
            + sum(
                signs[atoms[j]["element"]] * bonds[frozenset((i, j))]["order"]
                for j in adj[i]
            )
            if all(atoms[j]["element"] in signs for j in adj[i])
            and not any(bonds[frozenset((i, j))]["aromatic"] for j in adj[i])
            else None
        )
        reason = (
            "charged_phosphorus_unsupported"
            if atom["formal_charge"] or any(atoms[j]["formal_charge"] for j in adj[i])
            else "phosphine_no_degree_three_source_type"
            if kind == "phosphine_like_trivalent"
            else "p_oxo_template_connections_conflict"
            if kind
            in (
                "phosphate_ester_pentavalent",
                "phosphonate_pentavalent",
                "other_oxo_pentavalent",
                "non_oxo_pentavalent",
            )
            else "unsupported_phosphorus_oxidation_or_connectivity"
        )
        result.append(
            {
                "site": i,
                "environment": e,
                "classification": kind,
                "formal_charge": atom["formal_charge"],
                "bond_order_sum": sum(orders),
                "formal_oxidation_number": oxidation,
                "oxidation_convention": "P-H: H +I; P-Si: Si positive; P-C/O/N/S/halogen: P positive; P-P: zero; description only",
                "source_supported_oxidation_state": "not_established",
                "neighbors": [
                    {
                        "id": j,
                        "element": atoms[j]["element"],
                        "formal_charge": atoms[j]["formal_charge"],
                    }
                    for j in sorted(adj[i])
                ],
                "selected_type": None,
                "execution_authorized": False,
                "reason": reason,
                "evidence": {
                    "p": "pcff.frc:163 general P connections=4; template:1557 general P",
                    "p=": "pcff.frc:164 phosphazene P connections=5; template:1562 oxo P has 4 neighbors and bond-order sum 5",
                    "interpretation": "connections means neighbor count elsewhere (c= has 3, n= has 2); no global reinterpretation as bond-order sum",
                    "rlb": "VERSION/elib only; no typing or oxidation-state specification",
                    "converter": "msi2lmp consumes supplied types; it does not resolve chemical typing",
                },
            }
        )
    return result


def _assessment(system, distribution, profile):
    require(profile in PROFILES, "Unregistered mixed Si/P profile")
    distribution.validate_integrity()
    system.validate()
    graph = chemical_graph(system)
    atoms, adj, _, env = environments(graph)
    elements = {a["element"] for a in atoms.values()}
    domain = PROFILES[profile]
    diagnostics = []
    if graph["has_box"] or graph["representation"] != "atomistic":
        diagnostics.append({"stage": "domain", "reason": "finite_atomistic_only"})
    valid = (
        {"Si", "P"} <= elements
        if domain == "mixed_sip_audit"
        else "P" in elements
        if domain == "phosphorus_audit"
        else "Si" in elements
    )
    if domain == "organosilicon_polymer":
        valid = valid and "C" in elements
    if domain == "siloxane":
        valid = valid and any(
            a["element"] == "O" and env[i]["neighbor_elements"].get("Si") == 2
            for i, a in atoms.items()
        )
    if not valid:
        diagnostics.append({"stage": "domain", "reason": "outside_selected_profile"})
    for i, a in atoms.items():
        if (
            a["formal_charge"]
            or a["metadata"]["isotope"]
            or a["metadata"]["radical_electrons"]
        ):
            diagnostics.append(
                {
                    "stage": "domain",
                    "sites": [i],
                    "reason": "charged_isotopic_or_radical_unsupported",
                }
            )
        if a["element"] == "Si":
            if a["atomic_number"] != 14 or env[i]["bond_orders"] != [1, 1, 1, 1]:
                diagnostics.append(
                    {
                        "stage": "typing",
                        "sites": [i],
                        "reason": "unsupported_silicon_valence_or_identity",
                    }
                )
            elif not set(env[i]["neighbor_elements"]) <= {"C", "H", "O"}:
                diagnostics.append(
                    {
                        "stage": "typing",
                        "sites": [i],
                        "reason": "mixed_si_neighbor_type_not_independently_qualified",
                        "environment": env[i],
                    }
                )
    p_sites = phosphorus_environments(graph)
    diagnostics.extend(
        {"stage": "typing", "sites": [p["site"]], "reason": p["reason"]}
        for p in p_sites
    )
    typing = inspection = None
    if not diagnostics:
        typing = type_pcff_atoms(
            system, distribution.source, profile="island_pcff_source_graph_v1"
        )
        diagnostics.extend(
            {"stage": "typing", **d} for d in typing.payload["diagnostics"]
        )
        if typing.complete:
            inspection = inspect_pcff_operational_support(
                system, typing, resolution_policy=MSI_POLICY, special_pairs=SPECIAL
            )
            if not inspection["charges_complete"]:
                diagnostics.append(
                    {"stage": "charges", "reason": "source_charge_incomplete"}
                )
            for q in inspection["interaction_queries"]:
                r = q["resolution"]
                if (
                    r["status"] not in ("assigned", "not_applicable")
                    or not q["dependencies_complete"]
                ):
                    diagnostics.append(
                        {
                            "stage": "parameters",
                            "reason": "required_source_row_or_dependency_unresolved",
                            "interaction": q,
                        }
                    )
                # A generic wildcard cannot authorize a new Si/P environment.
                if r["status"] == "assigned" and any(
                    atoms[i]["element"] in ("Si", "P") for i in q["sites"]
                ):
                    for row in r.get("selected", []):
                        raw = inspection["source_rows"][row["record_id"]]
                        if any(t.startswith("*") for t in raw["types"]):
                            diagnostics.append(
                                {
                                    "stage": "parameters",
                                    "reason": "exact_si_p_source_row_required",
                                    "interaction": q,
                                }
                            )
            diagnostics.extend(
                {"stage": "model", **d} for d in inspection["model_diagnostics"]
            )
    # No broader model is authorized merely because a coverage query is complete.
    # The independently qualified J41 silane is the sole retained executable.
    silane = (
        len(atoms) == 5
        and sum(a["element"] == "Si" for a in atoms.values()) == 1
        and sum(a["element"] == "H" for a in atoms.values()) == 4
        and all(
            len(adj[i]) == (4 if a["element"] == "Si" else 1) for i, a in atoms.items()
        )
    )
    if not diagnostics and not silane:
        diagnostics.append(
            {
                "stage": "independent_reference",
                "reason": "source_complete_but_fixture_not_independently_qualified",
            }
        )
    return {
        "schema": "island_pcff_mixed_sip_assessment_v1",
        "profile": profile,
        "domain": domain,
        "policy": POLICY,
        "distribution_hashes": distribution.hashes,
        "graph": graph,
        "graph_identity": identity(graph),
        "atom_order": sorted(atoms),
        "coordinates": [
            [i, list(map(float, system.coordinates.get(i)))] for i in sorted(atoms)
        ],
        "phosphorus_environments": p_sites,
        "typing": typing.payload if typing else None,
        "inspection": inspection,
        "charge_origin": POLICY["charge_origin"],
        "complete": not diagnostics,
        "execution_authorized": not diagnostics,
        "diagnostics": diagnostics,
        "periodic_status": "not_authorized",
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFMixedSiPAssessment:
    json_text: str
    distribution: PCFFSiPDistribution

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @boundary
    def validate_integrity(self, system):
        p = self.payload
        require(
            pack(p) == pack(_assessment(system, self.distribution, p["profile"])),
            "Contradictory mixed Si/P assessment, source, graph or atom ordering",
        )


@boundary
def assess_pcff_mixed_sip(system, distribution, *, profile):
    return PCFFMixedSiPAssessment(
        pack(_assessment(system, distribution, profile)), distribution
    )


@boundary
def load_pcff_mixed_sip_assessment(path, distribution, *, system):
    result = PCFFMixedSiPAssessment(Path(path).read_text(), distribution)
    result.validate_integrity(system)
    return result


@boundary
def prepare_pcff_mixed_sip(system, assessment):
    assessment.validate_integrity(system)
    require(
        assessment.payload["execution_authorized"],
        "Incomplete or unqualified mixed Si/P source: "
        + str(assessment.payload["diagnostics"]),
    )
    baseline = assess_pcff_sip(
        system, assessment.distribution.source, profile="island_pcff_silane_source_v1"
    )
    return prepare_pcff_sip(system, baseline)
