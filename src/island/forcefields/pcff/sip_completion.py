"""C2 source-completion decisions; an absent row never becomes a coefficient.

Only the independently verified finite silane has an execution profile. Other
domains are audits, and retain C1's strict exact/equivalent source requirements.
"""

import re
from dataclasses import dataclass
from pathlib import Path

from island.charge_references.records import pack, unpack

from .charges import identity
from .mixed_sip import PROFILES, assess_pcff_mixed_sip, prepare_pcff_mixed_sip
from .source import FLAGS, boundary, require

SILANE_PROFILE = "island_pcff_finite_silane_c2_v1"
EXECUTABLE_PROFILES = (SILANE_PROFILE,)
DOMAINS = tuple(PROFILES.values())
POLICY = {
    "version": "island_pcff_sip_completion_c2_v1",
    "missing_source_row": "reject; no zero, interpolation or invented coefficient",
    "phosphorus": "all source candidates must have consistent connectivity and chemical evidence",
    "new_execution": "requires complete source and independent numerical qualification",
    "charge_origin": "PCFF/source_native_bond_increments",
    "periodic": "not_authorized",
}


def p_source_evidence(distribution):
    """Keep all P declarations and complete template blocks, including conflicts."""
    atom_rows = []
    for section in distribution.source.inventory["sections"]:
        if section["name"] == "atom_types":
            atom_rows.extend(
                {**r, **r["data"]}
                for r in section["records"]
                if r["data"]["element"] == "P"
            )
    text = distribution.templates.decode("ascii")
    blocks = []
    for m in re.finditer(r"(?m)^type:\s*(p=?)\s*\n[\s\S]*?^end_type[ \t]*$", text):
        blocks.append(
            {"type": m[1], "line": text[: m.start()].count("\n") + 1, "raw": m[0]}
        )
    return {"atom_rows": atom_rows, "template_blocks": blocks}


def _decision(system, distribution, domain):
    require(domain in DOMAINS, "Unknown Si/P audit domain")
    profile = next(k for k, v in PROFILES.items() if v == domain)
    legacy = assess_pcff_mixed_sip(system, distribution, profile=profile)
    p = legacy.payload
    source = p_source_evidence(distribution)
    phosphorus = []
    graph_atoms = {a["id"]: a for a in p["graph"]["sites"]}
    for state in p["phosphorus_environments"]:
        neighbors = state["neighbors"]
        candidates = []
        for row in source["atom_rows"]:
            degree_ok = row["connections"] == len(neighbors)
            # A broad (>P) template has no oxidation/connectivity qualification.
            # The p= oxo template contradicts its five-connection declaration.
            reasons = []
            if not degree_ok:
                reasons.append("source_connections_mismatch")
            if state["formal_charge"] or any(n["formal_charge"] for n in neighbors):
                reasons.append("charged_environment_unqualified")
            reasons.append(
                "general_template_has_no_state_specific_disambiguation"
                if row["type"] == "p"
                else "phosphazene_definition_does_not_establish_this_environment"
            )
            candidates.append(
                {
                    "type": row["type"],
                    "source_line": row["line"],
                    "source_raw": row["raw"],
                    "source_connections": row["connections"],
                    "degree_matches": degree_ok,
                    "template_lines": [
                        b["line"]
                        for b in source["template_blocks"]
                        if b["type"] == row["type"]
                    ],
                    "rejection_reasons": reasons,
                    "authorized": False,
                }
            )
        # Persist explicit H/O substitution and the oxygen's second shell.
        oxygen = []
        for n in neighbors:
            if n["element"] == "O":
                attached = []
                for bond in p["graph"]["bonds"]:
                    if n["id"] in bond["sites"]:
                        j = next(i for i in bond["sites"] if i != n["id"])
                        attached.append(
                            {
                                "site": j,
                                "element": graph_atoms[j]["element"],
                                "order": bond["order"],
                            }
                        )
                oxygen.append(
                    {
                        "site": n["id"],
                        "neighbors": sorted(attached, key=lambda a: a["site"]),
                    }
                )
        phosphorus.append(
            {
                "chemical_state": state,
                "degree": len(neighbors),
                "hydrogen_neighbors": sum(n["element"] == "H" for n in neighbors),
                "oxygen_substitution": oxygen,
                "candidates": candidates,
                "selected_type": None,
                "status": "typing_blocked",
            }
        )
    inspection = p["inspection"]
    blockers = []
    if inspection:
        for q in inspection["interaction_queries"]:
            if (
                q["resolution"]["status"] not in ("assigned", "not_applicable")
                or not q["dependencies_complete"]
            ):
                blockers.append(q)
    return {
        "schema": "island_pcff_sip_completion_c2_v1",
        "domain": domain,
        "policy": POLICY,
        "c1_assessment": p,
        "c1_identity": legacy.identity,
        "source_evidence": source,
        "phosphorus_decisions": phosphorus,
        "blocking_interactions": blockers,
        "source_complete": p["complete"],
        "execution_profile": SILANE_PROFILE if p["execution_authorized"] else None,
        "execution_authorized": p["execution_authorized"],
        "ewald_status": "ewald_unvalidated",
        "pressure_status": "pressure_unavailable",
        "long_npt_status": "blocked",
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFSiPCompletion:
    json_text: str
    distribution: object

    @property
    def payload(self):
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @boundary
    def validate_integrity(self, system):
        require(
            pack(self.payload)
            == pack(_decision(system, self.distribution, self.payload["domain"])),
            "Contradictory C2 source, graph, ordering, charge or interaction decision",
        )


@boundary
def assess_pcff_sip_completion(system, distribution, *, domain):
    return PCFFSiPCompletion(
        pack(_decision(system, distribution, domain)), distribution
    )


@boundary
def load_pcff_sip_completion(path, distribution, *, system):
    result = PCFFSiPCompletion(Path(path).read_text(), distribution)
    result.validate_integrity(system)
    return result


@boundary
def prepare_pcff_sip_completion(system, decision, *, profile):
    decision.validate_integrity(system)
    require(profile in EXECUTABLE_PROFILES, "Unqualified C2 execution profile")
    require(
        decision.payload["execution_authorized"]
        and decision.payload["execution_profile"] == profile,
        "Source-incomplete or typing-ambiguous Si/P domain; no C2 execution authorized",
    )
    legacy = assess_pcff_mixed_sip(
        system,
        decision.distribution,
        profile=decision.payload["c1_assessment"]["profile"],
    )
    return prepare_pcff_mixed_sip(system, legacy)
