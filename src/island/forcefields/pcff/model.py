"""Versioned LAMMPS-compatible definition, separate from historical assignments.

Model specification is separate from evaluator construction and readiness claims.
"""

from collections import Counter, deque
from copy import deepcopy
from dataclasses import dataclass
from math import isfinite
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.workflows.storage import publish

from .charges import identity
from .class2 import PCFFClass2Result
from .source import FLAGS, boundary, require
from .validation_cache import record_identity, validated_record

SCHEMA = "island_pcff_class2_model_v1"
PROFILE = {
    "name": "island_lammps_pcff_acyclic_cho_v1",
    "lammps_revision": "e891a3e10973c1a729e391a0aefaa02fd70f8c0f",
    "frc_sha256": "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c",
    "typing_profile": "island_pcff_acyclic_cho_v1",
    "units": {
        "energy": "kJ/mol",
        "distance": "angstrom",
        "angle": "radian",
        "charge": "elementary_charge",
    },
    "equations": {
        "quartic_bond": "K2*(r-r0)^2+K3*(r-r0)^3+K4*(r-r0)^4",
        "quartic_angle": "K2*(theta-theta0)^2+K3*(theta-theta0)^3+K4*(theta-theta0)^4",
        "bond-bond": "K*(r12-r01)*(r23-r02)",
        "bond-angle": "(Kleft*(r12-r01)+Kright*(r23-r02))*(theta123-theta0)",
        "torsion_3": "sum(n=1..3,Vn*(1-cos(n*phi-phase_n)))",
        "middle_bond-torsion_3": "(r23-r0)*sum(n=1..3,Fn*cos(n*phi))",
        "end_bond-torsion_3": "(r12-r01)*sum(Fleft_n*cos(n*phi))+(r34-r02)*sum(Fright_n*cos(n*phi))",
        "angle-torsion_3": "(theta123-theta01)*sum(Fleft_n*cos(n*phi))+(theta234-theta02)*sum(Fright_n*cos(n*phi))",
        "angle-angle-torsion_1": "K*(theta123-theta01)*(theta234-theta02)*cos(phi)",
        "bond-bond_1_3": "K*(r12-r01)*(r34-r02)",
        "angle-angle": "K*(theta123-theta01)*(theta324-theta02)",
        "lj_9_6": "epsilon_ij*(2*(rmin_ij/r)^9-3*(rmin_ij/r)^6)",
        "coulomb": "coulomb_constant*qi*qj/r",
    },
    "torsion_geometry": "v1=x1-x2; t=(x3-x2)/|x3-x2|; v3=x4-x3; a=v1-(v1.t)t; b=v3-(v3.t)t; phi=atan2((t cross a).b,a.b); cis=0,trans=pi",
    "reversal": "full proper reversal preserves signed phi; reverse left/right coefficient and equilibrium blocks; AA swaps outer 1/4 only, center 2/shared arm 3 fixed",
    "angle_angle_multiplicity": "three couplings per unordered neighbor triple, including four triples at tetrahedral carbon",
    "aat_units": "kJ/(mol*radian^2); equilibrium angles radians",
    "frc_discrepancies": [
        "torsion_3 FRC plus-cos comment; executable minus-cos adopted with unchanged coefficients/phases",
        "AAT FRC linear-Phi comment; executable cos(phi) adopted with unchanged scalar",
    ],
    "bb13_policy": "msi2lmp_non_cp_zero_v1: initialize K=0; only cp-containing supplied-type propers trigger lookup; no other missing-term fallback",
    "bb13_evidence": "tools/msi2lmp/src/GetParameters.c: initialization and cp_type/quo_cp conditional",
    "torsion_torsion": "excluded by named operational model scope; not a physical assertion inferred from an empty section",
    "wilson": "excluded: supported saturated explicit-H CHO has no degree-three center; angle-angle retained",
    "mixing": {
        "rmin": "((ri^6+rj^6)/2)^(1/6)",
        "epsilon": "2*sqrt(ei*ej)*ri^3*rj^3/(ri^6+rj^6)",
        "distance_meaning": "minimum, not 12-6 sigma",
    },
    "coulomb_constant": 332.06371 * 4.184,
    "coulomb_constant_units": "kJ*angstrom/(mol*e^2)",
    "coulomb_evidence": "src/update.cpp units real qqr2e=332.06371 kcal*angstrom/(mol*e^2)",
    "boundary": "finite nonperiodic unconstrained, all finite pairs, no switching, shift, tail correction or reciprocal electrostatics",
    "special_pairs": "required caller choice for shortest bond-path distances 1,2,3; all other pairs weight 1; not a recovered source fact",
}


@boundary
def special_pair_policy(*, lj, coulomb):
    """Explicit independent 1-2/1-3/1-4 weights, each in [0,1]. No default."""
    for values in (lj, coulomb):
        require(
            type(values) in (list, tuple) and len(values) == 3,
            "Three special-pair weights required",
        )
        require(
            all(
                type(v) in (int, float) and isfinite(v) and 0 <= v <= 1 for v in values
            ),
            "Invalid special-pair weights",
        )
    return {
        "schema": "island_pcff_special_pairs_v1",
        "origin": "explicit_caller_model_choice",
        "lj": list(lj),
        "coulomb": list(coulomb),
        "distance": "shortest_authoritative_bond_path",
    }


def bb13_policy_applies(
    family, supplied_types, dependencies, *, resolution_policy=None
):
    """Existing H4 converter compatibility rule; not a source-derived zero.

    Keep this one predicate shared by model construction and diagnostic
    adjudication. Do not broaden it to other families or missing dependencies.
    """
    from .fallbacks import MSI_POLICY

    return (
        resolution_policy != MSI_POLICY
        and family == "bond-bond_1_3"
        and "cp" not in supplied_types
        and all(d["status"] == "assigned" for d in dependencies)
    )


def definition(assignment, policy):
    from .expanded import PROFILE_NAME

    fallback = assignment["schema"] == "island_pcff_source_class2_assignment_v2"
    expanded = (
        fallback or assignment["schema"] == "island_pcff_source_class2_assignment_v1"
    )
    profile = deepcopy(PROFILE)
    if expanded:
        profile.update(
            name="island_lammps_pcff_source_graph_v1",
            typing_profile=PROFILE_NAME,
            wilson="degree-three centers; K*(mean(three signed asin out-of-plane angles)-chi0)^2; zero chi0 only until nonzero permutation semantics verified",
        )
        profile["equations"]["wilson_out_of_plane"] = (
            "K*(mean(asin(u1.(u2 cross u3)/sin(theta23)), cyclic)-chi0)^2"
        )
    if fallback:
        from .fallbacks import policy_evidence

        profile.update(
            name="island_lammps_pcff_source_fallbacks_v"
            + assignment["resolution_policy"][-1],
            typing_profile=assignment["charge_record"]["automatic_typing"]["profile"][
                "name"
            ],
            resolution_policy=policy_evidence(assignment["resolution_policy"]),
        )
        profile["equations"].update(
            quadratic_bond="K2*(r-r0)^2",
            quadratic_angle="K2*(theta-theta0)^2",
            torsion_1="Kphi*(1+cos(n*phi-phase))",
        )
        from .fallbacks import MSI_POLICY

        if assignment["resolution_policy"] == MSI_POLICY:
            profile["bb13_policy"] = "source_row_required_v1: no compatibility zero"
            profile["bb13_evidence"] = (
                "Explicit strict-source model choice; historical converter non-cp initialization does not supply a source row"
            )
    require(
        policy == special_pair_policy(lj=policy["lj"], coulomb=policy["coulomb"]),
        "Contradictory special-pair policy",
    )
    require(
        assignment["source"]["sha256"] == PROFILE["frc_sha256"], "Model source mismatch"
    )
    auto = assignment["charge_record"]["automatic_typing"]
    require(
        auto["profile"]["name"] == profile["typing_profile"],
        "Model typing profile mismatch",
    )
    graph = auto["graph"]
    terms = []
    diagnostics = []
    for a in assignment["assignments"]:
        family = a["family"]
        if family == "nonbond(9-6)":
            if a["status"] != "assigned":
                diagnostics.append(
                    {"id": a["id"], "reason": "nonbonded row unresolved"}
                )
            continue
        if family == "wilson_out_of_plane":
            if a["status"] == "not_applicable":
                continue
            if (
                not expanded
                or a["status"] != "assigned"
                or a["normalized_values"][1] != 0
            ):
                diagnostics.append(
                    {
                        "id": a["id"],
                        "reason": "Wilson missing/ambiguous or nonzero equilibrium permutation semantics unresolved",
                    }
                )
                continue
        term = {
            "assignment_id": a["id"],
            "family": family,
            "sites": a["sites"],
            "equilibria": [d["equilibrium_value"] for d in a["dependencies"]],
            "dependencies": a["dependencies"],
            "raw_status": a["status"],
        }
        if bb13_policy_applies(
            family,
            a["supplied_types"],
            a["dependencies"],
            resolution_policy=assignment.get("resolution_policy"),
        ):
            term.update(
                coefficients=[0.0],
                origin="policy_derived_zero",
                source_rows=[],
                policy=PROFILE["bb13_policy"],
            )
        elif a["status"] == "assigned":
            term.update(
                coefficients=a["normalized_values"],
                origin="source_row",
                source_rows=a["selected"],
            )
        else:
            diagnostics.append({"id": a["id"], "reason": a.get("reason", a["status"])})
            continue
        terms.append(term)
    adjacency = {s["id"]: set() for s in graph["sites"]}
    for b in graph["bonds"]:
        i, j = b["sites"]
        adjacency[i].add(j)
        adjacency[j].add(i)
    pairs = []
    for i in sorted(adjacency):
        distances = {i: 0}
        queue = deque([i])
        while queue:
            j = queue.popleft()
            if distances[j] == 3:
                continue
            for k in sorted(adjacency[j]):
                if k not in distances:
                    distances[k] = distances[j] + 1
                    queue.append(k)
        for j, d in sorted(distances.items()):
            if j > i:
                pairs.append(
                    {
                        "sites": [i, j],
                        "path_length": d,
                        "lj_weight": policy["lj"][d - 1],
                        "coulomb_weight": policy["coulomb"][d - 1],
                    }
                )
    nonbonded = []
    charges = assignment["charge_record"]["native_charge_record"]["partial_charges"]
    for a in assignment["assignments"]:
        if a["family"] == "nonbond(9-6)" and a["status"] == "assigned":
            nonbonded.append(
                {
                    "site": a["sites"][0],
                    "rmin": a["normalized_values"][0],
                    "epsilon": a["normalized_values"][1],
                    "charge": charges[a["sites"][0]],
                    "source_rows": a["selected"],
                }
            )
    return {
        "schema": "island_pcff_source_model_v2"
        if fallback
        else "island_pcff_source_model_v1"
        if expanded
        else SCHEMA,
        "compatibility_profile": profile,
        "source": assignment["source"],
        "assignment_identity": identity(assignment),
        "typing_identity": assignment["typing_identity"],
        "charge_identity": assignment["charge_identity"],
        "graph_identity": assignment["graph_identity"],
        "special_pairs": deepcopy(policy),
        "special_pair_inventory": pairs,
        "terms": terms,
        "nonbonded": nonbonded,
        "diagnostics": diagnostics,
        "raw_parameter_coverage_complete": assignment["parameter_coverage_complete"],
        "model_definition_complete": not diagnostics,
        "term_origins": dict(Counter(t["origin"] for t in terms)),
        "numerical_verification": "not_attested_by_definition",
        "simulation_snapshot_available": False,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFModelSpecification:
    """Compact signed model referencing an external, unmodified H3 assignment."""

    json_text: str
    assignment: PCFFClass2Result

    @boundary
    @validated_record
    def validate_integrity(self, system=None):
        require(type(self.assignment) is PCFFClass2Result, "H3 assignment required")
        self.assignment.validate_integrity(system)
        data = unpack(self.json_text)
        expected = definition(unpack(self.assignment.json_text), data["special_pairs"])
        require(pack(data) == pack(expected), "Contradictory PCFF model specification")

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return record_identity(self)

    def numerical_terms(self):
        """Validate once; return owned bonded kernels. Nonbonded records stay separate."""
        from .terms import Class2Term, FallbackClass2Term, SourceClass2Term

        p = self.payload
        require(p["model_definition_complete"], "Incomplete model definition")
        kernel = (
            FallbackClass2Term
            if p["schema"] == "island_pcff_source_model_v2"
            else SourceClass2Term
            if p["schema"] == "island_pcff_source_model_v1"
            else Class2Term
        )
        return tuple(
            kernel(
                t["family"],
                tuple(t["sites"]),
                tuple(t["coefficients"]),
                tuple(t["equilibria"]),
            )
            for t in p["terms"]
        )


@boundary
def define_pcff_model(assignment, *, special_pairs):
    require(type(assignment) is PCFFClass2Result, "H3 assignment required")
    p = assignment.payload
    result = PCFFModelSpecification(pack(definition(p, special_pairs)), assignment)
    return result


@boundary
def save_pcff_model(result, path):
    require(type(result) is PCFFModelSpecification, "PCFF model specification required")
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_model(path, assignment, *, system=None):
    result = PCFFModelSpecification(Path(path).read_text(), assignment)
    result.validate_integrity(system)
    return result
