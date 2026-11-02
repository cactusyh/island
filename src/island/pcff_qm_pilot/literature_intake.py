"""Traceable literature cells, with no authority to assign or execute a force field.

PDF packages are deliberately absent from this offline, ordinary-Python module.
The reviewed cell ledger is a pinned input, not a corpus that may be re-signed
into new source authority. Printed zero cells and missing cells are distinct.
"""

import json
import math
from dataclasses import dataclass
from decimal import Decimal

from island.pcff_completion import canonical, digest, require, signed

SCHEMA = "island.literature-parameter.v1"
STATUSES = (
    "exact_native_pcff_match",
    "exact_author_primary_external_match",
    "exact_literature_external_match",
    "compatible_but_non_native",
    "different_resolution",
    "different_force_function",
    "different_convention",
    "ambiguous",
    "not_applicable",
)
ENERGY = {
    "quartic_bond": "E=sum(n=2..4) k_n (r-r0)^n",
    "quartic_angle": "E=sum(n=2..4) k_n (theta-theta0)^n; angular displacements in radians",
    "torsion_3": "E=sum(n=1..3) k_n[1-cos(n*phi)]; trans=pi",
    "bond-bond": "E=K(r_ij-r0_ij)(r_jk-r0_jk)",
    "bond-angle": "E=[K_L(r_ij-r0_ij)+K_R(r_jk-r0_jk)](theta_ijk-theta0_ijk)",
    "angle-angle": "E=K(theta_ijk-theta0_ijk)(theta_kjl-theta0_kjl); j center, k shared",
    "middle_bond-torsion_3": "E=(r_jk-r0_jk)sum(n=1..3) K_n cos(n*phi_ijkl)",
    "end_bond-torsion_3": "E=sum(n=1..3)[K_Ln(r_ij-r0_ij)+K_Rn(r_kl-r0_kl)]cos(n*phi_ijkl)",
    "angle-torsion_3": "E=sum(n=1..3)[K_Ln(theta_ijk-theta0_ijk)+K_Rn(theta_jkl-theta0_jkl)]cos(n*phi_ijkl)",
    "angle-angle-torsion_1": "E=K(theta_ijk-theta0_ijk)(theta_jkl-theta0_jkl)cos(phi_ijkl)",
    "bond-bond_1_3": "E=K(r_ij-r0_ij)(r_kl-r0_kl)",
    "nonbond(9-6)": "E=epsilon[2(r_min/r)^9-3(r_min/r)^6]",
    "LJ12_6": "E=4*epsilon[(sigma/r)^12-(sigma/r)^6]",
    "OPLS_nonbond": "E=4*epsilon[(sigma/r)^12-(sigma/r)^6]+q_i*q_j/(4*pi*epsilon0*r)",
    "exp6": "E=A*exp(-B*r)-C/r^6; Coulomb separately",
    "smith_torsion": "E=-0.5*sum(n=1..3) k_n[1+cos(n*phi)]; Si-O-Si-O additionally multiplied by sin(theta)^2/sin(theta0)^2 (Table4 footnote d)",
    "RB_torsion": "E=sum(n=0..5) C_n cos(phi)^n; source Ryckaert-Bellemans convention",
    "harmonic_angle": "E=0.5*k*(theta-theta0)^2",
    "harmonic_angle_lammps": "E=K*(theta-theta0)^2",
    "harmonic_bond_lammps": "E=K*(r-r0)^2",
    "rigid_bond": "r=r0 constrained; no finite stretching coefficient supplied",
    "charge": "E_Coulomb=sum(i<j) q_i*q_j/(4*pi*epsilon0*r_ij); exclusions and scaling belong to source model",
    "bond_increment": "q_i+=delta_ij; q_j-=delta_ij; Coulomb energy uses resulting source charges",
    "wilson_out_of_plane": "E=K*chi^2; source Wilson out-of-plane definition",
    "class2_improper": "E=K*((chi1+chi2+chi3)/3-chi0)^2",
    "class2_AA_triplet": "E=M1*dt_ijk*dt_kjl+M2*dt_ijk*dt_ijl+M3*dt_ijl*dt_kjl; j central (LAMMPS ordering)",
    "class2_torsion": "E=sum(n=1..3) K_n[1-cos(n*phi-phi_n)]",
    "opls_torsion": "E=0.5*[V1(1+cos(phi))+V2(1-cos(2phi))+V3(1+cos(3phi))+V4(1-cos(4phi))]",
    "harmonic_torsion": "E=K[1+d*cos(n*phi)]",
    "charmm_torsion": "E=K[1+cos(n*phi-delta)]; fourth field is 1-4 weight",
}
# This is an analytic chain rule, not independently adjustable force offsets.
FORCE = (
    "F_a=-sum_q (partial E/partial q)*grad_a(q), differentiating the same scalar E; "
    "grad_i(r_ij)=(x_i-x_j)/r_ij; grad_i(theta_ijk)="
    "(cos(theta)*u_ji-u_jk)/(r_ji*sin(theta)); grad_j=-grad_i-grad_k; "
    "for torsions set b1=x_j-x_i,b2=x_k-x_j,b3=x_l-x_k,n1=b1 cross b2,n2=b2 cross b3; "
    "phi=atan2(h,b),h=(n1 cross n2) dot (b2/|b2|),b=n1 dot n2; "
    "grad_a(phi)=(b*grad_a(h)-h*grad_a(b))/(h*h+b*b), with product/cross-product derivatives; "
    "d cos(n*phi)/d phi=-n*sin(n*phi); d(dt^n)/dt=n*dt^(n-1); "
    "d exp(-B*r)/dr=-B*exp(-B*r); rigid constraints do not supply a finite potential force"
)


def lammps_contract(cell):
    """Name every real-unit coefficient from the declared style, never pad zeros."""
    command, sub, model = cell["command"], cell["substyle"], cell["model"]
    class2 = model == "COMPASS" or (
        model == "HUANG" and command in ("bond_coeff", "angle_coeff")
    )
    if command == "pair_coeff":
        return (
            "nonbond(9-6)" if model == "COMPASS" else "LJ12_6",
            ["epsilon", "r_min" if model == "COMPASS" else "sigma"],
            ["kcal/mol", "angstrom"],
        )
    if command == "bond_coeff":
        return (
            (
                "quartic_bond",
                ["r0", "k2", "k3", "k4"],
                ["angstrom"] + [f"kcal/(mol*angstrom^{n})" for n in (2, 3, 4)],
            )
            if class2
            else (
                "harmonic_bond_lammps",
                ["k", "r0"],
                ["kcal/(mol*angstrom^2)", "angstrom"],
            )
        )
    if command == "angle_coeff":
        if sub == "bb":
            return (
                "bond-bond",
                ["k", "r0L", "r0R"],
                ["kcal/(mol*angstrom^2)", "angstrom", "angstrom"],
            )
        if sub == "ba":
            return (
                "bond-angle",
                ["kL", "kR", "r0L", "r0R"],
                ["kcal/(mol*angstrom*radian)"] * 2 + ["angstrom"] * 2,
            )
        return (
            (
                "quartic_angle",
                ["theta0", "k2", "k3", "k4"],
                ["degree"] + [f"kcal/(mol*radian^{n})" for n in (2, 3, 4)],
            )
            if class2
            else (
                "harmonic_angle_lammps",
                ["k", "theta0"],
                ["kcal/(mol*radian^2)", "degree"],
            )
        )
    if command == "dihedral_coeff" and class2:
        if sub == "mbt":
            return (
                "middle_bond-torsion_3",
                ["k1", "k2", "k3", "r0"],
                ["kcal/(mol*angstrom)"] * 3 + ["angstrom"],
            )
        if sub in ("ebt", "at"):
            angular = sub == "at"
            return (
                "angle-torsion_3" if angular else "end_bond-torsion_3",
                [f"k{end}{n}" for end in ("L", "R") for n in (1, 2, 3)]
                + (["theta0L", "theta0R"] if angular else ["r0L", "r0R"]),
                ["kcal/(mol*radian)" if angular else "kcal/(mol*angstrom)"] * 6
                + ["degree" if angular else "angstrom"] * 2,
            )
        if sub == "aat":
            return (
                "angle-angle-torsion_1",
                ["k", "theta0L", "theta0R"],
                ["kcal/(mol*radian^2)", "degree", "degree"],
            )
        if sub == "bb13":
            return (
                "bond-bond_1_3",
                ["k", "r0L", "r0R"],
                ["kcal/(mol*angstrom^2)", "angstrom", "angstrom"],
            )
        return (
            "class2_torsion",
            ["k1", "phi1", "k2", "phi2", "k3", "phi3"],
            ["kcal/mol", "degree"] * 3,
        )
    if command == "dihedral_coeff":
        if model in ("OPLS", "HUANG"):
            return "opls_torsion", ["v1", "v2", "v3", "v4"], ["kcal/mol"] * 4
        if model == "UA":
            return (
                "harmonic_torsion",
                ["k", "d", "n"],
                ["kcal/mol", "dimensionless", "dimensionless"],
            )
        return (
            "charmm_torsion",
            ["k", "n", "delta", "weight"],
            ["kcal/mol", "dimensionless", "degree", "dimensionless"],
        )
    if command == "improper_coeff" and sub == "aa":
        return (
            "class2_AA_triplet",
            ["m1", "m2", "m3", "theta01", "theta02", "theta03"],
            ["kcal/(mol*radian^2)"] * 3 + ["degree"] * 3,
        )
    if command == "improper_coeff" and class2:
        return "class2_improper", ["k", "chi0"], ["kcal/(mol*radian^2)", "degree"]
    return "unresolved", cell["fields"], cell["raw_units"]


def normalize(value, unit):
    if value is None:
        return {"value": None, "factor": None, "unit": unit}
    factor, target = Decimal(1), unit
    if unit == "degree":
        factor, target = Decimal(str(math.pi / 180)), "radian"
    elif unit == "nm":
        factor, target = Decimal(10), "angstrom"
    elif "kJ" in unit:
        factor, target = Decimal(1) / Decimal("4.184"), unit.replace("kJ", "kcal")
    return {
        "value": str(
            Decimal(value) / Decimal("4.184")
            if "kJ" in unit
            else Decimal(value) * factor
        ),
        "factor": str(factor),
        "unit": target,
    }


def value_origin(cell):
    src, family = cell["source"], cell["family"]
    if src == "8b0e8b74303b":
        return "converted_external_archive; underlying native row authority not established"
    if src == "6e86aa640854":
        return "external_benchmark_transcription; not native PCFF authority"
    if src == "a1415925e044" and family == "exp6":
        if cell["ordered_sites"] in (
            ["Si", "O"],
            ["Si", "C"],
            ["O", "C"],
            ["O", "H"],
            ["C", "H"],
        ):
            return "derived_by_source: Aij=sqrt(Aii*Ajj), Bij=(Bii+Bjj)/2, Cij=sqrt(Cii*Cjj); printed rounded values retained"
        return "fitted_Si_or_transferred_O_C_H_by_source; Table4 source-ref column"
    if src in ("5372747c08f6", "1e95027a1202"):
        if family in ("rigid_bond", "harmonic_angle"):
            return "external_literature_transferred_by_source; PolCA Tables3–4 source column"
        if family == "RB_torsion" and cell["ordered_sites"] in (
            ["C", "C", "Si", "OB"],
            ["C", "OC", "Si", "OB"],
            ["OB", "Si", "OB", "Si"],
        ):
            return "derived_by_source_environment_reuse; TableS72 comment; no PCFF transfer"
    return "external_reported_parameter; fitted_or_transferred_by_original_source, not fitted in J71"


def energy_contract(family, cell):
    if family == "torsion_3" and cell["source"] == "65bda3a5b14c":
        return "E=sum(n=1..3) V_n[1-cos(n*phi-phi0_n)]; source Eq1b, Appendix row phases not supplied"
    return ENERGY.get(family, "extraction_unresolved: no usable potential equation")


def phase_contract(family):
    if family in (
        "angle-angle-torsion_1",
        "middle_bond-torsion_3",
        "end_bond-torsion_3",
        "angle-torsion_3",
    ):
        return "external source cos(n*phi) with zero cross-term phase; AAT cos(phi). Not native FRC literal-Phi authority."
    if family in ("class2_torsion", "torsion_3"):
        return "trans=pi; E=sum K_n[1-cos(n*phi-phi_n)]; printed phi_n degrees for LAMMPS class2, zero for Sun Appendix torsion_3"
    if family == "RB_torsion":
        return "source eq11 powers of cos(phi); no conversion of cis/trans reference into native FRC; source RB convention retained"
    if family in (
        "opls_torsion",
        "smith_torsion",
        "harmonic_torsion",
        "charmm_torsion",
    ):
        return ENERGY[family]
    return "not_applicable"


def orientation_contract(family):
    if family in ("angle-angle", "class2_AA_triplet"):
        return "j central; shared neighbor and each angle pair explicit; no unrestricted four-site reversal; triplet M order retained"
    if family == "bond-angle":
        return (
            "reverse i-j-k exchanges left/right coefficients and source bond references"
        )
    if family in ("end_bond-torsion_3", "angle-torsion_3"):
        return "reverse i-j-k-l exchanges left/right coefficient triplets and equilibrium references"
    if "torsion" in family or family == "bond-bond_1_3":
        return "reverse consecutive four-site chain; keep source dihedral reference and exchange endpoint references where present"
    if family == "bond_increment":
        return "ordered i-j increment: add delta to i, subtract from j; reversal changes sign"
    if family == "charge":
        return "one named source site; no inferred charge normalization"
    if family in (
        "quartic_angle",
        "harmonic_angle",
        "harmonic_angle_lammps",
        "bond-bond",
    ):
        return "reverse i-j-k preserves center j; exchange endpoint reference data"
    return (
        "symmetric source pair or explicitly named site; raw source site order retained"
    )


def make_record(cell, source):
    cell = json.loads(json.dumps(cell))
    family, fields, units = cell["family"], cell["fields"], cell["raw_units"]
    if family == "LAMMPS_command":
        family, fields, units = lammps_contract(cell)
    length_ok = len(fields) == len(cell["raw_values"])
    if not length_ok:
        fields, units, family = cell["fields"], cell["raw_units"], "unresolved"
    model = cell.get("model", source.get("force_field", "multiple_benchmark"))
    resolution = (
        "united_atom"
        if model in ("PolCA", "UA")
        else source.get("resolution", "all_atom")
    )
    if resolution == "mixed":
        resolution = "all_atom"
    issue = cell.get("extraction_issue")
    if family == "unresolved":
        issue = "command style/field count unresolved"
    missing = [k for k, v in zip(fields, cell["raw_values"], strict=True) if v is None]
    if family == "torsion_3" and cell["source"] == "65bda3a5b14c":
        missing += ["source_equation_phase_authority_not_printed_in_Appendix_row"]
    wildcard = any("*" in s for s in cell["ordered_sites"])
    dependencies = (
        "source equilibrium bond/angle rows must be bound with exact ordered environments; no global equilibrium propagation is inferred"
        if family
        in (
            "bond-bond",
            "bond-angle",
            "angle-angle",
            "class2_AA_triplet",
            "bond-bond_1_3",
            "middle_bond-torsion_3",
            "end_bond-torsion_3",
            "angle-torsion_3",
            "angle-angle-torsion_1",
        )
        else "none for this scalar term; whole-model nonbonded/exclusion contract remains separate"
    )
    return signed(
        {
            "schema": SCHEMA,
            "record_id": cell["cell_id"],
            "source_file_sha256": source["sha256"],
            "citation": {
                k: source.get(k) for k in ("title", "authors", "journal_year", "doi")
            },
            "location": cell["location"],
            "source_force_field_family": model,
            "source_model_resolution": resolution,
            "ordered_sites": cell["ordered_sites"],
            "interaction_family": family,
            "fields": fields,
            "raw_value": cell["raw_values"],
            "raw_unit": units,
            "normalized": [
                normalize(v, u) for v, u in zip(cell["raw_values"], units, strict=True)
            ],
            "scalar_energy_equation": energy_contract(family, cell),
            "analytic_force_equation": FORCE,
            "equilibrium_values": {
                k: v
                for k, v in zip(fields, cell["raw_values"], strict=True)
                if k.startswith(("r0", "theta0", "chi0"))
            },
            "phase_convention": (
                "source Eq1b contains phi0_n; Appendix row supplies V_n only; phase authority unresolved"
                if family == "torsion_3" and cell["source"] == "65bda3a5b14c"
                else phase_contract(family)
            ),
            "multiplicity": "source Fourier n=1..3, per declared term"
            if "torsion_3" in family
            else "one per source-declared interaction; no extra prefactor inferred",
            "orientation_reversal_rule": orientation_contract(family),
            "cross_term_dependencies": dependencies,
            "provenance_category": "literature_external",
            "value_origin": value_origin(cell),
            "extraction_confidence": "extraction_unresolved"
            if issue
            else "direct_printed_numeric_cells",
            "extraction_issue": issue,
            "missing_convention_fields": missing
            + (
                ["exact_source_equilibrium_dependencies"]
                if dependencies.startswith("source")
                else []
            )
            + (["wildcard_exact_environment"] if wildcard else []),
            "explicit_printed_zero_fields": [
                k
                for k, v in zip(fields, cell["raw_values"], strict=True)
                if v is not None and Decimal(v) == 0
            ],
            "zero_filled": False,
            "wildcard_resolved": False,
            "transferability_status": "not_established; source domain only",
            "native_pcff_authority": False,
            "execution_eligible": False,
        }
    )


@dataclass(frozen=True)
class LiteratureRecord:
    """Immutable serialized cell, independently checked against its reviewed source."""

    json_text: str

    @property
    def payload(self):
        return json.loads(self.json_text)

    def validate(self, cell, source):
        p = self.payload
        require(
            p == make_record(cell, source),
            "literature_source_location_equation_units_mutation",
        )
        return p


def element(site):
    for prefix, label in (("si", "Si"), ("o", "O"), ("h", "H"), ("c", "C")):
        if site.lower().startswith(prefix):
            return label
    return None


def candidate_sites(cell):
    sites = cell["ordered_sites"]
    if cell["source"] == "8b0e8b74303b" and cell.get("model") == "COMPASS":
        lookup = {"6": "C", "13": "H", "39": "O", "44": "Si"}
        return [lookup.get(s) for s in sites]
    if cell.get("family") == "angle-angle" and "/" in "-".join(sites):
        left, right = "-".join(sites).split("/")
        a, b = left.split("-"), right.split("-")
        # ijk / kjl, source Appendix footnote: j center, k shared.
        if len(a) == len(b) == 3 and a[1] == b[1] and a[2] == b[0]:
            sites = [a[0], a[1], a[2], b[2]]
    return [element(s) for s in sites]


def compatibility(
    record, target, elements, *, same_environment=False, complete_contract=False
):
    """Element resemblance yields a candidate, never exact PCFF source authority."""
    if record["interaction_family"] != target["family"] and (
        record["interaction_family"] != "class2_AA_triplet"
        or target["family"] != "angle-angle"
    ):
        return "not_applicable"
    want = [element(s) for s in target["ordered_types"]]
    if None in elements or len(elements) != len(want):
        return "not_applicable"
    if target["family"] == "angle-angle":
        if record["interaction_family"] == "class2_AA_triplet":
            # A LAMMPS improper AA command contains all three angle pairs.
            # This is only a non-native candidate; selecting a component also
            # requires its M/reference-order contract.
            match = elements[1] == want[1] and sorted(
                elements[::2] + elements[3:]
            ) == sorted(want[::2] + want[3:])
        else:
            # The individual FRC/Sun AA row fixes j(center) AND k(shared).
            match = elements in (want, [want[3], want[1], want[2], want[0]])
    else:
        match = elements in (want, want[::-1])
    if not match:
        return "not_applicable"
    if record["source_model_resolution"] != "all_atom":
        return "different_resolution"
    if record["source_force_field_family"] not in ("PCFF", "PCFF_related"):
        return (
            "compatible_but_non_native"
            if record["source_force_field_family"] == "COMPASS"
            else "different_force_function"
        )
    if (
        record["extraction_confidence"] == "extraction_unresolved"
        or record["missing_convention_fields"]
    ):
        return "ambiguous"
    if target["family"] == "angle-angle-torsion_1":
        return (
            "different_convention"  # native literal Phi authority remains unresolved.
        )
    if not same_environment or not complete_contract:
        return "ambiguous"
    return "exact_literature_external_match"


def coherent_bundle_gate(matrix, mappings, records):
    """Necessary source gates only; even success cannot authorize registration."""
    ids = {r["id"] for r in matrix}
    require(len(ids) == len(matrix), "duplicate_target_instance")
    require(
        len({r["record_id"] for r in records}) == len(records),
        "duplicate_literature_record",
    )
    require(
        all(not r["zero_filled"] and not r["wildcard_resolved"] for r in records),
        "zero_fill_wildcard_rejected",
    )
    exact = {
        m["target_id"]
        for m in mappings
        if m["status"].startswith("exact_") and m.get("coherent_contract") is True
    }
    return ids == exact and bool(ids)


def converter_gate(bundle, *, matrix=None, mappings=None, records=None):
    require(
        bundle.get("parameter_completeness") == "complete"
        and bundle.get("coherent_contract") is True,
        "partial_bundle_converter_forbidden",
    )
    require(
        bundle.get("resolved_instances") == 378
        and bundle.get("remaining_row_ids") == [],
        "partial_bundle_instances_forbidden",
    )
    require(
        matrix is not None and mappings is not None and records is not None,
        "complete_bundle_requires_source_mapping_evidence",
    )
    require(
        len(matrix) == 378 and coherent_bundle_gate(matrix, mappings, records),
        "partial_bundle_source_mapping_rejected",
    )
    require(
        not bundle.get("zero_fill") and not bundle.get("wildcard_fallback"),
        "permissive_converter_forbidden",
    )
    return True


def record_object(cell, source):
    return LiteratureRecord(canonical(make_record(cell, source)))


def semantic_declarations(record):
    expected = {
        "decision": "source_conflict",
        "parameter_completeness": "incomplete",
        "accepted_coefficients": 0,
        "automatic_assignment": False,
        "new_QM": 0,
        "QM_fits": 0,
        "MD_steps": 0,
        "reserved_validation_activated": False,
        "reserved_validation_consumed": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "scientific_validation": "not_established",
        "chemical_transfer": "not_established",
        "original_author_artifact_status": "not_yet_supplied",
        "author_confirmations": "28 user-attested",
    }
    require(
        record.get("identity")
        == digest({k: v for k, v in record.items() if k != "identity"}),
        "literature_outer_identity",
    )
    for key, value in expected.items():
        require(
            type(record.get(key)) is type(value) and record[key] == value,
            "J71_scientific_declaration",
            field=key,
        )
    return expected
