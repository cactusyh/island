"""Source-consistent Xiang Qeq fragment mapping, with unresolved MD policies.

This additive v2 audit does not change the v1 literature supplement. Figure 7
supports an explicit E/M partition, not arbitrary terminal-oxygen inference or
execution authority. No numerical adapter or native OPLS defaults are used.
"""

from collections import Counter
from copy import deepcopy
from dataclasses import dataclass
from decimal import Decimal

from island.charge_references.records import pack, unpack
from island.forcefields.external_pdms import LIMITS, digest, pdms_final_graph
from island.forcefields.pdms_literature import INPUT_ERRORS, pinned_bytes
from island.workflows.bundle import system_data, system_from

from .pdms_literature import (
    SOURCE_FILES,
    OPLSPDMSLiteratureError,
    OPLSPDMSLiteratureSupplement,
)

SCHEMA = "island_external_opls_pdms_conventions_v2"
DATA_ID = "island_external_opls_compatible_pdms_xiang2025_conventions_v2"
FIGURE7_SHA = "e30c927bb5bc7ff76608865b5b7ddd95995517303c9501911a9bd09cf44c8b7d"
QEQ = {
    "O": ["-0.6", "-0.578"],
    "Si": ["0.9", "0.7"],
    "C": ["-0.342", "-0.386"],
    "H": ["0.064", "0.083"],
}
MISSING = {
    "units": "author input units; LJ/bond/angle/torsion coefficient units and equations",
    "mixing": "author pair_modify mix or explicit unlike-type pair coefficients",
    "special_pairs": "author special_bonds LJ/Coulomb 1-2, 1-3, 1-4 factors",
    "exclusions": "author exclusions, special-pair generation and overrides",
    "mass": "author mass records and isotope/version policy",
    "improper_applicability": "author improper/inversion topology and coefficients or explicit inapplicability",
    "author_input": "version-pinned original Moltemplate and generated LAMMPS inputs",
    "author_fragment_mapping": "author-script confirmation of the explicit Figure 7 E/M partition",
}


def source_evidence(supplement):
    """Pinned visual transcription plus machine-reconstructable text citations."""
    sources = supplement.payload["sources"]
    main, si = SOURCE_FILES
    result = {
        "doi": "10.1021/acs.jpcb.4c08471",
        "source_hashes": {n: s["sha256"] for n, s in sources.items()},
        "figure7": {
            "source_file": main,
            "page": 7,
            "printed_page": 1870,
            "location": "Figure 7; E/M fragments and Qeq Charge columns",
            "image_resource": "I20.png",
            "image_sha256": FIGURE7_SHA,
            "transcription_method": "visual_v1; numerical Qeq columns cross-checked against Table S4",
            "E_fragments": ["-OSi(CH3)3", "-Si(CH3)3"],
            "M_fragment": "-Si(CH3)2O-",
            "column_names": ["M Monomer", "E Monomer"],
            "Qeq_Charge": QEQ,
            "scope": "source-consistent partition; exact author atom mapping still unverified",
        },
        "citations": [
            {
                "source_file": main,
                "page": 7,
                "line": 15,
                "raw_text": "Qeqchargesareemployedfortherestofthesimulations.The",
            },
            {
                "source_file": si,
                "page": 5,
                "line": 44,
                "raw_text": "using the charge equilibration method5 (Qeq) in Materials Studio 7.0. These charges were",
            },
            {
                "source_file": si,
                "page": 6,
                "line": 7,
                "raw_text": "pair_style lj/cut/coul/long 11 11",
            },
            {
                "source_file": main,
                "page": 2,
                "line": 107,
                "raw_text": "distanceof1nmisappliedtoCOMPASS,UA,anddreiding",
            },
            {
                "source_file": main,
                "page": 2,
                "line": 108,
                "raw_text": "(HUANGmodeluses1.2andOPLSuses1.1)fornonbonded",
            },
        ],
        "unit_evidence": {
            "length": "main p2 OPLS cutoff 1.1 nm; Table S4 cutoff 11",
            "reported_total_energy": "SI p8 Figure S1(a) raster axis: Total Energy [kcal/mol]",
            "conclusion": "consistent with angstrom/kcal; not a complete coefficient-unit declaration",
            "lammps_units": None,
            "coefficient_units": None,
            "conversion_performed": False,
            "figure_s1": {
                "source_file": si,
                "page": 8,
                "image_resource": "Im1.jpg",
                "image_sha256": "d13eb071fe6d2534e5cee877a658749b5f13c67faba7b7f880961676bcdaeabb",
                "visual_axis_transcription": "Total Energy [kcal/mol]",
            },
        },
    }
    return deepcopy({**result, "evidence_checksum": digest(result)})


def _blocked(stage, **context):
    raise OPLSPDMSLiteratureError(
        "typing_blocked", stage, [{"stage": stage, **context}]
    )


def _partition(system, graph, e_oxygen_site_id):
    """Recognize only a single linear Me3Si-O-(SiMe2-O)n-SiMe3 graph."""
    if len(graph["components"]) != 1:
        _blocked("single_linear_molecule_required")
    sites = {a["id"]: a for a in graph["sites"]}
    adjacent = {i: system.topology.neighbors(i) for i in sites}
    elements = {i: a["element"] for i, a in sites.items()}
    terminal, oxygens = [], []
    for i, element in elements.items():
        neighbors = Counter(elements[j] for j in adjacent[i])
        allowed = {
            "Si": ({"C": 3, "O": 1}, {"C": 2, "O": 2}),
            "O": ({"Si": 2},),
            "C": ({"Si": 1, "H": 3},),
            "H": ({"C": 1},),
        }
        if neighbors not in allowed[element]:
            _blocked(
                "unsupported_local_environment",
                site=i,
                element=element,
                neighbor_elements=dict(neighbors),
            )
        if element == "Si" and neighbors == {"C": 3, "O": 1}:
            terminal.append(i)
        if element == "O":
            oxygens.append(i)
    # Connected, degree-constrained siloxane backbone + two ends excludes cycles.
    if len(terminal) != 2:
        _blocked("two_trimethylsilyl_ends_required", terminal_silicon=terminal)
    candidates = sorted(
        {o for si in terminal for o in adjacent[si] if elements[o] == "O"}
    )
    if e_oxygen_site_id is None:
        if len(oxygens) == 1:
            selected = oxygens[0]
            origin = "unique_HMDS_oxygen"
        else:
            raise OPLSPDMSLiteratureError(
                "typing_ambiguous",
                "explicit_E_oxygen_required",
                [
                    {
                        "stage": "explicit_E_oxygen_required",
                        "candidate_site_ids": candidates,
                    }
                ],
            )
    else:
        if type(e_oxygen_site_id) is not int or e_oxygen_site_id not in candidates:
            _blocked(
                "invalid_E_oxygen_site",
                supplied=e_oxygen_site_id,
                candidate_site_ids=candidates,
            )
        selected = e_oxygen_site_id
        origin = "caller_explicit_graph_bound_E_oxygen"
    end_sites = set(terminal) | {selected}
    for si in terminal:
        carbons = {j for j in adjacent[si] if elements[j] == "C"}
        end_sites.update(carbons)
        end_sites.update(h for c in carbons for h in adjacent[c] if elements[h] == "H")
    rows = [
        {
            "site_id": i,
            "element": elements[i],
            "fragment": "E" if i in end_sites else "M",
            "literal_pair_label": elements[i],
            "literal_bonded_label": elements[i].lower(),
            "charge": QEQ[elements[i]][int(i in end_sites)],
        }
        for i in system.topology.sites
    ]
    # Source decimals are summed verbatim. No correction of any atom is possible.
    total = sum((Decimal(row["charge"]) for row in rows), Decimal(0))
    if total != 0:
        _blocked("source_partition_charge_inconsistent", exact_charge_sum=str(total))
    return {
        "rows": rows,
        "e_oxygen_site_id": selected,
        "selection_origin": origin,
        "terminal_silicon": sorted(terminal),
        "exact_charge_sum": str(total),
        "charge_origin": "Figure 7 Qeq E/M columns, independently matched to Table S4",
        "unit": "elementary_charge",
        "normalization_or_redistribution": False,
        "author_atom_mapping_verified": False,
    }


def _document(system, supplement, e_oxygen_site_id):
    if type(supplement) is not OPLSPDMSLiteratureSupplement:
        raise OPLSPDMSLiteratureError(
            "source_incompatible", "pinned_v1_supplement_required"
        )
    supplement = OPLSPDMSLiteratureSupplement(
        supplement.json_text, supplement.selection
    )
    for element, values in QEQ.items():
        if (
            supplement.lookup("partial_charge", [element])[0]["original_values"]
            != values
        ):
            raise OPLSPDMSLiteratureError(
                "source_incompatible", "figure7_table_s4_mismatch"
            )
    graph, partition, diagnostics = None, None, []
    try:
        graph = pdms_final_graph(system)
    except INPUT_ERRORS as error:
        diagnostics.append({"stage": "final_graph_invalid", "reason": str(error)})
        status = "typing_blocked"
    else:
        try:
            partition = _partition(system, graph, e_oxygen_site_id)
            status = "parameter_incomplete"
        except OPLSPDMSLiteratureError as error:
            status, diagnostics = error.status, error.diagnostics
    p = {
        "schema": SCHEMA,
        "data_identity": DATA_ID,
        "v1_supplement_identity": supplement.identity,
        "source_family": "external:opls_compatible_pdms",
        "system": system_data(system),
        "graph": graph,
        "graph_identity": None if graph is None else digest(graph),
        "atom_order": list(system.topology.sites),
        "requested_E_oxygen_site_id": e_oxygen_site_id,
        "source_evidence": source_evidence(supplement),
        "charge_partition": partition,
        "exact_charge_vector": None
        if partition is None
        else [r["charge"] for r in partition["rows"]],
        "diagnostics": diagnostics,
        "missing_dependencies": MISSING,
        "status": status,
        "stages": {
            "source_columns_crosschecked": True,
            "source_consistent_charge_partition": partition is not None,
            "author_atom_mapping_verified": False,
            "coefficient_units_complete": False,
            "mixing_complete": False,
            "special_pairs_complete": False,
            "executable": False,
        },
        "profile_identity": None,
        "parameter_assignment": None,
        "execution_authorized": False,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    return {**p, "identity": digest(p)}


class OPLSPDMSConventionsBlockedError(OPLSPDMSLiteratureError):
    def __init__(self, result):
        self.result = result
        p = result.payload
        diagnostics = p["diagnostics"] + [
            {"stage": "policy_evidence_missing", "dependency": key, "required": value}
            for key, value in p["missing_dependencies"].items()
        ]
        super().__init__(
            p["status"], "external_pdms_conventions_incomplete", diagnostics
        )


@dataclass(frozen=True)
class OPLSPDMSConventionAudit:
    """Replay-validated graph-bound evidence; self-signing cannot grant MD."""

    json_text: str
    supplement: OPLSPDMSLiteratureSupplement

    def __post_init__(self):
        try:
            p = self.payload
            expected = _document(
                system_from(p["system"]),
                self.supplement,
                p["requested_E_oxygen_site_id"],
            )
            if pack(p) != pack(expected):
                raise ValueError("record_or_authority_changed")
        except INPUT_ERRORS as error:
            raise OPLSPDMSLiteratureError(
                "source_incompatible", "convention_record_invalid"
            ) from error

    @property
    def payload(self):
        return unpack(self.json_text)

    def validate_integrity(self, system, *, references):
        self.supplement.verify_sources(references)
        p = self.payload
        try:
            current = _document(
                system, self.supplement, p["requested_E_oxygen_site_id"]
            )
            # Coordinates can change; IDs, order, chemical state and memberships cannot.
            for key in (
                "graph",
                "atom_order",
                "charge_partition",
                "status",
                "diagnostics",
            ):
                if pack(current[key]) != pack(p[key]):
                    raise ValueError("stale_graph_or_charge_partition")
            for key in ("sites", "bonds", "representation", "box", "metadata"):
                if pack(current["system"][key]) != pack(p["system"][key]):
                    raise ValueError("stale_input_chemical_state")
        except INPUT_ERRORS as error:
            raise OPLSPDMSLiteratureError(
                "source_incompatible", "stale_graph_or_charge_partition"
            ) from error

    def require_executable(self):
        raise OPLSPDMSConventionsBlockedError(self)


def assess_opls_pdms_conventions(
    system, supplement, *, references, e_oxygen_site_id=None
):
    """Verify explicit local sources and map only source-consistent E/M fragments."""
    try:
        supplement.verify_sources(references)
        return OPLSPDMSConventionAudit(
            pack(_document(system, supplement, e_oxygen_site_id)), supplement
        )
    except OPLSPDMSLiteratureError:
        raise
    except INPUT_ERRORS as error:
        raise OPLSPDMSLiteratureError(
            "typing_blocked", "final_graph_invalid"
        ) from error


def load_opls_pdms_conventions(path, expected_sha256, *, supplement, references):
    """Explicit hash-bound import; no runtime discovery, defaults or downloads."""
    try:
        raw = pinned_bytes(path, expected_sha256)
        result = OPLSPDMSConventionAudit(raw.decode("utf-8"), supplement)
        result.supplement.verify_sources(references)
        return result
    except OPLSPDMSLiteratureError:
        raise
    except INPUT_ERRORS as error:
        status = (
            "asset_unavailable"
            if isinstance(error, FileNotFoundError)
            else "source_incompatible"
        )
        raise OPLSPDMSLiteratureError(
            status, "explicit_conventions_record_invalid"
        ) from error
