"""Pinned author-input audit, isolated from the Qeq paper and native OPLS.

No network, numerical backend, profile registration or default-policy transfer.
Callers explicitly supply the downloaded author directory and local engine docs.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from decimal import Decimal
from hashlib import sha256
from itertools import combinations
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.forcefields.external_pdms import LIMITS, digest, require
from island.forcefields.pdms_literature import INPUT_ERRORS, pinned_bytes

from .pdms_conventions import source_evidence
from .pdms_literature import OPLSPDMSLiteratureError, OPLSPDMSLiteratureSupplement

COMMIT = "0053fda8449797e7e2b0790426cfc52dcd023e3e"
REPOSITORY = "https://github.com/COMO-lab/PDMS_benchmark"
CANONICAL_REPOSITORY = "https://github.com/WuGroup-XJTLU/PDMS_benchmark"
DATA_ID = "island_external_opls_pdms_author_input_audit_v3"
PINS = {
    "README.md": "79f2c4c56fe32f44b344b1d0c6c13c628fd1ef9b64941cda80e8d62d8866191a",
    "param/OPLS/opls_system.in": "dfd22dac270ad11bf50001573e7d581ea924e287036e34aa85cded67913d6cdd",
    "param/OPLS/opls_system.in.init": "4ae3754ef4f428f60eff9f1a33f80a040ec17b45e41ac1e48bba1be3d3ee2592",
    "param/OPLS/opls_system.in.settings": "d22a19757398bdf9fb04536f4acf74a7ad3dabce686aaba65d58675bba0016e5",
    "param/OPLS/opls_system.data": "8b03e6e8f3c51525f011d9d64e02d69d81faab78cd17ba66676b208d1c4b9c24",
}
TYPES = {98: "H", 122: "O", 868: "Si", 871: "C"}
DFT = {
    "O": ["-1.261", "-1.255"],
    "Si": ["2.121", "1.874"],
    "C": ["-1.102", "-1.083"],
    "H": ["0.224", "0.2225"],
}
SECTIONS = (
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
    "Atoms",
    "Velocities",
    "Bonds",
    "Angles",
    "Dihedrals",
    "Impropers",
)
COMMANDS = {
    "units": ["real"],
    "atom_style": ["full"],
    "bond_style": ["harmonic"],
    "angle_style": ["harmonic"],
    "dihedral_style": ["opls"],
    "improper_style": ["harmonic"],
    "pair_style": ["lj/charmm/coul/long", "9.0", "11.0"],
    "pair_modify": ["mix", "geometric"],
    "special_bonds": ["lj/coul", "0.0", "0.0", "0.5"],
    "kspace_style": ["pppm", "0.0001"],
}
DOC_PHRASES = {
    "units.rst": [
        "For style *real*, these are the units:",
        "* energy = kcal/mol",
        "* distance = Angstroms",
    ],
    "bond_harmonic.rst": ["E = K (r - r_0)^2"],
    "angle_harmonic.rst": ["E = K (", "degrees", "radians"],
    "dihedral_opls.rst": ["usual 1/2 factor is not included"],
    "pair_charmm.rst": ["epsilon_14", "sigma_14"],
    "special_bonds.rst": ["1-2", "1-3", "1-4"],
    "pair_modify.rst": ["geometric"],
}


def _cite(name, number, raw):
    return {
        "source_file": name,
        "source_sha256": PINS[name],
        "line": number,
        "raw_text": raw,
    }


def _sections(text):
    result, current = {}, None
    for n, raw in enumerate(text.splitlines(), 1):
        value = raw.split("#", 1)[0].strip()
        if not value:
            continue
        if value in SECTIONS:
            require(value not in result, "duplicate_section")
            current = result[value] = {}
        elif current is not None:
            tokens = value.split()
            key = int(tokens[0])
            require(key not in current, "duplicate_source_id")
            current[key] = (n, raw, tokens[1:])
    return result


def _conversion(kind, values):
    v = list(map(Decimal, values))
    if kind == "nonbonded":
        return {
            "epsilon_kJ_per_mol": str(v[0] * Decimal("4.184")),
            "sigma_nm": str(v[1] * Decimal("0.1")),
        }
    if kind == "bond":
        return {
            "openmm_k_kJ_per_mol_nm2": str(v[0] * Decimal("836.8")),
            "r0_nm": str(v[1] * Decimal("0.1")),
        }
    if kind == "angle":
        return {
            "openmm_k_kJ_per_mol_rad2": str(v[0] * Decimal("8.368")),
            "theta0_radian_expression": values[1] + "*pi/180",
        }
    return {
        "opls_K_kJ_per_mol": [str(x * Decimal("4.184")) for x in v],
        "equation": "sum_n K[n]/2*(1+(-1)^(n+1)*cos(n*phi)); n=1..4",
    }


def _engine_evidence(root):
    result = {}
    for name, phrases in DOC_PHRASES.items():
        raw = (Path(root) / name).read_bytes()
        lines = raw.decode().splitlines()
        citations = []
        for phrase in phrases:
            found = [(n, line) for n, line in enumerate(lines, 1) if phrase in line]
            require(
                bool(found), "engine_equation_evidence_missing:" + name + ":" + phrase
            )
            citations.append({"line": found[0][0], "raw_text": found[0][1]})
        result[name] = {"sha256": sha256(raw).hexdigest(), "citations": citations}
    return result


def _audit(sources, supplement, docs):
    require(
        type(supplement) is OPLSPDMSLiteratureSupplement, "pinned_supplement_required"
    )
    supplement = OPLSPDMSLiteratureSupplement(
        supplement.json_text, supplement.selection
    )
    text = {
        name: pinned_bytes(Path(sources) / name, pin).decode()
        for name, pin in PINS.items()
    }
    prefix = "param/OPLS/opls_system"
    main = [
        (n, raw, raw.split("#", 1)[0].strip())
        for n, raw in enumerate(text[prefix + ".in"].splitlines(), 1)
    ]
    require(
        [v for _, _, v in main if v]
        == [
            'include "opls_system.in.init"',
            "read_data ./opls_system.data",
            'include "opls_system.in.settings"',
        ],
        "changed_include_order",
    )
    init = {}
    for n, raw in enumerate(text[prefix + ".in.init"].splitlines(), 1):
        tokens = raw.split("#", 1)[0].split()
        if tokens:
            require(tokens[0] not in init, "duplicate_initialization_command")
            init[tokens[0]] = (tokens[1:], _cite(prefix + ".in.init", n, raw))
    require({k: v[0] for k, v in init.items()} == COMMANDS, "unexpected_author_policy")
    section = _sections(text[prefix + ".data"])
    atoms = section["Atoms"]
    require(len(atoms) == 32448, "incomplete_author_topology")
    require(
        {int(r[2][1]) for r in atoms.values()} == set(TYPES),
        "unrecognized_author_atom_type",
    )
    adjacency = {i: set() for i in atoms}
    for _, _, row in section["Bonds"].values():
        a, b = map(int, row[1:])
        require(a in atoms and b in atoms and a != b, "invalid_source_bond")
        require(b not in adjacency[a], "duplicate_source_bond")
        adjacency[a].add(b)
        adjacency[b].add(a)
        require(atoms[a][2][0] == atoms[b][2][0], "cross_molecule_bond")
    members = defaultdict(list)
    charges = defaultdict(Counter)
    environments = defaultdict(Counter)
    for i, (_, _, row) in atoms.items():
        molecule, atom_type, q = int(row[0]), int(row[1]), row[2]
        members[molecule].append(i)
        charges[atom_type][q] += 1
        environments[atom_type][
            tuple(sorted(TYPES[int(atoms[j][2][1])] for j in adjacency[i]))
        ] += 1
    allowed = {
        98: {("C",)},
        122: {("Si", "Si")},
        868: {("C", "C", "O", "O"), ("C", "C", "C", "O")},
        871: {("H", "H", "H", "Si")},
    }
    require(
        all(set(environments[t]) <= allowed[t] for t in TYPES),
        "unsupported_author_valence_or_environment",
    )
    expected_angles = {
        min((a, b, c), (c, b, a))
        for b in atoms
        for a, c in combinations(sorted(adjacency[b]), 2)
    }
    expected_dihedrals = set()
    for _, _, row in section["Bonds"].values():
        b, c = map(int, row[1:])
        for a in adjacency[b] - {c}:
            for d in adjacency[c] - {b, a}:
                expected_dihedrals.add(min((a, b, c, d), (d, c, b, a)))
    for name, expected in (
        ("Angles", expected_angles),
        ("Dihedrals", expected_dihedrals),
    ):
        listed = [tuple(map(int, row[1:])) for _, _, row in section[name].values()]
        canonical = [min(path, path[::-1]) for path in listed]
        require(
            len(canonical) == len(set(canonical)) and set(canonical) == expected,
            "incomplete_or_duplicate_author_interaction_set:" + name,
        )
    # Inspect the literal author's DFT branch; no Qeq substitution or charge edit.
    partitions = []
    for molecule, ids in sorted(members.items()):
        pending, reached = [ids[0]], set()
        while pending:
            i = pending.pop()
            if i not in reached:
                reached.add(i)
                pending.extend(adjacency[i] - reached)
        require(reached == set(ids), "stale_molecule_membership")
        end_si = [
            i
            for i in ids
            if int(atoms[i][2][1]) == 868
            and sum(int(atoms[j][2][1]) == 122 for j in adjacency[i]) == 1
        ]
        require(len(end_si) == 2, "unsupported_author_end_groups")
        end_o = [
            i
            for i in ids
            if int(atoms[i][2][1]) == 122
            and atoms[i][2][2] == format(float(DFT["O"][1]), ".16e")
        ]
        require(
            len(end_o) == 1 and bool(adjacency[end_o[0]] & set(end_si)),
            "ambiguous_author_end_oxygen",
        )
        end = set(end_si) | set(end_o)
        for si in end_si:
            carbons = {j for j in adjacency[si] if int(atoms[j][2][1]) == 871}
            end.update(carbons)
            end.update(
                h for c in carbons for h in adjacency[c] if int(atoms[h][2][1]) == 98
            )
        require(len(end) == 27, "incomplete_author_end_partition")
        for i in ids:
            e = TYPES[int(atoms[i][2][1])]
            require(
                atoms[i][2][2] == format(float(DFT[e][int(i in end)]), ".16e"),
                "not_figure7_DFT_write_data_value",
            )
        partitions.append(
            {
                "molecule_id": molecule,
                "atoms": len(ids),
                "terminal_silicon": sorted(end_si),
                "E_oxygen_site_id": end_o[0],
                "exact_serialized_charge_sum": str(
                    sum(Decimal(atoms[i][2][2]) for i in ids)
                ),
            }
        )
    settings = {}
    for n, raw in enumerate(text[prefix + ".in.settings"].splitlines(), 1):
        tokens = raw.split("#", 1)[0].split()
        if tokens:
            kind = tokens.pop(0)
            require(
                kind
                in {
                    "pair_coeff",
                    "bond_coeff",
                    "angle_coeff",
                    "dihedral_coeff",
                    "improper_coeff",
                },
                "unresolved_settings_command",
            )
            key = (kind, int(tokens[0]))
            require(key not in settings, "conflicting_duplicate_coefficient")
            values = tokens[2:] if kind == "pair_coeff" else tokens[1:]
            if kind == "pair_coeff":
                require(tokens[0] == tokens[1], "unexpected_unlike_pair_override")
            require(
                all(Decimal(v).is_finite() for v in values), "nonfinite_coefficient"
            )
            settings[key] = (values, _cite(prefix + ".in.settings", n, raw))
    candidates = {}
    for row in supplement.payload["rows"]:
        if row["kind"] != "partial_charge":
            labels = tuple(row["literal_labels"])
            candidates[row["kind"], min(labels, labels[::-1])] = row
    selected = []
    for family, kind, command, coeff_section in (
        ("Atoms", "nonbonded", "pair_coeff", "Pair Coeffs"),
        ("Bonds", "bond", "bond_coeff", "Bond Coeffs"),
        ("Angles", "angle", "angle_coeff", "Angle Coeffs"),
        ("Dihedrals", "proper_torsion", "dihedral_coeff", "Dihedral Coeffs"),
    ):
        uses = defaultdict(list)
        for _, _, row in section[family].values():
            type_id = int(row[1] if family == "Atoms" else row[0])
            path = (
                (TYPES[type_id],) * 2
                if family == "Atoms"
                else tuple(TYPES[int(atoms[int(i)][2][1])].lower() for i in row[1:])
            )
            uses[type_id].append(min(path, path[::-1]))
        for type_id, paths in sorted(uses.items()):
            labels = set(paths)
            require(len(labels) == 1, "chemically_ambiguous_numeric_type")
            labels = next(iter(labels))
            si_row = candidates[kind, labels]
            values, citation = settings[command, type_id]
            number, raw, original = section[coeff_section][type_id]
            require(
                list(map(Decimal, values))
                == list(map(Decimal, original[: len(values)])),
                "data_settings_conflict",
            )
            if family == "Atoms":
                require(
                    list(map(Decimal, original)) == list(map(Decimal, values * 2)),
                    "unresolved_special14_pair_coefficients",
                )
            require(
                list(map(Decimal, values))
                == list(map(Decimal, si_row["original_values"])),
                "paper_author_coefficient_conflict",
            )
            row = {
                "family": kind,
                "author_type_id": type_id,
                "labels": list(labels),
                "interaction_count": len(paths),
                "original_values": values,
                "data_evidence": _cite(prefix + ".data", number, raw),
                "selected_settings_evidence": citation,
                "paper_evidence": si_row["source_evidence"],
                "converted_values": _conversion(kind, values),
                "conversion_version": "lammps_real_to_nm_kJ_equation_preserving_v1",
                "execution_authorized": False,
            }
            selected.append({**row, "evidence_checksum": digest(row)})
    require(
        len(selected) == 17 and len(section.get("Impropers", {})) == 0,
        "incomplete_active_interactions",
    )
    p = {
        "schema": "island_opls_pdms_author_input_audit_v3",
        "data_identity": DATA_ID,
        "repository": REPOSITORY,
        "canonical_repository": CANONICAL_REPOSITORY,
        "repository_id": 902160948,
        "commit": COMMIT,
        "source_hashes": PINS,
        "v1_supplement_identity": supplement.identity,
        "license": {
            "explicit_repository_license": None,
            "redistribution_status": "not_established",
        },
        "include_order": [
            _cite(prefix + ".in", n, raw) for n, raw, value in main if value
        ],
        "active_commands": {k: v[1] for k, v in init.items()},
        "author_variant": {
            "units": "real",
            "LJ_mixing": "geometric_sigma_and_epsilon",
            "LJ_scaling_by_bond_distance_1_2_3": ["0.0", "0.0", "0.5"],
            "Coulomb_scaling_by_bond_distance_1_2_3": ["0.0", "0.0", "0.5"],
            "pair_style": ["lj/charmm/coul/long", "9.0", "11.0"],
            "charge_method": "DFT/NPA Figure 7; binary64 write_data literals",
            "extra_neigh_exclusions_in_include_closure": [],
            "improper_interactions_in_author_graph": 0,
        },
        "engine_documentation": _engine_evidence(docs),
        "selected_coefficient_rows": selected,
        "masses": {
            str(t): {
                "element": TYPES[t],
                "original_mass": section["Masses"][t][2][0],
                "evidence": _cite(prefix + ".data", *section["Masses"][t][:2]),
            }
            for t in TYPES
        },
        "source_graph": {
            "atom_count": len(atoms),
            "molecule_count": len(members),
            "interaction_counts": {
                s: len(section.get(s, {}))
                for s in ("Bonds", "Angles", "Dihedrals", "Impropers")
            },
            "atom_order_sha256": digest(list(atoms)),
            "ordered_charge_vector_sha256": digest(
                [row[2][2] for row in atoms.values()]
            ),
            "exact_serialized_total_charge": str(
                sum(Decimal(row[2][2]) for row in atoms.values())
            ),
            "molecules": partitions,
            "type_charge_histograms": {
                str(t): dict(c) for t, c in sorted(charges.items())
            },
            "type_environment_histograms": {
                str(t): [
                    {"neighbors": list(k), "count": v} for k, v in sorted(c.items())
                ]
                for t, c in sorted(environments.items())
            },
        },
        "paper_comparison": {
            "source_evidence": source_evidence(supplement),
            "bonded_and_self_LJ_rows_match": True,
            "pair_style_matches": False,
            "charge_branch_matches_table_s4": False,
            "paper_Qeq_model_reproduced": False,
            "diagnostics": [
                {
                    "stage": "pair_style_conflict",
                    "paper": "lj/cut/coul/long 11 11",
                    "author_active": "lj/charmm/coul/long 9.0 11.0",
                },
                {
                    "stage": "charge_branch_conflict",
                    "paper_table_s4": "Qeq",
                    "author_data": "DFT/NPA",
                },
            ],
        },
        "status": "scientifically_unqualified",
        "stages": {
            "author_source_identity_verified": True,
            "author_unit_policy_resolved": True,
            "author_LJ_and_special_pairs_resolved": True,
            "author_DFT_E_M_partition_verified": True,
            "active_author_coefficient_coverage_complete": True,
            "paper_Qeq_execution_qualified": False,
            "numerically_verified": False,
            "workflow_verified": False,
            "scientifically_validated": False,
        },
        "missing_dependencies": [
            "paper-versus-author model selection and provenance",
            "Qeq input matching Table S4 is not present in this OPLS source",
            "explicit asset use/redistribution provenance",
            "original Moltemplate version and typing-rule provenance for new final graphs",
            "independent backend energy/force and workflow qualification",
        ],
        "profile_identity": None,
        "execution_authorized": False,
        "execution_profiles": [],
        "backend_executions": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }
    return {**p, "identity": digest(p)}


@dataclass(frozen=True)
class AuthorOPLSAudit:
    json_text: str
    sources: Path
    supplement: OPLSPDMSLiteratureSupplement
    engine_docs: Path

    def __post_init__(self):
        try:
            require(
                pack(self.payload)
                == pack(_audit(self.sources, self.supplement, self.engine_docs)),
                "author_record_or_authority_changed",
            )
        except INPUT_ERRORS as error:
            raise OPLSPDMSLiteratureError(
                "source_incompatible", "author_input_record_invalid"
            ) from error

    @property
    def payload(self):
        return unpack(self.json_text)

    def validate_integrity(self):
        self.__post_init__()

    def require_executable(self):
        error = OPLSPDMSLiteratureError(
            "scientifically_unqualified",
            "author_variant_not_qualified",
            [
                {"stage": "qualification_missing", "dependency": d}
                for d in self.payload["missing_dependencies"]
            ],
        )
        error.result = self
        raise error


def audit_author_opls(sources, supplement, *, references, engine_docs):
    try:
        supplement.verify_sources(references)
        return AuthorOPLSAudit(
            pack(_audit(sources, supplement, engine_docs)),
            Path(sources),
            supplement,
            Path(engine_docs),
        )
    except OPLSPDMSLiteratureError:
        raise
    except INPUT_ERRORS as error:
        status = (
            "asset_unavailable"
            if isinstance(error, FileNotFoundError)
            else "source_incompatible"
        )
        raise OPLSPDMSLiteratureError(
            status,
            "author_input_audit_failed",
            [{"stage": "author_input_audit_failed", "reason": str(error)}],
        ) from error


def load_author_opls_audit(
    path, expected_sha256, *, sources, supplement, references, engine_docs
):
    """Reconstruct an explicit export; its independent pin is never self-reported."""
    try:
        raw = pinned_bytes(path, expected_sha256)
        supplement.verify_sources(references)
        return AuthorOPLSAudit(
            raw.decode(), Path(sources), supplement, Path(engine_docs)
        )
    except OPLSPDMSLiteratureError:
        raise
    except INPUT_ERRORS as error:
        status = (
            "asset_unavailable"
            if isinstance(error, FileNotFoundError)
            else "source_incompatible"
        )
        raise OPLSPDMSLiteratureError(status, "author_input_export_invalid") from error
