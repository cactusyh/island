"""Source-bound, local-research Xiang author DFT/NPA PDMS assignment.

This is an external OPLS-compatible model, not native OPLS/PCFF, not the
paper's Qeq branch, and not a grant to redistribute the supplied assets.
Historical audits and their execution decisions are not changed.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from io import BytesIO
from itertools import combinations
from pathlib import Path
from types import SimpleNamespace
from zipfile import BadZipFile, ZipFile

from island.charge_references.records import pack, unpack
from island.exceptions import IslandError
from island.forcefields.external_pdms import LIMITS, digest, pdms_final_graph
from island.forcefields.opls_external.graph import chemical_metadata
from island.forcefields.oplsaa.pdms_author import (
    COMMANDS,
    DFT,
    PINS,
    TYPES,
    _sections,
)
from island.forcefields.oplsaa.pdms_conventions import _partition
from island.forcefields.oplsaa.pdms_literature import OPLSPDMSLiteratureError

FAMILY = "external_pdms_opls"
PROFILE = "island_external_opls_pdms_xiang2025_dft_local_v1"
SCHEMA = "island_pdms_opls_assignment_v1"
DOI = "10.1021/acs.jpcb.4c08471"
KINDS = {"bond": "Bonds", "angle": "Angles", "proper_torsion": "Dihedrals"}


class PDMSOPLSError(IslandError):
    def __init__(self, status, stage, *, result=None, **details):
        self.status = status
        self.diagnostics = [{"stage": stage, **details}]
        self.result = result
        super().__init__(f"{status}: {stage}")


def check(condition, stage, status="parameter_incomplete", **details):
    if not condition:
        raise PDMSOPLSError(status, stage, **details)


@dataclass(frozen=True)
class PDMSOPLSOptions:
    archive_path: object = None
    expected_sha256: str | None = None
    e_oxygen_site_ids: tuple[int, ...] = ()
    use_scope: str = "local_research"


def _citation(name, line, raw):
    return {
        "source_file": name,
        "sha256": PINS[name],
        "line": line,
        "raw_text": raw,
        "doi": DOI,
    }


@lru_cache(maxsize=2)
def _source(raw):
    """Cache only parsing; actual bytes are independently hashed on every load."""
    with ZipFile(BytesIO(raw)) as archive:
        names = archive.namelist()
        content = {}
        for name, pin in PINS.items():
            matches = [n for n in names if n == name or n.endswith("/" + name)]
            check(len(matches) == 1, "missing_or_duplicate_source_file", file=name)
            data = archive.read(matches[0])
            check(sha256(data).hexdigest() == pin, "source_hash_mismatch", file=name)
            content[name] = data.decode()
    prefix = "param/OPLS/opls_system"
    init = {}
    for line in content[prefix + ".in.init"].splitlines():
        tokens = line.split("#", 1)[0].split()
        if tokens:
            check(tokens[0] not in init, "duplicate_source_policy")
            init[tokens[0]] = tokens[1:]
    check(init == COMMANDS, "source_policy_changed")
    data = _sections(content[prefix + ".data"])
    atoms = data["Atoms"]
    hist = defaultdict(Counter)
    for _, _, row in atoms.values():
        hist[int(row[1])][row[2]] += 1
    settings = {}
    for n, line in enumerate(content[prefix + ".in.settings"].splitlines(), 1):
        tokens = line.split("#", 1)[0].split()
        if tokens:
            command, tid = tokens[0], int(tokens[1])
            key = command, tid
            check(key not in settings, "duplicate_source_coefficient")
            start = 3 if command == "pair_coeff" else 2
            if start == 3:
                check(tokens[1] == tokens[2], "unlike_pair_override_unsupported")
            settings[key] = (
                tokens[start:],
                _citation(prefix + ".in.settings", n, line),
            )
    source_rows = []
    for kind, section, command, coeffs in (
        ("nonbonded", "Atoms", "pair_coeff", "Pair Coeffs"),
        ("bond", "Bonds", "bond_coeff", "Bond Coeffs"),
        ("angle", "Angles", "angle_coeff", "Angle Coeffs"),
        ("proper_torsion", "Dihedrals", "dihedral_coeff", "Dihedral Coeffs"),
    ):
        uses = defaultdict(set)
        for _, _, row in data[section].values():
            tid = int(row[1] if kind == "nonbonded" else row[0])
            labels = (
                (TYPES[tid],)
                if kind == "nonbonded"
                else tuple(TYPES[int(atoms[int(i)][2][1])] for i in row[1:])
            )
            uses[tid].add(min(labels, labels[::-1]))
        for tid, candidates in sorted(uses.items()):
            check(len(candidates) == 1, "ambiguous_source_labels")
            values, evidence = settings[command, tid]
            n, text, original = data[coeffs][tid]
            check(
                [Fraction(v) for v in values]
                == [Fraction(v) for v in original[: len(values)]],
                "conflicting_source_rows",
            )
            if kind == "nonbonded":
                check(
                    list(map(Fraction, original)) == list(map(Fraction, values * 2)),
                    "special14_coefficients_differ",
                )
            row = {
                "family": kind,
                "author_type_id": tid,
                "labels": list(next(iter(candidates))),
                "values": values,
                "units": "LAMMPS real: kcal/mol, angstrom, degree equilibrium, radian displacement",
                "settings_evidence": evidence,
                "data_evidence": _citation(prefix + ".data", n, text),
                "conversion_version": "lammps_real_to_nm_kJ_equation_preserving_v1",
            }
            source_rows.append({**row, "evidence_checksum": digest(row)})
    check(len(source_rows) == 17, "incomplete_source_rows")
    return {
        "rows": source_rows,
        "masses": {TYPES[t]: data["Masses"][t][2][0] for t in TYPES},
        "charge_literals": {TYPES[t]: list(hist[t]) for t in TYPES},
        "source_hashes": PINS,
    }


def load_source(options):
    check(type(options) is PDMSOPLSOptions, "explicit_options_required")
    check(
        options.use_scope == "local_research",
        "unsupported_use_scope",
        "asset_unavailable",
    )
    check(
        options.archive_path is not None,
        "explicit_archive_required",
        "asset_unavailable",
    )
    check(
        isinstance(options.expected_sha256, str) and len(options.expected_sha256) == 64,
        "independent_archive_sha256_required",
        "asset_unavailable",
    )
    try:
        raw = Path(options.archive_path).read_bytes()
        check(
            sha256(raw).hexdigest() == options.expected_sha256,
            "archive_hash_mismatch",
            "asset_unavailable",
        )
        return unpack(pack(_source(raw)))
    except PDMSOPLSError:
        raise
    except (
        OSError,
        ValueError,
        KeyError,
        BadZipFile,
        UnicodeError,
        RuntimeError,
    ) as error:
        raise PDMSOPLSError(
            "asset_unavailable", "source_read_failed", reason=str(error)
        ) from error


def final_graph(system):
    try:
        graph = pdms_final_graph(system)
        graph["atom_order"] = list(system.topology.sites)
        graph["polymer_provenance"] = chemical_metadata(system.metadata)
        check(
            not system.topology.impropers,
            "unsupported_explicit_impropers",
            "typing_blocked",
        )
        check(
            system.metadata.get("production_validated", False) is False
            and system.metadata.get("simulation_readiness", "not_established")
            == "not_established",
            "promoted_input_readiness",
            "typing_blocked",
        )
        return graph
    except PDMSOPLSError:
        raise
    except (ValueError, KeyError, AttributeError, TypeError, IslandError) as error:
        raise PDMSOPLSError(
            "typing_blocked", "final_graph_invalid", reason=str(error)
        ) from error


def _document(system, options):
    graph = final_graph(system)
    source = load_source(options)
    anchors = options.e_oxygen_site_ids
    check(
        type(anchors) is tuple and all(type(i) is int for i in anchors),
        "explicit_integer_oxygen_anchors_required",
        "typing_blocked",
    )
    check(
        not anchors or len(anchors) == len(graph["components"]),
        "one_oxygen_anchor_per_component_required",
        "typing_blocked",
    )
    mapping = []
    for molecule, component in enumerate(graph["components"], 1):
        subset = {
            **graph,
            "components": [component],
            "sites": [a for a in graph["sites"] if a["id"] in component],
        }
        view = SimpleNamespace(
            topology=SimpleNamespace(
                sites={i: system.topology.sites[i] for i in component},
                neighbors=system.topology.neighbors,
            )
        )
        try:
            partition = _partition(
                view, subset, anchors[molecule - 1] if anchors else None
            )
        except OPLSPDMSLiteratureError as error:
            raise PDMSOPLSError(
                error.status, error.diagnostics[0]["stage"], context=error.diagnostics
            ) from error
        for r in partition["rows"]:
            e, sid = r["element"], r["site_id"]
            literal = format(float(DFT[e][int(r["fragment"] == "E")]), ".16e")
            check(
                literal in source["charge_literals"][e], "source_charge_literal_missing"
            )
            check(
                Fraction(str(system.topology.sites[sid].mass))
                == Fraction(source["masses"][e]),
                "source_mass_mismatch",
                "typing_blocked",
                site_id=sid,
            )
            mapping.append(
                {
                    "site_id": sid,
                    "element": e,
                    "molecule_id": molecule,
                    "fragment": r["fragment"],
                    "charge_e": literal,
                    "mass_dalton": source["masses"][e],
                    "author_type_id": next(t for t in TYPES if TYPES[t] == e),
                }
            )
    by_id = {r["site_id"]: r for r in mapping}
    mapping = [by_id[i] for i in graph["atom_order"]]
    adjacent = {i: set(system.topology.neighbors(i)) for i in by_id}
    bonds = sorted(min((a, b), (b, a)) for a, b in system.topology.bonds)
    angles = sorted(
        min((a, b, c), (c, b, a))
        for b in by_id
        for a, c in combinations(sorted(adjacent[b]), 2)
    )
    propers = sorted(
        {
            min((a, b, c, d), (d, c, b, a))
            for b, c in bonds
            for a in adjacent[b] - {c}
            for d in adjacent[c] - {b, a}
        }
    )
    rows = {(r["family"], tuple(r["labels"])): r for r in source["rows"]}
    assignments, missing = {}, []
    for kind, paths in (
        ("bond", bonds),
        ("angle", angles),
        ("proper_torsion", propers),
        ("nonbonded", [(i,) for i in by_id]),
    ):
        selected = []
        for path in paths:
            labels = tuple(by_id[i]["element"] for i in path)
            row = rows.get((kind, min(labels, labels[::-1])))
            if row is None:
                missing.append(
                    {"family": kind, "sites": list(path), "labels": list(labels)}
                )
            else:
                selected.append(
                    {"sites": list(path), "row_checksum": row["evidence_checksum"]}
                )
        assignments[kind] = selected
    check(not missing, "missing_interaction_rows", missing=missing)
    payload = {
        "schema": SCHEMA,
        "profile_identity": PROFILE,
        "family": FAMILY,
        "graph": graph,
        "graph_identity": digest(graph),
        "mapping": mapping,
        "archive_sha256": options.expected_sha256,
        "source": source,
        "e_oxygen_site_ids": list(anchors),
        "interactions": assignments,
        "charge_policy": {
            "origin": "Figure 7 DFT/NPA; exact author serialized literals",
            "normalization": False,
            "redistribution": False,
            "sum_e": str(sum((Fraction(r["charge_e"]) for r in mapping), Fraction())),
        },
        "model": {
            "bond": "harmonic K*(r-r0)^2",
            "angle": "harmonic K*(theta-theta0)^2",
            "proper": "OPLS four term with 1/2",
            "improper": "none in declared domain",
            "pair": "lj/charmm/coul/long 9 11",
            "mixing": "geometric_sigma_epsilon",
            "special_lj": [0, 0, 0.5],
            "special_coulomb": [0, 0, 0.5],
            "tail_correction": False,
            "paper_Qeq_model": False,
        },
        "use_scope": "user_requested_local_research",
        "redistribution_permission": "not_established",
        "scientific_status": "scientifically_unqualified",
        "limitations": LIMITS,
        "native_opls_coverage": False,
        "native_pcff_coverage": False,
        "status": "prepared",
        "diagnostics": [],
    }
    return {**payload, "identity": digest(payload)}


@dataclass(frozen=True)
class PDMSOPLSPreparation:
    json_text: str
    options: PDMSOPLSOptions

    @property
    def payload(self):
        return unpack(self.json_text)

    def require_executable(self, system):
        # Reconstruct, including raw source hashes; a re-signed forged JSON fails.
        expected = prepare_pdms_opls(system, self.options)
        check(
            expected.json_text == self.json_text,
            "stale_or_tampered_preparation",
            "typing_blocked",
        )
        p = self.payload
        if p["status"] != "prepared":
            error = PDMSOPLSError(
                p["status"], p["diagnostics"][0]["stage"], result=self
            )
            error.diagnostics = p["diagnostics"]
            raise error

    validate_integrity = require_executable


def prepare_pdms_opls(system, options):
    try:
        payload = _document(system, options)
    except PDMSOPLSError as error:
        payload = {
            "schema": SCHEMA,
            "status": error.status,
            "diagnostics": error.diagnostics,
            "profile_identity": None,
            "assignment": None,
            "limitations": LIMITS,
        }
    return PDMSOPLSPreparation(pack(payload), options)


def save_preparation(prepared, path):
    raw = prepared.json_text.encode()
    with Path(path).open("xb") as f:
        f.write(raw)
    return sha256(raw).hexdigest()


def load_preparation(path, expected_sha256, system, options):
    try:
        raw = Path(path).read_bytes()
        check(sha256(raw).hexdigest() == expected_sha256, "preparation_hash_mismatch")
        result = PDMSOPLSPreparation(raw.decode(), options)
        result.require_executable(system)
        return result
    except PDMSOPLSError:
        raise
    except (OSError, UnicodeError) as error:
        raise PDMSOPLSError(
            "asset_unavailable", "preparation_read_failed", reason=str(error)
        ) from error


def export_lammps(system, prepared, directory):
    """Deterministic complete real-unit export; no converter or zero defaults."""
    from island.pdms import box_lengths

    prepared.require_executable(system)
    p = prepared.payload
    lengths = box_lengths(system)
    root = Path(directory)
    check(not root.exists(), "output_directory_exists")
    rows = p["source"]["rows"]
    by_hash = {r["evidence_checksum"]: r for r in rows}
    by_id = {r["site_id"]: n for n, r in enumerate(p["mapping"], 1)}
    atom_types = {e: n for n, e in enumerate(("H", "O", "Si", "C"), 1)}
    types = {}
    lines = [f"ISLAND external PDMS {p['identity']}", "", f"{len(by_id)} atoms"]
    for kind, section in KINDS.items():
        lines.append(f"{len(p['interactions'][kind])} {section.lower()}")
        chosen = sorted(
            {r["row_checksum"] for r in p["interactions"][kind]},
            key=lambda h: by_hash[h]["author_type_id"],
        )
        types[kind] = {h: n for n, h in enumerate(chosen, 1)}
    lines += ["0 impropers", "", "4 atom types"]
    for kind, name in (
        ("bond", "bond"),
        ("angle", "angle"),
        ("proper_torsion", "dihedral"),
    ):
        lines.append(f"{len(types[kind])} {name} types")
    lines += [
        "",
        *[f"0 {length:.17g} {axis}lo {axis}hi" for length, axis in zip(lengths, "xyz")],
        "",
        "Masses",
        "",
    ]
    lines += [f"{tid} {p['source']['masses'][e]}" for e, tid in atom_types.items()]
    lines += ["", "Atoms # full", ""]
    for atom in p["mapping"]:
        xyz = system.coordinates.get(atom["site_id"])
        lines.append(
            f"{by_id[atom['site_id']]} {atom['molecule_id']} {atom_types[atom['element']]} "
            f"{atom['charge_e']} " + " ".join(format(float(x), ".17g") for x in xyz)
        )
    for kind, section in KINDS.items():
        lines += ["", section, ""]
        for n, row in enumerate(p["interactions"][kind], 1):
            lines.append(
                f"{n} {types[kind][row['row_checksum']]} "
                + " ".join(str(by_id[i]) for i in row["sites"])
            )
    settings = []
    for row in rows:
        kind = row["family"]
        if kind == "nonbonded":
            tid = atom_types[row["labels"][0]]
            settings.append(f"pair_coeff {tid} {tid} " + " ".join(row["values"]))
        elif row["evidence_checksum"] in types[kind]:
            command = "dihedral" if kind == "proper_torsion" else kind
            settings.append(
                f"{command}_coeff {types[kind][row['evidence_checksum']]} "
                + " ".join(row["values"])
            )
    init = (
        "units real\natom_style full\nboundary p p p\nbond_style harmonic\n"
        "angle_style harmonic\ndihedral_style opls\npair_style lj/charmm/coul/long 9 11\n"
        "pair_modify mix geometric\nspecial_bonds lj/coul 0 0 0.5\n"
        "kspace_style pppm 0.0001\nread_data system.data\ninclude coefficients.in\n"
    )
    root.mkdir(parents=True)
    files = {
        "system.data": "\n".join(lines) + "\n",
        "coefficients.in": "\n".join(settings) + "\n",
        "system.in": init,
    }
    for name, text in files.items():
        (root / name).write_text(text)
    return {n: sha256(t.encode()).hexdigest() for n, t in files.items()}
