"""Lazy ParmEd adapter for a deliberately restricted resolved Amber prmtop subset."""

import hashlib
from collections import defaultdict, deque
from itertools import pairwise
from math import isclose, isfinite, sqrt
from pathlib import Path

from island.core import AtomSite, MolecularSystem
from island.exceptions import (
    AmberImportError,
    InvalidAmberMappingError,
    UnsupportedAmberFeatureError,
)
from island.forcefields.amber.models import (
    ImportedAmberResult,
    ImportedSelection,
    PeriodicImproperParameter,
    _digest,
)
from island.forcefields.charges import ProvidedChargeEngine
from island.forcefields.nonbonded import NonbondedPolicy
from island.forcefields.parameters.models import (
    HarmonicAngleParameter,
    HarmonicBondParameter,
    LennardJonesParameter,
    PeriodicTorsionTerm,
    ProperTorsionParameter,
)
from island.forcefields.typing.signatures import graph_signature

ADAPTER_VERSION = "1"
KCAL_TO_KJ = 4.184
ANGSTROM_TO_NM = 0.1
RMIN_TO_SIGMA = 2 / (2 ** (1 / 6))


def _unsupported(message: str) -> None:
    raise UnsupportedAmberFeatureError(message)


def _finite(value: object, label: str, *, minimum: float | None = None) -> float:
    if (not isinstance(value, (int, float)) or isinstance(value, bool)
        or not isfinite(value) or (minimum is not None and value < minimum)):
        _unsupported(f"{label} must be finite and >= {minimum}")
    return float(value)


def _canonical(sites: tuple[int, ...]) -> tuple[int, ...]:
    return min(sites, sites[::-1])


def _local_pair_distances(topology: object) -> dict[tuple[int, int], int]:
    """Discover only graph pairs within three bonds, never an atom-pair matrix."""
    distances = {}
    for origin in topology.sites:
        visited = {origin: 0}
        pending = deque([origin])
        while pending:
            current = pending.popleft()
            if visited[current] == 3:
                continue
            for neighbor in topology.neighbors(current):
                if neighbor not in visited:
                    visited[neighbor] = visited[current] + 1
                    pending.append(neighbor)
        for other, distance in visited.items():
            if origin < other:
                distances[(origin, other)] = distance
    return distances


def _source_exclusions(parm: object, mapping: dict[int, int]) -> set[tuple[int, int]]:
    counts = parm.parm_data["NUMBER_EXCLUDED_ATOMS"]
    values = parm.parm_data["EXCLUDED_ATOMS_LIST"]
    if (len(counts) != len(mapping)
        or any(not isinstance(count, int) or count < 0 for count in counts)
        or sum(counts) != len(values)):
        _unsupported("Malformed Amber exclusion counts/list")
    pairs = set()
    cursor = 0
    for source_index, count in enumerate(counts):
        for encoded in values[cursor:cursor + count]:
            if encoded == 0:
                continue  # Amber placeholder for an empty list
            other = int(encoded) - 1
            if other not in mapping or other == source_index:
                _unsupported("Invalid source exclusion atom index")
            pairs.add(tuple(sorted((mapping[source_index], mapping[other]))))
        cursor += count
    return pairs


def _validate_lj_coefficients(parm: object, by_index: dict[int, LennardJonesParameter]) -> None:
    """Compare source A/B coefficients by numerical type, never atom-name strings."""
    ntypes = parm.ptr("NTYPES")
    if set(by_index) != set(range(1, ntypes + 1)):
        _unsupported("Unused or missing numerical LJ types cannot be audited")
    matrix = parm.parm_data["NONBONDED_PARM_INDEX"]
    a_values = parm.parm_data["LENNARD_JONES_ACOEF"]
    b_values = parm.parm_data["LENNARD_JONES_BCOEF"]
    if len(matrix) != ntypes * ntypes:
        _unsupported("Invalid Amber NONBONDED_PARM_INDEX matrix")
    for i, first in by_index.items():
        for j, second in by_index.items():
            position = int(matrix[(i - 1) * ntypes + (j - 1)])
            if position <= 0 or position > len(a_values) or position > len(b_values):
                _unsupported(f"Unsupported pair coefficient index for LJ types {i}, {j}")
            sigma = (first.sigma + second.sigma) / 2 / ANGSTROM_TO_NM
            epsilon = sqrt(first.epsilon * second.epsilon) / KCAL_TO_KJ
            expected_a = 4 * epsilon * sigma ** 12
            expected_b = 4 * epsilon * sigma ** 6
            actual_a, actual_b = a_values[position - 1], b_values[position - 1]
            if not (isclose(actual_a, expected_a, rel_tol=2e-5, abs_tol=1e-8)
                    and isclose(actual_b, expected_b, rel_tol=2e-5, abs_tol=1e-8)):
                _unsupported(
                    f"Amber LJ type pair ({i}, {j}) has non-Lorentz–Berthelot "
                    "coefficients or a pair override"
                )


def import_amber_prmtop(
    system: MolecularSystem,
    prmtop: str | Path,
    source_atom_to_site: dict[int, int],
    *,
    source: str,
    force_field: str = "unknown",
    charge_method: str = "unknown",
    charge_tolerance: float = 1e-4,
) -> ImportedAmberResult:
    """Import a fixed-charge Amber prmtop with explicit **zero-based** atom mapping.

    The input system is never changed. ParmEd is imported only when called.
    This adapter requires an exact graph and a representable, global Amber-style
    1–2/1–3 exclusion and 1–4 scaling policy.
    """
    try:
        import parmed as pmd
    except ImportError as error:
        raise AmberImportError(
            "ParmEd is required for Amber import; install island[amber]"
        ) from error

    if not isinstance(system, MolecularSystem) or system.representation != "atomistic":
        raise InvalidAmberMappingError("An existing atomistic MolecularSystem is required")
    system.topology.validate_bond_graph()
    if not isinstance(source, str) or not source.strip():
        raise AmberImportError("Source provenance must be a non-empty string")
    for label, value in (("force_field", force_field), ("charge_method", charge_method)):
        if not isinstance(value, str) or not value.strip():
            raise AmberImportError(f"{label} declaration must be non-empty or 'unknown'")
    if not isinstance(source_atom_to_site, dict) or any(
        not isinstance(k, int) or isinstance(k, bool)
        or not isinstance(v, int) or isinstance(v, bool)
        for k, v in source_atom_to_site.items()
    ):
        raise InvalidAmberMappingError("Mapping must have integer source indices and stable IDs")
    ids = set(system.topology.sites)
    mapping = dict(source_atom_to_site)
    if (set(mapping) != set(range(len(ids))) or set(mapping.values()) != ids
        or len(mapping) != len(ids)):
        raise InvalidAmberMappingError(
            "Mapping must be a bijection from zero-based source atom indices "
            "to every stable site ID"
        )
    path = Path(prmtop)
    source_bytes = path.read_bytes()
    try:
        parm = pmd.load_file(str(path))
    except Exception as error:
        raise AmberImportError(f"ParmEd could not read {path}: {error}") from error
    if type(parm) is not pmd.amber.AmberParm:
        _unsupported(f"Only standard AmberParm is supported, got {type(parm).__name__}")
    for pointer in ("IFPERT", "IFCAP", "NUMEXTRA"):
        if parm.pointers.get(pointer, 0):
            _unsupported(f"Unsupported Amber topology pointer: {pointer}")
    if parm.pointers.get("NPHB", 0):
        _unsupported("10-12 hydrogen-bond terms are unsupported")
    if len(parm.atoms) != len(ids):
        raise InvalidAmberMappingError("Source and authoritative atom counts differ")
    for name in (
        "cmaps", "urey_bradleys", "rb_torsions", "adjusts", "trigonal_angles",
        "out_of_plane_bends", "pi_torsions", "stretch_bends", "torsion_torsions",
    ):
        if getattr(parm, name, ()):
            _unsupported(f"Unsupported source interaction family: {name}")
    for flag in parm.parm_data:
        if flag.startswith(("CMAP", "POLARIZABILITY", "LENNARD_JONES_14_")):
            _unsupported(f"Unsupported Amber flag: {flag}")
    if any(abs(float(x)) > 1e-12 for x in parm.parm_data.get("HBOND_ACOEF", [])):
        _unsupported("10-12 hydrogen-bond coefficients are unsupported")
    if any(abs(float(x)) > 1e-12 for x in parm.parm_data.get("HBOND_BCOEF", [])):
        _unsupported("10-12 hydrogen-bond coefficients are unsupported")

    by_index: dict[int, LennardJonesParameter] = {}
    atom_types: dict[int, str] = {}
    type_indices: dict[int, int] = {}
    site_assignments: dict[int, ImportedSelection] = {}
    charges: dict[int, float] = {}
    library_name, library_version = "amber_prmtop_resolved", "1"
    for index, atom in enumerate(parm.atoms):
        site_id = mapping[index]
        site = system.topology.sites[site_id]
        if not isinstance(site, AtomSite) or int(atom.atomic_number) != (
            site.atomic_number if site.atomic_number is not None else _element_number(site.element)
        ):
            raise InvalidAmberMappingError(
                f"Source atom {index} element does not match stable site {site_id}"
            )
        if int(atom.atomic_number) <= 0 or atom.mass <= 0:
            _unsupported(f"Virtual/extra site at source atom {index} is unsupported")
        nb_idx = int(atom.nb_idx)
        if nb_idx <= 0:
            _unsupported(f"Missing LJ numerical type for source atom {index}")
        # Numerical LJ index, not textual atom.type, is the resolved LJ identity.
        label = f"amber_nb_{nb_idx}"
        epsilon = _finite(atom.epsilon, f"epsilon at atom {index}", minimum=0)
        rmin_half = _finite(atom.rmin, f"Rmin/2 at atom {index}", minimum=0)
        if rmin_half <= 0:
            _unsupported(f"Zero Rmin/2 at source atom {index} cannot define sigma")
        for name in ("epsilon_14", "rmin_14"):
            value = getattr(atom, name, None)
            if value is not None:
                reference = epsilon if name == "epsilon_14" else rmin_half
                if not isclose(value, reference, rel_tol=1e-7, abs_tol=1e-9):
                    _unsupported(f"Distinct {name} for atom {index} is unsupported")
        lj = LennardJonesParameter(
            parameter_id=f"source_lj_type_{nb_idx}", atom_type=label,
            epsilon=epsilon * KCAL_TO_KJ,
            sigma=rmin_half * RMIN_TO_SIGMA * ANGSTROM_TO_NM,
            source=source, library_name=library_name, library_version=library_version,
        )
        if nb_idx in by_index and by_index[nb_idx] != lj:
            _unsupported(f"Numerical LJ type {nb_idx} has inconsistent atom parameters")
        by_index[nb_idx] = lj
        atom_types[site_id] = label
        type_indices[site_id] = nb_idx
        charges[site_id] = _finite(atom.charge, f"charge at atom {index}")
        site_assignments[site_id] = ImportedSelection(
            (site_id,), (label,), lj, f"ATOM {index + 1}", (nb_idx,),
        )
    _validate_lj_coefficients(parm, by_index)

    source_bonds = {
        tuple(sorted((mapping[b.atom1.idx], mapping[b.atom2.idx]))) for b in parm.bonds
    }
    if len(source_bonds) != len(parm.bonds) or source_bonds != set(system.topology.bonds):
        raise InvalidAmberMappingError("Mapped source bonds differ from authoritative connectivity")

    def make_selection(sites: tuple[int, ...], parameter: object, reference: str) -> ImportedSelection:
        return ImportedSelection(
            sites, tuple(atom_types[i] for i in sites), parameter, reference,
            tuple(type_indices[i] for i in sites),
        )

    bonds = {}
    for row, bond in enumerate(parm.bonds, 1):
        if bond.type is None:
            _unsupported(f"Bond row {row} has no parameter")
        sites = tuple(sorted((mapping[bond.atom1.idx], mapping[bond.atom2.idx])))
        parameter = HarmonicBondParameter(
            f"bond_row_{row}", tuple(atom_types[i] for i in sites),
            _finite(bond.type.k, f"bond k row {row}", minimum=0) * 2 * KCAL_TO_KJ * 100,
            _finite(bond.type.req, f"bond r0 row {row}", minimum=0) * ANGSTROM_TO_NM,
            source, library_name, library_version,
        )
        bonds[sites] = make_selection(sites, parameter, f"BOND {row}")

    angles = {}
    for row, angle in enumerate(parm.angles, 1):
        if angle.type is None:
            _unsupported(f"Angle row {row} has no parameter")
        sites = _canonical(tuple(mapping[a.idx] for a in (
            angle.atom1, angle.atom2, angle.atom3,
        )))
        if sites in angles:
            _unsupported(f"Duplicate angle interaction {sites}")
        parameter = HarmonicAngleParameter(
            f"angle_row_{row}", tuple(atom_types[i] for i in sites),
            _finite(angle.type.k, f"angle k row {row}", minimum=0) * 2 * KCAL_TO_KJ,
            _finite(angle.type.theteq, f"angle theta0 row {row}", minimum=0),
            source, library_name, library_version,
        )
        angles[sites] = make_selection(sites, parameter, f"ANGLE {row}")
    expected_angles = {
        _canonical((a, c, b)) for c in ids
        for a in system.topology.neighbors(c)
        for b in system.topology.neighbors(c) if a < b
    }
    if set(angles) != expected_angles:
        _unsupported("Source angle inventory does not cover authoritative bond angles")

    proper_terms: dict[tuple[int, ...], list[PeriodicTorsionTerm]] = defaultdict(list)
    improper_terms: dict[tuple[int, ...], list[PeriodicTorsionTerm]] = defaultdict(list)
    proper_rows: dict[tuple[int, ...], list[int]] = defaultdict(list)
    improper_rows: dict[tuple[int, ...], list[int]] = defaultdict(list)
    active_pairs: dict[tuple[int, int], tuple[float, float]] = {}
    for row, dihedral in enumerate(parm.dihedrals, 1):
        if dihedral.type is None or isinstance(dihedral.type, (list, tuple)):
            _unsupported(f"Torsion row {row} has an unsupported parameter form")
        sites = tuple(mapping[a.idx] for a in (
            dihedral.atom1, dihedral.atom2, dihedral.atom3, dihedral.atom4,
        ))
        if len(set(sites)) != 4:
            _unsupported(f"Torsion row {row} does not have four distinct atoms")
        periodicity = dihedral.type.per
        if (not isinstance(periodicity, (int, float)) or isinstance(periodicity, bool)
            or not isfinite(periodicity) or periodicity <= 0
            or int(periodicity) != periodicity):
            _unsupported(f"Torsion row {row} has unsupported periodicity")
        phase = _finite(dihedral.type.phase, f"torsion phase row {row}")
        term = PeriodicTorsionTerm(
            _finite(dihedral.type.phi_k, f"torsion k row {row}", minimum=0) * KCAL_TO_KJ,
            int(periodicity), phase % 360,
        )
        if dihedral.improper:
            if not all(i in system.topology.neighbors(sites[2]) for i in (
                sites[0], sites[1], sites[3],
            )):
                _unsupported(f"Improper row {row} does not have atom 3 central")
            improper_terms[sites].append(term)
            improper_rows[sites].append(row)
            if not dihedral.ignore_end:
                _unsupported(f"Improper row {row} unexpectedly creates a 1-4 pair")
            continue
        key = _canonical(sites)
        if not all(b in system.topology.neighbors(a) for a, b in pairwise(sites)):
            _unsupported(f"Proper torsion row {row} is not a bonded path")
        proper_terms[key].append(term)
        proper_rows[key].append(row)
        if not dihedral.ignore_end:
            pair = tuple(sorted((sites[0], sites[3])))
            scee = _finite(dihedral.type.scee, f"SCEE row {row}", minimum=0)
            scnb = _finite(dihedral.type.scnb, f"SCNB row {row}", minimum=0)
            if scee <= 0 or scnb <= 0:
                _unsupported(f"Zero SCEE/SCNB denominator in torsion row {row}")
            if pair in active_pairs:
                _unsupported(f"Duplicate effective 1-4 pair {pair}; ignore_end is inconsistent")
            active_pairs[pair] = (1 / scnb, 1 / scee)
    expected_propers = {
        _canonical((a, b, c, d))
        for b, c in system.topology.bonds
        for a in system.topology.neighbors(b) - {c}
        for d in system.topology.neighbors(c) - {b}
        if len({a, b, c, d}) == 4
    }
    if set(proper_terms) != expected_propers:
        _unsupported("Source proper torsions do not cover all authoritative 4-site paths")
    if system.topology.impropers:
        required = {
            (i.site1, i.site2, i.site3, i.site4)
            for i in system.topology.impropers
        }
        if set(improper_terms) != required:
            _unsupported("Authoritative improper inventory differs from source ordering")

    proper = {}
    for key, terms in proper_terms.items():
        rows = proper_rows[key]
        parameter = ProperTorsionParameter(
            f"proper_rows_{'_'.join(map(str, rows))}",
            tuple(atom_types[i] for i in key), tuple(terms), source,
            library_name, library_version,
        )
        proper[key] = make_selection(key, parameter, f"DIHEDRAL rows {rows}")
    impropers = {}
    for key, terms in improper_terms.items():
        rows = improper_rows[key]
        parameter = PeriodicImproperParameter(
            f"improper_rows_{'_'.join(map(str, rows))}",
            tuple(atom_types[i] for i in key), tuple(terms), source,
            library_name, library_version,
        )
        impropers[key] = make_selection(key, parameter, f"IMPROPER rows {rows}")

    excluded = _source_exclusions(parm, mapping)
    expected_excluded = set()
    expected_14 = set()
    for pair, distance in _local_pair_distances(system.topology).items():
        expected_excluded.add(pair)
        if distance == 3:
            expected_14.add(pair)
    if excluded != expected_excluded:
        _unsupported(
            "Source exclusions do not match supported shortest-path 1-2/1-3/1-4 policy; "
            f"missing={sorted(expected_excluded - excluded)}, "
            f"extra={sorted(excluded - expected_excluded)}"
        )
    if set(active_pairs) != expected_14:
        _unsupported(
            "Effective source 1-4 pairs do not match shortest-path policy; "
            f"missing={sorted(expected_14 - set(active_pairs))}, "
            f"extra={sorted(set(active_pairs) - expected_14)}"
        )
    factors = set(active_pairs.values())
    if len(factors) > 1:
        _unsupported("Nonuniform SCEE/SCNB cannot be represented by one global policy")
    lj14, coul14 = next(iter(factors), (1.0, 1.0))
    if not 0 <= lj14 <= 1 or not 0 <= coul14 <= 1:
        _unsupported("Source 1-4 scaling lies outside supported [0, 1] range")
    policy = NonbondedPolicy(
        name="amber_import_verified_lb", version="1", mixing_rule="lorentz_berthelot",
        lj_scale_12=0, coulomb_scale_12=0, lj_scale_13=0, coulomb_scale_13=0,
        lj_scale_14=lj14, coulomb_scale_14=coul14, source=source,
    )
    charge_result = ProvidedChargeEngine().assign(
        system, charges, source=f"ParmEd-decoded charges from {source}",
        tolerance=charge_tolerance,
    )
    result = ImportedAmberResult(
        graph_signature=graph_signature(system.topology),
        source_sha256=hashlib.sha256(source_bytes).hexdigest(), source=source,
        parser_version=pmd.__version__, adapter_version=ADAPTER_VERSION,
        mapping=mapping, mapping_signature=_digest(sorted(mapping.items())),
        atom_types=atom_types, source_type_indices=type_indices,
        site_assignments=site_assignments, bond_assignments=bonds,
        angle_assignments=angles, proper_torsion_assignments=proper,
        improper_assignments=impropers, source_exclusions=tuple(sorted(excluded)),
        source_14_pairs=tuple(sorted(active_pairs)), charge_result=charge_result,
        nonbonded_policy=policy,
        provenance={
            "force_field": force_field, "force_field_status": "externally_declared" if force_field != "unknown" else "unknown",
            "charge_method": charge_method, "charge_method_status": "externally_declared" if charge_method != "unknown" else "unknown",
            "source_format": "Amber prmtop", "source_indices": "zero_based",
            "source_atom_type_names": [str(atom.type) for atom in parm.atoms],
            "source_charge": "ParmEd-decoded elementary charge; no raw Amber scale applied",
            "conversion": "kcal/mol→kJ/mol 4.184; Å→nm 0.1; harmonic k×2; Rmin/2→sigma 2/2^(1/6); phases degrees",
        },
        result_signature="",
    )
    from dataclasses import replace

    signed = replace(result, result_signature=result.content_signature())
    signed.validate_integrity(system)
    return signed


def _element_number(symbol: str) -> int:
    """Minimal element lookup without requiring RDKit at the import boundary."""
    common = {
        "H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "P": 15,
        "S": 16, "Cl": 17, "Br": 35, "I": 53,
    }
    if symbol not in common:
        _unsupported(f"Element {symbol!r} requires an authoritative atomic_number")
    return common[symbol]
