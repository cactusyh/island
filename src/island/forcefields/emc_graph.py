"""Assign an EMC field to the complete, already-connected molecular graph."""

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from dataclasses import field as dataclass_field
from math import fsum, isfinite
from time import monotonic

from island.exceptions import ForceFieldError
from island.forcefields.emc import FIELDS, _build_script, _field_sources


@dataclass(frozen=True)
class EMCGraphAssignment:
    field: str
    field_sha256: dict[str, str]
    emc_version: str
    data: str
    parameters: str
    # LAMMPS atom ID -> source site IDs. For UA the first site is the
    # representative atom; subsequent H sites are folded into its potential.
    site_groups: dict[int, tuple[int, ...]]
    added_hydrogen_ids: tuple[int, ...]
    coordinate_source: str
    emc_groups: dict[str, str]
    input_smiles: str
    input_smiles_site_ids: tuple[int, ...]
    automatic_ring_depth: bool = False
    ring_scan_depth: int | str = 8
    serialization: str = "input"
    assignment_attempts: list[dict] = dataclass_field(default_factory=list)
    timings_seconds: dict[str, float] = dataclass_field(default_factory=dict)

    @property
    def site_to_atom_id(self):
        return {
            site: atom for atom, sites in self.site_groups.items() for site in sites
        }


def _section(text, name):
    return text.partition("\n" + name + "\n\n")[2].split("\n\n", 1)[0]


def _field_elements(directory, field):
    """Use the provider's type-to-element declaration (some libraries use nonstandard masses)."""
    if field.startswith("pcff"):
        text = (directory / f"cff/{field}/{field}.frc").read_text()
        section = re.search(r"(?ms)^#atom_types[^\n]*\n(.*?)(?:^#|\Z)", text)
        if section:
            return {
                r[2]: r[4]
                for line in section[1].splitlines()
                if (r := line.split()) and r[0][0].isdigit()
            }
    else:
        text = (directory / (field + ".prm")).read_text()
        section = re.search(r"(?ms)^ITEM\s+MASS\s*\n(.*?)^ITEM\s+END", text)
        if section:
            return {
                r[0]: r[2]
                for line in section[1].splitlines()
                if (r := line.split()) and not r[0].startswith("#")
            }
    raise ForceFieldError(f"Missing element declarations in {field}")


class _GraphMismatch(ForceFieldError):
    """EMC interpreted a whole-group SMILES differently from the input."""


def _map_output(source, orders, data, psf, elements, ua):
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    masses = {
        int(r[0]): float(r[1])
        for r in (l.split() for l in _section(data, "Masses").splitlines())
    }
    rows = {
        int(r[0]): r for r in (l.split() for l in _section(data, "Atoms").splitlines())
    }
    header = re.search(r"(?m)^\s*(\d+) !NATOM[^\n]*\n", psf)
    if not header or int(header[1]) != len(rows):
        raise _GraphMismatch("PSF and LAMMPS atom counts disagree")
    groups = defaultdict(list)
    for line in psf[header.end() :].splitlines()[: len(rows)]:
        fields = line.split()
        groups[fields[1]].append(int(fields[0]))
    if groups.keys() != orders.keys() or sorted(
        i for ids in groups.values() for i in ids
    ) != sorted(rows):
        raise _GraphMismatch("EMC component membership differs from the input graph")
    if len({row[1] for row in rows.values()}) != len(groups):
        raise _GraphMismatch("EMC molecule count differs from the input graph")
    if any(row[-1] not in elements for row in rows.values()):
        raise _GraphMismatch("EMC emitted atom types without an element declaration")
    edges = set()
    neighbors = defaultdict(list)
    for line in _section(data, "Bonds").splitlines():
        a, b = map(int, line.split()[2:4])
        edge = tuple(sorted((a, b)))
        if edge in edges or a == b or a not in rows or b not in rows:
            raise _GraphMismatch("EMC emitted a duplicate or invalid bond")
        edges.add(edge)
        neighbors[a].append(b)
        neighbors[b].append(a)

    mol = source.mol
    mapped = {}
    for name, native in groups.items():
        order = orders[name]
        particles = native
        if ua:
            particles = [i for i in native if elements[rows[i][-1]] != "H"]
        if len(particles) != len(order):
            raise _GraphMismatch(f"EMC changed the size of input component {name}")
        for i, idx in zip(particles, order, strict=True):
            mapped[i] = [idx]
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetSymbol() != elements.get(rows[i][-1]):
                raise _GraphMismatch(
                    f"EMC element mismatch at input site {source.sites[idx]}"
                )
            if ua:
                hydrogens = [
                    h.GetIdx() for h in atom.GetNeighbors() if h.GetAtomicNum() == 1
                ]
                explicit = sorted(
                    h for h in neighbors[i] if elements[rows[h][-1]] == "H"
                )
                if len(explicit) > len(hydrogens):
                    raise _GraphMismatch("Extra united-atom hydrogen")
                for j, h in enumerate(explicit):
                    mapped[h] = [hydrogens[j]]
                folded = hydrogens[len(explicit) :]
                if folded:
                    expected = table.GetAtomicWeight(atom.GetAtomicNum()) + len(
                        folded
                    ) * table.GetAtomicWeight(1)
                    if abs(masses[int(rows[i][2])] - expected) > 0.15:
                        raise _GraphMismatch(
                            "United-atom mass disagrees with folded hydrogens"
                        )
                mapped[i].extend(folded)
        if len({rows[i][1] for i in native}) != 1:
            raise _GraphMismatch("EMC split an input component")
        if any(i not in mapped for i in native):
            raise _GraphMismatch(f"Unmapped EMC particle in component {name}")
        expected_charge = sum(
            mol.GetAtomWithIdx(idx).GetFormalCharge()
            for i in native
            for idx in mapped[i]
        )
        if (
            abs(fsum(float(rows[i][3]) for i in native) - expected_charge)
            > len(native) * 0.00005 + 1e-9
        ):
            raise _GraphMismatch(
                "EMC partial charges disagree with input formal charge"
            )

    if sorted(idx for group in mapped.values() for idx in group) != list(
        range(mol.GetNumAtoms())
    ):
        raise _GraphMismatch("EMC site map is incomplete")
    inverse = {idx: i for i, group in mapped.items() for idx in group}
    expected_edges = {
        tuple(sorted((inverse[b.GetBeginAtomIdx()], inverse[b.GetEndAtomIdx()])))
        for b in mol.GetBonds()
    }
    if edges != {edge for edge in expected_edges if edge[0] != edge[1]}:
        raise _GraphMismatch("EMC bond endpoints disagree with the mapped input graph")
    return rows, mapped


def assign_emc_graph(
    graph,
    family,
    *,
    seed=2026,
    emc_root=None,
    field_dir=None,
    cache_dir=None,
    serialization="retry",
    backend_timeout=None,
):
    """Type the actual bonds and termini, preserving source IDs and 3D coordinates.

    Accepts MolecularSystem, ChemicalGraph, FinalChemicalGraph or an RDKit Mol.
    RDKit IDs are zero-based. Absent H sites receive fresh IDs. The all-atom
    geometry carries the requested stereochemistry independently of EMC packing.
    """
    from rdkit import Chem

    from island.forcefields._emc_graph_input import emc_groups, graph_input
    from island.forcefields._emc_isotopes import split_isotope_types

    if family not in FIELDS:
        raise ValueError(f"Unknown EMC field {family!r}")
    if type(seed) is not int or not 0 < seed < 2**31:
        raise ValueError("seed must be a positive signed 32-bit integer")
    traversals = (
        "input",
        "canonical",
        "reverse",
        *(f"permuted-{i}" for i in range(1, 6)),
    )
    if serialization not in ("retry", *traversals):
        raise ValueError("Unknown EMC serialization strategy")
    if backend_timeout is not None and (
        type(backend_timeout) not in (int, float)
        or not isfinite(backend_timeout)
        or backend_timeout <= 0
    ):
        raise ValueError("backend_timeout must be finite and positive")
    started = monotonic()
    field = FIELDS[family]
    try:
        source = graph_input(graph, seed)
    except Exception as exc:
        exc.assignment_stage = "input_geometry"
        exc.elapsed_seconds = monotonic() - started
        raise
    geometry_seconds = monotonic() - started
    root, binary, directory, hashes = _field_sources(field, emc_root, field_dir)
    cache_request = None
    if cache_dir is not None:
        from pathlib import Path

        from island.forcefields._assignment_cache import (
            read_cache,
            request,
            source_binding,
        )

        cache_path = Path(cache_dir).resolve()
        if any(
            cache_path.is_relative_to(p.resolve())
            for p in (directory, root / "scripts", root / "lib", root / "templates")
        ):
            raise ValueError(
                "Assignment cache must be outside parameter and runtime source directories"
            )
        binding = source_binding(root, binary, directory)
        cache_request = request(source, field, seed, binding)
        cache_request["serialization"] = serialization
        cached = read_cache(cache_dir, cache_request)
        if cached is not None:
            return cached
    ua = field.endswith("-ua")
    elements = _field_elements(directory, field)
    selected = traversals if serialization == "retry" else (serialization,)
    attempts = []
    backend_started = monotonic()
    runtime_options = (
        {"timeout": backend_timeout} if backend_timeout is not None else {}
    )
    ring_depth = max(
        8, max((len(r) for r in source.mol.GetRingInfo().AtomRings()), default=0)
    )
    for traversal in selected:
        if ring_depth > 8:
            runtime_options["ring_depth"] = ring_depth
        auto_depth = False
        attempt_started = monotonic()
        try:
            for kekule in (False, True):
                groups, orders = emc_groups(source.mol, ua, kekule, traversal=traversal)
                definitions = "\n".join(
                    f"{name} {smiles}" for name, smiles in groups.items()
                )
                clusters = "\n".join(f"{name} {name} 1" for name in groups)
                script = f"""ITEM OPTIONS
        field {field}
        number true
        density 0.03
        system_charge false
        ITEM END
        ITEM GROUPS
        {definitions}
        ITEM END
        ITEM CLUSTERS
        {clusters}
        ITEM END
        """
                try:
                    try:
                        data, params, psf = _build_script(
                            script,
                            field,
                            root,
                            binary,
                            directory,
                            seed,
                            auto_depth,
                            psf=True,
                            **runtime_options,
                        )
                    except ForceFieldError as exc:
                        if auto_depth or "increase scanning depth?" not in str(exc):
                            raise
                        auto_depth = True
                        runtime_options.pop("ring_depth", None)
                        data, params, psf = _build_script(
                            script,
                            field,
                            root,
                            binary,
                            directory,
                            seed,
                            True,
                            psf=True,
                            **runtime_options,
                        )
                except ForceFieldError as exc:
                    exc.source_sites = [
                        source.sites[orders[group][int(site)]]
                        for group, site in re.findall(
                            r"\{group, site\} = \{(\w+), (\d+)\}", str(exc)
                        )
                        if group in orders and int(site) < len(orders[group])
                    ]
                    raise
                if not re.search(r"(?m)^pair_coeff\s", params) or re.search(
                    r"(?i)(?<!\S)[+-]?(?:nan|inf)(?!\S)", data + "\n" + params
                ):
                    raise ForceFieldError("Missing or non-finite EMC coefficients")
                try:
                    rows, mapped = _map_output(source, orders, data, psf, elements, ua)
                    break
                except _GraphMismatch:
                    if kekule or not any(
                        a.GetIsAromatic() for a in source.mol.GetAtoms()
                    ):
                        raise
        except ForceFieldError as exc:
            attempts.append(
                {
                    "traversal": traversal,
                    "passed": False,
                    "seconds": monotonic() - attempt_started,
                    "diagnostic": str(exc),
                }
            )
            exc.assignment_stage = "backend"
            exc.assignment_attempts = attempts
            exc.elapsed_seconds = monotonic() - started
            retryable = re.search(
                r"Inconsistent connectivity in rings|Missing rules|Missing force field parameters",
                str(exc),
                re.IGNORECASE,
            )
            retryable = retryable or (
                backend_timeout is not None and "EMC backend timed out" in str(exc)
            )
            if (
                traversal == selected[-1]
                or not retryable
                or not any(a.GetIsAromatic() for a in source.mol.GetAtoms())
            ):
                raise
        else:
            attempts.append(
                {
                    "traversal": traversal,
                    "passed": True,
                    "seconds": monotonic() - attempt_started,
                }
            )
            break

    atoms = _section(data, "Atoms")
    isotopes = []
    for i, indices in mapped.items():
        idx = indices[0]
        row = rows[i]
        row[4:7] = [format(float(x), ".16g") for x in source.xyz[idx]]
        atom = source.mol.GetAtomWithIdx(idx)
        if atom.GetIsotope():
            isotopes.append((row, atom.GetAtomicNum(), atom.GetIsotope()))
    if isotopes:
        data, params = split_isotope_types(data, params, list(rows.values()), isotopes)
    else:
        data = data.replace(
            atoms + "\n\n",
            "\n".join(" ".join(row) for row in rows.values()) + "\n\n",
            1,
        )
    for axis, (lo, hi) in zip("xyz", source.bounds, strict=True):
        data, count = re.subn(
            rf"(?m)^[^\n]*\s{axis}lo {axis}hi$",
            f"{lo:.16g} {hi:.16g}  {axis}lo {axis}hi",
            data,
        )
        if count != 1 or not (isfinite(lo) and isfinite(hi)):
            raise ForceFieldError("Invalid LAMMPS box")
    version = re.match(r"LAMMPS output created by (EMC [^\n]+)", data)
    if not version:
        raise ForceFieldError("Missing EMC version")
    site_groups = {
        i: tuple(source.sites[j] for j in indices) for i, indices in mapped.items()
    }
    input_smiles = Chem.MolToSmiles(source.mol, canonical=False)
    input_traversal = json.loads(source.mol.GetProp("_smilesAtomOutputOrder"))
    result = EMCGraphAssignment(
        field,
        hashes,
        version[1],
        data,
        params,
        site_groups,
        source.generated,
        source.coordinate_source,
        groups,
        input_smiles,
        tuple(source.sites[i] for i in input_traversal),
        automatic_ring_depth=auto_depth,
        ring_scan_depth="auto" if auto_depth else ring_depth,
        serialization=traversal,
        assignment_attempts=attempts,
        timings_seconds={
            "input_geometry": geometry_seconds,
            "backend": monotonic() - backend_started,
        },
    )
    if cache_request is not None:
        from island.forcefields._assignment_cache import source_binding, write_cache

        if source_binding(root, binary, directory) != cache_request["binding"]:
            raise ForceFieldError(
                "Parameter sources or implementation changed during assignment"
            )
        write_cache(cache_dir, cache_request, result)
    return result
