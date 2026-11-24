"""Assign force-field parameters with the RadonPy and EMC reference engines."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from island.forcefields.emc import FIELDS, EMCAssignment, assign_emc
from island.forcefields.emc_graph import EMCGraphAssignment, assign_emc_graph
from island.forcefields.emc_md import EMCMDResult, run_emc_md
from island.forcefields.radonpy import (
    FAMILIES,
    QM_CHARGES,
    RadonPyAssignment,
    _assign,
    assign_radonpy,
)

if TYPE_CHECKING:
    from rdkit.Chem import Mol

    from island.core import MolecularSystem
    from island.forcefields.assignment_core.schema import ChemicalGraph
    from island.graph.final import FinalChemicalGraph

__all__ = [
    "EMCAssignment",
    "EMCGraphAssignment",
    "EMCMDResult",
    "RadonPyAssignment",
    "assign_emc",
    "assign_emc_graph",
    "assign_polymer_parameters",
    "assign_radonpy",
    "run_emc_md",
]


def assign_polymer_parameters(
    psmiles: str | MolecularSystem | ChemicalGraph | FinalChemicalGraph | Mol,
    forcefield: str,
    *,
    dp: int | None = None,
    seed: int = 2026,
    charge: str | None = None,
    emc_root: str | Path | None = None,
    head_end_group: str | None = None,
    tail_end_group: str | None = None,
    atom_order: str | None = None,
    **backend_options,
) -> RadonPyAssignment | EMCAssignment | EMCGraphAssignment:
    """Assign a repeat unit at the given DP, or a complete graph for PCFF/OPLS."""
    if forcefield in FIELDS:
        if (
            atom_order is not None
            or charge is not None
            or head_end_group is not None
            or tail_end_group is not None
        ):
            raise ValueError(
                "EMC uses its field-specific atom order, charge method and end groups"
            )
        if not isinstance(psmiles, str):
            if dp is not None:
                raise ValueError(
                    "An already-built graph has a fixed connectivity; omit dp"
                )
            return assign_emc_graph(
                psmiles, forcefield, seed=seed, emc_root=emc_root, **backend_options
            )
        return assign_emc(
            psmiles, forcefield, dp=dp, seed=seed, emc_root=emc_root, **backend_options
        )
    if forcefield not in FAMILIES:
        raise ValueError(f"Unknown force field {forcefield!r}")
    if not isinstance(psmiles, str):
        raise TypeError("Use assign_radonpy to type an existing GAFF-family molecule")
    if emc_root is not None:
        raise ValueError("emc_root only applies to EMC force fields")
    from island.builders import build_linear_polymer

    # RadonPy also represents its two linkers as tritium. Real isotopes in
    # star-labelled repeat units keep their chemical meaning.
    if (
        isinstance(psmiles, str)
        and psmiles.count("*") < 2
        and psmiles.count("*") + psmiles.count("[3H]") == 2
    ):
        psmiles = psmiles.replace("[3H]", "*")
    end_groups = {
        "head": "*[H]" if head_end_group is None else head_end_group,
        "tail": "*[H]" if tail_end_group is None else tail_end_group,
    }
    if atom_order not in (None, "island", "radonpy"):
        raise ValueError("atom_order must be island or radonpy")
    if atom_order == "radonpy" and any(s != "*[H]" for s in end_groups.values()):
        raise ValueError("RadonPy polymer atom order requires hydrogen caps")
    stereo = any(
        isinstance(s, str) and any(c in s for c in "/\\@")
        for s in (psmiles, *end_groups.values())
    )

    # The topology builder records cap atom indices, but cannot store E/Z.
    # Removing only directional bonds keeps those indices stable for restoration.
    def topology_cap(s):
        return s.translate(str.maketrans("", "", "/\\")) if isinstance(s, str) else s

    embed = charge in QM_CHARGES
    system = build_linear_polymer(
        psmiles,
        dp=dp,
        generate_3d=embed and not stereo,
        random_seed=seed,
        head_end_group=topology_cap(end_groups["head"]),
        tail_end_group=topology_cap(end_groups["tail"]),
    )
    if stereo:
        from ._repeat_stereo import stereo_chain

        mol, sites = stereo_chain(system, psmiles, dp, seed, embed, end_groups)
    elif atom_order == "radonpy":
        from island.chemistry import to_rdkit

        graph = to_rdkit(system)
        mol = graph.mol
        sites = tuple(graph.rdkit_index_to_site_id[i] for i in range(mol.GetNumAtoms()))
    else:
        return assign_radonpy(system, forcefield, charge=charge, **backend_options)
    if atom_order == "radonpy":
        from ._repeat_stereo import radonpy_order

        mol, sites = radonpy_order(mol, sites, psmiles, dp)
    return _assign(mol, sites, forcefield, charge=charge, **backend_options)
