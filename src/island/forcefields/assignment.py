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
    emc_mode: str = "graph",
    **backend_options,
) -> RadonPyAssignment | EMCAssignment | EMCGraphAssignment:
    """Assign a finite capped chain or a complete graph with the same semantics.

    EMC defaults to complete-graph typing and returns EMCGraphAssignment.
    emc_mode="repeat" explicitly selects historical EMC polymerization (including
    nchains, density, stereo_policy and its UA methyl-cap convention).
    """
    if emc_mode not in ("graph", "repeat"):
        raise ValueError("emc_mode must be graph or repeat")
    if forcefield in FIELDS:
        if atom_order is not None or charge is not None:
            raise ValueError("EMC uses its field-specific atom order and charge method")
        if emc_mode == "repeat":
            if (
                not isinstance(psmiles, str)
                or head_end_group is not None
                or tail_end_group is not None
            ):
                raise ValueError(
                    "EMC repeat mode requires PSMILES with backend end groups"
                )
            return assign_emc(
                psmiles,
                forcefield,
                dp=dp,
                seed=seed,
                emc_root=emc_root,
                **backend_options,
            )
        legacy = {"nchains", "density", "stereo_policy"} & backend_options.keys()
        if legacy:
            raise ValueError(
                "Use emc_mode='repeat' for EMC packing/legacy stereo options, or supply a prebuilt packed graph: "
                + ", ".join(sorted(legacy))
            )
        if isinstance(psmiles, str):
            from island.builders import build_linear_polymer

            psmiles = build_linear_polymer(
                psmiles,
                dp=dp,
                generate_3d=False,
                random_seed=seed,
                head_end_group=head_end_group,
                tail_end_group=tail_end_group,
            )
        elif any(v is not None for v in (dp, head_end_group, tail_end_group)):
            raise ValueError(
                "An already-built graph has fixed connectivity and end groups; omit dp and caps"
            )
        return assign_emc_graph(
            psmiles, forcefield, seed=seed, emc_root=emc_root, **backend_options
        )
    if emc_mode != "graph":
        raise ValueError("emc_mode only applies to EMC force fields")
    if forcefield not in FAMILIES:
        raise ValueError(f"Unknown force field {forcefield!r}")
    if not isinstance(psmiles, str):
        if any(
            v is not None
            for v in (dp, atom_order, head_end_group, tail_end_group, emc_root)
        ):
            raise ValueError(
                "A built GAFF graph fixes its size, atom order and end groups"
            )
        from island.forcefields._emc_graph_input import graph_input

        source = graph_input(psmiles, seed)
        from rdkit import Chem

        conformer = Chem.Conformer(source.mol.GetNumAtoms())
        for i, position in enumerate(source.xyz):
            conformer.SetAtomPosition(i, position)
        source.mol.RemoveAllConformers()
        source.mol.AddConformer(conformer)
        return _assign(
            source.mol,
            tuple(source.sites),
            forcefield,
            charge=charge,
            **backend_options,
        )
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
    system = build_linear_polymer(
        psmiles,
        dp=dp,
        generate_3d=charge in QM_CHARGES,
        random_seed=seed,
        head_end_group=end_groups["head"],
        tail_end_group=end_groups["tail"],
    )
    if atom_order != "radonpy":
        return assign_radonpy(system, forcefield, charge=charge, **backend_options)
    from island.chemistry import to_rdkit

    graph = to_rdkit(system)
    mol = graph.mol
    sites = tuple(graph.rdkit_index_to_site_id[i] for i in range(mol.GetNumAtoms()))
    if atom_order == "radonpy":
        from ._repeat_stereo import radonpy_order

        mol, sites = radonpy_order(mol, sites, psmiles, dp)
    return _assign(mol, sites, forcefield, charge=charge, **backend_options)
