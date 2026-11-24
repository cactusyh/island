"""Assign force-field parameters with the RadonPy and EMC reference engines."""

from pathlib import Path

from island.forcefields.emc import EMCAssignment, assign_emc
from island.forcefields.radonpy import RadonPyAssignment, assign_radonpy

__all__ = [
    "EMCAssignment",
    "RadonPyAssignment",
    "assign_emc",
    "assign_polymer_parameters",
    "assign_radonpy",
]


def assign_polymer_parameters(
    psmiles: str,
    forcefield: str,
    *,
    dp: int,
    seed: int = 2026,
    charge: str | None = None,
    emc_root: str | Path | None = None,
) -> RadonPyAssignment | EMCAssignment:
    """Assign parameters to one linear H-capped homopolymer with a reference backend.

    RadonPy uses depiction coordinates for graph-only assignment and a 3D
    conformer for quantum charges. EMC returns its own LAMMPS atom ordering.
    """
    if forcefield in ("opls", "oplsaa", "pcff"):
        if charge is not None:
            raise ValueError("EMC uses the charge policy of its selected force field")
        return assign_emc(psmiles, forcefield, dp=dp, seed=seed, emc_root=emc_root)
    if forcefield not in ("gaff", "gaff2"):
        raise ValueError(f"Unknown force field {forcefield!r}")
    if emc_root is not None:
        raise ValueError("emc_root only applies to OPLS-AA and PCFF")
    from island.builders import build_linear_polymer

    system = build_linear_polymer(
        psmiles,
        dp=dp,
        generate_3d=charge in ("RESP", "ESP", "Mulliken", "Lowdin"),
        head_end_group="*[H]",
        tail_end_group="*[H]",
        random_seed=seed,
    )
    return assign_radonpy(system, forcefield, charge=charge)
