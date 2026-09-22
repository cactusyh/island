"""Application of graph-level reaction transformations."""

from island.core.system import MolecularSystem
from island.exceptions import ReactionError, TopologyError
from island.reactions.template import BreakBond, ChangeBondOrder, FormBond


def apply_transformation(
    system: MolecularSystem,
    transformation: FormBond | BreakBond | ChangeBondOrder,
    *,
    rebuild_derived: bool = True,
) -> None:
    """Apply one topology-only transformation to a system in place."""
    try:
        if isinstance(transformation, FormBond):
            system.topology.add_bond(
                transformation.site1, transformation.site2, transformation.order
            )
        elif isinstance(transformation, BreakBond):
            system.topology.remove_bond(transformation.site1, transformation.site2)
        elif isinstance(transformation, ChangeBondOrder):
            system.topology.remove_bond(transformation.site1, transformation.site2)
            system.topology.add_bond(
                transformation.site1, transformation.site2, transformation.order
            )
        else:
            raise ReactionError(f"Unsupported transformation: {type(transformation)!r}")
    except TopologyError as error:
        raise ReactionError(str(error)) from error
    if rebuild_derived:
        system.topology.rebuild_derived_interactions()
