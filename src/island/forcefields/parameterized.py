"""Container for force-field assignments kept apart from molecular structure."""

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, cast

from island.core.system import MolecularSystem

if TYPE_CHECKING:
    from island.forcefields.charges.models import ChargeAssignment
    from island.forcefields.nonbonded.models import NonbondedPolicy


@dataclass
class ParameterizedSystem:
    """A molecular system paired with force-field assignments.

    Per-site assignments remain external to ``AtomSite``. Future backends may
    store atom types, partial charges, and nonbonded parameters in
    ``site_assignments`` without mutating chemical identity.
    """

    system: MolecularSystem
    backend_name: str | None = None
    site_assignments: dict[int, Any] = field(default_factory=dict)
    interaction_assignments: dict[str, dict[tuple[int, ...], Any]] = field(
        default_factory=dict
    )
    metadata: dict[str, Any] = field(default_factory=dict)
    charge_assignments: dict[int, "ChargeAssignment"] = field(default_factory=dict)
    nonbonded_policy: "NonbondedPolicy | None" = None
    aggregate_signature: str | None = None

    @classmethod
    def from_assignment(
        cls, system: MolecularSystem, assignment_result: object
    ) -> "ParameterizedSystem":
        """Create an owned system/assignment snapshot from a complete result.

        The input system is copied. Parameter selections are immutable and their
        mappings are copied, so later mutations of caller-owned objects cannot
        silently alter this validated snapshot.
        """
        from island.forcefields.parameters.models import ParameterAssignmentResult

        if not isinstance(assignment_result, ParameterAssignmentResult):
            raise TypeError("assignment_result must be a ParameterAssignmentResult")
        return cast(
            "ParameterizedSystem", assignment_result.to_parameterized_system(system)
        )

    @classmethod
    def from_components(
        cls,
        system: MolecularSystem,
        parameter_result: object,
        charge_result: object,
        nonbonded_policy: object,
    ) -> "ParameterizedSystem":
        """Compose an owned snapshot from validated independent components."""
        from island.forcefields.composition import compose_parameterized_system

        return cast(
            "ParameterizedSystem",
            compose_parameterized_system(
                system, parameter_result, charge_result, nonbonded_policy
            ),
        )
