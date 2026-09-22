"""Container for force-field assignments kept apart from molecular structure."""

from dataclasses import dataclass, field
from typing import Any

from island.core.system import MolecularSystem


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
