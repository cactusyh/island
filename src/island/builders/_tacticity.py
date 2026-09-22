"""Internal tacticity preparation kept separate from polymer graph assembly."""

from collections.abc import Mapping

from island.builders.sequence import PolymerSequence
from island.chemistry.psmiles import RepeatUnit
from island.chemistry.stereochemistry import (
    StereochemicalSequence,
    generate_stereochemical_sequence,
    require_controllable_stereocenter,
)
from island.exceptions import UnsupportedStereochemistryError


def prepare_stereochemical_sequence(
    library: Mapping[str, RepeatUnit],
    polymer_sequence: PolymerSequence,
    *,
    tacticity: str | None,
    stereo_seed: int,
    atactic_fraction: float,
) -> StereochemicalSequence | None:
    """Validate Phase 3.6A scope and generate a homopolymer stereo sequence."""
    if tacticity is None:
        return None
    used_types = set(polymer_sequence.identities)
    if len(used_types) != 1:
        raise UnsupportedStereochemistryError(
            "Phase 3.6A tacticity control is limited to a single repeat-unit "
            "type; mixed repeat-unit sequences are unsupported"
        )
    repeat_type = next(iter(used_types))
    center = require_controllable_stereocenter(library[repeat_type])
    if center.cip_label not in {"R", "S"}:
        raise UnsupportedStereochemistryError(
            "The controllable repeat-unit center lacks an assigned R/S state"
        )
    return generate_stereochemical_sequence(
        polymer_sequence.number_of_repeat_units,
        tacticity,
        reference_state=center.cip_label,
        stereo_seed=stereo_seed,
        atactic_fraction=atactic_fraction,
    )
