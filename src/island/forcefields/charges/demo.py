"""Synthetic charge-table fixtures; not scientific charge models."""

from island.forcefields.charges.models import (
    AtomTypeChargeEntry,
    AtomTypeChargeTable,
)
from island.forcefields.typing import island_demo_v1_ruleset
from island.forcefields.typing.signatures import ruleset_signature

_NAME = "island_demo_charges_v1"
_VERSION = "1.0.0"
_SOURCE = "SYNTHETIC SOFTWARE-TEST CHARGES — NOT FOR SCIENTIFIC SIMULATION."


def island_demo_charges_v1() -> AtomTypeChargeTable:
    """Return an explicit zero-charge PE table for workflow demonstrations."""
    ruleset = island_demo_v1_ruleset()
    common = {
        "source": _SOURCE,
        "table_name": _NAME,
        "table_version": _VERSION,
    }
    return AtomTypeChargeTable(
        name=_NAME,
        version=_VERSION,
        required_ruleset_name=ruleset.name,
        required_ruleset_version=ruleset.version,
        required_ruleset_signature=ruleset_signature(ruleset),
        entries=(
            AtomTypeChargeEntry(
                "demo_charge_c_aliphatic",
                "demo_c_aliphatic",
                0.0,
                **common,
            ),
            AtomTypeChargeEntry(
                "demo_charge_h_on_carbon",
                "demo_h_on_carbon",
                0.0,
                **common,
            ),
        ),
        source=_SOURCE,
    )
