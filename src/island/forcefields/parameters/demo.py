"""Synthetic software-test parameters for Phase 4B demonstrations only."""

from island.forcefields.parameters.models import (
    HarmonicAngleParameter,
    HarmonicBondParameter,
    LennardJonesParameter,
    ParameterLibrary,
    PeriodicTorsionTerm,
    ProperTorsionParameter,
)
from island.forcefields.typing import island_demo_v1_ruleset
from island.forcefields.typing.signatures import ruleset_signature

_NAME = "island_demo_parameters_v1"
_VERSION = "1.0.0"
_SOURCE = "SYNTHETIC SOFTWARE-TEST PARAMETERS — NOT FOR SCIENTIFIC SIMULATION."


def island_demo_parameters_v1() -> ParameterLibrary:
    """Return synthetic PE-coverage records for testing assignment machinery."""
    carbon = "demo_c_aliphatic"
    hydrogen = "demo_h_on_carbon"
    common = {"source": _SOURCE, "library_name": _NAME, "library_version": _VERSION}
    ruleset = island_demo_v1_ruleset()
    return ParameterLibrary(
        name=_NAME,
        version=_VERSION,
        required_ruleset_name=ruleset.name,
        required_ruleset_version=ruleset.version,
        required_ruleset_signature=ruleset_signature(ruleset),
        records=(
            LennardJonesParameter("demo_lj_c", carbon, 0.30, 0.34, **common),
            LennardJonesParameter("demo_lj_hc", hydrogen, 0.10, 0.25, **common),
            HarmonicBondParameter(
                "demo_bond_cc", (carbon, carbon), 250000.0, 0.154, **common
            ),
            HarmonicBondParameter(
                "demo_bond_ch", (carbon, hydrogen), 280000.0, 0.109, **common
            ),
            HarmonicAngleParameter(
                "demo_angle_ccc",
                (carbon, carbon, carbon),
                300.0,
                112.0,
                **common,
            ),
            HarmonicAngleParameter(
                "demo_angle_cch",
                (carbon, carbon, hydrogen),
                250.0,
                110.0,
                **common,
            ),
            HarmonicAngleParameter(
                "demo_angle_hch",
                (hydrogen, carbon, hydrogen),
                220.0,
                109.5,
                **common,
            ),
            ProperTorsionParameter(
                "demo_torsion_cccc",
                (carbon, carbon, carbon, carbon),
                (
                    PeriodicTorsionTerm(0.80, 1, 0.0),
                    PeriodicTorsionTerm(0.35, 3, 180.0),
                ),
                **common,
            ),
            ProperTorsionParameter(
                "demo_torsion_hccc",
                (hydrogen, carbon, carbon, carbon),
                (PeriodicTorsionTerm(0.25, 3, 0.0),),
                **common,
            ),
            ProperTorsionParameter(
                "demo_torsion_hcch",
                (hydrogen, carbon, carbon, hydrogen),
                (PeriodicTorsionTerm(0.20, 3, 0.0),),
                **common,
            ),
        ),
        description=(
            "Synthetic exact-type records covering explicit-hydrogen polyethylene; "
            "not a production force field."
        ),
        source=_SOURCE,
    )
