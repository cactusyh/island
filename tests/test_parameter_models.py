from dataclasses import replace
from math import nan

import pytest

from island.exceptions import InvalidParameterDefinitionError
from island.forcefields.parameters import (
    HarmonicAngleParameter,
    HarmonicBondParameter,
    LennardJonesParameter,
    ParameterLibrary,
    PeriodicTorsionTerm,
    ProperTorsionParameter,
    island_demo_parameters_v1,
    parameter_library_signature,
)

COMMON = {
    "source": "synthetic test source",
    "library_name": "test_library",
    "library_version": "1",
}


def test_supported_parameter_models_expose_forms_units_and_provenance() -> None:
    lj = LennardJonesParameter("lj", "A", 0.2, 0.3, **COMMON)
    bond = HarmonicBondParameter("bond", ("A", "B"), 10.0, 0.15, **COMMON)
    angle = HarmonicAngleParameter("angle", ("A", "B", "C"), 20.0, 109.5, **COMMON)
    torsion = ProperTorsionParameter(
        "torsion",
        ("A", "B", "C", "D"),
        (
            PeriodicTorsionTerm(1.0, 1, 0.0),
            PeriodicTorsionTerm(0.5, 3, 180.0),
        ),
        **COMMON,
    )
    assert (lj.family, lj.functional_form) == ("site", "lj_12_6")
    assert (bond.force_constant_unit, bond.length_unit) == (
        "kJ/(mol*nm^2)",
        "nm",
    )
    assert (angle.force_constant_unit, angle.angle_unit) == (
        "kJ/(mol*rad^2)",
        "degree",
    )
    assert torsion.functional_form == "periodic_torsion"
    assert len(torsion.terms) == 2
    assert torsion.terms[1].phase_unit == "degree"
    assert torsion.source == COMMON["source"]


@pytest.mark.parametrize(
    "factory",
    [
        lambda: LennardJonesParameter("bad", "A", -0.1, 0.3, **COMMON),
        lambda: LennardJonesParameter("bad", "A", 0.2, nan, **COMMON),
        lambda: LennardJonesParameter(
            "bad", "A", 0.2, 0.3, epsilon_unit="kcal/mol", **COMMON
        ),
        lambda: HarmonicBondParameter("bad", ("A", "B"), -1.0, 0.15, **COMMON),
        lambda: HarmonicBondParameter(
            "bad", ("A", "B"), 1.0, 0.15, length_unit="angstrom", **COMMON
        ),
        lambda: HarmonicAngleParameter("bad", ("A", "B", "C"), 1.0, 181.0, **COMMON),
        lambda: HarmonicAngleParameter(
            "bad", ("A", "B", "C"), 1.0, 109.5, angle_unit="radian", **COMMON
        ),
        lambda: PeriodicTorsionTerm(1.0, 0, 0.0),
        lambda: PeriodicTorsionTerm(-0.1, 3, 0.0),
        lambda: PeriodicTorsionTerm(1.0, 3, 360.0),
        lambda: PeriodicTorsionTerm(1.0, 3, 0.0, phase_unit="radian"),
        lambda: ProperTorsionParameter("bad", ("A", "B", "C", "D"), (), **COMMON),
        lambda: HarmonicBondParameter("bad", ("A", "*"), 1.0, 0.1, **COMMON),
    ],
)
def test_invalid_values_units_and_unsupported_patterns_are_rejected(
    factory: object,
) -> None:
    with pytest.raises(InvalidParameterDefinitionError):
        factory()


def test_explicit_zero_lj_and_torsion_amplitudes_are_valid_records() -> None:
    lj = LennardJonesParameter("zero_lj", "A", 0.0, 0.3, **COMMON)
    term = PeriodicTorsionTerm(0.0, 3, 180.0)
    torsion = ProperTorsionParameter(
        "zero_torsion", ("A", "B", "C", "D"), (term,), **COMMON
    )
    assert lj.epsilon == 0.0
    assert torsion.terms == (term,)


def test_library_rejects_duplicates_identity_mismatch_and_unsupported_records() -> None:
    record = LennardJonesParameter("lj", "A", 0.2, 0.3, **COMMON)
    kwargs = {
        "name": "test_library",
        "version": "1",
        "required_ruleset_name": "rules",
        "required_ruleset_version": "1",
        "required_ruleset_signature": "abc",
        "description": "synthetic test library",
        "source": "synthetic test source",
    }
    with pytest.raises(InvalidParameterDefinitionError, match="Duplicate"):
        ParameterLibrary(records=(record, record), **kwargs)
    mismatch = replace(record, library_name="other")
    with pytest.raises(InvalidParameterDefinitionError, match="declares library"):
        ParameterLibrary(records=(mismatch,), **kwargs)
    with pytest.raises(InvalidParameterDefinitionError, match="unsupported"):
        ParameterLibrary(records=(object(),), **kwargs)  # type: ignore[arg-type]


def test_library_signature_covers_values_and_is_record_order_independent() -> None:
    library = island_demo_parameters_v1()
    reversed_library = replace(library, records=tuple(reversed(library.records)))
    assert parameter_library_signature(library) == parameter_library_signature(
        reversed_library
    )
    changed_record = replace(library.records[0], epsilon=0.31)
    changed_library = replace(library, records=(changed_record, *library.records[1:]))
    assert parameter_library_signature(library) != parameter_library_signature(
        changed_library
    )
