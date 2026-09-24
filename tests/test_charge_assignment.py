from copy import deepcopy
from dataclasses import replace
from math import inf, nan

import pytest

from island import AtomSite, Coordinates, MolecularSystem, Topology
from island.exceptions import (
    IncompleteChargeAssignmentError,
    InvalidChargeAssignmentResultError,
    InvalidChargeDefinitionError,
)
from island.forcefields import ProvidedChargeEngine


def system_with_atoms(
    atoms: tuple[tuple[int, str, int], ...],
    bonds: tuple[tuple[int, int], ...] = (),
) -> MolecularSystem:
    topology = Topology()
    coordinates = Coordinates()
    atomic_numbers = {"C": 6, "N": 7, "Na": 11, "Cl": 17}
    masses = {"C": 12.0, "N": 14.0, "Na": 23.0, "Cl": 35.5}
    for index, (site_id, element, formal_charge) in enumerate(atoms):
        topology.add_site(
            AtomSite(
                site_id,
                f"{element}{site_id}",
                masses[element],
                element=element,
                atomic_number=atomic_numbers[element],
                formal_charge=formal_charge,
            )
        )
        coordinates.set(site_id, [float(index), 0.0, 0.0])
    for site1, site2 in bonds:
        topology.add_bond(site1, site2, order=1)
    return MolecularSystem(topology, coordinates)


def test_provided_nonzero_charges_on_neutral_molecule() -> None:
    system = system_with_atoms(((10, "C", 0), (30, "C", 0)), ((10, 30),))
    result = ProvidedChargeEngine().assign(
        system,
        {10: 0.25, 30: -0.25},
        source="synthetic neutral test",
    )
    assert result.complete
    assert result.assignments[10].charge == 0.25
    assert result.assignments[30].charge == -0.25
    assert result.total_assigned_charge == 0.0
    assert result.total_formal_charge == 0.0
    assert result.unit == "elementary_charge"


def test_charged_molecule_remains_nonzero() -> None:
    system = system_with_atoms(((50, "N", 1),))
    result = ProvidedChargeEngine().assign(
        system, {50: 1.0}, source="synthetic cation test", target_charge=1
    )
    assert result.complete
    assert result.total_assigned_charge == 1.0
    assert result.component_diagnostics[0].expected_charge == 1.0


def test_disconnected_components_prevent_charge_error_cancellation() -> None:
    ions = system_with_atoms(((10, "Na", 1), (90, "Cl", -1)))
    correct = ProvidedChargeEngine().assign(
        ions, {10: 1.0, 90: -1.0}, source="synthetic ion pair"
    )
    assert correct.complete
    assert [item.expected_charge for item in correct.component_diagnostics] == [
        1.0,
        -1.0,
    ]

    with pytest.raises(IncompleteChargeAssignmentError) as captured:
        ProvidedChargeEngine().assign(
            ions, {10: 0.0, 90: 0.0}, source="bad cancelling ion pair"
        )
    result = captured.value.result
    assert result.total_assigned_charge == result.total_formal_charge == 0.0
    assert not result.complete
    assert sum(item.reason == "charge_mismatch" for item in result.diagnostics) == 2


def test_whole_system_tolerance_is_checked_after_each_component() -> None:
    system = system_with_atoms(((10, "C", 0), (30, "C", 0)))
    with pytest.raises(IncompleteChargeAssignmentError) as captured:
        ProvidedChargeEngine().assign(
            system,
            {10: 0.8e-6, 30: 0.8e-6},
            source="small component residuals",
            tolerance=1.0e-6,
        )
    result = captured.value.result
    assert all(item.within_tolerance for item in result.component_diagnostics)
    assert result.total_charge_residual == pytest.approx(1.6e-6)
    assert not result.total_within_tolerance
    assert not result.complete


@pytest.mark.parametrize(
    "charges",
    [
        {10: 0.0},
        {10: 0.0, 30: 0.0, 99: 0.0},
        {10: nan, 30: 0.0},
        {10: inf, 30: 0.0},
        {10: True, 30: 0.0},
    ],
)
def test_provided_charge_mapping_rejects_missing_extra_and_nonfinite_values(
    charges: dict[int, float],
) -> None:
    system = system_with_atoms(((10, "C", 0), (30, "C", 0)), ((10, 30),))
    with pytest.raises(InvalidChargeDefinitionError):
        ProvidedChargeEngine().assign(system, charges, source="invalid fixture")


def test_explicit_target_must_match_authoritative_formal_charge() -> None:
    system = system_with_atoms(((50, "N", 1),))
    with pytest.raises(InvalidChargeDefinitionError, match="target charge"):
        ProvidedChargeEngine().assign(
            system, {50: 1.0}, source="synthetic", target_charge=0.0
        )


def test_coordinate_changes_preserve_charge_compatibility() -> None:
    system = system_with_atoms(((10, "C", 0), (30, "C", 0)), ((10, 30),))
    charges = {10: 0.1, 30: -0.1}
    engine = ProvidedChargeEngine()
    result = engine.assign(system, charges, source="coordinate-independent")
    moved = system.copy()
    moved.coordinates.translate([100.0, -20.0, 3.0])
    assert result.is_compatible_with(moved)
    assert result.input_signature == engine.input_signature(
        moved, charges, source="coordinate-independent"
    )


def test_charge_copying_metadata_isolation_and_malformed_result_rejection() -> None:
    system = system_with_atoms(((10, "C", 0), (30, "C", 0)), ((10, 30),))
    result = ProvidedChargeEngine().assign(
        system, {10: 0.1, 30: -0.1}, source="copy test"
    )
    copied = deepcopy(result)
    copied.validate_integrity(system)
    with pytest.raises(TypeError):
        copied.assignments[99] = copied.assignments[10]  # type: ignore[index]

    metadata = {"nested": {"values": [1]}}
    reconstructed = replace(result, metadata=metadata)
    metadata["nested"]["values"].append(2)
    assert reconstructed.metadata["nested"]["values"] == [1]

    assignments = dict(result.assignments)
    assignments[10] = replace(assignments[10], charge=0.2)
    malformed = replace(result, assignments=assignments, metadata=dict(result.metadata))
    before = system.to_dict()
    assert not malformed.is_compatible_with(system)
    with pytest.raises(InvalidChargeAssignmentResultError):
        malformed.validate_integrity(system)
    assert system.to_dict() == before


def test_recomputed_output_signature_cannot_hide_changed_provided_charge_input() -> None:
    from island.forcefields.charges.signatures import charge_result_signature

    system = system_with_atoms(((10, "C", 0), (30, "C", 0)), ((10, 30),))
    result = ProvidedChargeEngine().assign(
        system, {10: 0.1, 30: -0.1}, source="signature test"
    )
    assignments = dict(result.assignments)
    assignments[10] = replace(assignments[10], charge=0.2)
    assignments[30] = replace(assignments[30], charge=-0.2)
    changed = replace(result, assignments=assignments)
    changed = replace(changed, result_signature=charge_result_signature(changed))
    with pytest.raises(
        InvalidChargeAssignmentResultError, match="provided charge input signature"
    ):
        changed.validate_integrity(system)


def test_graph_or_charge_values_change_reuse_signatures() -> None:
    system = system_with_atoms(((10, "C", 0), (30, "C", 0)), ((10, 30),))
    engine = ProvidedChargeEngine()
    first = engine.assign(system, {10: 0.1, 30: -0.1}, source="signature test")
    second = engine.assign(system, {10: 0.2, 30: -0.2}, source="signature test")
    assert first.input_signature != second.input_signature
    assert first.result_signature != second.result_signature
    changed_input = engine.input_signature(
        system, {10: 0.2, 30: -0.2}, source="signature test"
    )
    assert not first.is_input_compatible_with(system, input_signature=changed_input)

    changed_graph = system.copy()
    changed_graph.topology.sites[10].formal_charge = 1
    assert not first.is_input_compatible_with(changed_graph)


def test_charge_signatures_ignore_site_and_bond_insertion_order() -> None:
    first = system_with_atoms(((10, "C", 0), (30, "C", 0)), ((10, 30),))
    second = system_with_atoms(((30, "C", 0), (10, "C", 0)), ((30, 10),))
    engine = ProvidedChargeEngine()
    first_result = engine.assign(first, {10: 0.1, 30: -0.1}, source="order test")
    second_result = engine.assign(second, {30: -0.1, 10: 0.1}, source="order test")
    assert first_result.input_signature == second_result.input_signature
    assert first_result.result_signature == second_result.result_signature
