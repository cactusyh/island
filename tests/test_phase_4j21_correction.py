"""Corrective J21 controls for density, persistence and periodic preparation."""

from dataclasses import replace

import pytest
from test_phase_4j21_packing import _chain, _pack_two

from island.exceptions import ValidationError
from island.forcefields import (
    AmberToolsOptions,
    ForceFieldRequest,
    OPLSOptions,
    PCFFOptions,
    PreparedForceFieldError,
    prepare_forcefield,
)
from island.graph.final import pack
from island.packing import (
    PeriodicPackingConfig,
    _digest,
    load_periodic_packing_plan,
    pack_multichain_periodic,
    save_periodic_packing_plan,
)


def _units():
    s0, g0, _ = _chain()
    s1, g1, _ = _chain(seed=2027)
    return [("chain-0", s0, g0), ("chain-1", s1, g1)]


def test_both_density_and_explicit_box_are_allowed_and_recorded():
    units = _units()
    mass = sum(
        site.mass for _, system, _ in units for site in system.topology.sites.values()
    )
    box = (40.0, 40.0, 40.0)
    density = mass / (6.02214076e23 * 64000.0e0 * 1e-24)
    result = pack_multichain_periodic(
        units,
        config=PeriodicPackingConfig(
            target_density=density,
            box_lengths=box,
            density_tolerance=1e-12,
        ),
        provenance="density correction",
        evidence=("controlled fixture",),
    )
    assert result.plan.payload["schema"] == "island_periodic_packing_plan_v2"
    assert result.plan.payload["target_density"] == density
    assert result.plan.payload["calculated_density"] == pytest.approx(density)
    assert result.plan.payload["density_tolerance"] == 1e-12
    assert result.system.box.lengths == box


def test_inconsistent_density_and_explicit_box_reject():
    with pytest.raises(ValidationError, match="density"):
        pack_multichain_periodic(
            _units(),
            config=PeriodicPackingConfig(
                target_density=0.9,
                box_lengths=(40.0, 40.0, 40.0),
                density_tolerance=1e-8,
            ),
        )


def test_density_tolerance_is_identity_bound():
    result = _pack_two()
    changed = replace(
        result.config, density_tolerance=result.config.density_tolerance * 2
    )
    assert changed.identity != result.config.identity
    payload = result.plan.payload
    payload["config"]["density_tolerance"] *= 2
    payload["identity"] = _digest(
        {key: value for key, value in payload.items() if key != "identity"}
    )
    with pytest.raises(ValidationError):
        result.plan.__class__(pack(payload)).validate_integrity(_units())


def test_expected_identity_rejects_resigned_plan(tmp_path):
    result = _pack_two()
    payload = result.plan.payload
    payload["provenance"] = "resigned"
    payload["identity"] = _digest(
        {key: value for key, value in payload.items() if key != "identity"}
    )
    path = tmp_path / "resigned.json"
    path.write_text(pack(payload))
    with pytest.raises(ValidationError, match="expected_identity"):
        load_periodic_packing_plan(path, expected_identity=result.plan.identity)


def test_transactional_save_preserves_existing_and_cleans_failed_publication(
    tmp_path, monkeypatch
):
    result = _pack_two()
    destination = tmp_path / "plan.json"
    destination.write_text("previous-valid-plan")
    with pytest.raises(ValidationError, match="destination exists"):
        save_periodic_packing_plan(result.plan, destination)
    assert destination.read_text() == "previous-valid-plan"

    destination.unlink()
    from island import packing

    def fail_replace(*args):
        raise OSError("injected publication failure")

    monkeypatch.setattr(packing.os, "replace", fail_replace)
    with pytest.raises(OSError, match="publication failure"):
        save_periodic_packing_plan(result.plan, destination)
    assert not destination.exists()
    assert not list(tmp_path.glob(".plan.json.*"))


@pytest.mark.parametrize(
    "family,options",
    [
        ("pcff", PCFFOptions("unused.frc", (0, 0, 1), (0, 0, 1))),
        ("oplsaa", OPLSOptions("unused.xml")),
        ("gaff", AmberToolsOptions("gaff", "provided", {i: 0.0 for i in range(1, 17)})),
        (
            "gaff2",
            AmberToolsOptions("gaff2", "provided", {i: 0.0 for i in range(1, 17)}),
        ),
    ],
)
def test_periodic_forcefield_preparation_rejects_before_backend(
    family, options, monkeypatch
):
    from island.forcefields import preparation

    result = _pack_two()
    called = []
    monkeypatch.setattr(preparation, "_description", lambda *args: called.append(args))
    with pytest.raises(
        PreparedForceFieldError, match="does not support periodic final graphs"
    ):
        prepare_forcefield(
            result.system,
            ForceFieldRequest(family, options, final_graph=result.graph),
        )
    assert called == []
