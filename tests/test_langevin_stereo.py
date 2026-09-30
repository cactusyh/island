"""Stereo is validated on all trials, not inferred from stored chiral tags."""

import builtins
from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("rdkit")
from test_dynamics import Harmonic
from test_langevin import options

from island.chemistry import from_smiles
from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import run_langevin
from island.exceptions import DynamicsUnavailableError, InvalidDynamicsResultError


def chiral():
    return from_smiles("F[C@H](Cl)Br", random_seed=2026)


def test_preserved_stereo_and_rejected_unsaved_trial():
    system = chiral()
    velocities = {site: (0.2, 0.1, -0.1) for site in system.topology.sites}
    good = run_langevin(
        system,
        Harmonic(0),
        velocities,
        options(5, timestep_fs=0.1, recording_interval=5),
    )
    assert good.completed and good.stereochemistry == "passed"
    good.to_system(system).validate()
    # One force-free step would reflect all x coordinates across the yz plane.
    reflected = {
        site: (-2 * system.coordinates.get(site)[0] / 0.001, 0, 0)
        for site in system.topology.sites
    }
    failed = run_langevin(
        system,
        Harmonic(0),
        reflected,
        options(5, timestep_fs=1, friction_per_ps=0, recording_interval=5),
    )
    assert failed.termination_reason == "stereochemistry_changed"
    assert failed.completed_steps == 0 and failed.attempted_step == 1
    assert failed.evaluations == 2 and failed.failed_trial_evaluations == 0
    assert failed.final_state.coordinates == failed.initial_state.coordinates
    assert failed.final_state.velocities == failed.initial_state.velocities
    failed.to_system(system, allow_incomplete=True).validate()


@pytest.mark.parametrize("scale", [(-1, 1, 1), (1, 1, 0)])
def test_reconstructed_invalid_stereo_rejected(scale):
    system = chiral()
    original = deepcopy(system.to_dict())
    result = run_langevin(
        system,
        Harmonic(0),
        {site: (0.1, 0.1, 0.1) for site in system.topology.sites},
        options(1),
    )
    last = result.final_state
    coordinates = {
        site: tuple(np.array(xyz) * scale) for site, xyz in last.coordinates.items()
    }
    checksum = coordinate_hash(coordinates)
    altered = replace(
        last,
        coordinates=coordinates,
        coordinate_fingerprint=checksum,
        evaluation=replace(last.evaluation, coordinate_fingerprint=checksum),
    )
    bad = replace(
        result,
        frames=(result.initial_state, altered),
        final_evaluation=replace(
            result.final_evaluation, coordinate_fingerprint=checksum
        ),
    )
    bad.validate_integrity()  # numerical consistency cannot establish coordinate stereo
    with pytest.raises(InvalidDynamicsResultError):
        bad.to_system(system, allow_incomplete=True)
    assert system.to_dict() == original


def test_required_rdkit_cannot_be_silently_skipped(monkeypatch):
    system = chiral()
    velocities = {site: (0.1, 0, 0) for site in system.topology.sites}
    result = run_langevin(system, Harmonic(0), velocities, options(1))
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "rdkit":
            raise ImportError("blocked dependency")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(DynamicsUnavailableError):
        run_langevin(system, Harmonic(0), velocities, options(1))
    with pytest.raises(DynamicsUnavailableError):
        result.to_system(system)
