"""Boundary loading is structural; applying/resuming also verifies coordinate stereo."""

import pytest

pytest.importorskip("rdkit")
from test_dynamics import Harmonic
from test_dynamics_stereo import chiral
from test_langevin import options

from island.core.coordinate_provenance import coordinate_hash
from island.dynamics import (
    DynamicsCheckpoint,
    DynamicsSegmentOptions,
    create_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
)
from island.dynamics import _checkpoint_data as data
from island.exceptions import DynamicsInputError


def test_assigned_checkpoint_stereo_and_inverted_boundary():
    system = chiral()
    velocities = {s: (0.0, 0.0, 0.0) for s in system.topology.sites}
    first = run_dynamics_segment(
        system, Harmonic(0), velocities, options(2, temperature_kelvin=0)
    )
    cp = create_dynamics_checkpoint(first)
    assert resume_dynamics(
        cp, system, Harmonic(0), DynamicsSegmentOptions(1, 3, 2)
    ).completed
    p = cp.payload
    for _, x in p["state"]["coordinates"]:
        x[0] *= -1
    fp = coordinate_hash(data.inventory(p["state"]["coordinates"]))
    p["state"]["coordinate_fingerprint"] = fp
    p["state"]["evaluation"]["coordinate_fingerprint"] = fp
    p["final_evaluation"]["coordinate_fingerprint"] = fp
    altered = DynamicsCheckpoint(data.canonical(p), data.checksum(p))
    with pytest.raises(DynamicsInputError):
        resume_dynamics(altered, system, Harmonic(0), DynamicsSegmentOptions(1, 3, 2))
