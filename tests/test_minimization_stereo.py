"""Coordinate stereo checks cannot be satisfied by retained chiral tags alone."""

import builtins
from dataclasses import replace

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("rdkit")
from test_minimization import Harmonic

from island import Coordinates
from island.chemistry import from_smiles
from island.exceptions import MinimizationInputError, MinimizationUnavailableError
from island.minimization import MinimizationOptions, minimize_geometry
from island.minimization.models import coordinate_hash


class Target(Harmonic):
    def __init__(self, target):
        super().__init__()
        self.target = target

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        shifted = {
            site: np.array(xyz) - self.target[site] for site, xyz in coordinates.items()
        }
        result = super().evaluate(shifted, coordinate_unit=coordinate_unit)
        return replace(result, coordinate_fingerprint=coordinate_hash(coordinates))


def test_stereo_preserved_and_inverted_trial_rejected():
    system = from_smiles("F[C@H](Cl)Br", random_seed=2026)
    xyz = {
        site: np.array(system.coordinates.get(site)) for site in system.topology.sites
    }
    translated = {site: value + 0.02 for site, value in xyz.items()}
    good = minimize_geometry(
        system, Target(translated), MinimizationOptions(force_tolerance=1e-6)
    )
    assert good.converged and good.stereochemistry == "passed"
    reflected = {site: value * (-1, 1, 1) for site, value in xyz.items()}
    bad = minimize_geometry(
        system, Target(reflected), MinimizationOptions(force_tolerance=1e-6)
    )
    assert not bad.converged
    assert bad.termination_reason == "stereochemistry_changed"
    assert (
        bad.stereochemistry == "passed"
    )  # diagnostic coordinates themselves remain valid
    with pytest.raises(MinimizationInputError):
        bad.to_system(system)
    bad.to_system(system, allow_unconverged=True)
    system.coordinates = Coordinates(reflected)
    with pytest.raises(MinimizationInputError, match="stereochemistry expected"):
        minimize_geometry(system, Target(reflected))
    system.coordinates = Coordinates(
        {site: (value[0], value[1], 0) for site, value in xyz.items()}
    )
    with pytest.raises(MinimizationInputError, match="planar/degenerate"):
        minimize_geometry(system, Target(xyz))


def test_assigned_stereo_cannot_skip_missing_dependency(monkeypatch):
    system = from_smiles("F[C@H](Cl)Br", random_seed=2026)
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "rdkit":
            raise ImportError("blocked")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(MinimizationUnavailableError, match="RDKit"):
        minimize_geometry(system, Harmonic())


@pytest.mark.parametrize("mode", ["inverted", "degenerate"])
@pytest.mark.parametrize("diagnostic", [False, True])
def test_reconstructed_stereo_claim_not_trusted_at_application(mode, diagnostic):
    from copy import deepcopy

    from island.exceptions import InvalidMinimizationResultError

    system = from_smiles("F[C@H](Cl)Br", random_seed=2026)
    before = deepcopy(system.to_dict())
    target = {
        site: np.array(system.coordinates.get(site)) + 0.02
        for site in system.topology.sites
    }
    result = minimize_geometry(system, Target(target))
    coordinates = {
        site: tuple(np.array(xyz) * ((-1, 1, 1) if mode == "inverted" else (1, 1, 0)))
        for site, xyz in result.coordinates.items()
    }
    history = [
        *result.history[:-1],
        replace(
            result.history[-1], coordinate_fingerprint=coordinate_hash(coordinates)
        ),
    ]
    bad = replace(
        result,
        coordinates=coordinates,
        history=history,
        final_evaluation=replace(
            result.final_evaluation, coordinate_fingerprint=coordinate_hash(coordinates)
        ),
        converged=not diagnostic,
        termination_reason="energy_stagnation" if diagnostic else "force_converged",
    )
    bad.validate_integrity()  # content alone cannot establish geometric stereo
    with pytest.raises(
        InvalidMinimizationResultError, match="stereochemistry|degenerate"
    ):
        bad.to_system(system, allow_unconverged=diagnostic)
    assert system.to_dict() == before


def test_missing_rdkit_at_application(monkeypatch):
    system = from_smiles("F[C@H](Cl)Br", random_seed=2026)
    target = {
        site: np.array(system.coordinates.get(site)) for site in system.topology.sites
    }
    result = minimize_geometry(system, Target(target))
    original = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name == "rdkit":
            raise ImportError("blocked")
        return original(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", blocked)
    with pytest.raises(MinimizationUnavailableError, match="RDKit"):
        result.to_system(system)
