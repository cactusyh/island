"""NoCutoff OPLS single points from independently resolved ISLAND records."""

from collections.abc import Mapping
from copy import deepcopy
from itertools import combinations
from math import radians, sqrt

import numpy as np

from island.charge_references.records import pack
from island.exceptions import EvaluationError, EvaluationInputError
from island.forcefields.oplsaa.models import chemical_graph
from island.forcefields.oplsaa.parameters import OPLSParameterizationResult

from .models import EvaluationResult, fingerprint
from .openmm import OpenMMBoundPotential, _openmm

COMPONENTS = ("bond", "angle", "rb_proper", "coulomb", "lj")
SETTINGS = {
    "method": "NoCutoff",
    "periodic": False,
    "constraints": False,
    "switching": False,
    "long_range_correction": False,
    "precision": "double",
    "integration_steps": 0,
    "mixing": "geometric",
    "force_sign": "-dE/dR",
    "implementation": "island_opls_singlepoint_v1",
}


def build_model(data, mm):
    """Only ISLAND numerical records enter this construction path."""
    model = mm.System()
    ids = sorted(data["sites"])
    index = {s: i for i, s in enumerate(ids)}
    for s in ids:
        model.addParticle(data["sites"][s]["mass_dalton"])
    bond, angle, rb = (
        mm.HarmonicBondForce(),
        mm.HarmonicAngleForce(),
        mm.RBTorsionForce(),
    )
    for row in data["bonds"]:
        r = row["converted"]
        bond.addBond(
            *(index[s] for s in row["sites"]),
            r["length_angstrom"] * 0.1,
            r["k_kj_mol_angstrom2"] * 100,
        )
    for row in data["angles"]:
        r = row["converted"]
        angle.addAngle(
            *(index[s] for s in row["sites"]),
            radians(r["angle_degrees"]),
            r["k_kj_mol_radian2"],
        )
    for row in data["rb"]:
        rb.addTorsion(
            *(index[s] for s in row["sites"]), *row["converted"]["coefficients_kj_mol"]
        )
    coulomb = mm.NonbondedForce()
    coulomb.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
    coulomb.setUseDispersionCorrection(False)
    lj = mm.CustomNonbondedForce(
        "4*sqrt(epsilon1*epsilon2)*(x^12-x^6);x=sqrt(sigma1*sigma2)/r"
    )
    lj.addPerParticleParameter("sigma")
    lj.addPerParticleParameter("epsilon")
    lj.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
    lj.setUseLongRangeCorrection(False)
    lj.setUseSwitchingFunction(False)
    scaled = mm.CustomBondForce("4*epsilon*((sigma/r)^12-(sigma/r)^6)")
    scaled.addPerBondParameter("sigma")
    scaled.addPerBondParameter("epsilon")
    for s in ids:
        row = data["sites"][s]
        coulomb.addParticle(row["charge_e"], 1.0, 0.0)
        lj.addParticle([row["sigma_angstrom"] * 0.1, row["epsilon_kj_mol"]])
    for pair in data["pairs"]:
        a, b = pair["sites"]
        i, j = index[a], index[b]
        r, t = data["sites"][a], data["sites"][b]
        coulomb.addException(
            i, j, r["charge_e"] * t["charge_e"] * pair["coulomb_scale"], 1.0, 0.0
        )
        lj.addExclusion(i, j)
        if pair["lj_scale"]:
            scaled.addBond(
                i,
                j,
                [
                    0.1 * sqrt(r["sigma_angstrom"] * t["sigma_angstrom"]),
                    pair["lj_scale"] * sqrt(r["epsilon_kj_mol"] * t["epsilon_kj_mol"]),
                ],
            )
    for group, force in enumerate((bond, angle, rb, coulomb, lj)):
        force.setForceGroup(group)
        model.addForce(force)
    scaled.setForceGroup(4)
    model.addForce(scaled)
    return model


class OPLSSinglePointEvaluator(OpenMMBoundPotential):
    """Owned Reference-platform OPLS single points; no dynamics/session API."""

    def __init__(self, system, parameters, source):
        if not isinstance(parameters, OPLSParameterizationResult):
            raise EvaluationInputError("Require validated OPLSParameterizationResult")
        parameters.validate_integrity(system, source)
        self._system = deepcopy(system)
        self._parameters = deepcopy(parameters)
        self._source = deepcopy(source)
        self._data = parameters.payload["resolved"]
        self._ids = sorted(system.topology.sites)
        self._index = {s: i for i, s in enumerate(self._ids)}
        self._parameter_fingerprint = parameters.identity
        self._model_fingerprint = fingerprint(
            {"parameters": parameters.identity, "settings": SETTINGS}
        )
        mm, _ = _openmm()
        self._xml = mm.XmlSerializer.serialize(build_model(self._data, mm))

    @property
    def parameter_fingerprint(self):
        return self._parameter_fingerprint

    @property
    def model_fingerprint(self):
        return self._model_fingerprint

    def validate_system(self, system):
        if chemical_graph(system) != chemical_graph(self._system):
            raise EvaluationInputError("OPLS bound graph/stereo/mass changed; rebind")

    @classmethod
    def from_parameterized_system(cls, snapshot, source):
        try:
            result = OPLSParameterizationResult(
                pack(snapshot.metadata["opls_parameters_v1"])
            )
            result.validate_integrity(snapshot.system, source)
            expected = result.to_parameterized_system(snapshot.system, source)
            for attr in (
                "backend_name",
                "site_assignments",
                "interaction_assignments",
                "charge_assignments",
                "nonbonded_policy",
                "aggregate_signature",
                "metadata",
            ):
                if getattr(expected, attr) != getattr(snapshot, attr):
                    raise EvaluationInputError(f"Changed snapshot content: {attr}")
            return cls(snapshot.system, result, source)
        except EvaluationError:
            raise
        except Exception as error:
            raise EvaluationInputError(f"Invalid OPLS snapshot: {error}") from error

    def evaluate_fresh(self, coordinates=None, *, coordinate_unit="angstrom"):
        """Independent verification: evaluate always constructs a new Context."""
        return self.evaluate(coordinates, coordinate_unit=coordinate_unit)

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        if coordinate_unit != "angstrom":
            raise EvaluationInputError("Expected angstrom coordinates")
        try:
            if coordinates is None:
                coordinates = {s: self._system.coordinates.get(s) for s in self._ids}
            if (
                not isinstance(coordinates, Mapping)
                or set(coordinates) != set(self._ids)
                or any(type(s) is not int for s in coordinates)
            ):
                raise ValueError("Coordinate coverage mismatch")
            if any(
                isinstance(v, (bool, np.bool_))
                or not isinstance(v, (int, float, np.integer, np.floating))
                for s in self._ids
                for v in coordinates[s]
            ):
                raise ValueError("Coordinates must be nonboolean real numbers")
            xyz = np.array([coordinates[s] for s in self._ids], dtype=float)
            if xyz.shape != (len(self._ids), 3) or not np.isfinite(xyz).all():
                raise ValueError("Require finite N x 3 coordinates")
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                for i, j in combinations(range(len(xyz)), 2):
                    if np.linalg.norm(xyz[i] - xyz[j]) < 1e-10:
                        raise ValueError("Coincident atoms")
                for row in self._data["angles"]:
                    a, b, c = (xyz[self._index[s]] for s in row["sites"])
                    if np.linalg.norm(np.cross(a - b, c - b)) <= 1e-12 * np.linalg.norm(
                        a - b
                    ) * np.linalg.norm(c - b):
                        raise ValueError("Collinear angular geometry")
        except Exception as error:
            raise EvaluationInputError(f"Invalid geometry: {error}") from error
        mm, unit = _openmm()
        context = integrator = model = None
        try:
            model = mm.XmlSerializer.deserialize(self._xml)
            integrator = mm.VerletIntegrator(0.001)
            context = mm.Context(
                model, integrator, mm.Platform.getPlatformByName("Reference")
            )
            context.setPositions(xyz * 0.1 * unit.nanometer)
            state = context.getState(getEnergy=True, getForces=True)
            energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
            force = (
                state.getForces(asNumpy=True).value_in_unit(
                    unit.kilojoule_per_mole / unit.nanometer
                )
                * 0.1
            )
            components = {
                name: context.getState(getEnergy=True, groups=1 << i)
                .getPotentialEnergy()
                .value_in_unit(unit.kilojoule_per_mole)
                for i, name in enumerate(COMPONENTS)
            }
            if (
                not np.isfinite(force).all()
                or not np.isfinite(energy)
                or not all(np.isfinite(v) for v in components.values())
            ):
                raise EvaluationError("Nonfinite OPLS backend output")
            coordinate_id = fingerprint(
                {
                    "unit": "angstrom",
                    "sites": list(zip(self._ids, xyz.tolist(), strict=True)),
                }
            )
            evaluation_id = fingerprint(
                {
                    "model": self.model_fingerprint,
                    "coordinates": coordinate_id,
                    "backend": mm.__version__,
                    "platform": "Reference",
                }
            )
            return EvaluationResult(
                energy,
                components,
                dict(zip(self._ids, map(tuple, force.tolist()), strict=True)),
                coordinate_id,
                self.parameter_fingerprint,
                self.model_fingerprint,
                evaluation_id,
                "OpenMM OPLS-AA",
                mm.__version__,
                "Reference",
                SETTINGS,
            )
        except EvaluationError:
            raise
        except Exception as error:
            raise EvaluationError(f"OPLS evaluation failed: {error}") from error
        finally:
            del context, integrator, model
