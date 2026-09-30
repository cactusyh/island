"""Explicit thermal sampling: local PCG64 streams, all 3N Cartesian degrees of freedom."""

from dataclasses import dataclass, field
from math import isclose
from types import MappingProxyType

import numpy as np

from island.core import AtomSite, MolecularSystem
from island.evaluation.models import fingerprint
from island.exceptions import (
    DynamicsInputError,
    InvalidDynamicsResultError,
    IslandError,
)
from island.minimization.integrity import number
from island.minimization.models import system_identity

from .integrity import require
from .models import kinetic_energy, owned_vectors, velocity_hash

# Exact SI defining constants: N_A * k_B / 1000, kJ/(mol*K).
GAS_CONSTANT = 0.00831446261815324
RNG_ALGORITHM = "numpy.random.Generator(PCG64)/standard_normal"


def validate_temperature_seed(temperature, seed):
    if not number(temperature, nonnegative=True):
        raise DynamicsInputError("temperature_kelvin must be finite and nonnegative")
    if type(seed) is not int or not 0 <= seed < 2**128:
        raise DynamicsInputError("Seed must be an integer in [0, 2**128)")


def _rng(seed):
    return np.random.Generator(np.random.PCG64(seed))


def _normal_increments(rng, shape):
    # Private test seam, not a replay/file API. C order = ascending site, xyz.
    return rng.standard_normal(shape)


def mass_inventory(system):
    if not isinstance(system, MolecularSystem):
        raise DynamicsInputError("Expected MolecularSystem")
    system.validate()
    if (
        system.representation != "atomistic"
        or system.box is not None
        or not system.topology.sites
    ):
        raise DynamicsInputError("Require nonempty nonperiodic atomistic system")
    if any(type(site) is not int for site in system.topology.sites):
        raise DynamicsInputError("Stable IDs must be integers")
    masses = {}
    for site in sorted(system.topology.sites):
        atom = system.topology.sites[site]
        if type(atom) is not AtomSite or not number(atom.mass) or atom.mass <= 0:
            raise DynamicsInputError(
                "Require positive finite atomistic masses in dalton"
            )
        masses[site] = float(atom.mass)
    return masses


def mass_hash(masses):
    return fingerprint({"unit": "dalton", "sites": sorted(masses.items())})


def initialization_hash(result):
    return fingerprint(
        {
            "temperature_kelvin": result.temperature_kelvin,
            "velocity_seed": result.velocity_seed,
            "rng": result.rng_algorithm,
            "numpy": result.numpy_version,
            "masses": result.mass_fingerprint,
            "system": result.system_fingerprint,
            "velocities": result.velocity_fingerprint,
            "dof": result.degrees_of_freedom,
            "policy": "3N_no_removal_no_rescaling_v1",
        }
    )


@dataclass(frozen=True)
class VelocityInitialization:
    velocities: object
    masses: object
    temperature_kelvin: float
    velocity_seed: int
    numpy_version: str
    system_fingerprint: str
    mass_fingerprint: str
    velocity_fingerprint: str
    kinetic_energy: float
    degrees_of_freedom: int
    instantaneous_temperature_kelvin: float
    initialization_fingerprint: str
    rng_algorithm: str = field(default=RNG_ALGORITHM, init=False)
    velocity_unit: str = field(default="angstrom/ps", init=False)
    mass_unit: str = field(default="dalton", init=False)
    energy_unit: str = field(default="kJ/mol", init=False)
    temperature_unit: str = field(default="kelvin", init=False)

    def __post_init__(self):
        try:
            object.__setattr__(self, "velocities", owned_vectors(self.velocities))
            object.__setattr__(self, "masses", MappingProxyType(dict(self.masses)))
        except (TypeError, ValueError) as error:
            raise InvalidDynamicsResultError(str(error)) from error

    def validate_integrity(self):
        try:
            validate_temperature_seed(self.temperature_kelvin, self.velocity_seed)
            require(
                self.rng_algorithm == RNG_ALGORITHM
                and type(self.numpy_version) is str
                and bool(self.numpy_version),
                "Invalid RNG identity",
            )
            require(
                (
                    self.velocity_unit,
                    self.mass_unit,
                    self.energy_unit,
                    self.temperature_unit,
                )
                == ("angstrom/ps", "dalton", "kJ/mol", "kelvin"),
                "Invalid thermal units",
            )
            require(
                isinstance(self.masses, MappingProxyType)
                and isinstance(self.velocities, MappingProxyType),
                "Thermal state must be owned immutable mappings",
            )
            owned_vectors(self.velocities)
            require(
                set(self.masses) == set(self.velocities)
                and all(
                    type(s) is int and number(m) and m > 0
                    for s, m in self.masses.items()
                ),
                "Invalid mass coverage",
            )
            require(
                type(self.degrees_of_freedom) is int
                and self.degrees_of_freedom == 3 * len(self.masses),
                "DOF must be 3N",
            )
            require(
                number(self.kinetic_energy, nonnegative=True)
                and isclose(
                    self.kinetic_energy,
                    kinetic_energy(self.masses, self.velocities),
                    rel_tol=1e-12,
                    abs_tol=1e-10,
                ),
                "Invalid kinetic energy",
            )
            require(
                number(self.instantaneous_temperature_kelvin, nonnegative=True)
                and isclose(
                    self.instantaneous_temperature_kelvin,
                    2 * self.kinetic_energy / (self.degrees_of_freedom * GAS_CONSTANT),
                    rel_tol=1e-12,
                    abs_tol=1e-10,
                ),
                "Invalid instantaneous temperature",
            )
            require(
                self.temperature_kelvin != 0 or self.kinetic_energy == 0,
                "Zero-temperature sample must be zero",
            )
            require(
                type(self.system_fingerprint) is str and bool(self.system_fingerprint),
                "Missing system identity",
            )
            require(
                self.mass_fingerprint == mass_hash(self.masses)
                and self.velocity_fingerprint == velocity_hash(self.velocities)
                and self.initialization_fingerprint == initialization_hash(self),
                "Thermal fingerprint mismatch",
            )
        except (
            IslandError,
            AttributeError,
            TypeError,
            ValueError,
            KeyError,
            OverflowError,
            FloatingPointError,
        ) as error:
            raise InvalidDynamicsResultError(
                f"Malformed thermal initialization: {error}"
            ) from error

    def __deepcopy__(self, memo):
        self.validate_integrity()
        return self


def initialize_velocities(system, *, temperature_kelvin, seed):
    """Sample independent Maxwell-Boltzmann velocities; T=0 gives zero velocities.

    No COM/rotation removal or exact-temperature rescaling. A local PCG64 stream
    consumes 3N normals in ascending stable-ID/xyz order, including at T=0.
    """
    from dataclasses import replace

    try:
        validate_temperature_seed(temperature_kelvin, seed)
        masses = mass_inventory(system)
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            scale = np.sqrt(
                100
                * GAS_CONSTANT
                * temperature_kelvin
                / np.array(list(masses.values()))
            )
            values = _normal_increments(_rng(seed), (len(masses), 3)) * scale[:, None]
            speeds = owned_vectors(dict(zip(masses, values, strict=True)))
            kinetic = kinetic_energy(masses, speeds)
        result = VelocityInitialization(
            speeds,
            masses,
            float(temperature_kelvin),
            seed,
            np.__version__,
            system_identity(system),
            mass_hash(masses),
            velocity_hash(speeds),
            kinetic,
            3 * len(masses),
            2 * kinetic / (3 * len(masses) * GAS_CONSTANT),
            "pending",
        )
        result = replace(result, initialization_fingerprint=initialization_hash(result))
        result.validate_integrity()
        return result
    except (
        IslandError,
        AttributeError,
        TypeError,
        ValueError,
        KeyError,
        OverflowError,
        FloatingPointError,
    ) as error:
        raise DynamicsInputError(
            f"Invalid thermal initialization input: {error}"
        ) from error
