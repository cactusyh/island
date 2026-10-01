"""Explicit bounded short-chain configuration, independent of optional engines."""

from dataclasses import dataclass, field, fields
from math import isfinite
from types import MappingProxyType

from island.dynamics import LangevinOptions
from island.exceptions import WorkflowError
from island.forcefields.ambertools.models import AmberToolsOptions
from island.forcefields.ambertools.policy import DEFAULT_MAX_ATOMS
from island.minimization import MinimizationOptions


@dataclass(frozen=True)
class WorkflowConfig:
    psmiles: str
    dp: int
    output_directory: str
    force_field: str
    charge_method: str
    provided_charges: object = None
    charge_tolerance: float = 1e-4
    ambertools_timeout_seconds: float = 600.0
    amberhome: str | None = None
    tacticity: str | None = None
    atactic_fraction: float = 0.5
    template_seed: int = 2026
    assembly_seed: int = 2026
    stereo_seed: int = 2026
    velocity_seed: int = 78123
    thermostat_seed: int = 99181
    minimization: MinimizationOptions = field(default_factory=MinimizationOptions)
    temperature_kelvin: float = 300.0
    friction_per_ps: float = 5.0
    timestep_fs: float = 0.1
    total_steps: int = 100
    segment_steps: int = 20
    recording_interval: int = 10
    max_evaluations_per_segment: int = 22
    max_frames_per_segment: int = 3
    max_segments: int = 100
    max_artifact_bytes: int = 100_000_000
    schema: str = "island_single_chain_config_v2"
    max_atoms: int = DEFAULT_MAX_ATOMS
    charge_source: str | None = None

    def __post_init__(self):
        try:
            if self.schema not in ("island_single_chain_config_v1", "island_single_chain_config_v2"):
                raise ValueError("Unsupported workflow config schema")
            if self.schema == "island_single_chain_config_v1" and (self.max_atoms != 100 or self.charge_source is not None):
                raise ValueError("Historical v1 configs retain max_atoms=100 and no charge_source field")
            for name in ("psmiles", "output_directory"):
                if type(getattr(self, name)) is not str or not getattr(self, name):
                    raise ValueError(f"{name} must be a nonempty string")
            for name in (
                "dp",
                "total_steps",
                "segment_steps",
                "max_segments",
                "max_artifact_bytes",
            ):
                if type(getattr(self, name)) is not int or getattr(self, name) <= 0:
                    raise ValueError(f"{name} must be a positive integer")
            for name in ("template_seed", "assembly_seed", "stereo_seed"):
                if (
                    type(getattr(self, name)) is not int
                    or not 0 <= getattr(self, name) < 2**31 - 100
                ):
                    raise ValueError(f"{name} must fit the builder's signed seed range")
            if (
                type(self.velocity_seed) is not int
                or not 0 <= self.velocity_seed < 2**128
            ):
                raise ValueError("Invalid velocity seed")
            if self.tacticity not in (None, "isotactic", "syndiotactic", "atactic"):
                raise ValueError("Unsupported tacticity")
            if (
                type(self.atactic_fraction) not in (int, float)
                or not isfinite(self.atactic_fraction)
                or not 0 <= self.atactic_fraction <= 1
            ):
                raise ValueError("Invalid atactic fraction")
            if type(self.minimization) is not MinimizationOptions:
                raise ValueError("Expected MinimizationOptions")
            self.minimization.__post_init__()
            if self.amberhome is not None and (
                type(self.amberhome) is not str or not self.amberhome
            ):
                raise ValueError("amberhome must be a path string or null")
            self.langevin(self.total_steps)
            options = self.amber_options()
            if options.provided_charges is not None:
                charges = dict(options.provided_charges)
                if any(
                    type(k) is not int or type(v) not in (int, float) or not isfinite(v)
                    for k, v in charges.items()
                ):
                    raise ValueError(
                        "Provided charges require integer IDs and finite values"
                    )
                object.__setattr__(self, "provided_charges", MappingProxyType(charges))
        except Exception as error:
            raise WorkflowError(f"Invalid workflow configuration: {error}") from error

    def amber_options(self, work_root=None):
        return AmberToolsOptions(
            self.force_field,
            self.charge_method,
            self.provided_charges,
            self.ambertools_timeout_seconds,
            self.charge_tolerance,
            work_root,
            self.amberhome,
            True,
            max_atoms=self.max_atoms, charge_source=self.charge_source,
        )

    def langevin(self, steps):
        return LangevinOptions(
            self.timestep_fs,
            self.temperature_kelvin,
            self.friction_per_ps,
            self.thermostat_seed,
            steps,
            self.max_evaluations_per_segment,
            self.max_frames_per_segment,
            self.recording_interval,
        )

    def to_dict(self):
        result = {f.name: getattr(self, f.name) for f in fields(self)}
        result["minimization"] = {
            f.name: getattr(self.minimization, f.name)
            for f in fields(self.minimization)
        }
        if self.provided_charges is not None:
            result["provided_charges"] = {
                str(k): v for k, v in self.provided_charges.items()
            }
        if self.schema == "island_single_chain_config_v1":
            result.pop("max_atoms")
            result.pop("charge_source")
        return result

    @classmethod
    def from_dict(cls, value):
        try:
            value = dict(value)
            if value.get("schema") == "island_single_chain_config_v1" and any(k in value for k in ("max_atoms", "charge_source")):
                raise ValueError("Size/source fields require workflow config v2")
            if "minimization" in value:
                value["minimization"] = MinimizationOptions(**value["minimization"])
            if value.get("provided_charges") is not None:
                charges = value["provided_charges"]
                if any(type(k) is not str or str(int(k)) != k for k in charges):
                    raise ValueError("Noncanonical charge ID")
                value["provided_charges"] = {int(k): v for k, v in charges.items()}
            return cls(**value)
        except Exception as error:
            raise WorkflowError(f"Malformed workflow configuration: {error}") from error
