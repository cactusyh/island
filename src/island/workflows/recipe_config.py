"""Validated JSON inputs for periodic polymer recipe runs."""

from dataclasses import asdict, dataclass, field
from math import isfinite

from island.forcefields.emc import FIELDS
from island.forcefields.radonpy import FAMILIES, QM_CHARGES

from .recipe_equilibration import RecipeEquilibration


def positive(value, name, *, integer=False):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not isfinite(value)
        or value <= 0
        or (integer and type(value) is not int)
    ):
        raise ValueError(
            f"{name} must be a positive {'integer' if integer else 'number'}"
        )


@dataclass(frozen=True)
class RecipeForceField:
    name: str
    charge: str | None = None
    profile: str | None = None
    field_dir: str | None = None

    def __post_init__(self):
        if self.name not in (*FIELDS, *FAMILIES):
            raise ValueError(f"Unknown force field {self.name}")
        if self.name in FAMILIES:
            if self.charge not in ("gasteiger", "zero", *QM_CHARGES):
                raise ValueError("GAFF recipes require an explicit charge method")
            if self.profile or self.field_dir:
                raise ValueError("EMC profiles and field_dir apply to PCFF/OPLS")
        elif self.charge is not None:
            raise ValueError(
                "PCFF/OPLS use charges from their selected parameter library"
            )
        expected = {
            "pcff-carboxamide-v1": "pcff",
            "pcff-polymer-v2": "pcff",
            "opls-epoxy-amine-v1": "opls/2024/opls-aa",
            "opls-polymer-v1": "opls/2024/opls-aa",
            "opls-polymer-v2": "opls/2024/opls-aa",
            "opls-polymer-v3": "opls/2024/opls-aa",
        }
        if self.profile and (
            self.field_dir or expected.get(self.profile) != FIELDS.get(self.name)
        ):
            raise ValueError(
                "Profile must match the field and requires its own installed directory"
            )


@dataclass(frozen=True)
class RecipeMD:
    temperature_k: float = 300.0
    pressure_atm: float = 1.0
    timestep_fs: float = 0.5
    nvt_steps: int = 1000
    npt_steps: int = 4000
    production_steps: int = 4000
    production_ensemble: str = "npt"
    sample_every: int = 100
    segment_steps: int = 2000

    def __post_init__(self):
        for name in ("temperature_k", "timestep_fs"):
            positive(getattr(self, name), name)
        if type(self.pressure_atm) not in (int, float) or not isfinite(
            self.pressure_atm
        ):
            raise ValueError("pressure_atm must be finite")
        for name in (
            "nvt_steps",
            "npt_steps",
            "production_steps",
            "sample_every",
            "segment_steps",
        ):
            positive(getattr(self, name), name, integer=True)
        if self.production_ensemble not in ("nvt", "npt"):
            raise ValueError("production_ensemble must be nvt or npt")


@dataclass(frozen=True)
class PolymerRecipe:
    architecture: str
    forcefield: RecipeForceField
    polymer: dict | None = None
    components: list[dict] | None = field(default=None, kw_only=True)
    cure: dict | None = None
    chains: int = 1
    density_g_cm3: float = 0.5
    seed: int = 2026
    md: RecipeMD = field(default_factory=RecipeMD)
    equilibration: RecipeEquilibration | None = field(default=None, kw_only=True)
    schema: int = 1

    def __post_init__(self):
        if (
            self.equilibration
            and self.equilibration.steps_per_segment % self.md.sample_every
        ):
            raise ValueError("Equilibration segment must align with sample_every")
        if self.schema != 1 or type(self.schema) is not int:
            raise ValueError("Unsupported recipe schema")
        if self.architecture not in (
            "homopolymer",
            "block",
            "branched",
            "copolymer",
            "mixture",
            "crosslinked",
        ):
            raise ValueError("Unknown polymer architecture")
        positive(self.chains, "chains", integer=True)
        positive(self.seed, "seed", integer=True)
        positive(self.density_g_cm3, "density_g_cm3")
        if self.seed > 900_000_000:
            raise ValueError("seed exceeds the LAMMPS velocity seed range")
        if self.architecture == "mixture":
            if (
                self.polymer is not None
                or self.cure is not None
                or self.components is None
                or self.chains != 1
            ):
                raise ValueError(
                    "Mixtures use components and per-component chain counts"
                )
            from .recipe_chains import expand

            expand(self)
            return
        if self.components is not None:
            raise ValueError("Components require mixture architecture")
        if self.architecture == "crosslinked":
            if self.polymer is not None or self.cure is None or self.chains != 1:
                raise ValueError("Crosslinked recipes specify component counts in cure")
            if self.forcefield.name not in ("gaff", "gaff2") and FIELDS.get(
                self.forcefield.name
            ) not in ("pcff", "opls/2024/opls-aa"):
                raise ValueError(
                    "Periodic reaction templates currently support PCFF and 2024 OPLS-AA"
                )
            if {"seed", "packing_density"} & self.cure.keys():
                raise ValueError("Specify seed and density_g_cm3 at recipe level")
            cfg = self.cure_config()
            from island.epoxy.chemistry import DGEBA

            if self.forcefield.name in ("gaff", "gaff2"):
                if self.forcefield.charge != "gasteiger":
                    raise ValueError("Periodic GAFF cure requires Gasteiger charges")
            elif (
                cfg.epoxy_smiles != DGEBA
                or cfg.hardener_smiles != cfg.template.hardener
            ):
                raise ValueError("Custom reaction monomers require GAFF/GAFF2")
            if (
                cfg.petmp
                and self.forcefield.name not in ("gaff", "gaff2")
                and FIELDS.get(self.forcefield.name) != "pcff"
            ):
                raise ValueError("PETMP requires the PCFF reaction library")
        elif self.polymer is None or self.cure is not None:
            raise ValueError("Finite-chain recipes require polymer and omit cure")
        else:
            self.polymer_spec()

    @classmethod
    def from_dict(cls, value):
        value = dict(value)
        value["forcefield"] = RecipeForceField(**value["forcefield"])
        value["md"] = RecipeMD(**value.get("md", {}))
        if value.get("equilibration") is not None:
            value["equilibration"] = RecipeEquilibration(**value["equilibration"])
        return cls(**value)

    def polymer_spec(self):
        from island.forcefields.assignment_core.schema import Branch, PolymerSpec

        data = dict(self.polymer)
        repeats = data.pop("repeat_units")
        dp = data.pop("dp", None)
        blocks = data.pop("blocks", None)
        sequence = data.pop("sequence", None)
        if sum(v is not None for v in (dp, blocks, sequence)) != 1:
            raise ValueError("Specify exactly one of dp, blocks or sequence")
        if blocks is not None:
            sequence = []
            for label, count in blocks:
                positive(count, "block length", integer=True)
                sequence.extend([label] * count)
        if sequence is not None:
            dp = len(sequence)
        branches = tuple(Branch(**b) for b in data.pop("branches", []))
        spec = PolymerSpec(
            tuple(repeats.items()),
            dp,
            sequence=tuple(sequence or ()),
            branches=branches,
            stereo_seed=self.seed,
            **data,
        )
        if spec.crosslinks or spec.periodic_topology:
            raise ValueError("Use a crosslinked recipe for periodic reaction chemistry")
        if self.architecture == "homopolymer" and (len(repeats) != 1 or branches):
            raise ValueError("Homopolymer requires one repeat unit and no branches")
        if self.architecture == "block" and (len(set(spec.sequence)) < 2 or branches):
            raise ValueError(
                "Block recipe requires at least two repeat identities and no branches"
            )
        if self.architecture == "copolymer" and (len(repeats) < 2 or branches):
            raise ValueError(
                "Copolymer requires multiple repeat identities and no branches"
            )
        if self.architecture == "branched" and not branches:
            raise ValueError("Branched recipe requires attachment sites")
        return spec

    def cure_config(self):
        from island.epoxy.cure import CureConfig

        return CureConfig(
            **self.cure, seed=self.seed, packing_density=self.density_g_cm3
        )

    def payload(self):
        return asdict(self)
