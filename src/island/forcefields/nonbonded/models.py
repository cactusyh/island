"""Explicit immutable LJ mixing and bonded-neighbor scaling policies."""

from dataclasses import dataclass, field
from math import isfinite, sqrt
from typing import Literal

from island.core import Topology
from island.exceptions import NonbondedPolicyError
from island.forcefields.parameters.models import LennardJonesParameter

LJMixingRule = Literal["lorentz_berthelot", "geometric"]


@dataclass(frozen=True)
class MixedLJParameters:
    """One on-demand mixed LJ 12-6 type pair."""

    atom_types: tuple[str, str]
    epsilon: float
    sigma: float
    epsilon_unit: str
    sigma_unit: str
    mixing_rule: LJMixingRule
    policy_name: str
    policy_version: str


@dataclass(frozen=True)
class PairScaling:
    """Independent LJ and Coulomb weights for one stable unordered pair."""

    site_ids: tuple[int, int]
    shortest_bond_distance: int | None
    relationship: str
    lj_scale: float
    coulomb_scale: float


@dataclass(frozen=True)
class NonbondedPolicy:
    """Versioned explicit LJ 12-6 mixing and shortest-path scaling policy."""

    name: str
    version: str
    mixing_rule: LJMixingRule
    lj_scale_12: float
    coulomb_scale_12: float
    lj_scale_13: float
    coulomb_scale_13: float
    lj_scale_14: float
    coulomb_scale_14: float
    source: str
    functional_form: str = field(default="lj_12_6", init=False)
    epsilon_unit: str = field(default="kJ/mol", init=False)
    sigma_unit: str = field(default="nm", init=False)

    def __post_init__(self) -> None:
        for name in ("name", "version", "source"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise NonbondedPolicyError(f"Policy {name} must be non-empty")
        if self.mixing_rule not in {"lorentz_berthelot", "geometric"}:
            raise NonbondedPolicyError(
                f"Unsupported LJ mixing rule: {self.mixing_rule!r}"
            )
        for name in (
            "lj_scale_12",
            "coulomb_scale_12",
            "lj_scale_13",
            "coulomb_scale_13",
            "lj_scale_14",
            "coulomb_scale_14",
        ):
            value = getattr(self, name)
            if (
                not isinstance(value, (int, float))
                or isinstance(value, bool)
                or not isfinite(value)
                or not 0 <= value <= 1
            ):
                raise NonbondedPolicyError(
                    f"{name} must be a finite scale factor in [0, 1]"
                )
            object.__setattr__(self, name, float(value))

    @property
    def mixing_formula(self) -> str:
        if self.mixing_rule == "lorentz_berthelot":
            return "sigma_ij=(sigma_i+sigma_j)/2; epsilon_ij=sqrt(epsilon_i*epsilon_j)"
        return "sigma_ij=sqrt(sigma_i*sigma_j); epsilon_ij=sqrt(epsilon_i*epsilon_j)"

    def mix_lj(
        self,
        first: LennardJonesParameter,
        second: LennardJonesParameter,
    ) -> MixedLJParameters:
        """Mix two LJ 12-6 type records without materializing pair matrices."""
        if not isinstance(first, LennardJonesParameter) or not isinstance(
            second, LennardJonesParameter
        ):
            raise NonbondedPolicyError(
                "LJ mixing requires LennardJonesParameter records"
            )
        if (
            first.functional_form != self.functional_form
            or second.functional_form != self.functional_form
        ):
            raise NonbondedPolicyError("Only LJ 12-6 mixing is supported")
        if (
            first.epsilon_unit != self.epsilon_unit
            or second.epsilon_unit != self.epsilon_unit
            or first.sigma_unit != self.sigma_unit
            or second.sigma_unit != self.sigma_unit
        ):
            raise NonbondedPolicyError("LJ record units do not match policy units")
        epsilon = sqrt(first.epsilon * second.epsilon)
        sigma = (
            (first.sigma + second.sigma) / 2
            if self.mixing_rule == "lorentz_berthelot"
            else sqrt(first.sigma * second.sigma)
        )
        return MixedLJParameters(
            atom_types=tuple(sorted((first.atom_type, second.atom_type))),
            epsilon=epsilon,
            sigma=sigma,
            epsilon_unit=self.epsilon_unit,
            sigma_unit=self.sigma_unit,
            mixing_rule=self.mixing_rule,
            policy_name=self.name,
            policy_version=self.version,
        )

    def mixed_parameters_for_types(
        self, parameter_result: object, first_type: str, second_type: str
    ) -> MixedLJParameters:
        """Select unique assigned type records and mix the requested pair."""
        records: dict[str, set[LennardJonesParameter]] = {
            first_type: set(),
            second_type: set(),
        }
        for selection in getattr(parameter_result, "site_assignments", {}).values():
            parameter = selection.parameter
            if (
                isinstance(parameter, LennardJonesParameter)
                and selection.atom_types[0] in records
            ):
                records[selection.atom_types[0]].add(parameter)
        missing_or_ambiguous = {
            atom_type: len(items)
            for atom_type, items in records.items()
            if len(items) != 1
        }
        if missing_or_ambiguous:
            raise NonbondedPolicyError(
                "Requested atom types do not have one unique assigned LJ record: "
                f"{missing_or_ambiguous}"
            )
        return self.mix_lj(
            next(iter(records[first_type])), next(iter(records[second_type]))
        )

    def scaling_for_pair(
        self, topology: Topology, site1: int, site2: int
    ) -> PairScaling:
        """Query shortest bond-path scaling for one stable unordered site pair."""
        if site1 == site2:
            raise NonbondedPolicyError("A nonbonded pair requires two distinct sites")
        if site1 not in topology.sites or site2 not in topology.sites:
            raise NonbondedPolicyError("Nonbonded pair references an unknown site")
        distance = shortest_bond_distance(topology, site1, site2, maximum=3)
        pair = tuple(sorted((site1, site2)))
        if distance == 1:
            return PairScaling(
                pair, distance, "1-2", self.lj_scale_12, self.coulomb_scale_12
            )
        if distance == 2:
            return PairScaling(
                pair, distance, "1-3", self.lj_scale_13, self.coulomb_scale_13
            )
        if distance == 3:
            return PairScaling(
                pair, distance, "1-4", self.lj_scale_14, self.coulomb_scale_14
            )
        return PairScaling(pair, None, "full", 1.0, 1.0)

    def local_pair_scalings(
        self, topology: Topology
    ) -> dict[tuple[int, int], PairScaling]:
        """Return only local 1-2/1-3/1-4 exceptions with non-full weights."""
        topology.validate_bond_graph()
        result = {}
        for site1 in sorted(topology.sites):
            distances = _local_distances(topology, site1, maximum=3)
            for site2 in distances:
                if site1 >= site2:
                    continue
                scaling = self.scaling_for_pair(topology, site1, site2)
                if scaling.lj_scale != 1.0 or scaling.coulomb_scale != 1.0:
                    result[(site1, site2)] = scaling
        return result


def shortest_bond_distance(
    topology: Topology, site1: int, site2: int, *, maximum: int = 3
) -> int | None:
    """Return the shortest distance up to ``maximum`` bonds, otherwise ``None``."""
    return _local_distances(topology, site1, maximum=maximum).get(site2)


def _local_distances(
    topology: Topology, source: int, *, maximum: int
) -> dict[int, int]:
    distances = {source: 0}
    pending = [source]
    while pending:
        current = pending.pop(0)
        if distances[current] >= maximum:
            continue
        for neighbor in sorted(topology.neighbors(current)):
            candidate = distances[current] + 1
            if neighbor not in distances or candidate < distances[neighbor]:
                distances[neighbor] = candidate
                pending.append(neighbor)
    distances.pop(source, None)
    return distances
