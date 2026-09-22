"""Chemical topology represented using stable site IDs."""

from dataclasses import dataclass
from itertools import combinations

from island.core.site import Site
from island.exceptions import TopologyError, ValidationError


@dataclass(frozen=True)
class Bond:
    site1: int
    site2: int
    order: float | None = None
    aromatic: bool = False

    def __post_init__(self) -> None:
        if self.site1 == self.site2:
            raise TopologyError("A bond cannot connect a site to itself")

    @property
    def key(self) -> tuple[int, int]:
        return tuple(sorted((self.site1, self.site2)))


@dataclass(frozen=True)
class Angle:
    site1: int
    site2: int
    site3: int

    @property
    def key(self) -> tuple[int, int, int]:
        reverse = (self.site3, self.site2, self.site1)
        return min((self.site1, self.site2, self.site3), reverse)


@dataclass(frozen=True)
class Dihedral:
    site1: int
    site2: int
    site3: int
    site4: int

    @property
    def key(self) -> tuple[int, int, int, int]:
        forward = (self.site1, self.site2, self.site3, self.site4)
        return min(forward, tuple(reversed(forward)))


@dataclass(frozen=True)
class Improper:
    site1: int
    site2: int
    site3: int
    site4: int


class Topology:
    """Molecular graph and graph-derived interactions."""

    def __init__(self) -> None:
        self.sites: dict[int, Site] = {}
        self.bonds: dict[tuple[int, int], Bond] = {}
        self.angles: dict[tuple[int, int, int], Angle] = {}
        self.dihedrals: dict[tuple[int, int, int, int], Dihedral] = {}
        self.impropers: list[Improper] = []

    def add_site(self, site: Site) -> None:
        if site.id in self.sites:
            raise TopologyError(f"Site id already exists: {site.id}")
        self.sites[site.id] = site

    def remove_site(self, site_id: int) -> Site:
        self._require_site(site_id)
        site = self.sites[site_id]
        for bond_key in [key for key in self.bonds if site_id in key]:
            del self.bonds[bond_key]
        del self.sites[site_id]
        self.rebuild_derived_interactions()
        self.impropers = [
            improper
            for improper in self.impropers
            if site_id not in self._members(improper)
        ]
        return site

    def get_site(self, site_id: int) -> Site:
        self._require_site(site_id)
        return self.sites[site_id]

    def add_bond(
        self,
        site1: int,
        site2: int,
        order: float | None = None,
        *,
        aromatic: bool = False,
    ) -> Bond:
        self._require_site(site1)
        self._require_site(site2)
        bond = Bond(site1, site2, order, aromatic)
        if bond.key in self.bonds:
            raise TopologyError(f"Duplicate bond: {bond.key}")
        self.bonds[bond.key] = bond
        return bond

    def remove_bond(self, site1: int, site2: int) -> Bond:
        key = tuple(sorted((site1, site2)))
        try:
            return self.bonds.pop(key)
        except KeyError as error:
            raise TopologyError(f"Bond does not exist: {key}") from error

    def neighbors(self, site_id: int) -> set[int]:
        self._require_site(site_id)
        return {
            b.site2 if b.site1 == site_id else b.site1
            for b in self.bonds.values()
            if site_id in b.key
        }

    def degree(self, site_id: int) -> int:
        return len(self.neighbors(site_id))

    def connected_components(self) -> list[set[int]]:
        remaining = set(self.sites)
        components: list[set[int]] = []
        while remaining:
            component: set[int] = set()
            pending = [next(iter(remaining))]
            while pending:
                current = pending.pop()
                if current in component:
                    continue
                component.add(current)
                pending.extend(self.neighbors(current) - component)
            remaining -= component
            components.append(component)
        return components

    def rebuild_derived_interactions(self) -> None:
        """Fully regenerate angles and proper dihedrals from the bond graph."""
        self.angles.clear()
        self.dihedrals.clear()
        for center in self.sites:
            for left, right in combinations(sorted(self.neighbors(center)), 2):
                angle = Angle(left, center, right)
                self.angles[angle.key] = angle
        for bond in self.bonds.values():
            for left in self.neighbors(bond.site1) - {bond.site2}:
                for right in self.neighbors(bond.site2) - {bond.site1}:
                    dihedral = Dihedral(left, bond.site1, bond.site2, right)
                    self.dihedrals[dihedral.key] = dihedral

    def validate(self) -> None:
        for bond in self.bonds.values():
            self._validate_members(bond)
        for interaction in [
            *self.angles.values(),
            *self.dihedrals.values(),
            *self.impropers,
        ]:
            self._validate_members(interaction)

    def _require_site(self, site_id: int) -> None:
        if site_id not in self.sites:
            raise TopologyError(f"Unknown site id: {site_id}")

    @staticmethod
    def _members(interaction: object) -> tuple[int, ...]:
        return tuple(
            getattr(interaction, field)
            for field in interaction.__dataclass_fields__
            if field.startswith("site")
        )

    def _validate_members(self, interaction: object) -> None:
        unknown = set(self._members(interaction)) - set(self.sites)
        if unknown:
            raise ValidationError(
                f"Interaction references unknown site IDs: {sorted(unknown)}"
            )
