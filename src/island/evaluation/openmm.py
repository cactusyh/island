"""Restricted Amber single-point adapter constructed from ISLAND records only."""

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import asdict, replace
from itertools import combinations
from math import radians

import numpy as np

from island.core import AtomSite, MolecularSystem
from island.exceptions import (
    EvaluationError,
    EvaluationInputError,
    EvaluationUnavailableError,
    IslandError,
    UnsupportedEvaluationError,
)
from island.forcefields.amber import ImportedAmberResult
from island.forcefields.ambertools.models import AmberToolsPreparationResult
from island.forcefields.parameterized import ParameterizedSystem
from island.forcefields.typing.signatures import graph_signature

from .models import EvaluationResult, fingerprint

COMPONENTS = ("bond", "angle", "proper_torsion", "periodic_improper", "nonbonded")
SETTINGS = {
    "nonbonded_method": "NoCutoff",
    "periodic": False,
    "constraints": False,
    "implicit_solvent": False,
    "switching": False,
    "dispersion_correction": False,
    "hydrogen_mass_repartitioning": False,
    "restraints": False,
    "integration_steps": 0,
    "force_sign": "-dE/dR",
}


def _openmm():
    try:
        import openmm
        from openmm import unit
    except ImportError as error:
        raise EvaluationUnavailableError(
            "OpenMM single-point evaluation requires the optional dependency: "
            "pip install 'island[evaluation]'"
        ) from error
    return openmm, unit


def _validate(system: MolecularSystem, imported: ImportedAmberResult) -> None:
    if not isinstance(system, MolecularSystem) or not isinstance(
        imported, ImportedAmberResult
    ):
        raise EvaluationInputError("Bind a MolecularSystem and an ImportedAmberResult")
    if system.representation != "atomistic" or not system.topology.sites:
        raise UnsupportedEvaluationError(
            "Only nonempty atomistic systems are supported"
        )
    if system.box is not None and any(system.box.periodic):
        raise UnsupportedEvaluationError(
            "Periodic evaluation is unsupported; use a nonperiodic system"
        )
    try:
        system.validate()
        imported.validate_integrity(system)
        for site in system.topology.sites.values():
            if (
                type(site) is not AtomSite
                or (site.atomic_number is None and not site.element)
                or (
                    site.atomic_number is not None
                    and (type(site.atomic_number) is not int or site.atomic_number <= 0)
                )
                or not np.isfinite(site.mass)
                or site.mass <= 0
            ):
                raise UnsupportedEvaluationError(
                    "Virtual sites and nonphysical atom masses are unsupported"
                )
        for name in (
            "site_assignments",
            "bond_assignments",
            "angle_assignments",
            "proper_torsion_assignments",
            "improper_assignments",
        ):
            for selection in getattr(imported, name).values():
                parameter = selection.parameter
                # Reconstruct to validate numeric values/units even if a caller
                # bypassed frozen fields and recomputed the result's hash.
                rebuilt = replace(parameter)
                if asdict(rebuilt) != asdict(parameter):
                    raise EvaluationInputError(
                        "Unsupported parameter functional form or family"
                    )
                for term in getattr(parameter, "terms", ()):
                    replace(term)
        policy = imported.nonbonded_policy
        rebuilt_policy = replace(policy)
        if asdict(rebuilt_policy) != asdict(policy):
            raise UnsupportedEvaluationError(
                "Unsupported nonbonded policy units or functional form"
            )
        if policy.mixing_rule != "lorentz_berthelot":
            raise UnsupportedEvaluationError(
                "Only Lorentz-Berthelot LJ mixing is supported"
            )
        if any(
            getattr(policy, field) != 0
            for field in (
                "lj_scale_12",
                "coulomb_scale_12",
                "lj_scale_13",
                "coulomb_scale_13",
            )
        ):
            raise UnsupportedEvaluationError(
                "Amber 1-2 and 1-3 exclusions must have zero scales"
            )
        for value in (system, imported, policy):
            if getattr(value, "virtual_sites", None) or getattr(
                value, "pair_overrides", None
            ):
                raise UnsupportedEvaluationError(
                    "Virtual sites and pair overrides are unsupported"
                )
        if len(set(imported.source_exclusions)) != len(
            imported.source_exclusions
        ) or len(set(imported.source_14_pairs)) != len(imported.source_14_pairs):
            raise EvaluationInputError("Duplicate source nonbonded pairs")
    except EvaluationError:
        raise
    except (IslandError, TypeError, ValueError, AttributeError, KeyError) as error:
        raise EvaluationInputError(
            f"Invalid authoritative system or Amber parameters: {error}"
        ) from error


def _graph(system: MolecularSystem) -> str:
    return fingerprint(
        {
            "chemical_graph": graph_signature(system.topology),
            "sites": [
                (
                    i,
                    system.topology.sites[i].element,
                    system.topology.sites[i].mass,
                    system.topology.sites[i].metadata.get("cip_label"),
                    system.topology.sites[i].metadata.get("chiral_tag"),
                )
                for i in sorted(system.topology.sites)
            ],
            "impropers": sorted(
                (i.site1, i.site2, i.site3, i.site4) for i in system.topology.impropers
            ),
        }
    )


class OpenMMSinglePointEvaluator:
    """Owned bound Amber potential. No file parsing, optimization or MD steps.

    ``evaluate`` accepts an explicitly angstrom-labelled stable-ID mapping, or
    uses the binding coordinates. ``evaluate_system`` additionally checks the
    authoritative graph of a new system. Backend contexts are local to each
    call, so frames never share mutable OpenMM state.
    """

    def __init__(
        self,
        system: MolecularSystem,
        imported_result: ImportedAmberResult,
        *,
        platform: str = "Reference",
        precision: str | None = None,
        periodic: bool = False,
    ) -> None:
        if periodic:
            raise UnsupportedEvaluationError(
                "Only nonperiodic NoCutoff evaluation is supported"
            )
        _validate(system, imported_result)
        self._system = deepcopy(system)
        self._imported = deepcopy(imported_result)
        self._graph_fingerprint = _graph(self._system)
        self._ids = tuple(sorted(self._system.topology.sites))
        self._index = {site: i for i, site in enumerate(self._ids)}
        self._parameter_fingerprint = imported_result.content_signature()
        self._model_fingerprint = fingerprint(
            {
                "schema": "island_openmm_singlepoint_v1",
                "graph": self._graph_fingerprint,
                "parameters": self._parameter_fingerprint,
                "settings": SETTINGS,
            }
        )
        mm, _unit = _openmm()
        try:
            selected = mm.Platform.getPlatformByName(platform)
            properties = set(selected.getPropertyNames())
            self._properties = {}
            if precision is not None:
                if "Precision" not in properties:
                    raise UnsupportedEvaluationError(
                        f"{platform} does not expose a Precision setting"
                    )
                if precision not in ("single", "mixed", "double"):
                    raise EvaluationInputError(
                        "precision must be single, mixed or double"
                    )
                self._properties["Precision"] = precision
            self._platform = selected.getName()
            # Keep no mutable OpenMM System exposed or shared with callers.
            self._xml = mm.XmlSerializer.serialize(self._build(mm))
        except EvaluationError:
            raise
        except Exception as error:
            raise EvaluationUnavailableError(
                f"Cannot construct OpenMM {platform} backend: {error}"
            ) from error
        self._coordinates(None, "angstrom")

    @classmethod
    def from_preparation(
        cls,
        preparation_system: MolecularSystem,
        result: AmberToolsPreparationResult,
        **settings,
    ) -> "OpenMMSinglePointEvaluator":
        if not isinstance(result, AmberToolsPreparationResult):
            raise EvaluationInputError("Expected an AmberToolsPreparationResult")
        result.validate_integrity(preparation_system)
        return cls(preparation_system, result.imported_result, **settings)

    @classmethod
    def from_parameterized_system(
        cls,
        snapshot: ParameterizedSystem,
        imported_result: ImportedAmberResult,
        **settings,
    ) -> "OpenMMSinglePointEvaluator":
        """Require the signed import and compare actual mutable snapshot contents."""
        if not isinstance(snapshot, ParameterizedSystem):
            raise EvaluationInputError("Expected an Amber-backed ParameterizedSystem")
        _validate(snapshot.system, imported_result)
        expected = imported_result.to_parameterized_system(snapshot.system)
        for field in (
            "backend_name",
            "site_assignments",
            "interaction_assignments",
            "charge_assignments",
            "nonbonded_policy",
            "aggregate_signature",
        ):
            if getattr(snapshot, field) != getattr(expected, field):
                raise EvaluationInputError(
                    f"Snapshot {field} contradicts validated Amber import"
                )
        return cls(snapshot.system, imported_result, **settings)

    @property
    def model_fingerprint(self) -> str:
        return self._model_fingerprint

    @property
    def parameter_fingerprint(self) -> str:
        return self._parameter_fingerprint

    def _build(self, mm):
        resolved = self._imported
        system = mm.System()
        forces = [
            mm.HarmonicBondForce(),
            mm.HarmonicAngleForce(),
            mm.PeriodicTorsionForce(),
            mm.PeriodicTorsionForce(),
            mm.NonbondedForce(),
        ]
        for group, force in enumerate(forces):
            force.setForceGroup(group)
            force.setName(COMPONENTS[group])
            system.addForce(force)
        bond, angle, proper, improper, nonbonded = forces
        nonbonded.setNonbondedMethod(mm.NonbondedForce.NoCutoff)
        nonbonded.setUseSwitchingFunction(False)
        nonbonded.setUseDispersionCorrection(False)
        for site in self._ids:
            system.addParticle(self._system.topology.sites[site].mass)
            lj = resolved.site_assignments[site].parameter
            nonbonded.addParticle(
                resolved.charge_result.assignments[site].charge, lj.sigma, lj.epsilon
            )
        for sites, selection in sorted(resolved.bond_assignments.items()):
            p = selection.parameter
            bond.addBond(
                *(self._index[i] for i in sites), p.equilibrium_length, p.force_constant
            )
        for sites, selection in sorted(resolved.angle_assignments.items()):
            p = selection.parameter
            angle.addAngle(
                *(self._index[i] for i in sites),
                radians(p.equilibrium_angle),
                p.force_constant,
            )
        for force, assignments in (
            (proper, resolved.proper_torsion_assignments),
            (improper, resolved.improper_assignments),
        ):
            for sites, selection in sorted(assignments.items()):
                # Never canonicalize the ordered improper or move its center.
                for term in selection.parameter.terms:
                    force.addTorsion(
                        *(self._index[i] for i in sites),
                        term.periodicity,
                        radians(term.phase),
                        term.force_constant,
                    )
        policy = resolved.nonbonded_policy
        pairs14 = set(resolved.source_14_pairs)
        for left, right in sorted(set(resolved.source_exclusions)):
            lj = policy.mix_lj(
                resolved.site_assignments[left].parameter,
                resolved.site_assignments[right].parameter,
            )
            scaled = (left, right) in pairs14
            charge_product = (
                resolved.charge_result.assignments[left].charge
                * resolved.charge_result.assignments[right].charge
                * policy.coulomb_scale_14
                if scaled
                else 0
            )
            epsilon = lj.epsilon * policy.lj_scale_14 if scaled else 0
            nonbonded.addException(
                self._index[left], self._index[right], charge_product, lj.sigma, epsilon
            )
        return system

    def _coordinates(self, coordinates, coordinate_unit):
        if not isinstance(coordinate_unit, str) or coordinate_unit != "angstrom":
            raise EvaluationInputError("Public coordinates must be labelled 'angstrom'")
        if coordinates is None:
            coordinates = {
                site: self._system.coordinates.get(site) for site in self._ids
            }
        if not isinstance(coordinates, Mapping):
            raise EvaluationInputError(
                "Coordinates must map stable site IDs to angstrom 3-vectors"
            )
        if any(type(site) is not int for site in coordinates) or set(
            coordinates
        ) != set(self._ids):
            raise EvaluationInputError("Coordinates require exact stable-site coverage")
        try:
            raw = np.asarray([coordinates[site] for site in self._ids])
            if raw.dtype.kind not in "iuf":
                raise ValueError("components must be real numbers")
            xyz = np.array(raw, dtype=float, copy=True)
        except (TypeError, ValueError, OverflowError) as error:
            raise EvaluationInputError(
                f"Malformed coordinate mapping: {error}"
            ) from error
        if xyz.shape != (len(self._ids), 3) or not np.isfinite(xyz).all():
            raise EvaluationInputError("Coordinates must be a finite N x 3 array")
        return xyz

    def _nonsingular(self, xyz):
        # Conservatively reject coincident sites (even an excluded pair) and
        # undefined angular/torsion geometries instead of accepting backend NaNs.
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            for i, j in combinations(range(len(xyz)), 2):
                if np.linalg.norm(xyz[i] - xyz[j]) < 1e-10:
                    raise EvaluationInputError(
                        f"Singular coincident sites {self._ids[i]}, {self._ids[j]}"
                    )
            triples = set(self._imported.angle_assignments)
            for assignments in (
                self._imported.proper_torsion_assignments,
                self._imported.improper_assignments,
            ):
                for a, b, c, d in assignments:
                    triples.update(((a, b, c), (b, c, d)))
            for a, b, c in triples:
                first, second = (
                    xyz[self._index[a]] - xyz[self._index[b]],
                    xyz[self._index[c]] - xyz[self._index[b]],
                )
                if np.linalg.norm(np.cross(first, second)) <= 1e-12 * np.linalg.norm(
                    first
                ) * np.linalg.norm(second):
                    raise EvaluationInputError(
                        f"Singular collinear angle/torsion at sites {(a, b, c)}"
                    )

    def evaluate_system(self, system: MolecularSystem) -> EvaluationResult:
        _validate(system, self._imported)
        if _graph(system) != self._graph_fingerprint:
            raise EvaluationInputError("Authoritative graph changed; rebind parameters")
        return self.evaluate({site: system.coordinates.get(site) for site in self._ids})

    def evaluate(
        self, coordinates=None, *, coordinate_unit="angstrom"
    ) -> EvaluationResult:
        xyz = self._coordinates(coordinates, coordinate_unit)
        mm, unit = _openmm()
        try:
            self._nonsingular(xyz)
            system = mm.XmlSerializer.deserialize(self._xml)
            integrator = mm.VerletIntegrator(
                0.001
            )  # Context requirement; never stepped.
            context = mm.Context(
                system,
                integrator,
                mm.Platform.getPlatformByName(self._platform),
                self._properties,
            )
            try:
                context.setPositions(xyz * 0.1 * unit.nanometer)
                state = context.getState(getEnergy=True, getForces=True)
                energy = state.getPotentialEnergy().value_in_unit(
                    unit.kilojoule_per_mole
                )
                forces = (
                    np.asarray(
                        state.getForces(asNumpy=True).value_in_unit(
                            unit.kilojoule_per_mole / unit.nanometer
                        ),
                        dtype=float,
                    )
                    * 0.1
                )
                components = {
                    name: context.getState(getEnergy=True, groups=1 << group)
                    .getPotentialEnergy()
                    .value_in_unit(unit.kilojoule_per_mole)
                    for group, name in enumerate(COMPONENTS)
                }
                platform = context.getPlatform()
                settings = {
                    **SETTINGS,
                    **{
                        "platform_" + key: platform.getPropertyValue(context, key)
                        for key in platform.getPropertyNames()
                    },
                    "precision": (
                        "double"
                        if self._platform == "Reference"
                        else self._properties.get("Precision", "platform default")
                    ),
                }
            finally:
                del context
                del integrator
            if (
                forces.shape != xyz.shape
                or not np.isfinite(forces).all()
                or not np.isfinite(energy)
                or not all(np.isfinite(v) for v in components.values())
            ):
                raise EvaluationError(
                    "OpenMM returned nonfinite energy/forces; check singular or extreme geometry"
                )
            coordinate_hash = fingerprint(
                {
                    "unit": "angstrom",
                    "sites": list(zip(self._ids, xyz.tolist(), strict=True)),
                }
            )
            evaluation_hash = fingerprint(
                {
                    "model": self._model_fingerprint,
                    "coordinates": coordinate_hash,
                    "backend": mm.__version__,
                    "platform": self._platform,
                    "settings": settings,
                }
            )
            return EvaluationResult(
                energy,
                components,
                dict(zip(self._ids, map(tuple, forces.tolist()), strict=True)),
                coordinate_hash,
                self._parameter_fingerprint,
                self._model_fingerprint,
                evaluation_hash,
                "OpenMM",
                mm.__version__,
                self._platform,
                settings,
            )
        except EvaluationError:
            raise
        except Exception as error:
            raise EvaluationError(
                f"OpenMM single-point evaluation failed: {error}"
            ) from error
