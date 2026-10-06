"""Finite Reference-platform PCFF single points from validated H4 model records.

All forces are differentiated by OpenMM. No FRC parsing occurs during evaluation.
"""

from collections.abc import Mapping
from copy import deepcopy
from itertools import combinations
from math import sqrt

import numpy as np

from island.exceptions import EvaluationError, EvaluationInputError

from .models import EvaluationResult, fingerprint
from .openmm import OpenMMBoundPotential, _openmm, _OpenMMResources
from .pcff_identity import (
    SETTINGS,  # noqa: F401 -- historical import path
    pcff_evaluation_identity,
)

COMPONENTS = (
    "quartic_bond",
    "quartic_angle",
    "bond-bond",
    "bond-angle",
    "torsion_3",
    "middle_bond-torsion_3",
    "end_bond-torsion_3",
    "angle-torsion_3",
    "angle-angle-torsion_1",
    "bond-bond_1_3",
    "angle-angle",
    "lj_9_6",
    "coulomb",
)


def expression(family):
    """Lengths in expressions are angstroms; angular geometry is already radians."""
    r12 = "(10*distance(p1,p2))"
    r23 = "(10*distance(p2,p3))"
    r34 = "(10*distance(p3,p4))"
    a123 = "angle(p1,p2,p3)"
    a234 = "angle(p2,p3,p4)"
    phi = "dihedral(p1,p2,p3,p4)"
    if family in ("quadratic_bond", "quadratic_angle"):
        x = r12 if family == "quadratic_bond" else a123
        return f"c1*({x}-c0)^2"
    if family == "torsion_1":
        return f"c0*(1+cos(c1*{phi}-c2))"
    if family == "wilson_out_of_plane":
        # Signed scalar triple product of the three outward center-to-arm vectors.
        # Length factors cancel; x/y/z here are OpenMM nm coordinates.
        triple = (
            "((x1-x2)*((y3-y2)*(z4-z2)-(z3-z2)*(y4-y2))"
            "+(y1-y2)*((z3-z2)*(x4-x2)-(x3-x2)*(z4-z2))"
            "+(z1-z2)*((x3-x2)*(y4-y2)-(y3-y2)*(x4-x2)))"
        )
        norm = "(distance(p1,p2)*distance(p3,p2)*distance(p4,p2))"
        angles = ("angle(p3,p2,p4)", "angle(p4,p2,p1)", "angle(p1,p2,p3)")
        chi = (
            "(" + "+".join(f"asin({triple}/({norm}*sin({a})))" for a in angles) + ")/3"
        )
        return f"c0*({chi}-c1)^2"
    if family in ("quartic_bond", "quartic_angle"):
        x = r12 if family == "quartic_bond" else a123
        return "+".join(f"c{n - 1}*({x}-c0)^{n}" for n in (2, 3, 4))
    if family == "bond-bond":
        return f"c0*({r12}-e0)*({r23}-e1)"
    if family == "bond-angle":
        return f"(c0*({r12}-e0)+c1*({r23}-e1))*({a123}-e2)"
    if family == "bond-bond_1_3":
        return f"c0*({r12}-e0)*({r34}-e1)"
    if family == "angle-angle":
        return f"c0*({a123}-e0)*(angle(p3,p2,p4)-e1)"
    if family == "torsion_3":
        return "+".join(
            f"c{2 * n - 2}*(1-cos({n}*{phi}-c{2 * n - 1}))" for n in (1, 2, 3)
        )
    left = "(" + "+".join(f"c{n - 1}*cos({n}*{phi})" for n in (1, 2, 3)) + ")"
    right = "(" + "+".join(f"c{n + 2}*cos({n}*{phi})" for n in (1, 2, 3)) + ")"
    if family == "middle_bond-torsion_3":
        return f"({r23}-e0)*{left}"
    if family == "end_bond-torsion_3":
        return f"({r12}-e0)*{left}+({r34}-e1)*{right}"
    if family == "angle-torsion_3":
        return f"({a123}-e0)*{left}+({a234}-e1)*{right}"
    if family == "angle-angle-torsion_1":
        return f"c0*({a123}-e0)*({a234}-e1)*cos({phi})"
    raise EvaluationInputError("Unsupported PCFF term: " + str(family))


def build_model(data, masses, mm):
    """Compile only owned, validated ISLAND numerical data; no upstream model."""
    model = mm.System()
    ids = sorted(masses)
    index = {sid: i for i, sid in enumerate(ids)}
    for sid in ids:
        model.addParticle(masses[sid])
    families = list(enumerate(COMPONENTS[:-2]))
    if data.get("schema") in (
        "island_pcff_source_model_v1",
        "island_pcff_source_model_v2",
    ):
        families.append((13, "wilson_out_of_plane"))
    if data.get("schema") == "island_pcff_source_model_v2":
        families.extend(
            enumerate(("quadratic_bond", "quadratic_angle", "torsion_1"), 14)
        )
    for group, family in families:
        rows = [t for t in data["terms"] if t["family"] == family]
        if not rows:
            continue
        n = len(rows[0]["sites"])
        nc = len(rows[0]["coefficients"])
        ne = len(rows[0]["equilibria"])
        force = mm.CustomCompoundBondForce(n, expression(family))
        for k in range(nc):
            force.addPerBondParameter(f"c{k}")
        for k in range(ne):
            force.addPerBondParameter(f"e{k}")
        for row in rows:
            force.addBond(
                [index[s] for s in row["sites"]],
                row["coefficients"] + row["equilibria"],
            )
        force.setForceGroup(group)
        model.addForce(force)
    lj = mm.CustomBondForce("weight*epsilon*(2*(rmin/(10*r))^9-3*(rmin/(10*r))^6)")
    for name in ("rmin", "epsilon", "weight"):
        lj.addPerBondParameter(name)
    coul = mm.CustomBondForce("138.935456264*weighted_charge_product/r")
    coul.addPerBondParameter("weighted_charge_product")
    sites = {r["site"]: r for r in data["nonbonded"]}
    special = {tuple(p["sites"]): p for p in data["special_pair_inventory"]}
    for i, j in combinations(ids, 2):
        a, b = sites[i], sites[j]
        p = special.get((i, j), {"lj_weight": 1, "coulomb_weight": 1})
        mean = (a["rmin"] ** 6 + b["rmin"] ** 6) / 2
        eps = sqrt(a["epsilon"] * b["epsilon"]) * a["rmin"] ** 3 * b["rmin"] ** 3 / mean
        lj.addBond(index[i], index[j], [mean ** (1 / 6), eps, p["lj_weight"]])
        coul.addBond(
            index[i], index[j], [a["charge"] * b["charge"] * p["coulomb_weight"]]
        )
    for group, force in ((11, lj), (12, coul)):
        force.setForceGroup(group)
        model.addForce(force)
    return model


class PCFFSinglePointEvaluator(OpenMMBoundPotential):
    """Owned H4 potential; evaluate() is fresh, open_session() opts into reuse."""

    def __init__(self, system, specification):
        from island.charge_references.records import unpack
        from island.forcefields.pcff.automatic import chemical_graph

        try:
            (
                self._settings,
                self._parameter_fingerprint,
                self._model_fingerprint,
            ) = pcff_evaluation_identity(specification, system=system)
            data = unpack(specification.json_text)
            self._graph = deepcopy(chemical_graph(system))
            self._system = deepcopy(system)
            self._data = data
            self._ids = sorted(system.topology.sites)
            self._index = {s: i for i, s in enumerate(self._ids)}
            self._components = COMPONENTS
            if data.get("schema") in (
                "island_pcff_source_model_v1",
                "island_pcff_source_model_v2",
            ):
                self._components = COMPONENTS + ("wilson_out_of_plane",)
            if data.get("schema") == "island_pcff_source_model_v2":
                self._components += ("quadratic_bond", "quadratic_angle", "torsion_1")
            self._angles = [
                r["sites"]
                for r in data["terms"]
                if r["family"] in ("quartic_angle", "quadratic_angle")
            ]
        except EvaluationError:
            raise
        except Exception as error:
            raise EvaluationInputError(f"Invalid PCFF binding: {error}") from error
        mm, _ = _openmm()
        try:
            masses = {s: system.topology.sites[s].mass for s in self._ids}
            self._xml = mm.XmlSerializer.serialize(build_model(data, masses, mm))
        except Exception as error:
            raise EvaluationError(f"PCFF compilation failed: {error}") from error
        self._platform, self._properties = "Reference", {}

    @property
    def parameter_fingerprint(self):
        return self._parameter_fingerprint

    @property
    def model_fingerprint(self):
        return self._model_fingerprint

    def validate_system(self, system):
        from island.forcefields.pcff.automatic import chemical_graph

        try:
            if chemical_graph(system) != self._graph:
                raise EvaluationInputError(
                    "PCFF bound graph/stereo/mass changed; rebind"
                )
        except EvaluationError:
            raise
        except Exception as error:
            raise EvaluationInputError(f"Invalid PCFF bound system: {error}") from error

    def _coordinates(self, coordinates, unit):
        try:
            if unit != "angstrom":
                raise ValueError("Expected explicit angstrom coordinates")
            if coordinates is None:
                coordinates = {s: self._system.coordinates.get(s) for s in self._ids}
            if (
                not isinstance(coordinates, Mapping)
                or set(coordinates) != set(self._ids)
                or any(type(s) is not int for s in coordinates)
            ):
                raise ValueError("Exact integer stable-site coverage required")
            if any(
                isinstance(v, (bool, np.bool_))
                or not isinstance(v, (int, float, np.integer, np.floating))
                for s in self._ids
                for v in coordinates[s]
            ):
                raise ValueError("Nonboolean real coordinates required")
            xyz = np.array([coordinates[s] for s in self._ids], dtype=float)
            if xyz.shape != (len(self._ids), 3) or not np.isfinite(xyz).all():
                raise ValueError("Finite N x 3 coordinates required")
            with np.errstate(over="raise", invalid="raise", divide="raise"):
                for i, j in combinations(range(len(xyz)), 2):
                    if np.linalg.norm(xyz[i] - xyz[j]) < 1e-10:
                        raise ValueError("Coincident sites")
                for ids in self._angles:
                    a, b, c = (xyz[self._index[s]] for s in ids)
                    if np.linalg.norm(np.cross(a - b, c - b)) <= 1e-12 * np.linalg.norm(
                        a - b
                    ) * np.linalg.norm(c - b):
                        raise ValueError("Singular angular geometry")
            return xyz
        except Exception as error:
            raise EvaluationInputError(f"Invalid PCFF coordinates: {error}") from error

    def open_session(self):
        from .pcff_session import PCFFEvaluationSession

        return PCFFEvaluationSession(self)

    def _evaluate_context(self, xyz, context, mm, unit):
        context.setPositions(xyz * 0.1 * unit.nanometer)
        state = context.getState(getEnergy=True, getForces=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        forces = (
            state.getForces(asNumpy=True).value_in_unit(
                unit.kilojoule_per_mole / unit.nanometer
            )
            * 0.1
        )
        components = {
            name: context.getState(getEnergy=True, groups=1 << i)
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
            for i, name in enumerate(self._components)
        }
        if (
            not np.isfinite(energy)
            or not np.isfinite(forces).all()
            or not all(np.isfinite(v) for v in components.values())
        ):
            raise EvaluationError("Nonfinite PCFF backend output")
        coords = fingerprint(
            {
                "unit": "angstrom",
                "sites": list(zip(self._ids, xyz.tolist(), strict=True)),
            }
        )
        calc = fingerprint(
            {
                "model": self.model_fingerprint,
                "coordinates": coords,
                "backend": mm.__version__,
                "platform": "Reference",
            }
        )
        return EvaluationResult(
            energy,
            components,
            dict(zip(self._ids, map(tuple, forces.tolist()), strict=True)),
            coords,
            self.parameter_fingerprint,
            self.model_fingerprint,
            calc,
            "OpenMM PCFF Class II",
            mm.__version__,
            "Reference",
            self._settings,
        )

    def evaluate_fresh(self, coordinates=None, *, coordinate_unit="angstrom"):
        return self.evaluate(coordinates, coordinate_unit=coordinate_unit)

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        xyz = self._coordinates(coordinates, coordinate_unit)
        mm, unit = _openmm()
        resources = None
        try:
            resources = _OpenMMResources(self, mm)
            context = resources.context
            result = self._evaluate_context(xyz, context, mm, unit)
        except BaseException as error:
            if resources is not None:
                _cleanup(resources, error)
            if not isinstance(error, Exception) or isinstance(error, EvaluationError):
                raise
            raise EvaluationError(f"PCFF evaluation failed: {error}") from error
        else:
            _cleanup(resources)
            return result


def _cleanup(resources, original=None):
    try:
        resources.close()
    except BaseException as error:
        if original is not None:
            original.add_note(f"PCFF cleanup also failed: {error}")
        elif isinstance(error, Exception):
            raise EvaluationError(f"PCFF cleanup failed: {error}") from error
        else:
            raise
