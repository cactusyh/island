"""Finite external OPLS evaluation; no native-source identity or periodic bypass."""

from copy import deepcopy

from island.evaluation.models import EvaluationResult, fingerprint
from island.evaluation.openmm import OpenMMBoundPotential, _openmm
from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.forcefields.external_pdms import require

from .assignment import ExternalOPLSPreparation


def build_model(data, mm):
    # The existing native builder implements exactly the geometric finite LJ,
    # Coulomb, harmonic bond/angle and RB equations. Source selection is separate.
    from island.evaluation.oplsaa import build_model as build_common_forms

    common = deepcopy(data)
    if data["torsion_form"] == "OPLSTorsionForce":
        common["rb"] = []
    model = build_common_forms(common, mm)
    if data["torsion_form"] == "OPLSTorsionForce":
        index = {i: n for n, i in enumerate(sorted(data["sites"]))}
        force = mm.CustomTorsionForce(
            "0.5*(k1*(1+cos(theta))+k2*(1-cos(2*theta))+k3*(1+cos(3*theta))+k4*(1-cos(4*theta)))"
        )
        for n in range(1, 5):
            force.addPerTorsionParameter("k" + str(n))
        for row in data["rb"]:
            force.addTorsion(
                *(index[i] for i in row["sites"]),
                row["converted"]["coefficients_kj_mol"],
            )
        force.setForceGroup(2)
        model.addForce(force)
    return model


class ExternalOPLSEvaluator(OpenMMBoundPotential):
    # Geometry/resource lifetime code is shared, not the native source contract.
    _coordinates = OPLSSinglePointEvaluator._coordinates

    def __init__(self, system, preparation):
        require(
            type(preparation) is ExternalOPLSPreparation,
            "external_preparation_required",
        )
        preparation.require_executable(system)
        self._preparation = preparation
        self._system = deepcopy(system)
        self._data = preparation.payload["assignment"]
        self._ids = sorted(system.topology.sites)
        self._index = {i: n for n, i in enumerate(self._ids)}
        self._parameter_fingerprint = preparation.identity
        self._model_fingerprint = fingerprint(
            {
                "external_opls_preparation": preparation.identity,
                "backend": "finite_external_opls_reference_v1",
            }
        )
        mm, _ = _openmm()
        self._xml = mm.XmlSerializer.serialize(build_model(self._data, mm))
        self._platform, self._properties = "Reference", {}

    @property
    def parameter_fingerprint(self):
        return self._parameter_fingerprint

    @property
    def model_fingerprint(self):
        return self._model_fingerprint

    def validate_system(self, system):
        self._preparation.require_executable(system)

    def evaluate(self, coordinates=None, *, coordinate_unit="angstrom"):
        self._preparation.require_executable(self._system)
        return OPLSSinglePointEvaluator.evaluate(
            self, coordinates, coordinate_unit=coordinate_unit
        )

    def evaluate_fresh(self, coordinates=None, *, coordinate_unit="angstrom"):
        return self.evaluate(coordinates, coordinate_unit=coordinate_unit)

    def _evaluate_context(self, xyz, context, mm, unit):
        # Reject changed assets before every numerical evaluation, including after
        # evaluator construction. Coordinates do not change chemical identity.
        self._preparation.require_executable(self._system)
        r = OPLSSinglePointEvaluator._evaluate_context(self, xyz, context, mm, unit)
        components = dict(r.energy_components)
        components["proper_torsion"] = components.pop("rb_proper")
        return EvaluationResult(
            r.potential_energy,
            components,
            r.forces,
            r.coordinate_fingerprint,
            r.parameter_fingerprint,
            r.model_fingerprint,
            r.evaluation_fingerprint,
            "OpenMM external OPLS-compatible",
            r.backend_version,
            r.platform,
            {
                **dict(r.settings),
                "implementation": "island_external_opls_evaluation_v1",
                "native_opls_coverage": False,
            },
        )
