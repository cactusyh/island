"""Source-bound immutable QM plans; unfilled settings never become QM results."""

import json
from dataclasses import dataclass

from . import _verified, canonical, digest, equation, require, signed
from .confirmation import reported_confirmation

FIXTURE_GRAPHS = {
    "disiloxane": "7aa04b55c75e48d32abd159a51e6f6e39fe51244521118c94df49484ba1a738f",
    "hmds": "d7052150c506aace667eee55c0106c3aa97c8ca5be9a29bc5ecf028e28bde5dd",
    "dp2": "507d05b898ec84c598764267cc8ddb494d77741800d909030c212b801113f611",
    "dp4": "45177e7114115a0214b0e3e855b7290de5fdfa2f3cf435271efbd27ecc2612d4",
    "dp6": "a80f0424786eaa596b543c58fdbe67c32272f1bc4cbb594507df1d3a8d84c6a3",
}

SETTINGS = (
    "qm_method",
    "basis_set",
    "dispersion_treatment",
    "convergence_criteria",
    "deformation_range",
    "expected_uncertainty",
    "held_out_thresholds",
)


def targets(system, assessment, fixture):
    p = _verified(system, assessment)
    require(
        FIXTURE_GRAPHS.get(fixture) == p["graph_identity"], "fixture_graph_mismatch"
    )
    grouped = {}
    for r in p["missing_terms"]:
        if "torsion" not in r["family"]:
            grouped.setdefault(r["query"], []).append(r)
    result = []
    for query, sites in sorted(grouped.items()):
        row = min(sites, key=lambda r: r["id"])
        result.append(
            signed(
                {
                    "schema": "island.pcff-qm-target.v1",
                    "family": "pcff",
                    "status": "qm_fit_required",
                    "query": query,
                    "interaction_type": row["family"],
                    "atom_type_environment": query.split(":", 1)[1].split(","),
                    "molecular_fixture": fixture,
                    "target_site": row,
                    "applicable_sites": sorted(r["id"] for r in sites),
                    **{
                        k: p[k]
                        for k in (
                            "graph_identity",
                            "typing_identity",
                            "charge_identity",
                            "atom_order",
                            "charge_vector",
                            "source_hashes",
                            "domain",
                            "candidate_profile",
                        )
                    },
                    "charge_model": p["charge_origin"],
                    "author_confirmation_identity": reported_confirmation()["identity"],
                    "functional_form": equation(assessment, row["family"]),
                    "required_deformation": (
                        "independent signed paired angle displacements about verified equilibria"
                        if row["family"] == "angle-angle"
                        else "independent signed first/third bond displacements about verified equilibria"
                    ),
                    "observables": [
                        "energies",
                        "forces",
                        "Hessians_or_frequencies_when_available",
                    ],
                    "units": {
                        "coordinates": "angstrom",
                        "energy": "kcal/mol",
                        "forces": "kcal/mol/angstrom",
                        "hessian": "kcal/mol/angstrom^2",
                        "frequencies": "cm^-1",
                        "charge": "e",
                        "angles": "radian",
                    },
                    **{k: None for k in SETTINGS},
                    "planning_status": "user_QM_settings_required",
                    "identifiability_requirements": [
                        "joint Class-II residual or independently verified target isolation",
                        "full rank weighted Jacobian for all fitted coefficients",
                        "more independent training observations than coefficients",
                        "opposite-sign and independent deformations; no duplicate geometries",
                        "fixed verified equilibria, typing and charge evidence",
                    ],
                    "held_out_validation_requirements": [
                        "disjoint geometries; declared error thresholds before fitting",
                        "energies and available forces/Hessians with uncertainties",
                        "validate only the exact declared graph domain",
                    ],
                    "uncertainty": {
                        "status": "reference_data_not_supplied",
                        "value": None,
                        "components": ["QM", "model", "statistical", "identifiability"],
                    },
                    "execution_authority": False,
                }
            )
        )
    return tuple(QMTarget(canonical(r)) for r in result)


@dataclass(frozen=True)
class QMTarget:
    json_text: str

    @property
    def payload(self):
        p = json.loads(self.json_text)
        require(
            p.get("identity")
            == digest({k: v for k, v in p.items() if k != "identity"}),
            "qm_target_checksum_mutation",
        )
        return p

    def validate(self, system, assessment):
        p = self.payload
        expected = targets(system, assessment, p["molecular_fixture"])
        require(any(t.payload == p for t in expected), "qm_target_binding_mutation")
        return p


def configure(target, settings):
    """A new immutable plan revision; never mutate or supply hidden QM defaults."""
    p = target.payload
    require(set(settings) == set(SETTINGS), "explicit_qm_settings_required")
    require(
        all(v is not None and v != "" and v != {} for v in settings.values()),
        "explicit_qm_settings_required",
    )
    for key in (
        "deformation_range",
        "expected_uncertainty",
        "held_out_thresholds",
        "convergence_criteria",
    ):
        require(
            isinstance(settings[key], dict), "structured_qm_setting_required", field=key
        )
    return signed(
        {
            "schema": "island.pcff-qm-plan.v1",
            "target_identity": p["identity"],
            "settings": settings,
            "reference_data_status": "not_supplied",
            "execution_authority": False,
        }
    )
