"""Offline PCFF evaluation identity; no OpenMM import or numerical construction.

Accept a native validated specification, not an unchecked settings override.
The tuple contains owned settings, parameter fingerprint and model fingerprint.
Historical H5 and expanded J1 fingerprint formulas are preserved exactly.
"""

from copy import deepcopy

from island.exceptions import EvaluationError, EvaluationInputError

from .models import fingerprint

SETTINGS = {
    "implementation": "island_pcff_singlepoint_v1",
    "compatibility_profile": "island_lammps_pcff_acyclic_cho_v1",
    "method": "NoCutoff",
    "periodic": False,
    "constraints": False,
    "switching": False,
    "tail_correction": False,
    "precision": "double",
    "integration_steps": 0,
    "force_sign": "-dE/dR",
    "mixing": "sixth_power",
    "coulomb_constant_nm": 138.935456264,
}


def pcff_evaluation_identity(specification, *, system=None):
    """Validate a PCFF model and derive its existing evaluation identities offline.

    Optional system validation uses the native graph/mass/stereo contract and
    permits coordinate replacement. No force-field typing or force evaluation
    occurs. Each call returns a new settings dictionary.
    """
    from island.charge_references.records import unpack
    from island.forcefields.pcff.charges import identity
    from island.forcefields.pcff.model import PCFFModelSpecification

    try:
        if type(specification) is not PCFFModelSpecification:
            raise EvaluationInputError("Validated PCFFModelSpecification required")
        specification.validate_integrity(system)
        data = unpack(specification.json_text)
        if not data["model_definition_complete"]:
            raise EvaluationInputError("Incomplete PCFF model definition")
        settings = deepcopy(SETTINGS)
        if data["schema"] == "island_pcff_source_model_v1":
            settings.update(
                implementation="island_pcff_source_singlepoint_v1",
                compatibility_profile=data["compatibility_profile"]["name"],
            )
        if data["schema"] == "island_pcff_source_model_v2":
            settings.update(
                implementation="island_pcff_fallback_singlepoint_v1",
                compatibility_profile=data["compatibility_profile"]["name"],
            )
        if data["schema"] == "island_pcff_typed_graph_model_v1":
            settings.update(
                implementation="island_pcff_typed_graph_singlepoint_v1",
                compatibility_profile=data["compatibility_profile"]["name"],
            )
        parameter = identity(data)
        model = fingerprint({"specification": parameter, "settings": settings})
        return settings, parameter, model
    except EvaluationError:
        raise
    except Exception as error:
        raise EvaluationInputError(f"Invalid PCFF binding: {error}") from error
