"""Additive public preparation facade; historical native dispatch is untouched.

Use this explicit namespace for extensions. Native requests still delegate to
their original facade. No plugin discovery, monkeypatching or source fallback.
"""

from dataclasses import dataclass

from . import preparation as native
from .opls_external import (
    FAMILY,
    ExternalOPLSError,
    ExternalOPLSOptions,
    ExternalOPLSPreparation,
    prepare_external_opls,
)


@dataclass(frozen=True)
class ExternalForceFieldRequest:
    family: str
    options: ExternalOPLSOptions

    def __post_init__(self):
        if self.family != FAMILY or type(self.options) is not ExternalOPLSOptions:
            raise native.ForceFieldRequestError(
                "Explicit external OPLS options required"
            )


def ForceFieldRequest(family, options, **kwargs):
    """Create an explicit native request or the independent external request."""
    if family == FAMILY:
        if kwargs:
            raise native.ForceFieldRequestError(
                "External OPLS binds its validated extension graph; native graph records cannot be substituted"
            )
        return ExternalForceFieldRequest(family, options)
    return native.ForceFieldRequest(family, options, **kwargs)


def prepare_forcefield(system, request):
    if type(request) is ExternalForceFieldRequest:
        ExternalForceFieldRequest(request.family, request.options)
        return prepare_external_opls(system, request.options)
    return native.prepare_forcefield(system, request)


def create_evaluator(system, prepared):
    if type(prepared) is ExternalOPLSPreparation:
        prepared.require_executable(system)
        from .opls_external.evaluator import ExternalOPLSEvaluator

        return ExternalOPLSEvaluator(system, prepared)
    return native.create_evaluator(system, prepared)


__all__ = [
    "ExternalForceFieldRequest",
    "ExternalOPLSError",
    "ExternalOPLSOptions",
    "ExternalOPLSPreparation",
    "ForceFieldRequest",
    "create_evaluator",
    "prepare_forcefield",
]
