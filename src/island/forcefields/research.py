"""Opt-in research models, preserving every historical preparation facade."""

from dataclasses import dataclass

from . import extensible as existing
from .pdms_opls import FAMILY, PDMSOPLSError, PDMSOPLSOptions, PDMSOPLSPreparation
from .pdms_opls import prepare_pdms_opls as _prepare


@dataclass(frozen=True)
class PDMSForceFieldRequest:
    family: str
    options: PDMSOPLSOptions

    def __post_init__(self):
        if self.family != FAMILY or type(self.options) is not PDMSOPLSOptions:
            raise PDMSOPLSError("typing_blocked", "explicit_pdms_options_required")


def ForceFieldRequest(family, options, **kwargs):
    if family == FAMILY:
        if kwargs:
            raise PDMSOPLSError("typing_blocked", "unexpected_pdms_request_fields")
        return PDMSForceFieldRequest(family, options)
    return existing.ForceFieldRequest(family, options, **kwargs)


def prepare_forcefield(system, request):
    if type(request) is PDMSForceFieldRequest:
        PDMSForceFieldRequest(request.family, request.options)
        return _prepare(system, request.options)
    return existing.prepare_forcefield(system, request)


def create_evaluator(system, prepared):
    if type(prepared) is PDMSOPLSPreparation:
        prepared.require_executable(system)
        from island.pdms import PDMSOPLSEvaluator

        return PDMSOPLSEvaluator(system, prepared)
    return existing.create_evaluator(system, prepared)


__all__ = [
    "ForceFieldRequest",
    "PDMSOPLSError",
    "PDMSOPLSOptions",
    "PDMSOPLSPreparation",
    "create_evaluator",
    "prepare_forcefield",
]
