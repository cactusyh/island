"""Optional pinned Foyer OPLS-AA typing, native charges and parameter records."""

from .adapter import type_atoms
from .models import (
    OPLSChargeResult,
    OPLSTypingResult,
    assign_native_charges,
    load_opls_result,
    save_opls_result,
)
from .parameters import OPLSParameterizationResult, parameterize_oplsaa
from .source import FoyerOPLSSource, load_oplsaa_source

__all__ = [
    "FoyerOPLSSource",
    "OPLSChargeResult",
    "OPLSParameterizationResult",
    "OPLSTypingResult",
    "assign_native_charges",
    "load_opls_result",
    "load_oplsaa_source",
    "parameterize_oplsaa",
    "save_opls_result",
    "type_atoms",
]
