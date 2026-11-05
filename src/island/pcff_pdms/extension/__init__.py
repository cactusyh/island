"""Offline PCFF-compatible PDMS extension contracts, separate from native PCFF.

This subpackage avoids changing the historical pcff_pdms/*.py identity glob.
Neither a target plan nor a synthetic numerical control registers a profile.
"""

from island.pcff_completion import canonical, digest, pin, require, signed
from island.pcff_qm.records import Record, record

PROFILE = "island_pcff_pdms_derived_v1"
FAMILY = "PCFF-compatible PDMS extension"
FLAGS = {
    "profile_origin": "qm_fitted_extension",
    "native_pcff_complete": False,
    "production_validated": False,
    "simulation_readiness": "not_established",
    "executable": False,
    "automatic_assignment_enabled": False,
    "accepted_coefficients": 0,
}
FIXTURES = ("hmds", "dp2", "dp4", "dp6")

__all__ = [
    "FAMILY",
    "FIXTURES",
    "FLAGS",
    "PROFILE",
    "Record",
    "canonical",
    "digest",
    "pin",
    "record",
    "require",
    "signed",
]
