"""Small backend contract; transfer and force-field compatibility are separate."""

from typing import Protocol

from island.charge_references.records import pack

from .core import FragmentChargeResult, boundary, charge_data


class FragmentChargeBackend(Protocol):
    def calculate(self, fragment) -> FragmentChargeResult: ...


class ProvidedFragmentChargeBackend:
    @boundary
    def __init__(
        self, charges, *, source="user-supplied; method provenance unverified"
    ):
        from copy import deepcopy

        self._charges = deepcopy(charges)
        self._source = source

    @boundary
    def calculate(self, fragment):
        from copy import deepcopy

        p = {
            "schema": "island_fragment_charges_v1",
            "fragment": fragment.payload,
            "raw_charges": deepcopy(self._charges),
            "backend": {
                "method": "provided",
                "source": self._source,
                "compatibility": "unverified",
            },
            "evidence": {},
            "data": {},
        }
        p["data"] = charge_data(p)
        result = FragmentChargeResult(pack(p))
        result.validate_integrity()
        return result
