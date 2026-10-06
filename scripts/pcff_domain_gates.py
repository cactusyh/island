"""Fail-closed J3 bounded gates, separate from full-source coverage."""

import math
from collections import Counter

REQUIRED = ("hydrogen_sulfide", "heavy_water")


def numerical_gate(rows, source_unchanged):
    if source_unchanged is not True or Counter(r.get("case") for r in rows) != Counter(
        REQUIRED
    ):
        return False
    for row in rows:
        if (
            row.get("passed") is not True
            or row.get("error")
            or not row.get("bundle_hashes")
            or not row.get("prepared_identity")
            or not row.get("model_identity")
        ):
            return False
        frames = row.get("numerical", [])
        if len(frames) != 2:
            return False
        for frame in frames:
            if frame.get("parameter_rows_and_inventory_match") is not True:
                return False
            for key in ("energy_error", "force_max_error"):
                value = frame.get(key)
                if (
                    type(value) not in (int, float)
                    or not math.isfinite(value)
                    or not 0 <= value <= 1e-5
                ):
                    return False
    return True
