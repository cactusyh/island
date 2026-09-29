"""Keep current coordinate labels aligned and superseded diagnostics historical."""

import hashlib
import json
from copy import deepcopy


def coordinate_hash(coordinates):
    """Canonical angstrom frame fingerprint shared by coordinate producers."""
    payload = {
        "unit": "angstrom",
        "sites": [(site, list(coordinates[site])) for site in sorted(coordinates)],
    }
    return hashlib.sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def previous_coordinate_source(metadata):
    polymer = metadata.get("polymer")
    if isinstance(polymer, dict):
        generation = polymer.get("coordinate_generation", {})
        return (
            metadata.get("coordinate_source")
            or (
                generation.get("coordinate_source")
                if isinstance(generation, dict)
                else None
            )
            or polymer.get("coordinates")
        )
    generation = metadata.get("coordinate_generation", {})
    return metadata.get("coordinate_source") or (
        generation.get("coordinate_source") if isinstance(generation, dict) else None
    )


def updated_coordinate_metadata(
    system, provenance, *, previous_fingerprint, minimization=None, dynamics=None
):
    """Return owned metadata; only coordinate provenance fields are replaced.

    History is a flat chronological list of superseded coordinate records. Raw
    prior records are retained, including legacy labels, rather than rewritten.
    Preparation and parameter provenance are deliberately outside this helper.
    """
    metadata = deepcopy(system.metadata)
    polymer = metadata.get("polymer")
    prior = {
        "coordinate_source": previous_coordinate_source(metadata),
        "coordinate_fingerprint": previous_fingerprint,
        "records": {
            key: deepcopy(metadata[key])
            for key in (
                "coordinate_source",
                "coordinate_generation",
                "local_minimization",
                "dynamics",
            )
            if key in metadata
        },
    }
    if isinstance(polymer, dict):
        prior["polymer_records"] = {
            key: deepcopy(polymer[key])
            for key in ("coordinates", "coordinate_generation")
            if key in polymer
        }
    metadata.setdefault("coordinate_history", []).append(prior)
    metadata["coordinate_source"] = provenance["coordinate_source"]
    metadata["coordinate_generation"] = deepcopy(provenance)
    if isinstance(polymer, dict):
        polymer["coordinates"] = provenance["coordinate_source"]
        polymer["coordinate_generation"] = deepcopy(provenance)
    metadata.pop("local_minimization", None)
    metadata.pop("dynamics", None)
    if minimization is not None:
        metadata["local_minimization"] = deepcopy(minimization)
    if dynamics is not None:
        metadata["dynamics"] = deepcopy(dynamics)
    return metadata
