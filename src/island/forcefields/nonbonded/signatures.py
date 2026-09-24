"""Deterministic nonbonded-policy content signatures."""

import hashlib
import json
from dataclasses import asdict

from island.forcefields.nonbonded.models import NonbondedPolicy


def nonbonded_policy_signature(policy: NonbondedPolicy) -> str:
    payload = asdict(policy)
    payload["mixing_formula"] = policy.mixing_formula
    payload["shortest_path_convention"] = "1-2/1-3/1-4; farther_or_disconnected_full"
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
