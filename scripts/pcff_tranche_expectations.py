"""Fail-closed software acceptance rules for the declared J13 negative controls."""

from collections import Counter
from pathlib import Path

from island.workflows import storage

CONTRACT_PATH = (
    Path(__file__).resolve().parents[1]
    / "docs/evidence/phase_4j13_negative_expectations_v1.json"
)
REQUIRED = {
    "main": {"protonated_amine", "peroxide", "sulfoxide"},
    "charged-sidechain": {"protonated_sidechain"},
}


def contract_for(declaration, raw, source_identity):
    contract = storage.read_json(CONTRACT_PATH)
    if contract["schema"] != "island_j13_negative_expectations_v1":
        raise ValueError("Unknown negative-control expectation contract")
    if (
        source_identity["sha256"] != contract["source_sha256"]
        or declaration["source"]["sha256"] != contract["source_sha256"]
    ):
        raise ValueError("Negative-control source mismatch")
    matches = [
        (name, c)
        for name, c in contract["control_sets"].items()
        if c["declaration_sha256"] == storage.checksum(raw)
    ]
    if len(matches) != 1:
        raise ValueError(
            "Declaration differs from the bound main/charged-sidechain control sets (missing, duplicate, empty or altered controls)"
        )
    name, selected = matches[0]
    required = REQUIRED[name]
    names = [c["name"] for c in declaration["control_graphs"]]
    if (
        Counter(names) != Counter(required)
        or set(selected["expectations"]) != required
        or Counter(selected["required_controls"]) != Counter(required)
    ):
        raise ValueError("Missing, duplicate, empty or unexpected required controls")
    stages = {
        "typing": (False, None, None),
        "charges": (True, False, False),
        "model": (True, True, False),
    }
    for expected in selected["expectations"].values():
        flags = tuple(
            expected[k]
            for k in ("typing_complete", "charges_complete", "model_complete")
        )
        wanted = stages.get(expected["stage"])
        if (
            expected["public_rejection"] is not True
            or wanted is None
            or any(a is not b for a, b in zip(flags, wanted, strict=True))
            or expected["exception"]
            != (
                "PreparedForceFieldError"
                if expected["stage"] == "model"
                else "PCFFError"
            )
        ):
            raise ValueError("Contradictory negative-control stage contract")
    return name, selected, storage.checksum(CONTRACT_PATH.read_bytes())


def assess(observed, expected):
    """Backend errors are evidence, never chemistry rejection by default."""
    failures = []
    for field in ("typing_complete", "charges_complete", "model_complete"):
        if observed.get(field) is not expected[field]:
            failures.append("Unexpected " + field)
    if observed.get("execution_error"):
        failures.append(
            "Unrelated execution failure at " + observed["execution_error"]["stage"]
        )
    if observed.get("public_preparation_complete"):
        failures.append("Unexpected successful public preparation")
    rejection = observed.get("public_exception", {})
    if rejection.get("type") != expected["exception"]:
        failures.append("Missing or wrong public exception type")
    if expected["stage"] == "model":
        diagnostics = observed.get("model_diagnostics", [])
        families = {r.get("id", "").split(":")[0] for r in diagnostics}
        if (
            not diagnostics
            or families != set(expected["missing_families"])
            or any(r.get("reason") != expected["model_reason"] for r in diagnostics)
        ):
            failures.append("Wrong missing-coupling model diagnostic")
        # Bind the public rejection to these exact native diagnostics, not a
        # permissive prefix match that could hide another preparation failure.
        message = expected["message_prefix"] + str(diagnostics)
    else:
        message = expected["message"]
    if rejection.get("message") != message:
        failures.append("Wrong public rejection reason/stage")
    if expected["stage"] == "typing":
        diagnostics = observed.get("typing_diagnostics", [])
        actual = Counter(
            (r.get("reason"), r.get("environment", {}).get("element"))
            for r in diagnostics
        )
        wanted = Counter(
            (r["reason"], r["element"]) for r in expected["typing_diagnostics"]
        )
        if actual != wanted:
            failures.append("Wrong typing chemistry diagnostic")
    else:
        if observed.get("typing_diagnostics") != []:
            failures.append("Unexpected typing diagnostics")
        independent = observed.get("independent_charge", {})
        if independent.get("complete") is not expected["charges_complete"]:
            failures.append("Independent charge outcome mismatch")
        if expected["stage"] == "charges":
            missing = independent.get("missing", [])
            if len(missing) != expected["missing_increment_count"] or any(
                r.get("reason") != expected["charge_reason"]
                or {s["path"]: s["types"] for s in r.get("searches", [])}
                != expected["missing_increment_searches"]
                or any(s.get("candidate_ids") != [] for s in r.get("searches", []))
                for r in missing
            ):
                failures.append("Wrong native missing-increment diagnostic")
            if observed.get("model_diagnostics") != expected["model_diagnostics"]:
                failures.append("Model was not blocked at native charges")
        elif independent.get("missing") != []:
            failures.append("Unexpected missing native charges")
    return failures
