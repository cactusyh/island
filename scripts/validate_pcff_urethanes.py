"""J15 declared typing/charge controls; incomplete polymers fail acceptance.

No numerical experiments run here. Independent numerical/workflow evidence is
required separately; a diagnostic matrix cannot authorize executable coverage.
"""

import argparse
from collections import Counter
from pathlib import Path

from inspect_pcff_tranche_controls import execute_control
from pcff_j5_reference import charge_raw, read_source
from pcff_j15_reference import labels_for

from island.forcefields.pcff import load_pcff_source
from island.workflows import storage
from island.workflows.bundle import system_from

DECLARATION = Path(__file__).parents[1] / "docs/evidence/phase_4j15_declaration.json"
CONTRACT = DECLARATION.with_name("phase_4j15_expectations.json")
CONTRACT_SHA256 = "0c9d761bdc338a07c6f66caa4757a2be5880c43455d05e6c83d82b21037c074b"


def contract():
    if storage.checksum(CONTRACT.read_bytes()) != CONTRACT_SHA256:
        raise ValueError("J15 expectation contract changed")
    return storage.decode(storage.read_json(CONTRACT))


def assess(result, expected):
    failures = []
    if result.get("execution_error"):
        failures.append("Execution error is not a chemical rejection")
    if expected["role"] == "negative_unchanged_v5":
        if result.get("public_preparation_complete"):
            failures.append("Unexpected public preparation success")
        for k, v in expected["expected"].items():
            if result.get(k) != v:
                failures.append("Historical negative stage/reason changed: " + k)
        return failures
    if (
        result.get("typing_complete") is not True
        or result.get("charges_complete") is not True
    ):
        failures.append("Expected complete typing and native charges")
    if result.get("independent_types_and_charges") is not True:
        failures.append("Independent labels, rows and charge vector not verified")
    counts = Counter(result.get("labels", {}))
    if any(counts[k] != n for k, n in expected["expected_motif_counts"].items()):
        failures.append("Wrong declared motif labels")
    if result.get("model_complete") is not expected["model_complete"]:
        failures.append("Wrong model stage")
    if result.get("model_diagnostics") != expected["model_diagnostics"]:
        failures.append("Wrong coupling/dependency diagnostic")
    if result.get("public_exception") != expected["public_exception"]:
        failures.append("Wrong public rejection reason")
    if bool(result.get("public_preparation_complete")) != expected["model_complete"]:
        failures.append("Unexpected public preparation outcome")
    return failures


def assess_set(cases, expected):
    if Counter(r.get("case") for r in cases) != Counter(expected.keys()):
        return ["Missing, duplicate, empty or unexpected required cases"]
    return [f"{r['case']}: {f}" for r in cases for f in assess(r, expected[r["case"]])]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    r = dict(  # noqa: C408
        schema="island_j15_acceptance_v1",
        cases=[],
        failures=[],
        expectations_passed=False,
        executable_polymer_gate=False,
        acceptance_passed=False,
        full_source_complete=False,
        production_validated=False,
        simulation_readiness="not_established",
    )
    try:
        e, d = contract(), storage.read_json(DECLARATION)
        source = load_pcff_source(a.source)
        if (
            storage.checksum(DECLARATION.read_bytes()) != e["declaration_sha256"]
            or source.identity["sha256"] != e["source_sha256"]
            or d["typing_profile"] != e["typing_profile"]
            or Counter(c["name"] for c in d["cases"]) != Counter(e["cases"].keys())
        ):
            raise ValueError("Declaration/source/profile/case mismatch")
        r.update(
            source=source.identity,
            declaration_sha256=e["declaration_sha256"],
            contract_sha256=CONTRACT_SHA256,
        )
        raw = read_source(source.raw)
        for c in d["cases"]:
            result = execute_control(c, d, source, a.source, a.output)
            if c["role"] == "positive" and not result.get("execution_error"):
                s = system_from(
                    storage.decode(
                        storage.read_json(a.output / (c["name"] + "-system.json"))
                    )
                )
                labels = labels_for(s)
                q = charge_raw(raw, s, labels)
                native = storage.decode(
                    storage.read_json(a.output / (c["name"] + "-inspection.json"))
                )["native_charge_record"]
                vector = native["native_charge_record"]["partial_charges"]
                result["independent_types_and_charges"] = (
                    native["automatic_typing"]["assignments"] == labels
                    and q == result.get("independent_charge")
                    and q["complete"]
                    and set(vector) == set(labels)
                    and all(
                        abs(vector[i] - float(q["partial_charges"][i]))
                        <= d["charge_atol_e"]
                        for i in vector
                    )
                )
            r["cases"].append(result)
        r["failures"] = assess_set(r["cases"], e["cases"])
        r["expectations_passed"] = not r["failures"]
        # This frozen declaration has no source-complete polymer. Never count
        # an expected incomplete outcome as the requested executable tranche.
        r["failures"].append(
            "Executable polymer gate unmet: required source couplings unresolved"
        )
    except Exception as exc:  # noqa: BLE001 -- retain evidence, fail closed
        r["failures"].append(type(exc).__name__ + ": " + str(exc))
    storage.publish(a.output / "report.json", storage.json_bytes(storage.encode(r)))
    print({k: v for k, v in r.items() if k not in ("cases", "source")})
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
