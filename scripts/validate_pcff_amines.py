"""J14 fail-closed typing/charge/model gate. Numerical/workflow gates are separate."""

import argparse
from collections import Counter
from pathlib import Path

from inspect_pcff_tranche_controls import execute_control
from pcff_j5_reference import charge_raw, read_source
from pcff_j14_reference import labels_for

from island.forcefields.pcff import load_pcff_source
from island.workflows import storage
from island.workflows.bundle import system_from

DECLARATION = Path(__file__).parents[1] / "docs/evidence/phase_4j14_declaration.json"
DECLARATION_SHA256 = "e3a8e3c40efc76203990cea0d4e8232abbb7493cb753296ec7a578583baa8c3d"
REQUIRED = {
    "primary",
    "secondary",
    "tertiary",
    "primary_protonated",
    "secondary_protonated",
    "tertiary_protonated",
    "quaternary",
    "amine_oligomer",
    "protonated_sidechain",
    "amide",
    "carbamate",
    "aromatic_amine",
    "pyridinium",
    "guanidinium",
    "radical",
}
NEGATIVE = {
    "amide": "model",
    "carbamate": "charges",
    "aromatic_amine": "charges",
    "pyridinium": "charges",
    "guanidinium": "charges",
    "radical": "typing",
}
MISSING = {
    "carbamate": [("hn", "n_2")] * 2,
    "aromatic_amine": [("hn2", "nn")] * 2,
    "pyridinium": [("cp", "nh+")] * 2 + [("hn", "nh+")],
    "guanidinium": [],
}


def assess(r):
    """Require the declared stage and actual chemistry reason, not any exception."""
    failures = []
    name = r.get("case")
    if name not in REQUIRED:
        return ["Unexpected case"]
    if r.get("execution_error"):
        failures.append("Execution error is not chemical rejection")
    stage = NEGATIVE.get(name)
    flags = (
        (False, None, None)
        if stage == "typing"
        else (True, False, False)
        if stage == "charges"
        else (True, True, False)
        if stage == "model"
        else (True, True, True)
    )
    if (
        tuple(
            r.get(k) for k in ("typing_complete", "charges_complete", "model_complete")
        )
        != flags
    ):
        failures.append("Wrong typing/charge/model stage")
    q = r.get("independent_charge", {})
    if stage is None:
        if r.get("public_preparation_complete") is not True or r.get(
            "public_exception"
        ):
            failures.append("Expected successful public preparation")
        if r.get("independent_types_and_charges") is not True:
            failures.append("Independent labels/charges were not verified")
    else:
        if r.get("public_preparation_complete"):
            failures.append("Unexpected public preparation success")
        exception = r.get("public_exception", {})
        diagnostics = r.get("model_diagnostics", [])
        message = (
            "Incomplete source graph typing"
            if stage == "typing"
            else "Complete typing and charges required"
            if stage == "charges"
            else "Incomplete PCFF model: " + str(diagnostics)
        )
        if exception != {
            "type": "PreparedForceFieldError" if stage == "model" else "PCFFError",
            "message": message,
        }:
            failures.append("Wrong public rejection reason")
        if stage == "typing":
            actual = Counter(
                (d["reason"], d["environment"]["element"])
                for d in r.get("typing_diagnostics", [])
            )
            if actual != Counter(
                {
                    ("isotope_radical_or_nonstandard_mass_rule_unresolved", "N"): 1,
                    ("chemical_rule_ambiguous", "H"): 1,
                }
            ):
                failures.append("Wrong radical diagnostic")
        elif stage == "charges":
            missing = q.get("missing", [])
            if Counter(
                tuple(sorted(m["searches"][0]["types"])) for m in missing
            ) != Counter(MISSING[name]) or any(
                m["reason"] != "source_parameter_missing" for m in missing
            ):
                failures.append("Wrong missing increment diagnostic")
            if q.get("complete") is not False or diagnostics != [
                {
                    "reason": "native_charge_incomplete; parameterization and model construction not attempted"
                }
            ]:
                failures.append("Wrong charge/model failure")
            if name == "guanidinium" and [
                (c["formal"], c["known_contribution_total"], c["residual"])
                for c in q.get("components", [])
            ] != [(1, "0.9999", "-0.0001")]:
                failures.append("Wrong guanidinium residual")
        else:
            reasons = {(x["id"].split(":")[0], x["reason"]) for x in diagnostics}
            missing_reason = (
                "no source row after declared exact/equivalent/wildcard searches"
            )
            wanted = {
                (f, missing_reason)
                for f in ("bond-bond", "bond-angle", "middle_bond-torsion_3")
            }
            wanted.add(
                (
                    "torsion_3",
                    "conflicting oriented candidates in msi2lmp match tier; file order is not physical authority",
                )
            )
            if reasons != wanted or q.get("complete") is not True:
                failures.append("Wrong amide coupling/ambiguity diagnostic")
    return failures


def assess_set(cases):
    if Counter(c.get("case") for c in cases) != Counter(REQUIRED):
        return ["Missing, duplicate, empty or unexpected required cases"]
    return [f"{r['case']}: {f}" for r in cases for f in assess(r)]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island_j14_amine_acceptance_v1",
        "acceptance_passed": False,
        "cases": [],
        "failures": [],
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        if storage.checksum(DECLARATION.read_bytes()) != DECLARATION_SHA256:
            raise ValueError("Amine declaration changed from the reviewed experiment")
        d = storage.read_json(DECLARATION)
        source = load_pcff_source(a.source)
        if (
            source.identity["sha256"] != d["source"]["sha256"]
            or d["typing_profile"] != "island_pcff_source_graph_v5"
            or Counter(c["name"] for c in d["cases"]) != Counter(REQUIRED)
        ):
            raise ValueError("Source/profile/required-case declaration mismatch")
        report.update(
            source=source.identity,
            declaration_sha256=storage.checksum(DECLARATION.read_bytes()),
        )
        for c in d["cases"]:
            r = execute_control(c, d, source, a.source, a.output)
            if c["name"] not in NEGATIVE and not r.get("execution_error"):
                s = system_from(
                    storage.decode(
                        storage.read_json(a.output / (c["name"] + "-system.json"))
                    )
                )
                expected = labels_for(s)
                independent = charge_raw(read_source(source.raw), s, expected)
                inspection = storage.decode(
                    storage.read_json(a.output / (c["name"] + "-inspection.json"))
                )
                actual = inspection["native_charge_record"]["automatic_typing"][
                    "assignments"
                ]
                r["independent_types_and_charges"] = (
                    actual == expected
                    and independent == r.get("independent_charge")
                    and independent["complete"]
                )
            r["acceptance_failures"] = assess(r)
            report["cases"].append(r)
        report["failures"] = assess_set(report["cases"])
        report["acceptance_passed"] = not report["failures"]
    except Exception as e:  # noqa: BLE001 -- failures retained, never accepted
        report["failures"].append(type(e).__name__ + ": " + str(e))
    storage.publish(
        a.output / "report.json", storage.json_bytes(storage.encode(report))
    )
    print(
        {
            "acceptance_passed": report["acceptance_passed"],
            "failures": report["failures"],
        }
    )
    return 0 if report["acceptance_passed"] and not a.require_full_source else 1


if __name__ == "__main__":
    raise SystemExit(main())
