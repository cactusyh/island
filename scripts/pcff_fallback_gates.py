"""Fail-closed J2 evidence checks; offline consistency, not authentication.

Legacy receipts record absolute errors after the original allclose assertions.
Their recorded errors must satisfy the declared absolute tolerance on their own;
missing per-component relative-error evidence is not invented. No scientific
rerun is made.
"""

import math
from pathlib import Path

from island.workflows.storage import checksum

REAL = "real_source_coefficients_on_declared_asymmetric_term_geometry"
SYNTHETIC = "synthetic nonzero-phase software fixture"
TERM_KEYS = {
    (family, REAL)
    for family in (
        "quadratic_bond",
        "quadratic_angle",
        "torsion_1",
        "wilson_out_of_plane",
    )
} | {("torsion_1", SYNTHETIC)}
BUNDLE_FILES = {"assignment.json", "model.json", "system.json", "manifest.json"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def case_map(rows):
    require(type(rows) is list, "Case list required")
    result = {}
    for row in rows:
        require(type(row) is dict and type(row.get("case")) is str, "Malformed case")
        name = row["case"]
        require(name not in result, "Duplicate case: " + name)
        result[name] = row
    return result


def mandatory_case(report):
    cases = case_map(report.get("cases"))
    require("dichlorine" in cases, "Missing mandatory dichlorine case")
    return cases["dichlorine"]


def successful_comparison(row):
    require(
        type(row) is dict and "error" not in row and "failure" not in row,
        "Failed numerical comparison",
    )
    require(row.get("passed", True) is True, "Failed numerical comparison")
    for field in ("energy_error", "force_max_error"):
        value = row.get(field)
        require(
            type(value) in (int, float) and math.isfinite(value) and 0 <= value <= 1e-5,
            "Missing/invalid/out-of-tolerance comparison: " + field,
        )


def verify_artifacts(case, root, source_path):
    """Native reconstruction validates source/graph/model/facade without Contexts."""
    import json

    from island.evaluation.pcff_identity import pcff_evaluation_identity
    from island.forcefields import PreparedForceFieldSources, load_prepared_forcefield
    from island.forcefields.pcff import load_pcff_source

    source = load_pcff_source(source_path)
    declaration = json.loads((root / "declaration.json").read_text())
    require(
        declaration["source_identity"] == source.identity, "Declaration/source mismatch"
    )
    bundle = root / "dichlorine/bundle"
    actual = {
        str(p.relative_to(bundle)): checksum(p.read_bytes())
        for p in bundle.rglob("*")
        if p.is_file()
    }
    require(
        set(actual) == BUNDLE_FILES and actual == case["bundle_hashes"],
        "Missing/altered bundle artifacts",
    )
    loaded = load_prepared_forcefield(
        bundle, sources=PreparedForceFieldSources(pcff_frc=source_path)
    )
    prepared = loaded.prepared
    require(prepared.identity == case["bundle_identity"], "Bundle identity mismatch")
    require(
        prepared.native_result.identity == case["model_identity"],
        "Native model identity mismatch",
    )
    system = loaded.system
    require(
        len(system.topology.sites) == 2
        and all(a.element == "Cl" for a in system.topology.sites.values())
        and len(system.topology.bonds) == 1
        and all(b.order == 1 for b in system.topology.bonds.values()),
        "Wrong dichlorine graph",
    )
    model_file = root / "dichlorine/model.json"
    require(
        model_file.read_bytes() == (bundle / "model.json").read_bytes(),
        "Case/bundle model mismatch",
    )
    require(
        (root / "dichlorine/parameters.json").read_bytes()
        == (bundle / "assignment.json").read_bytes(),
        "Case/bundle assignment mismatch",
    )
    # Both declared whole-system reference geometries must actually be retained.
    for frame in (0, 1):
        reference = root / f"dichlorine/reference{frame}"
        for filename in (
            "data.lmp",
            "in.lmp",
            "energy.txt",
            "forces.dump",
            "log.lammps",
        ):
            require(
                (reference / filename).is_file()
                and (reference / filename).stat().st_size > 0,
                "Missing whole-system reference artifact: " + str(reference / filename),
            )
    require(
        (root / "dichlorine/reference0/data.lmp").read_bytes()
        != (root / "dichlorine/reference1/data.lmp").read_bytes(),
        "Duplicate comparison geometries",
    )
    _, parameter, model = pcff_evaluation_identity(
        prepared.native_result, system=system
    )
    return {
        "source": source.identity,
        "bundle_identity": prepared.identity,
        "parameter_fingerprint": parameter,
        "model_fingerprint": model,
    }


def assignment_and_term_gate(report, root, source_path):
    """Return a durable diagnostic instead of promoting incomplete evidence."""
    try:
        require(
            report.get("source_unchanged") is True,
            "Source integrity check failed/missing",
        )
        require(
            "term_failure" not in report and "error" not in report,
            "Execution failure in report",
        )
        case = mandatory_case(report)
        require(case.get("smiles") == "ClCl", "Dichlorine identity mismatch")
        require(
            "error" not in case and "failure" not in case, "Dichlorine execution error"
        )
        require(
            case.get("status") == "complete"
            and all(
                case.get(k) is True
                for k in ("typed", "charges_complete", "model_complete")
            ),
            "Incomplete dichlorine preparation",
        )
        comparisons = case.get("numerical")
        require(
            type(comparisons) is list and len(comparisons) == 2,
            "Two whole-system comparisons required",
        )
        for index, comparison in enumerate(comparisons):
            successful_comparison(comparison)
            require(
                comparison.get("frame", index) == index,
                "Wrong/duplicate comparison frame",
            )
        for key in ("model_identity", "bundle_identity"):
            value = case.get(key)
            require(
                type(value) is str
                and len(value) == 64
                and all(c in "0123456789abcdef" for c in value),
                "Missing/invalid " + key,
            )
        require(
            type(case.get("bundle_hashes")) is dict
            and set(case["bundle_hashes"]) == BUNDLE_FILES,
            "Bundle publication incomplete",
        )
        terms = report.get("term_checks")
        require(
            type(terms) is list and len(terms) == len(TERM_KEYS),
            "Distinct term checks required",
        )
        keys = []
        for term in terms:
            successful_comparison(term)
            keys.append((term.get("family"), term.get("kind")))
            if term.get("kind") == REAL:
                require(
                    term.get("selection", {}).get("status") == "assigned"
                    and type(term.get("source_row")) is dict,
                    "Missing term source selection",
                )
                errors = term.get("finite_difference_errors")
                require(
                    type(errors) is list
                    and len(errors) == 3
                    and all(
                        type(x) in (int, float) and math.isfinite(x) and x >= 0
                        for x in errors
                    )
                    and errors[-1] <= 1e-5,
                    "Missing/failed finite-difference check",
                )
        require(set(keys) == TERM_KEYS, "Missing/duplicate declared term check")
        verified = verify_artifacts(case, Path(root), source_path)
        return {"passed": True, "verified": verified, "diagnostics": []}
    except Exception as exc:  # noqa: BLE001 -- fail closed, retain actionable gate reason
        return {"passed": False, "diagnostics": [f"{type(exc).__name__}: {exc}"]}


def corrected_case(baseline, checked):
    """Resolve by name and bind correction to the original declared scientific case."""
    original = mandatory_case(baseline)
    corrected = mandatory_case(checked)
    require(
        baseline.get("inventory") == checked.get("inventory"),
        "Corrected source inventory mismatch",
    )
    for key in ("case", "smiles", "model_identity"):
        require(
            key in original and original[key] == corrected.get(key),
            "Corrected case mismatch: " + key,
        )
    return corrected


def operational_gate(checked_gate, workflow, historical):
    if not checked_gate["passed"]:
        return False
    verified = checked_gate["verified"]
    return (
        workflow.get("passed") is True
        and "error" not in workflow
        and workflow.get("inputs_unchanged") is True
        and workflow.get("rng_equal") is True
        and workflow.get("parameter_fingerprint") == verified["parameter_fingerprint"]
        and workflow.get("model_fingerprint") == verified["model_fingerprint"]
        and historical.get("unchanged") is True
    )
