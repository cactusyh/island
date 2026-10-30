"""Independent claims. This policy never grants execution authority from JSON."""

from dataclasses import dataclass

from island.forcefields.opls_registry.api import require
from island.forcefields.selection.inventory import canonical, digest

SCHEMA = "island.parameter-completeness-policy.v1"
CONVENTIONS = (
    "chemical_domain",
    "charge_origin",
    "units",
    "functional_forms",
    "torsion_convention",
    "exclusions",
    "coulomb_14",
    "lj_14",
    "lj_mixing",
    "cutoff",
    "long_range",
    "provenance",
    "source_hashes",
    "evidence",
)
# Source absence for these changes the physical Hamiltonian, unlike precision,
# interpolation and neighbor-list choices, or a pressure estimator.
PDMS_MODEL_GAPS = ("dielectric", "periodic_image_policy")


def assess(
    required, assigned, conventions, *, typing_complete, charge_complete, source_valid
):
    """Exact inventory check, independent of all numerical/engine outcomes.

    An explicit not_applicable row must carry evidence. Dependency failures are
    propagated rather than replaced by zero coefficients. This function produces
    a report only; source adapters and integrity checks own parameter authority.
    """
    required = sorted(required, key=lambda r: r["id"])
    ids = [r["id"] for r in assigned]
    require(len(ids) == len(set(ids)), "source_incompatible", "duplicate_interactions")
    rows = {r["id"]: r for r in assigned}
    inventory = []
    for req in required:
        row = rows.get(req["id"])
        good = row is not None and row.get("status") == "assigned"
        if row and row.get("status") == "not_applicable":
            good = bool(row.get("reason") or row.get("evidence"))
        inventory.append(
            {**req, "assignment": row, "status": "complete" if good else "incomplete"}
        )
    complete_ids = {r["id"] for r in inventory if r["status"] == "complete"}
    # Fixed point: a coupling cannot be complete with a missing equilibrium row.
    while True:
        bad = {
            r["id"]
            for r in inventory
            if r["status"] == "complete"
            and any(d not in complete_ids for d in r.get("dependencies", []))
        }
        if not bad:
            break
        complete_ids -= bad
        for r in inventory:
            if r["id"] in bad:
                r["status"] = "incomplete"
                r["reason"] = "missing_dependency"
    missing = [r for r in inventory if r["status"] != "complete"]
    convention_gaps = [
        name
        for name in CONVENTIONS
        if not conventions.get(name, {}).get("value")
        or not conventions[name].get("evidence")
    ]
    extra_gaps = sorted(
        name
        for name, row in conventions.items()
        if row.get("required", False)
        and (row.get("value") is None or not row.get("evidence"))
    )
    complete = (
        bool(required)
        and typing_complete
        and charge_complete
        and source_valid
        and not missing
        and not convention_gaps
        and not extra_gaps
    )
    value = {
        "parameter_completeness": "complete" if complete else "incomplete",
        "parameter_complete": bool(complete),
        "typing_status": "complete" if typing_complete else "incomplete",
        "charge_status": "complete" if charge_complete else "incomplete",
        "interaction_inventory": inventory,
        "missing_terms": missing,
        "missing_conventions": sorted(set(convention_gaps + extra_gaps)),
        "source_integrity": "passed" if source_valid else "blocked",
        "deterministic_inventory_identity": digest(inventory),
    }
    # Also checks round-trip representation (no NaN, ambiguous map keys, etc.).
    require(
        canonical(value) == canonical(__import__("json").loads(canonical(value))),
        "source_incompatible",
        "non_deterministic_serialization",
    )
    return value


@dataclass(frozen=True)
class ClaimRequirements:
    parameter_complete: bool = False
    lammps_runtime: bool = False
    openmm_runtime: bool = False
    cross_engine_equivalence: bool = False
    scientific_validation: bool = False

    def failures(self, report):
        wanted = {
            "parameter_complete": ("parameter_completeness", "complete"),
            "lammps_runtime": ("lammps_runtime_status", "passed"),
            "openmm_runtime": ("openmm_runtime_status", "passed"),
            "cross_engine_equivalence": ("cross_engine_equivalence", "passed"),
            "scientific_validation": ("scientific_validation", "passed"),
        }
        return [
            {
                "claim": name,
                "field": field,
                "required": expected,
                "observed": report.get(field, "not_run"),
            }
            for name, (field, expected) in wanted.items()
            if getattr(self, name) and report.get(field) != expected
        ]

    def enforce(self, report):
        failures = self.failures(report)
        require(
            not failures,
            "claim_not_satisfied",
            "independent_claim_requirements",
            failures=failures,
        )
