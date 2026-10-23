"""One strict J38 matrix contract; diagnostic readability is not acceptance."""

import math

from island.analysis.npt import FLAGS, _body, _digest, _signed, frame_measurements
from scripts.phase_4j38_fixture import SPECS

BASE = "fb4dc8bed7a4e69e42c3186a22a5c4e5da39906b"
CASES = tuple(f"{family}/{key}" for family in ("oplsaa", "gaff2") for key in SPECS)
GATES = (
    "software_execution",
    "source_completeness",
    "numerical_reference",
    "finite_time_stability",
    "thermodynamic_equilibration",
    "scientific_validation",
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def validate_identity(value):
    require(value.get("identity") == _digest(_body(value)), "Changed record identity")


def validate_row(key, row):
    require(key in CASES and row.get("case") == key, "Case inventory/identity mismatch")
    require(
        row.get("family") == key.split("/")[0]
        and row.get("benchmark") == key.split("/")[1],
        "Family/benchmark mismatch",
    )
    require(
        all(type(row.get(k)) is type(v) and row[k] == v for k, v in FLAGS.items()),
        "Readiness changed",
    )
    require(set(row["gates"]) == set(GATES), "Gate inventory mismatch")
    require(
        all(
            g["status"]
            in {
                "passed",
                "failed",
                "blocked",
                "insufficient_sampling",
                "not_executed",
                "asset_unavailable",
            }
            and g.get("reason")
            for g in row["gates"].values()
        ),
        "Missing explicit gate diagnostic",
    )
    require(
        row["gates"]["thermodynamic_equilibration"]["status"] == "insufficient_sampling"
        and row["gates"]["scientific_validation"]["status"] == "blocked",
        "Short acceptance run cannot establish equilibrium/scientific validation",
    )
    executed = row["gates"]["software_execution"]["status"] == "passed"
    require(row.get("execution_met") is executed, "Contradictory execution claim")
    if executed:
        for binding in (
            "final_graph_identity",
            "packing_identity",
            "source_identity",
            "assignment_identity",
            "charge_identity",
            "typed_graph_identity",
            "workflow_identity",
            "trajectory_identity",
            "checkpoint_identity",
            "analysis_identity",
            "comparison_identity",
        ):
            v = row.get(binding)
            require(
                type(v) is str
                and len(v) == 64
                and all(c in "0123456789abcdef" for c in v),
                "Missing execution binding: " + binding,
            )
        require(
            row["frame_count"] >= 33 and row["density_blocks"] >= 4,
            "Missing frame/block evidence",
        )
        require(
            all(
                row.get(k) is True
                for k in (
                    "split_equal",
                    "checkpoint_byte_equal",
                    "rng_equal",
                    "relocation",
                    "fresh_process",
                    "box_volume_changed",
                )
            ),
            "Missing continuation evidence",
        )
        require(
            row.get("charge_origin")
            == ("RESP/Psi4" if key.startswith("gaff2/") else "OPLS-AA source"),
            "Wrong charge origin",
        )
        require(
            row["gates"]["source_completeness"]["status"] == "passed",
            "Executed without complete source",
        )
    else:
        require(
            not any(
                row.get(k)
                for k in (
                    "split_equal",
                    "checkpoint_byte_equal",
                    "rng_equal",
                    "fresh_process",
                )
            ),
            "Failed case has success evidence",
        )
        require(row.get("diagnostic"), "Failure diagnostic missing")
    expected = all(row["gates"][k]["status"] == "passed" for k in GATES[:4])
    require(row.get("software_acceptance_met") is expected, "Incorrect software gate")
    require(
        row.get("acceptance_met") is False,
        "Finite benchmark cannot prove all scientific gates",
    )
    validate_identity(row)


def validate_matrix(matrix):
    validate_identity(matrix)
    require(
        matrix.get("schema") == "island.j38-matrix.v1"
        and matrix["base_commit"] == BASE,
        "Wrong matrix schema/base",
    )
    require(set(matrix["cases"]) == set(CASES), "Incomplete benchmark/family inventory")
    for key, row in matrix["cases"].items():
        validate_row(key, row)
    require(matrix["acceptance_met"] is False, "Scientific acceptance unestablished")
    require(
        matrix["software_acceptance_met"]
        is all(r["software_acceptance_met"] for r in matrix["cases"].values()),
        "Matrix software gate mismatch",
    )
    require(
        all(
            type(matrix.get(k)) is type(v) and matrix[k] == v for k, v in FLAGS.items()
        ),
        "Matrix readiness changed",
    )
    return matrix


def gate(status, reason):
    return {"status": status, "reason": reason}


def base_row(key):
    family, benchmark = key.split("/")
    return dict(
        case=key,
        family=family,
        benchmark=benchmark,
        execution_met=False,
        software_acceptance_met=False,
        acceptance_met=False,
        gates={
            **{
                k: gate("not_executed", "No completed workflow evidence")
                for k in GATES[:4]
            },
            "thermodynamic_equilibration": gate(
                "insufficient_sampling",
                "One short acceptance trajectory; missing independent replicas and long-time convergence",
            ),
            "scientific_validation": gate(
                "blocked", "No experimental or production validation in this phase"
            ),
        },
        **FLAGS,
    )


def stable(frames, analysis):
    """Declared bounded software sanity criteria, not equilibrium criteria."""
    frames = frame_measurements(frames, analysis["total_mass_dalton"])
    npt = [f for f in frames if f["stage"] == "npt"]
    return (
        len(npt) >= 32
        and analysis["density"]["blocks"] >= 4
        and all(
            math.isfinite(f[k])
            for f in npt
            for k in (
                "density_g_cm3",
                "volume_nm3",
                "potential_energy_kj_mol",
                "kinetic_energy_kj_mol",
                "measured_temperature_kelvin",
            )
        )
        and all(0 < f["measured_temperature_kelvin"] < 2000 for f in npt)
        and max(f["density_g_cm3"] for f in npt) / min(f["density_g_cm3"] for f in npt)
        < 2
    )


def sign(row):
    return _signed(row)
