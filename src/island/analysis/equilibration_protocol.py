"""Additive, fail-closed J40-C1 protocol; no historical executor changes."""

from dataclasses import asdict, dataclass

from island.analysis.equilibration import EquilibrationPolicy, diagnose, require
from island.analysis.equilibration_run import HeatingSchedule
from island.analysis.npt import FLAGS, _body, _digest, _signed


@dataclass(frozen=True)
class CorrectiveProtocol:
    """Predeclared finite budget, not a promise of thermodynamic equilibrium.

    0.5 fs steps: 10 ps heating, 10 ps NVT, 100 ps NPT per evaluation,
    at most 1 ns NPT. Existing statistical bounds are not changed.
    """

    schema: str = "island.corrective-equilibration-protocol.v1"
    timestep_fs: float = 0.5
    heating: tuple = ((5000, 100.0), (10000, 200.0), (20000, 300.0))
    nvt_steps: int = 40000
    npt_window_steps: int = 200000
    maximum_windows: int = 10
    frame_interval: int = 100
    replicas: int = 3
    pressure_interval: int = 25
    target_temperature_kelvin: float = 300.0
    target_pressure_bar: float = 1.0
    minimization: str = "existing_J23_minimizer"
    burn_in: str = "fixed_first_quarter_of_each_cumulative_NPT_history"
    stopping_rule: str = "all_replica_and_reference_gates_pass_or_budget_exhausted"

    def __post_init__(self):
        # This version is intentionally frozen: a changed protocol needs a new version.
        expected = {
            "timestep_fs": 0.5,
            "heating": ((5000, 100.0), (10000, 200.0), (20000, 300.0)),
            "nvt_steps": 40000,
            "npt_window_steps": 200000,
            "maximum_windows": 10,
            "frame_interval": 100,
            "replicas": 3,
            "pressure_interval": 25,
            "target_temperature_kelvin": 300.0,
            "target_pressure_bar": 1.0,
            "schema": "island.corrective-equilibration-protocol.v1",
            "minimization": "existing_J23_minimizer",
            "burn_in": "fixed_first_quarter_of_each_cumulative_NPT_history",
            "stopping_rule": "all_replica_and_reference_gates_pass_or_budget_exhausted",
        }
        require(asdict(self) == expected, "Changed versioned corrective protocol")

    @property
    def payload(self):
        return _signed(
            dict(**asdict(self), statistical_policy=EquilibrationPolicy().payload)
        )


def prerequisites(reference, finite_difference, pressure):
    """No dynamics may start before reference and measurement gates are corrected."""
    failures = []
    if reference.get("passed") is not True:
        failures.append("cross_engine_discrepancy")
    if finite_difference.get("passed") is not True:
        failures.append("finite_difference_interval_not_converged")
    # A barostat target or a self-declared flag is not an observable implementation.
    if pressure != {
        "status": "measured_and_independently_validated",
        "includes_long_range_and_constraint_virial": True,
        "reference_agreement": True,
    }:
        failures.append("pressure_unavailable")
    return failures


def assess(workflows, *, reference, finite_difference, pressure):
    diagnostic = diagnose(workflows)  # retains unchanged J40 numerical criteria
    failures = prerequisites(reference, finite_difference, pressure)
    if len(workflows) < CorrectiveProtocol().replicas:
        failures.append("insufficient_independent_replicas_for_corrective_protocol")
    for i, replica in enumerate(diagnostic["replicas"]):
        for name, metric in replica["metrics"].items():
            if metric["status"] != "ok":
                failures.append(f"replica_{i}_{name}_{metric['status']}")
        if not replica["metrics"]["temperature"]["target_agreement"]:
            failures.append(f"replica_{i}_temperature_target_unmet")
    if any(x["status"] != "ok" for x in diagnostic["replica_agreement"].values()):
        failures.append("replica_agreement_insufficient_sampling")
    return _signed(
        dict(
            schema="island.corrective-equilibration-assessment.v1",
            protocol=CorrectiveProtocol().payload,
            analysis=diagnostic,
            reference_identity=_digest(reference),
            finite_difference_identity=_digest(finite_difference),
            pressure=pressure,
            failed_stages=failures,
            scientific_acceptance=not failures,
            **FLAGS,
        )
    )


def validate_assessment(record, workflows, **inputs):
    require(
        record.get("identity") == _digest(_body(record)), "Changed assessment identity"
    )
    require(record == assess(workflows, **inputs), "Stale corrective assessment")


def run_until_stable(run_window, *, reference, finite_difference, pressure):
    """Bounded scheduler, independent replicas supplied by a real engine adapter.

    The callback must persist cumulative real replica workflows after each window.
    No synthetic data or replay of another replica is permitted. This scheduler
    deliberately has no default pressure estimator.
    """
    failures = prerequisites(reference, finite_difference, pressure)
    require(not failures, "Dynamics blocked before execution: " + ", ".join(failures))
    protocol = CorrectiveProtocol()
    for window in range(1, protocol.maximum_windows + 1):
        workflows = run_window(protocol, window)
        require(
            len(workflows) == protocol.replicas, "Three independent replicas required"
        )
        for workflow in workflows:
            config = workflow["manifest"]["config"]
            require(
                workflow["manifest"]["heating"]
                == HeatingSchedule(
                    protocol.heating,
                    protocol.nvt_steps,
                    protocol.target_temperature_kelvin,
                ).payload,
                "Changed staged heating protocol",
            )
            for key, value in {
                "timestep_fs": protocol.timestep_fs,
                "temperature_kelvin": protocol.target_temperature_kelvin,
                "pressure_bar": protocol.target_pressure_bar,
                "pressure_interval": protocol.pressure_interval,
                "checkpoint_interval": protocol.frame_interval,
                "platform": "Reference",
                "thermostat_type": "LangevinMiddleIntegrator",
                "barostat_type": "MonteCarloBarostat",
            }.items():
                require(config[key] == value, "Changed protocol setting: " + key)
            require(
                config["npt_steps"] == window * protocol.npt_window_steps,
                "Window is not the declared cumulative NPT history",
            )
            require(
                config["nvt_steps"] == protocol.nvt_steps,
                "Staged protocol NVT mismatch",
            )
        result = assess(
            workflows,
            reference=reference,
            finite_difference=finite_difference,
            pressure=pressure,
        )
        if result["scientific_acceptance"]:
            return result
    return result  # insufficient samples remain blocked at the finite budget
