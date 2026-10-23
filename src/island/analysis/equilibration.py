"""Finite-run diagnostics; passing these criteria is not physical validation."""

from dataclasses import asdict, dataclass

import numpy as np

from island.analysis.npt import FLAGS, _body, _digest, _signed, frame_measurements
from island.exceptions import AnalysisError


def require(ok, message):
    if not ok:
        raise AnalysisError(message)


@dataclass(frozen=True)
class EquilibrationPolicy:
    burn_fraction: float = 0.25
    blocks: int = 4
    minimum_frames: int = 512
    minimum_ess: float = 50.0
    block_tau_ratio: float = 5.0
    maximum_relative_drift: float = 0.02
    maximum_temperature_error: float = 0.05
    replica_z: float = 2.0

    def __post_init__(self):
        require(0 <= self.burn_fraction < 0.5, "Invalid burn-in fraction")
        require(
            type(self.blocks) is int and self.blocks >= 4,
            "At least four blocks required",
        )
        require(
            type(self.minimum_frames) is int and self.minimum_frames >= self.blocks * 2,
            "Invalid sample floor",
        )
        for key in (
            "minimum_ess",
            "block_tau_ratio",
            "maximum_relative_drift",
            "maximum_temperature_error",
            "replica_z",
        ):
            value = getattr(self, key)
            require(
                type(value) in (int, float) and np.isfinite(value) and value > 0,
                "Invalid diagnostic threshold",
            )

    @property
    def payload(self):
        return _signed(dict(schema="island.equilibration-policy.v1", **asdict(self)))


def series(values, policy):
    """Biased ACF with initial-positive paired sequence and contiguous blocks.

    Standard error uses the larger of block and autocorrelation estimates.
    Burn-in is fixed in advance, never selected to manufacture convergence.
    """
    x = np.asarray(values, dtype=float)
    require(x.ndim == 1 and np.isfinite(x).all(), "Nonfinite diagnostic series")
    start = int(len(x) * policy.burn_fraction)
    x = x[start:]
    if len(x) < policy.blocks * 2:
        return {
            "status": "insufficient_sampling",
            "retained": len(x),
            "burn_in_frames": start,
        }
    centered = x - x.mean()
    variance = float(np.mean(centered**2))
    acf = [1.0]
    tau = 1.0
    if variance > 0:
        acf += [
            float(np.dot(centered[:-lag], centered[lag:]) / (len(x) * variance))
            for lag in range(1, len(x) // 2)
        ]
        for lag in range(1, len(acf) - 1, 2):
            pair = acf[lag] + acf[lag + 1]
            if pair <= 0:
                break
            tau += 2 * pair
    size = len(x) // policy.blocks
    means = np.array([chunk.mean() for chunk in np.array_split(x, policy.blocks)])
    block_sd = float(np.std(means, ddof=1))
    ess = float(len(x) / max(1, tau))
    se = max(block_sd / np.sqrt(policy.blocks), np.sqrt(variance / ess))
    scale = max(abs(float(x.mean())), np.finfo(float).tiny)
    drift = float(np.polyfit(np.arange(len(x)), x, 1)[0] * (len(x) - 1))
    independent = variance > 0 and size >= policy.block_tau_ratio * tau
    passed = (
        len(x) >= policy.minimum_frames
        and ess >= policy.minimum_ess
        and independent
        and abs(drift) / scale <= policy.maximum_relative_drift
    )
    return {
        "status": "ok" if passed else "insufficient_sampling",
        "burn_in_frames": start,
        "retained": len(x),
        "mean": float(x.mean()),
        "standard_deviation": float(x.std(ddof=1)),
        "block_means": means.tolist(),
        "block_mean_standard_deviation": block_sd,
        "standard_error": float(se),
        "autocorrelation": acf,
        "integrated_autocorrelation_frames": tau,
        "effective_sample_size": ess,
        "block_size_min": size,
        "block_independence": bool(independent),
        "drift": drift,
        "relative_drift": abs(drift) / scale,
        "half_mean_difference": float(
            x[len(x) // 2 :].mean() - x[: len(x) // 2].mean()
        ),
        "rolling_mean": np.convolve(x, np.ones(size) / size, mode="valid").tolist(),
    }


def diagnose(workflows, policy=None):
    policy = policy or EquilibrationPolicy()
    require(len(workflows) >= 2, "Independent replicas required")
    configs = [w["manifest"]["config"] for w in workflows]
    for key in ("seed", "velocity_seed", "barostat_seed"):
        require(
            len({c[key] for c in configs}) == len(configs),
            "Replica seeds must be distinct",
        )
    comparable = lambda c: {
        k: v
        for k, v in c.items()
        if k not in ("identity", "seed", "velocity_seed", "barostat_seed")
    }
    require(
        all(comparable(c) == comparable(configs[0]) for c in configs),
        "Replica configurations differ beyond seeds",
    )
    require(
        len({w["manifest"]["backend_identity"] for w in workflows}) == 1,
        "Replica backend mismatch",
    )
    replicas = []
    for w in workflows:
        frames = frame_measurements(w["frames"], w["manifest"]["total_mass_dalton"])
        npt = [f for f in frames if f["stage"] == "npt"]
        require(npt, "NPT did not execute")
        metrics = {
            name: series([f[key] for f in npt], policy)
            for name, key in (
                ("density", "density_g_cm3"),
                ("volume", "volume_nm3"),
                ("temperature", "measured_temperature_kelvin"),
            )
        }
        pressure = [f["measured_pressure_bar"] for f in npt]
        metrics["pressure"] = (
            {
                "status": "blocked",
                "reason": "Backend retains no measured pressure; target is not an observation",
            }
            if any(p is None for p in pressure)
            else series(pressure, policy)
        )
        temp = metrics["temperature"].get("mean")
        metrics["temperature"]["target_agreement"] = (
            temp is not None
            and abs(temp / configs[0]["temperature_kelvin"] - 1)
            <= policy.maximum_temperature_error
        )
        replicas.append(
            {"trajectory_identity": _digest(w["frames"]), "metrics": metrics}
        )
    agreement = {}
    for name in ("density", "volume", "temperature"):
        ms = [r["metrics"][name] for r in replicas]
        okay = all(m.get("status") == "ok" for m in ms)
        pairs = []
        for i, a in enumerate(ms):
            for b in ms[i + 1 :]:
                delta = abs(a.get("mean", 0) - b.get("mean", 0))
                uncertainty = float(
                    np.hypot(a.get("standard_error", 0), b.get("standard_error", 0))
                )
                pairs.append(
                    {"mean_difference": delta, "combined_standard_error": uncertainty}
                )
                okay &= delta <= policy.replica_z * uncertainty
        agreement[name] = {
            "status": "ok" if okay else "insufficient_sampling",
            "comparisons": pairs,
        }
    passed = (
        all(m["status"] == "ok" for r in replicas for m in r["metrics"].values())
        and all(r["metrics"]["temperature"]["target_agreement"] for r in replicas)
        and all(a["status"] == "ok" for a in agreement.values())
    )
    return _signed(
        dict(
            schema="island.equilibration-analysis.v1",
            policy=policy.payload,
            replicas=replicas,
            replica_agreement=agreement,
            equilibration_status="converged_under_declared_finite_run_criteria"
            if passed
            else "insufficient_sampling",
            acceptance_met=passed,
            scientific_validation_status="not_established",
            **FLAGS,
        )
    )


def validate_analysis(record, workflows):
    require(
        record.get("identity") == _digest(_body(record)), "Diagnostic identity changed"
    )
    raw = record["policy"]
    require(
        raw.get("identity") == _digest(_body(raw))
        and raw.get("schema") == "island.equilibration-policy.v1",
        "Policy identity changed",
    )
    policy = EquilibrationPolicy(
        **{k: v for k, v in raw.items() if k not in ("schema", "identity")}
    )
    require(
        record == diagnose(workflows, policy),
        "Stale or fabricated equilibration diagnostic",
    )
    return record
