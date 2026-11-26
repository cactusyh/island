"""Bounded, checkpointed NPT stationarity diagnostics for polymer recipes."""

import json
from dataclasses import dataclass, field

import numpy as np

from island.analysis.equilibration import EquilibrationPolicy, series
from island.forcefields.emc_graph import _section


@dataclass(frozen=True)
class RecipeEquilibration:
    max_segments: int = 10
    steps_per_segment: int = 20000
    consecutive_passes: int = 2
    policy: dict = field(default_factory=dict)
    maximum_bond_length_A: float = 3.5
    pressure_tolerance_atm: float = 1.0
    production_on_exhaustion: bool = False

    def __post_init__(self):
        for name in ("max_segments", "steps_per_segment", "consecutive_passes"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError(name + " must be a positive integer")
        if self.consecutive_passes > self.max_segments:
            raise ValueError("Consecutive passes exceed the segment budget")
        for name in ("maximum_bond_length_A", "pressure_tolerance_atm"):
            value = getattr(self, name)
            if type(value) not in (int, float) or not np.isfinite(value) or value <= 0:
                raise ValueError(name + " must be finite and positive")
        if type(self.production_on_exhaustion) is not bool:
            raise ValueError("production_on_exhaustion must be boolean")
        EquilibrationPolicy(**self.policy)


def frames(path):
    with path.open() as handle:
        for line in handle:
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError("Invalid equilibration trajectory")
            step = int(next(handle))
            if next(handle).strip() != "ITEM: NUMBER OF ATOMS":
                raise ValueError("Missing particle count")
            count = int(next(handle))
            if next(handle).strip() != "ITEM: BOX BOUNDS pp pp pp":
                raise ValueError("Equilibration requires an orthorhombic periodic box")
            bounds = np.array([next(handle).split() for _ in range(3)], float)
            columns = next(handle).split()[2:]
            if columns != ["id", "mol", "type", "xu", "yu", "zu", "vx", "vy", "vz"]:
                raise ValueError("Invalid equilibration trajectory columns")
            rows = np.array([next(handle).split() for _ in range(count)], float)
            if (
                not np.array_equal(rows[:, 0], np.arange(1, count + 1))
                or not np.isfinite(rows).all()
            ):
                raise ValueError("Invalid trajectory IDs or coordinates")
            yield step, rows[:, 3:6], bounds[:, 1] - bounds[:, 0]


def structure(xyz, box, masses, bonds):
    """Bonded-path unwrapping; winding networks have no finite-chain Rg."""
    neighbors = [[] for _ in xyz]
    lengths = []
    for i, j in bonds:
        neighbors[i].append(j)
        neighbors[j].append(i)
        dv = xyz[j] - xyz[i]
        dv -= np.rint(dv / box) * box
        lengths.append(float(np.linalg.norm(dv)))
    if not lengths:
        raise ValueError("Polymer diagnostics require bonds")
    seen, radii, closures = set(), [], []
    for root in range(len(xyz)):
        if root in seen:
            continue
        positions = {root: xyz[root]}
        stack = [root]
        while stack:
            i = stack.pop()
            for j in neighbors[i]:
                dv = xyz[j] - xyz[i]
                dv -= np.rint(dv / box) * box
                candidate = positions[i] + dv
                if j not in positions:
                    positions[j] = candidate
                    stack.append(j)
                else:
                    delta = (candidate - positions[j]) / box
                    if np.max(np.abs(delta)) > 1e-6:
                        if not np.allclose(delta, np.rint(delta), atol=1e-6):
                            raise ValueError("Inconsistent periodic graph closure")
                        closures.append(np.rint(delta).astype(int))
        ids = list(positions)
        seen.update(ids)
        coordinates = np.array([positions[i] for i in ids])
        center = np.average(coordinates, weights=masses[ids], axis=0)
        radii.append(
            float(
                np.average(
                    np.sum((coordinates - center) ** 2, axis=1), weights=masses[ids]
                )
            )
        )
    return {
        "mean_bond_length": float(np.mean(lengths)),
        "minimum_bond_length": min(lengths),
        "maximum_bond_length": max(lengths),
        "mean_rg2": float(np.mean(radii)) if not closures else None,
        "winding_rank": int(np.linalg.matrix_rank(closures)) if closures else 0,
        "components": len(radii),
    }


def observations(directory, md, steps):
    data = (directory / "system.data").read_text()
    mass_table = {
        int(r[0]): float(r[1])
        for line in _section(data, "Masses").splitlines()
        if (r := line.split())
    }
    atoms = sorted(
        [line.split() for line in _section(data, "Atoms").splitlines()],
        key=lambda r: int(r[0]),
    )
    masses = np.array([mass_table[int(r[2])] for r in atoms])
    bonds = [
        tuple(int(x) - 1 for x in line.split()[2:4])
        for line in _section(data, "Bonds").splitlines()
    ]
    thermo = np.loadtxt(directory / "thermo.dat", ndmin=2)
    if thermo.shape[1] != 6 or not np.isfinite(thermo).all():
        raise ValueError("Invalid equilibration thermodynamic samples")
    thermo = thermo[thermo[:, 0] != 0]
    expected = list(range(md.sample_every, steps + 1, md.sample_every))
    if (
        thermo.shape != (len(expected), 6)
        or not np.isfinite(thermo).all()
        or thermo[:, 0].tolist() != expected
    ):
        raise ValueError("Missing/nonfinite equilibration thermodynamic samples")
    configurations = {
        step: (xyz, box)
        for step, xyz, box in frames(directory / "trajectory.lammpstrj")
        if step > 0
    }
    if set(configurations) != set(expected):
        raise ValueError("Thermodynamic and structural samples do not align")
    rows = []
    for step, temperature, potential, pressure, density, volume in thermo:
        xyz, box = configurations[int(step)]
        if len(xyz) != len(masses) or not np.isclose(np.prod(box), volume, rtol=1e-6):
            raise ValueError("Trajectory and thermodynamic cell/atom mismatch")
        rows.append(
            {
                "step": int(step),
                "temperature": temperature,
                "potential": potential,
                "pressure": pressure,
                "density": density,
                **structure(xyz, box, masses, bonds),
            }
        )
    return rows


def assess(rows, cfg, temperature, pressure, *, network=False):
    policy = EquilibrationPolicy(**cfg.policy)
    keys = ["temperature", "potential", "density", "mean_bond_length"]
    if not network:
        keys.append("mean_rg2")
    if any(r.get(k) is None for r in rows for k in keys):
        return {
            "passed": False,
            "failures": ["finite_chain_winding_detected"],
            "samples": len(rows),
        }
    metrics = {key: series([r[key] for r in rows], policy) for key in keys}
    failures = [
        key + ":insufficient_sampling_or_drift"
        for key, value in metrics.items()
        if value["status"] != "ok"
    ]
    thermal = metrics["temperature"]
    if (
        "mean" not in thermal
        or abs(thermal["mean"] / temperature - 1) > policy.maximum_temperature_error
    ):
        failures.append("temperature_target")
    p = series([r["pressure"] for r in rows], policy)
    # Pressure's mean can be zero; use absolute drift/uncertainty, not relative drift.
    pressure_ok = (
        p.get("retained", 0) >= policy.minimum_frames
        and p.get("effective_sample_size", 0) >= policy.minimum_ess
        and p.get("block_independence", False)
        and abs(p["mean"] - pressure)
        <= cfg.pressure_tolerance_atm + policy.replica_z * p["standard_error"]
        and abs(p["drift"])
        <= cfg.pressure_tolerance_atm + policy.replica_z * p["standard_error"]
    )
    p["target_and_sampling_passed"] = bool(pressure_ok)
    metrics["pressure"] = p
    if not pressure_ok:
        failures.append("pressure_target_or_sampling")
    if any(
        r["minimum_bond_length"] < 0.2
        or r["maximum_bond_length"] > cfg.maximum_bond_length_A
        for r in rows
    ):
        failures.append("bond_geometry")
    if len({r["winding_rank"] for r in rows}) != 1:
        failures.append("periodic_winding_changed")
    if len({r["components"] for r in rows}) != 1:
        failures.append("component_count_changed")
    return {
        "passed": not failures,
        "samples": len(rows),
        "metrics": metrics,
        "failures": failures,
        "evidence_level": "single_run_stationarity",
        "physical_validation": False,
    }


def run(root, recipe, source, lammps, *, stop_after=None):
    from .recipe import checkpoint, digest, write_json
    from .recipe_md import run_segment

    cfg = recipe.equilibration
    history, outputs, streak = [], [], 0
    decision = None
    for index in range(cfg.max_segments):
        name = f"equilibration-{index + 1:04d}"
        parent_hash = digest(source / "checkpoint.json")
        with checkpoint(root, name) as (work, exists):
            if not exists:
                run_segment(
                    source,
                    work,
                    recipe.md,
                    cfg.steps_per_segment,
                    "npt",
                    recipe.seed,
                    lammps,
                    restart=True,
                )
                sampled = observations(work, recipe.md, cfg.steps_per_segment)
                write_json(work / "observations.json", sampled)
                decision = assess(
                    history + sampled,
                    cfg,
                    recipe.md.temperature_k,
                    recipe.md.pressure_atm,
                    network=recipe.architecture == "crosslinked",
                )
                write_json(
                    work / "diagnostic.json",
                    {"parent_checkpoint_sha256": parent_hash, "assessment": decision},
                )
            sampled = json.loads((work / "observations.json").read_text())
            record = json.loads((work / "diagnostic.json").read_text())
            if record["parent_checkpoint_sha256"] != parent_hash:
                raise ValueError("Equilibration parent checkpoint changed")
            decision = assess(
                history + sampled,
                cfg,
                recipe.md.temperature_k,
                recipe.md.pressure_atm,
                network=recipe.architecture == "crosslinked",
            )
            if decision != record["assessment"]:
                raise ValueError(
                    "Equilibration diagnostic differs from saved samples/policy"
                )
        history.extend(sampled)
        source = root / name
        outputs.append(json.loads((source / "result.json").read_text()))
        streak = streak + 1 if decision["passed"] else 0
        converged = streak >= cfg.consecutive_passes
        report = {
            "converged": converged,
            "segments": index + 1,
            "steps": (index + 1) * cfg.steps_per_segment,
            "consecutive_passes": streak,
            "assessment": decision,
            "stop_reason": "diagnostics_passed"
            if converged
            else "budget_exhausted"
            if index + 1 == cfg.max_segments
            else "continue",
        }
        write_json(root / "equilibration.json", report)
        if stop_after == "equilibration":
            return report, source, outputs, True
        if converged:
            break
    return report, source, outputs, False
