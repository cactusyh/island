"""Two-size, three-seed cure benchmark and functional-group-resolved aggregation."""

import json
from pathlib import Path

import numpy as np

from .cure import CureConfig, run_cure


def recipes():
    for n in [4, 8]:
        for seed in [41, 173, 2026]:
            yield n, n // 2, seed
    yield 4, 3, 41


def study(out, lib, *, lammps, packmol):
    out = Path(out)
    for n, p, seed in recipes():
        cfg = CureConfig(
            dgeba=n,
            pacm=p,
            seed=seed,
            pre_steps=4000,
            batch_steps=800,
            post_steps=4000,
            cutoff=5.5,
            batch_size=2,
            max_idle=6,
        )
        run_cure(out / f"d{n}-p{p}-s{seed}", lib, cfg, lammps=lammps, packmol=packmol)
    return summarize(out)


def summarize(directory):
    base = Path(directory)
    cases = []
    for n, p, seed in recipes():
        path = base / f"d{n}-p{p}-s{seed}"
        records = [
            json.loads(f.read_text()) for f in sorted(path.glob("stage-*/record.json"))
        ]
        if not records or records[-1]["phase"] != "post":
            raise ValueError(f"Incomplete cure: {path}")
        first = next(
            (
                r["statistics"]["epoxy_conversion"]
                for r in records
                if r["statistics"]["winding_rank"]
            ),
            None,
        )
        cases.append(
            {
                "dgeba": n,
                "pacm": p,
                "seed": seed,
                "config": json.loads((path / "config.json").read_text())["config"],
                "first_observed_winding": first,
                "reason": records[-1]["reason"],
                "statistics": records[-1]["statistics"],
                "history": [
                    {
                        **r["statistics"],
                        "phase": r["phase"],
                        "index": r["index"],
                        "box_before_selection_A": r.get("box_A"),
                        "candidate_count": r["candidate_count"],
                        "rejected": r["rejected"],
                        "thermo": r["thermo"],
                        "events": r["events"],
                    }
                    for r in records
                ],
            }
        )
    probability = []
    for n in [4, 8]:
        group = [c for c in cases if c["dgeba"] == n and n == 2 * c["pacm"]]
        for alpha in (0.25, 0.5, 0.625, 0.75, 0.875, 1.0):
            # First archived graph at or beyond the threshold. Report censored
            # runs separately; never treat a stalled run as a measured negative.
            samples = [
                next((s for s in c["history"] if s["epoxy_conversion"] >= alpha), None)
                for c in group
            ]
            observed = [s for s in samples if s is not None]
            probability.append(
                {
                    "dgeba": n,
                    "target": alpha,
                    "n_observed": len(observed),
                    "n_stalled": len(group) - len(observed),
                    "actual_conversions": [s["epoxy_conversion"] for s in observed],
                    "winding_probability": float(
                        np.mean([s["winding_rank"] > 0 for s in observed])
                    )
                    if observed
                    else None,
                }
            )
    return {"cases": cases, "probability": probability}


def property_report(study_dir, measurements, extended=None):
    """Aggregate paired volumes at identical T/P; keep individual rate/axis fits."""
    base, measurements = Path(study_dir), Path(measurements)
    pairs = []

    def read(path):
        return json.loads((path / "result.json").read_text())

    def volume_pair(n, p, seed, root):
        name = f"d{n}-p{p}-s{seed}"
        records = sorted((base / name).glob("stage-*/record.json"))
        stat = json.loads(records[-1].read_text())["statistics"]
        before = read(root / (name + "-uncured"))
        after = read(root / (name + "-cold"))
        if any(
            before[k] != after[k] for k in ("temperature_K", "pressure_atm", "steps")
        ):
            raise ValueError("Unmatched T/P or sampling protocol for volume shrinkage")
        return {
            "dgeba": n,
            "pacm": p,
            "seed": seed,
            "temperature_K": before["temperature_K"],
            "pressure_atm": before["pressure_atm"],
            "steps": before["steps"],
            "epoxy_conversion": stat["epoxy_conversion"],
            "density_before": before["mean_density"],
            "density_after": after["mean_density"],
            "shrinkage_percent": 100
            * (1 - after["mean_volume"] / before["mean_volume"]),
        }

    for n, p, seed in recipes():
        pairs.append(volume_pair(n, p, seed, measurements))
    groups = []
    for n in [4, 8]:
        group = [s for s in pairs if s["dgeba"] == n and n == 2 * s["pacm"]]
        values = {
            k: [s[k] for s in group] for k in ("density_after", "shrinkage_percent")
        }
        groups.append(
            dict(
                dgeba=n,
                pacm=n // 2,
                n_seeds=len(group),
                **{
                    k: {
                        "mean": float(np.mean(v)),
                        "sample_sd": float(np.std(v, ddof=1)),
                    }
                    for k, v in values.items()
                },
            )
        )
    fits = {}
    for name, root in [("d4-p2-s173", measurements), ("d8-p4-s173", extended)]:
        if root is None:
            continue
        paths = sorted(Path(root).glob("cool-*/result.json")) + sorted(
            Path(root).glob("pull-*/result.json")
        )
        if len(paths) != 8:
            raise ValueError(
                "Expected two cooling rates and two tensile rates across x, y, z"
            )
        fits[name] = [json.loads(path.read_text()) for path in paths]
    report = {"volume_pairs": pairs, "by_size": groups, "property_fits": fits}
    if extended is not None:
        report["extended_volume_pair"] = volume_pair(8, 4, 173, Path(extended))
    return report


def characterize(study_dir, out, lammps, extended=False):
    """Reproduce the matched 300 K volume panel and two cooling/pull rates."""
    from .properties import conditions, cooling, tension

    base, out = Path(study_dir), Path(out)
    out.mkdir(parents=True, exist_ok=False)

    def store(fun, dst, *args, **opts):
        r = fun(args[0], dst, lammps, *args[1:], **opts)
        (dst / "result.json").write_text(json.dumps(r, indent=2, allow_nan=False))

    for n, p, seed in recipes():
        name = f"d{n}-p{p}-s{seed}"
        stages = sorted((base / name).glob("stage-*"))
        if (
            not stages
            or json.loads((stages[-1] / "record.json").read_text())["phase"] != "post"
        ):
            raise ValueError(f"Incomplete cure: {name}")
        store(conditions, out / (name + "-uncured"), stages[0], 300, 12000)
        store(conditions, out / (name + "-cold"), stages[-1], 300, 12000)

    def sweeps(root, name, steps, rates, pulls):
        last = max((base / name).glob("stage-*"))
        hot, cold = root / (name + "-hot"), root / (name + "-cold")
        store(conditions, hot, last, 550, steps)
        if not cold.exists():
            store(conditions, cold, last, 300, steps)
        for rate in rates:
            store(cooling, root / f"cool-{rate}", hot, rate)
        for rate in pulls:
            for axis in "xyz":
                store(tension, root / f"pull-{rate}{axis}", cold, rate, axis)

    sweeps(out, "d4-p2-s173", 12000, [10, 20], [0.002, 0.005])
    if extended:
        sweeps(out / "extended", "d8-p4-s173", 60000, [2, 4], [0.0005, 0.001])
        store(
            conditions,
            out / "extended" / "d8-p4-s173-uncured",
            base / "d8-p4-s173" / "stage-0000",
            300,
            60000,
        )
    return property_report(base, out, out / "extended" if extended else None)
