"""Paired finite-size and observation-time statistics for a matched-conversion panel."""

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import t


def estimates(x):
    """Student-t interval across independent packing seeds."""
    v = np.asarray(x, float)
    if not len(v):
        return {"n": 0, "mean": None, "sample_sd": None, "ci95": None}
    mean = float(v.mean())
    sd = float(v.std(ddof=1)) if len(v) > 1 else None
    interval = None
    if sd is not None:
        delta = float(t.ppf(0.975, len(v) - 1)) * sd / np.sqrt(len(v))
        interval = [float(mean - delta), float(mean + delta)]
    return {"n": len(v), "mean": mean, "sample_sd": sd, "ci95": interval}


def summarize(records):
    groups = defaultdict(list)
    used = set()
    for case in records:
        n, p, seed = (case[k] for k in ("dgeba", "pacm", "seed"))
        if n != 2 * p or (n, seed) in used:
            raise ValueError("Expected unique stoichiometric replicates")
        used.add((n, seed))
        groups[n].append(case)
    sizes = []
    for n, cases in sorted(groups.items()):
        matched = [c for c in cases if c["statistics"]["epoxy_conversion"] == 0.75]
        observations = defaultdict(lambda: defaultdict(list))
        measured = [c for c in matched if "conditioning_failure" not in c]
        for c in measured:
            before = c["volumes"]["before"]
            after = c["volumes"]["after"]
            if set(before) != set(after) or "12000" not in before:
                raise ValueError("Unpaired observation windows")
            for steps in after:
                observations[steps]["density_after"].append(
                    after[steps]["density_g_cm3"]
                )
                observations[steps]["shrinkage_percent"].append(
                    100 * (1 - after[steps]["volume_A3"] / before[steps]["volume_A3"])
                )
        windows = {}
        for steps, values in observations.items():
            if len(values["density_after"]) != len(measured):
                raise ValueError("Incomplete time series")
            windows[steps] = {k: estimates(v) for k, v in values.items()}
            windows[steps]["duration_ps"] = int(steps) * 0.0005
            if steps != "12000":
                windows[steps]["change_from_6ps"] = {
                    k: estimates(np.array(v) - observations["12000"][k])
                    for k, v in values.items()
                }
        sizes.append(
            {
                "dgeba": n,
                "pacm": n // 2,
                "n_total": len(cases),
                "n_at_target": len(matched),
                "n_paired_measurements": len(measured),
                "seeds_with_paired_measurements": [c["seed"] for c in measured],
                "seeds_at_target": [c["seed"] for c in matched],
                "n_winding_at_target": sum(
                    c["statistics"]["winding_rank"] > 0 for c in matched
                ),
                "elapsed_worker_seconds": sum(c["elapsed_seconds"] for c in cases),
                "windows": windows,
            }
        )
    return {
        "target_conversion": 0.75,
        "temperature_K": 300,
        "pressure_atm": 1,
        "sizes": sizes,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_file", type=Path)
    args = parser.parse_args()
    print(json.dumps(summarize(json.loads(args.case_file.read_text())), indent=2))
