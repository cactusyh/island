"""Precision observables and bounded execution, reusing the J63 native capture.

This module neither chooses physical coefficients nor changes historical adapters.
Energy subtraction precedes conversion; decimal text and binary rounding are audited.
"""

import argparse
import fcntl
import json
import math
import os
import re
import resource
import time
from decimal import Decimal, localcontext
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, write
from island.pcff_qm.reference import CONVERSION

from .campaign import verify_completed
from .control_engine import molecule
from .hmds_engine import capture
from .hmds_paths import response


def read(path):
    return json.loads(Path(path).read_text())


def protocol_contrast(left, right, allowed):
    """Reject mixed geometry/method/settings disguised as a one-factor contrast."""
    keys = (
        "coordinates",
        "atom_order",
        "elements",
        "graph_identity",
        "charge",
        "multiplicity",
        "method",
        "basis",
        "grid",
        "thresholds",
        "scf_type",
    )
    changed = {k for k in keys if left[k] != right[k]}
    settings = {
        k
        for k in set(left["settings"]) | set(right["settings"])
        if left["settings"].get(k) != right["settings"].get(k)
    }
    expected = {
        "grid": {"dft_radial_points", "dft_spherical_points"},
        "basis": {"basis"},
        "thresholds": {"e_convergence", "d_convergence"},
    }
    require(
        changed <= set(allowed)
        and settings <= set().union(*(expected.get(k, set()) for k in allowed)),
        "confounded_precision_contrast",
        changed=sorted(changed),
        settings=sorted(settings),
    )
    return {"changed_fields": sorted(changed), "changed_settings": sorted(settings)}


def decimal_energy(directory):
    text = (Path(directory) / "psi4.out").read_text()
    values = re.findall(r"^\s+Total Energy =\s*([-+\d.]+)\s*$", text, re.MULTILINE)
    require(len(values) == 1, "unique_native_total_energy_required")
    value = Decimal(values[0])
    raw = read(Path(directory) / "result.json")["energy_hartree"]
    require(
        math.isfinite(raw) and abs(float(value) - raw) < 2e-10,
        "native_decimal_energy_mismatch",
    )
    return value, raw


def rounding_audit(energies, step):
    """Absolute bounds address serialization/arithmetic, not SCF/integration error."""
    require(len(energies) == 3 and step > 0, "curvature_triple_required")
    with localcontext() as context:
        context.prec = 50
        center, plus, minus = [e[0] for e in energies]
        ef = Decimal(str(CONVERSION["energy_factor"]))
        h2 = Decimal(str(step)) ** 2
        numerator = plus - 2 * center + minus
        curvature = numerator * ef / h2
        printed_bound = sum(
            Decimal(w) * Decimal(10) ** d.as_tuple().exponent / 2
            for w, (d, _) in zip((2, 1, 1), energies, strict=True)
        )
        binary = [v for _, v in energies]
        before = math.fsum([binary[1], binary[2], -2 * binary[0]])
        binary_bound = sum(
            w * math.ulp(v) for w, v in zip((2, 1, 1), binary, strict=True)
        )
        after_conversion = (
            binary[1] * float(ef) - 2 * binary[0] * float(ef) + binary[2] * float(ef)
        ) / step**2
        return {
            "native_decimal_energies_hartree": [str(d) for d, _ in energies],
            "second_difference_hartree": str(numerator),
            "decimal_energy_curvature": float(curvature),
            "binary_subtract_before_conversion_curvature": before * float(ef) / step**2,
            "legacy_convert_before_subtraction_curvature": after_conversion,
            "printed_decimal_rounding_bound": float(printed_bound * ef / h2),
            "binary_serialization_arithmetic_bound": binary_bound * float(ef) / step**2,
            "scope": "Rounding upper bounds only; neither integration nor SCF uncertainty.",
        }


def path_response(samples, geometries, energies, step):
    require(
        len(samples) == len(geometries) == len(energies) == 3, "path_triple_required"
    )
    for sample in samples:
        f = np.asarray(sample["forces"], float)
        require(f.shape == (27, 3) and np.isfinite(f).all(), "full_HMDS_force_required")
    r = response(*samples, *geometries, step)
    audit = rounding_audit(energies, step)
    r["legacy_normalized_energy_curvature"] = r["curvature_from_energies"]
    r["curvature_from_energies"] = audit["decimal_energy_curvature"]
    r["rounding"] = audit
    return r


def check_budget(plan, attempts, elapsed):
    require(
        len(attempts) < plan["limits"]["maximum_new_gradients"]
        and elapsed < plan["limits"]["campaign_wall_seconds"],
        "precision_campaign_budget_exhausted",
    )


def gradient(work, key):
    """One frozen job, serialized; outer ledger also counts initialization failures."""
    work = Path(work)
    plan = read(work / "scientific-plan.json")
    jobs = read(work / "jobs.json")["jobs"]
    require(key in jobs, "unplanned_precision_job")
    r = read(work / "requests" / (key + ".json"))
    require(
        digest(r) == jobs[key]["request_identity"], "frozen_precision_request_mutation"
    )
    out = work / "outputs" / key
    if out.exists():
        return verify_completed(out, digest(r))
    # Keep the lock throughout execution: no race between cap check and capture.
    with (work / "execution.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        attempts = read(work / "attempts.json")
        require(key in attempts["calls"], "counted_process_attempt_required")
        require(
            len(attempts["calls"]) <= plan["limits"]["maximum_new_gradients"],
            "precision_campaign_budget_exhausted",
        )
        import psi4

        psi4.set_num_threads(r["threads"])
        psi4.set_memory(f"{jobs[key]['memory_GiB']} GiB")
        psi4.core.IOManager.shared_object().set_default_path(os.environ["PSI_SCRATCH"])
        psi4.set_options(r["settings"])
        mol = molecule(psi4, r)
        write(
            work / "resources" / (key + ".json"),
            {
                "wrapper_sha256": pin(__file__),
                "capture_sha256": pin(capture.__code__.co_filename),
                "memory_bytes": psi4.get_memory(),
                "threads": psi4.get_num_threads(),
                "interpreter": os.sys.executable,
                "engine": psi4.__version__,
                "started_unix": time.time(),
            },
        )
        try:
            capture(
                work,
                key,
                r,
                lambda: psi4.gradient(r["method"], molecule=mol, return_wfn=True),
            )
        finally:
            usage = resource.getrusage(resource.RUSAGE_SELF)
            write(
                work / "resources" / (key + "-usage.json"),
                {
                    "maximum_RSS_KiB": usage.ru_maxrss,
                    "user_seconds": usage.ru_utime,
                    "system_seconds": usage.ru_stime,
                    "block_inputs": usage.ru_inblock,
                    "block_outputs": usage.ru_oublock,
                    "resource_scope": "Linux getrusage self, including Python/Psi4; not configured-memory usage",
                },
            )
            psi4.core.clean()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--key", required=True)
    args = parser.parse_args()
    gradient(args.work.resolve(), args.key)
