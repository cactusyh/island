"""Read-only acceptance of retained PE bundles; no preparation or propagation.

Declared before execution: atol=rtol=1e-12 for pair-distance Rg^2,
long-double direct tensor and independently resolved Re. Sources must remain
quiescent during this independent two-reader comparison.
"""

import argparse
import json
from pathlib import Path

import numpy as np

from island.analysis import analyze_workflow, export_analysis
from island.workflows import read_workflow_frames, storage, workflow_status
from island.workflows.chain import _load_bundle


def validate(directory, output):
    root = Path(directory)
    report = analyze_workflow(root)
    p = report.payload
    manifest = workflow_status(root)
    system, _, _ = _load_bundle(root, manifest)
    frames = read_workflow_frames(root)
    errors = {"Rg_squared_angstrom2": 0.0, "tensor_angstrom2": 0.0, "Re_angstrom": 0.0}
    for frame, row in zip(frames, p["frames"], strict=True):
        assert frame.coordinate_fingerprint == row["coordinate_fingerprint"]
        g = row["geometry"]
        ids = g["selected_ids"]
        xyz = np.array([frame.coordinates[s] for s in ids], dtype=np.longdouble)
        mass = np.array(
            [system.topology.sites[s].mass for s in ids], dtype=np.longdouble
        )
        pair_rg2 = (
            sum(
                mass[i] * mass[j] * sum((xyz[i] - xyz[j]) ** 2)
                for i in range(len(ids))
                for j in range(i)
            )
            / sum(mass) ** 2
        )
        center = sum(mass[i] * xyz[i] for i in range(len(ids))) / sum(mass)
        tensor = sum(
            mass[i] * np.outer(xyz[i] - center, xyz[i] - center)
            for i in range(len(ids))
        ) / sum(mass)
        np.testing.assert_allclose(
            g["radius_of_gyration"] ** 2, pair_rg2, atol=1e-12, rtol=1e-12
        )
        np.testing.assert_allclose(g["gyration_tensor"], tensor, atol=1e-12, rtol=1e-12)
        errors["Rg_squared_angstrom2"] = max(
            errors["Rg_squared_angstrom2"],
            float(abs(g["radius_of_gyration"] ** 2 - pair_rg2)),
        )
        errors["tensor_angstrom2"] = max(
            errors["tensor_angstrom2"],
            float(np.max(abs(np.array(g["gyration_tensor"]) - tensor))),
        )
        polymer = system.metadata["polymer"]
        pair = (polymer["head_site_id"], polymer["tail_site_id"])
        assert tuple(g["endpoint_ids"]) == pair
        a, b = [np.array(frame.coordinates[s], dtype=np.longdouble) for s in pair]
        re = np.sqrt(sum((a - b) ** 2))
        np.testing.assert_allclose(g["end_to_end_distance"], re, atol=1e-12, rtol=1e-12)
        errors["Re_angstrom"] = max(
            errors["Re_angstrom"], float(abs(g["end_to_end_distance"] - re))
        )
    assert (
        storage.checksum((root / "manifest.json").read_bytes())
        == p["source"]["manifest_file_sha256"]
    )
    export_analysis(report, output)
    return {
        "source": str(root),
        "manifest_checksum": p["source"]["manifest_checksum"],
        "evidence": p["source"]["evidence"],
        "samples": len(frames),
        "steps": [r["step"] for r in p["frames"]],
        "status": p["source"]["status"],
        "atol": 1e-12,
        "rtol": 1e-12,
        "maximum_absolute_errors": errors,
        "new_dynamics_or_parameterization": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    print(json.dumps(validate(args.workflow, args.output), indent=2, allow_nan=False))
