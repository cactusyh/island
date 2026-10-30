"""Standalone geometry audit; does not import the production geometry functions."""

import json
import math
import sys
import tarfile
from collections import Counter
from hashlib import sha256
from pathlib import Path

import numpy as np

root = Path(sys.argv[1])
archive = root / "geometry-input-bundle.tar.xz"
counts = Counter()
maxima = {}
directories = 0
with tarfile.open(archive, "r:xz") as tar:
    files = {m.name: tar.extractfile(m).read() for m in tar.getmembers() if m.isfile()}
    for name in sorted(files):
        if not name.endswith("/metadata.json"):
            continue
        directories += 1
        directory = name.rsplit("/", 1)[0]
        metadata = json.loads(files[name])
        samples = json.loads(files[directory + "/geometries.json"])
        order = metadata["atom_order"]
        coords = metadata["scan"]["coordinates"]
        for sample in samples:
            counts["geometries"] += 1
            xyz = np.array(sample["coordinates"])
            for j, c in enumerate(coords):
                v = xyz[[order.index(i) for i in c["sites"]]]
                kind = c["kind"]
                if kind == "bond":
                    d = v[1] - v[0]
                    value = math.sqrt(sum(float(z) ** 2 for z in d))
                elif kind == "angle":
                    a, b = v[0] - v[1], v[2] - v[1]
                    value = math.atan2(
                        float(np.linalg.norm(np.cross(a, b))), float(a @ b)
                    )
                else:
                    axis = v[2] - v[1]
                    axis = axis / math.sqrt(float(axis @ axis))
                    a, b = v[0] - v[1], v[3] - v[2]
                    a = a - float(a @ axis) * axis
                    b = b - float(b @ axis) * axis
                    value = math.atan2(
                        float(np.linalg.det(np.array([axis, a, b]))), float(a @ b)
                    )
                desired = sample["baseline_coordinates"][j] + sample["offsets"][j]
                residual = value - desired
                if kind == "torsion":
                    residual = math.remainder(residual, 2 * math.pi)
                counts[kind] += 1
                maxima[kind] = max(maxima.get(kind, 0), abs(residual))
                assert abs(residual) < 1e-10, (name, sample["id"], kind, residual)
record = {
    "schema": "island.j57-independent-geometry.v1",
    "status": "passed",
    "outside_production_geometry_functions": True,
    "input_archive_sha256": sha256(archive.read_bytes()).hexdigest(),
    "directories": directories,
    "counts": dict(counts),
    "max_absolute_residuals": maxima,
    "bond_units": "angstrom",
    "angular_units": "radian",
    "angle_formula": "atan2(norm(cross(a,b)),dot(a,b))",
    "torsion_formula": "atan2(det(axis,projected_a,projected_b),dot(projected_a,projected_b))",
    "QM_executed": False,
    "fitted_parameters_created": False,
    "scientific_validation": "not_established",
}
canonical = json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False)
record["identity"] = sha256(canonical.encode()).hexdigest()
Path(sys.argv[2]).write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
print(json.dumps(record, sort_keys=True))
