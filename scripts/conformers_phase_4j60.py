import json
from itertools import permutations, product
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, signed, write
from island.pcff_qm.geometry import measure

E = Path("docs/evidence/phase_4j60")
d = json.loads((E / "joint-target-design.json").read_text())
perms = []
for a, b in product(permutations((3, 4, 5)), permutations((6, 7, 8))):
    perms.append([0, 1, 2, *a, *b])
    perms.append([2, 1, 0, *b, *a])


def aligned_rms(a, b):
    a = np.array(a)
    b = np.array(b)
    a -= a.mean(axis=0)
    b -= b.mean(axis=0)
    vals = []
    for perm in perms:
        q = b[perm]
        u, _, v = np.linalg.svd(q.T @ a)
        r = u @ np.diag([1, 1, np.linalg.det(u @ v)]) @ v
        vals.append(float(np.sqrt(np.mean(np.sum((q @ r - a) ** 2, axis=1)))))
    return min(vals)


rows = []
for h in d["held_out"]:
    distances = [
        (t["id"], aligned_rms(h["coordinates"], t["coordinates"]))
        for t in d["selected_training"]
    ]
    nearest = min(distances, key=lambda x: x[1])
    rows.append(
        {
            "id": h["id"],
            "closest_training_id": nearest[0],
            "minimum_proper_rotation_and_graph_automorphism_RMS_A": nearest[1],
            "torsions_radian": [
                measure(
                    h["coordinates"],
                    list(range(1, 10)),
                    {"kind": "torsion", "sites": s},
                )
                for s in ([4, 1, 2, 3], [1, 2, 3, 7])
            ],
            "validation_not_used_for_design": True,
        }
    )
print(rows)
write(
    E / "conformer-separation.json",
    signed(
        {
            "schema": "island.j60-conformer-separation.v1",
            "design_sha256": pin(E / "joint-target-design.json"),
            "automorphisms": len(perms),
            "permitted_symmetry": "permutations of equivalent hydrogens on each Si and whole-end exchange; proper rotations/translations only",
            "held_out": rows,
            "held_out_mutual_minimum_RMS_A": aligned_rms(
                d["held_out"][0]["coordinates"], d["held_out"][1]["coordinates"]
            ),
            "limitation": "Distinct conformers of this single Si-H graph; no methyl-PDMS transfer test.",
        }
    ),
)
