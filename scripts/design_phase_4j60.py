import json
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, signed, write
from island.pcff_qm.geometry import _adjust
from island.pcff_qm_pilot.campaign import selected
from island.pcff_qm_pilot.hamiltonian import load_model, parts

E = "docs/evidence/phase_4j60/"
m = load_model(".")
center = np.array(
    json.loads(Path("docs/evidence/phase_4j59/center-effects.json").read_text())[
        "coordinates"
    ]
)
system, a, _ = selected(
    ".", "../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/frc_files", center
)
coords = [
    ("bond", [1, 4], 0.05),
    ("bond", [1, 2], 0.05),
    ("bond", [2, 3], 0.05),
    ("bond", [3, 7], 0.05),
    ("angle", [4, 1, 2], 0.12),
    ("angle", [4, 1, 5], 0.12),
    ("angle", [1, 2, 3], 0.15),
    ("angle", [2, 3, 7], 0.12),
    ("torsion", [4, 1, 2, 3], 0.7),
    ("torsion", [1, 2, 3, 7], 0.7),
]
rng = np.random.default_rng(4060)
candidates = []
matrices = []
names = parts(m, center)["columns"]
for k in range(100):
    x = center.copy()
    offsets = []
    for kind, sites, extent in coords:
        offset = float(rng.uniform(-extent, extent))
        _adjust(system, x, list(range(1, 10)), {"kind": kind, "sites": sites}, offset)
        offsets.append(offset)
    p = parts(m, x)
    b = np.vstack(
        [
            p["energy_jacobian"][None, :] / 0.05,
            p["force_jacobian"].reshape(27, -1) / 0.5,
        ]
    )
    candidates.append(
        {
            "id": f"candidate-{k:03d}",
            "coordinates": x.tolist(),
            "geometry_identity": digest(x.tolist()),
            "offsets": offsets,
            "split": "train",
        }
    )
    matrices.append(b)
# Greedy block-D optimal screening only of complete Jacobian, no QM outcomes.
scale = np.sqrt(sum((b * b).sum(axis=0) for b in matrices))
scaled = [b / scale for b in matrices]
selected_ids = []
information = np.eye(len(names)) * 1e-10
for k in range(18):
    scores = [
        np.linalg.slogdet(information + b.T @ b)[1]
        if i not in selected_ids
        else -np.inf
        for i, b in enumerate(scaled)
    ]
    i = int(np.argmax(scores))
    selected_ids.append(i)
    information += scaled[i].T @ scaled[i]
validation = []
for k, (first, second) in enumerate(((1.05, -0.85), (-1.0, 0.9))):
    x = center.copy()
    for sites, offset in (([4, 1, 2, 3], first), ([1, 2, 3, 7], second)):
        _adjust(
            system, x, list(range(1, 10)), {"kind": "torsion", "sites": sites}, offset
        )
    # Distinct torsional conformer, with separately declared bond/angle perturbations.
    for kind, sites, extent in coords[:8]:
        _adjust(
            system,
            x,
            list(range(1, 10)),
            {"kind": kind, "sites": sites},
            extent * (0.3 if k == 0 else -0.4),
        )
    validation.append(
        {
            "id": f"untouched-conformer-{k}",
            "coordinates": x.tolist(),
            "geometry_identity": digest(x.tolist()),
            "split": "validation_only",
            "torsion_offsets": [first, second],
            "pair_distance_matrix_max_change_A": float(
                np.max(
                    abs(
                        np.linalg.norm(x[:, None] - x[None, :], axis=2)
                        - np.linalg.norm(center[:, None] - center[None, :], axis=2)
                    )
                )
            ),
        }
    )


def metrics(bs):
    b = np.vstack(bs)
    norm = np.linalg.norm(b, axis=0)
    u = b / norm
    s = np.linalg.svd(u, compute_uv=False)
    corr = u.T @ u
    end = names.index("end_bond-torsion_3:hsi,sio,osi,sio:block1:n3")
    mid = names.index("middle_bond-torsion_3:hsi,sio,osi,sio:block0:n3")
    return {
        "rank": int(sum(s > s[0] * 1e-10)),
        "columns": len(names),
        "condition": float(s[0] / s[-1]),
        "scaled_singular_values": s.tolist(),
        "end_middle_n3_cosine": float(corr[end, mid]),
    }


old = json.loads(
    Path("docs/evidence/phase_4j59/revised-QZ-center-scans.json").read_text()
)
oldbs = []
for scan in old.values():
    for point in scan["points"]:
        p = parts(m, point["coordinates"])
        oldbs.append(
            np.vstack(
                [
                    p["energy_jacobian"][None, :] / 0.05,
                    p["force_jacobian"].reshape(27, -1) / 0.5,
                ]
            )
        )
write(
    E + "joint-target-design.json",
    signed(
        {
            "schema": "island.j60-joint-design.v1",
            "graph_identity": m["graph_identity"],
            "hamiltonian": "island_research_class2_sih_disiloxane_j60_v1",
            "model_identity": m["identity"],
            "columns": names,
            "independent_coefficient_count": len(names),
            "symmetry": "AA outer exchange only; path reversal with end block exchange; same type environment tied across sites",
            "coordinate_inventory": [
                {"kind": k, "sites": s, "maximum_offset": e} for k, s, e in coords
            ],
            "selection_seed": 4060,
            "candidate_count": len(candidates),
            "selection": "greedy block log-determinant using full energy/force Jacobian; no QM response used",
            "selected_training": [candidates[i] for i in selected_ids],
            "held_out": validation,
            "old_single_conformer_grid": metrics(oldbs),
            "selected_complete_design": metrics([matrices[i] for i in selected_ids]),
            "fitting_criteria_source": "plan/reference-plan.json",
            "energy_offsets": "one training reference intercept, never per-pair or per-conformer tuning",
            "regularization": "none",
            "bounds": "each variable +/-1000 natural kcal/mol,angstrom,radian powers",
            "held_out_leakage_policy": "held-outs excluded from Jacobian selection, parameter fitting and weights/bounds/threshold selection",
            "execution": "conditional on reference gate; no QM values fabricated",
        }
    ),
)
print(metrics(oldbs))
print(metrics([matrices[i] for i in selected_ids]))
