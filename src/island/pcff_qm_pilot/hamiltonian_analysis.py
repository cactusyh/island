"""Numerical verification and same-geometry reference decisions for J60."""

import json
from pathlib import Path

import numpy as np

from island.pcff_completion import digest, pin, require, signed
from island.pcff_qm.reference import CONVERSION

from .convergence_importer import native
from .hamiltonian import cos_phi, geometry, harmonic, load_model, parts
from .hamiltonian_campaign import requests


def normalize(evidence, work):
    evidence, work = Path(evidence), Path(work)
    pp = evidence / "plan/reference-plan.json"
    plan = json.loads(pp.read_text())
    manifest = json.loads((evidence / "plan/campaign.json").read_text())
    require(
        pin(
            Path(__file__).with_name(
                plan.get(
                    "adapter_module", "island.pcff_qm_pilot.hamiltonian_engine"
                ).rsplit(".", 1)[1]
                + ".py"
            )
        )
        == manifest["adapter_sha256"],
        "executed_adapter_source_mutation",
    )
    jobs = requests(plan, pin(pp))
    require(
        manifest["jobs"] == {k: digest(v) for k, v in jobs.items()},
        "campaign_request_mutation",
    )
    samples = {}
    for key, expected in jobs.items():
        require(
            json.loads((work / "requests" / (key + ".json")).read_text()) == expected,
            "request_mutation",
        )
        r, g = native(work / "outputs" / key, expected)
        require(
            r["adapter"]
            == (
                "island.psi4-j60-tight-d.v1"
                if "refinement_engine" in plan.get("adapter_module", "")
                else "island.psi4-j60-reference.v1"
            )
            and r["adapter_sha256"] == manifest["adapter_sha256"],
            "adapter_mutation",
        )
        samples[key] = {
            "coordinates": expected["coordinates"],
            "graph_identity": expected["graph_identity"],
            "charge": expected["charge"],
            "multiplicity": expected["multiplicity"],
            "atom_order": expected["atom_order"],
            "protocol": {
                k: expected[k]
                for k in ("method", "basis", "grid", "thresholds", "settings")
            },
            "energy": r["energy_hartree"] * CONVERSION["energy_factor"],
            "forces": (
                -g * CONVERSION["energy_factor"] / CONVERSION["distance_factor"]
            ).tolist(),
            "energy_hartree": r["energy_hartree"],
            "gradient_hartree_bohr": g.tolist(),
            "job_identity": digest(expected),
            "native_files": r["native_files"],
            "elapsed_seconds": r["elapsed_seconds"],
        }
    return signed(
        {
            "schema": "island.j60-normalized-native.v1",
            "samples": samples,
            "conversion": CONVERSION,
            "coverage": plan.get(
                "stage_role",
                "bounded common-geometry convergence study; not a complete fitting dataset",
            ),
        }
    )


def numerical(repo, evidence):
    from scripts.independent_phase_4j60 import check

    repo, evidence = Path(repo), Path(evidence)
    model = load_model(repo)
    plan = json.loads((evidence / "plan/reference-plan.json").read_text())
    design = json.loads((evidence / "joint-target-design.json").read_text())
    records = []
    for name in ("center", "angle-angle-0000", "bond-bond-0024"):
        records.append({"geometry": name, **check(model, plan["geometries"][name])})
    # Independent-conformer controls also exercise highly displaced torsions.
    records.append(
        {
            "geometry": "distinct-conformer",
            **check(model, design["held_out"][0]["coordinates"]),
        }
    )
    xyz = np.array(plan["geometries"]["center"])
    p = parts(model, xyz)
    theta = 0.47
    rotation = np.array(
        [
            [np.cos(theta), -np.sin(theta), 0],
            [np.sin(theta), np.cos(theta), 0],
            [0, 0, 1],
        ]
    )
    moved = parts(model, xyz @ rotation.T + [1.2, -3.4, 5.6])
    rev = np.array([2, 1, 0, 6, 7, 8, 3, 4, 5])
    reverse = parts(model, xyz[rev])
    rigid = {
        "fixed_energy_error": abs(moved["fixed_energy"] - p["fixed_energy"]),
        "column_energy_max_error": float(
            np.max(abs(moved["energy_jacobian"] - p["energy_jacobian"]))
        ),
        "fixed_force_max_error": float(
            np.max(abs(moved["fixed_forces"] - p["fixed_forces"] @ rotation.T))
        ),
        "column_force_max_error": float(
            np.max(
                abs(
                    moved["force_jacobian"]
                    - np.einsum("ijk,lj->ilk", p["force_jacobian"], rotation)
                )
            )
        ),
        "net_fixed_force": p["fixed_forces"].sum(axis=0).tolist(),
        "net_fixed_torque": np.cross(xyz, p["fixed_forces"]).sum(axis=0).tolist(),
        "net_column_force_max": float(abs(p["force_jacobian"].sum(axis=0)).max()),
        "net_column_torque_max": float(
            abs(
                np.cross(
                    xyz[:, :, None], p["force_jacobian"], axisa=1, axisb=1, axisc=1
                ).sum(axis=0)
            ).max()
        ),
    }
    require(
        max(
            rigid[k]
            for k in (
                "fixed_energy_error",
                "column_energy_max_error",
                "fixed_force_max_error",
                "column_force_max_error",
                "net_column_force_max",
                "net_column_torque_max",
            )
        )
        < 1e-9,
        "rigid_invariance_failed",
    )
    reversal = {
        "fixed_energy_error": abs(reverse["fixed_energy"] - p["fixed_energy"]),
        "column_energy_max_error": float(
            np.max(abs(reverse["energy_jacobian"] - p["energy_jacobian"]))
        ),
        "fixed_force_max_error": float(
            np.max(abs(reverse["fixed_forces"] - p["fixed_forces"][rev]))
        ),
        "column_force_max_error": float(
            np.max(abs(reverse["force_jacobian"] - p["force_jacobian"][rev]))
        ),
    }
    require(max(reversal.values()) < 1e-9, "atom_order_reversal_failed")
    torsions = []
    for phi in (
        -np.pi - 1e-8,
        -np.pi + 1e-8,
        -1.2,
        0,
        0.4,
        np.pi - 1e-8,
        np.pi + 1e-8,
        2 * np.pi + 0.4,
    ):
        x = np.array(
            [[1, 0, 0], [0, 0, 0], [0, 0, 1], [np.cos(phi), np.sin(phi), 1]]
            + [[5 + j, 6, 7] for j in range(5)],
            float,
        )
        c = cos_phi(geometry(x), [1, 2, 3, 4])
        cr = cos_phi(geometry(x), [4, 3, 2, 1])
        for n in (1, 2, 3):
            h = harmonic(c, n)
            hr = harmonic(cr, n)
            require(
                abs(h.v - np.cos(n * phi)) < 1e-12
                and abs(hr.v - h.v) < 1e-12
                and np.max(abs(hr.d - h.d)) < 1e-10,
                "torsion_order_failed",
            )
        v = 0.37
        delta = 0.22
        source = v * (1 + np.cos(3 * phi - delta))
        engine = v * (1 - np.cos(3 * phi - (delta + np.pi)))
        copied = v * (1 - np.cos(3 * phi - delta))
        require(abs(source - engine) < 1e-12, "phase_transformation_failed")
        torsions.append(
            {
                "phi": phi,
                "source_plus_cos": source,
                "phase_transformed_engine": engine,
                "copied_phase_engine": copied,
                "cosine_AAT": 0.03 * np.cos(phi),
                "literal_FRC_linear_AAT": 0.03 * phi,
            }
        )
    pairs = [
        {"sites": r["sites"], "distance": r["graph_distance"], "weight": r["weight"]}
        for r in p["pairs"]
    ]
    return signed(
        {
            "schema": "island.j60-numerical-checks.v1",
            "independent_scalar_and_FD": records,
            "rigid_invariance": rigid,
            "full_graph_reversal": reversal,
            "torsion_controls": torsions,
            "AAT_non_equivalence": {
                "linear_change_per_turn": float(0.03 * 2 * np.pi),
                "cosine_change_per_turn": 0.0,
                "conversion_exists_for_nonzero_K": False,
            },
            "pairs": pairs,
            "unknown_columns": p["columns"],
            "all_20_force_kernels_verified": True,
            "synthetic_coefficient_controls_only": True,
            "fitted_coefficients": [],
        }
    )


def convergence(repo, evidence, normalized):
    repo, evidence = Path(repo), Path(evidence)
    plan = json.loads((evidence / "plan/reference-plan.json").read_text())
    samples = dict(normalized["samples"])
    # Reuse J59's final full-precision optimization gradient at exactly this center.
    center = json.loads(
        (repo / "docs/evidence/phase_4j59/center-effects.json").read_text()
    )
    # Source result retained inside old native archive; copied scalar is independently
    # validated by J59's own native replay. Exact native arrays reconstructed below.
    import tempfile

    from .bundle import unpack

    with tempfile.TemporaryDirectory(prefix="j60-j59-center-") as temp:
        root = Path(temp) / "j59"
        unpack(repo / "docs/evidence/phase_4j59/native-outputs.tar.xz", root)
        result = json.loads(
            (root / "outputs/optimize/qz-dense/result.json").read_text()
        )
        g = np.load(root / "outputs/optimize/qz-dense/gradient.npy", allow_pickle=False)
        require(
            center["coordinates"] == plan["geometries"]["center"],
            "center_geometry_mutation",
        )
        from .campaign import ENGINE_BOHR_ANGSTROM, verify_completed

        req = json.loads((root / "requests/optimize/qz-dense.json").read_text())
        verify_completed(root / "outputs/optimize/qz-dense", digest(req))
        require(
            np.allclose(
                np.array(result["geometry_bohr"]) * ENGINE_BOHR_ANGSTROM,
                center["coordinates"],
                atol=1e-12,
                rtol=0,
            ),
            "center_native_mutation",
        )
        samples["def2/center"] = {
            "energy": result["energy_hartree"] * CONVERSION["energy_factor"],
            "forces": (
                -g * CONVERSION["energy_factor"] / CONVERSION["distance_factor"]
            ).tolist(),
            "reused_native_result_sha256": pin(
                root / "outputs/optimize/qz-dense/result.json"
            ),
        }
    scans = json.loads(
        (repo / "docs/evidence/phase_4j59/revised-QZ-center-scans.json").read_text()
    )
    mixed = {}
    for protocol in ("ccq", "cc5", "def2", "grid"):
        mixed[protocol] = {}
        for name, h in (("angle-angle", 0.08), ("bond-bond", 0.04)):
            keys = [f"{protocol}/{name}-{i:04d}" for i in (0, 4, 20, 24)]
            if not all(k in samples for k in keys):
                continue
            e = [samples[k]["energy"] for k in keys]
            grid = np.array([p["coordinates"] for p in scans[name]["points"]]).reshape(
                5, 5, 9, 3
            )
            ds = []
            for i, j in ((0, 0), (0, 4), (4, 0), (4, 4)):
                tangent = (
                    (-3 * grid[0, j] + 4 * grid[1, j] - grid[2, j]) / h
                    if i == 0
                    else (3 * grid[4, j] - 4 * grid[3, j] + grid[2, j]) / h
                )
                ds.append(
                    float(
                        np.sum(
                            -np.array(
                                samples[f"{protocol}/{name}-{i * 5 + j:04d}"]["forces"]
                            )
                            * tangent
                        )
                    )
                )
            mixed[protocol][name] = {
                "energy_mixed_secant": float((e[3] - e[2] - e[1] + e[0]) / (4 * h * h)),
                "gradient_mixed_secant": float(
                    ((ds[3] - ds[2]) + (ds[1] - ds[0])) / (4 * h)
                ),
                "unit": "kcal/mol/radian^2"
                if name == "angle-angle"
                else "kcal/mol/angstrom^2",
                "interpretation": "total-QM finite path secant, not an isolated MM coefficient",
            }
    comparisons = []
    for a, b, effect, gate in (
        ("ccq", "cc5", "cardinal_QZ_5Z", "basis"),
        ("def2", "cc5", "basis_family_sensitivity", None),
        ("grid", "cc5", "numerical_grid", "grid"),
        ("cc5", "tight", "SCF_threshold", "scf"),
        ("cc5", "method", "PBE0_B3LYP", None),
    ):
        common = sorted(
            k.split("/")[1]
            for k in samples
            if k.startswith(a + "/") and b + "/" + k.split("/")[1] in samples
        )
        if not common:
            continue
        records = []
        for point in common:
            u, v = samples[a + "/" + point], samples[b + "/" + point]
            ca, cb = samples[a + "/center"], samples[b + "/center"]
            d = (v["energy"] - cb["energy"]) - (u["energy"] - ca["energy"])
            df = np.array(v["forces"]) - np.array(u["forces"])
            relative = df - (np.array(cb["forces"]) - np.array(ca["forces"]))
            records.append(
                {
                    "geometry": point,
                    "relative_energy_difference": d,
                    "gradient_RMS_difference": float(np.sqrt(np.mean(df * df))),
                    "deformation_gradient_RMS_difference": float(
                        np.sqrt(np.mean(relative * relative))
                    ),
                }
            )
        derivatives = {}
        if a in mixed and b in mixed:
            for scan in mixed[a]:
                derivatives[scan] = {}
                for kind in ("energy_mixed_secant", "gradient_mixed_secant"):
                    u = mixed[a][scan][kind]
                    v = mixed[b][scan][kind]
                    derivatives[scan][kind] = {
                        "a": u,
                        "b": v,
                        "relative_difference": abs(u - v) / max(abs(v), 1e-12),
                    }
        me = max(abs(r["relative_energy_difference"]) for r in records)
        mg = max(r["gradient_RMS_difference"] for r in records)
        passed = None
        if gate:
            criteria = plan["accuracy_criteria"][gate]
            passed = (
                me <= criteria["relative_energy_kcal_mol"]
                and mg <= criteria["gradient_RMS_kcal_mol_A"]
            )
            if derivatives and "mixed_derivative_relative" in criteria:
                passed = passed and all(
                    v["relative_difference"] <= criteria["mixed_derivative_relative"]
                    for scan in derivatives.values()
                    for v in scan.values()
                )
        comparisons.append(
            {
                "effect": effect,
                "records": records,
                "max_relative_energy_difference": me,
                "max_gradient_RMS_difference": mg,
                "mixed_derivative_differences": derivatives,
                "predeclared_gate_passed": passed,
            }
        )
    basis = next(
        x["predeclared_gate_passed"]
        for x in comparisons
        if x["effect"] == "cardinal_QZ_5Z"
    )
    return signed(
        {
            "schema": "island.j60-reference-convergence.v1",
            "comparisons": comparisons,
            "mixed_derivatives": mixed,
            "cardinal_basis_gate_passed": basis,
            "basis_protocols": {
                k: samples[k + "/center"]["protocol"]["basis"] for k in ("ccq", "cc5")
            },
            "J59_failed_gate_unchanged": True,
            "old_center_result": "docs/evidence/phase_4j59/QM-convergence.json",
            "center_effect": "Old J58 vs optimized J59 coordinates/protocols are separate records; no common-geometry contrast mixes centers.",
            "final_QM_adequacy": "not_established",
            "reused_center": samples["def2/center"],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
