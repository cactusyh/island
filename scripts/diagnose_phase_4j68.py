"""Localize the checked model-span obstruction and freeze two scalar hypotheses."""

import argparse
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hmds_constrained import physical_problem
from island.pcff_qm_pilot.hmds_local_baseline import (
    derivative_check,
    extend,
    fixed_families,
    references,
    span_sensitivity,
)
from island.pcff_qm_pilot.hmds_precision import read
from scripts.analyze_phase_4j68 import ROOT, reconstruct


def localize(repo):
    root = repo / ROOT
    _plan, model, rows, data, _lo, _hi = reconstruct(repo)
    witness = read(root / "solves/bound-free-Clarabel.json")
    require(
        witness["classification"] == "checked_lower_bound_above_one"
        and read(root / "solves/bound-free-SCS.json")["classification"]
        == "checked_lower_bound_above_one",
        "J68_hypothesis_activation",
    )
    values = np.array(witness["solution"]["values"])
    projected = read(root / "model-span-projection.json")
    residual = np.array(
        witness["checks"]["subsets"]["full26"]["signed_Cartesian_force_errors"]
    )
    p = physical_problem(data)
    records = []
    for i, row in enumerate(rows):
        record = fixed_families(model, row["coordinates"])
        adjustable = {}
        for j, name in enumerate(data["names"]):
            family = name.split(":")[0] if not name.startswith("delta_K2") else name
            if family not in adjustable:
                adjustable[family] = np.zeros((27, 3))
            adjustable[family] += data["force"][i, :, :, j] * values[j]
        all_fixed = record["families"]
        fsum = sum(np.array(v["forces"]) for v in all_fixed.values())
        final = fsum + sum(adjustable.values())
        atoms = []
        for atom in (5, 2, 6):
            unit = residual[i, atom - 1] / max(
                np.linalg.norm(residual[i, atom - 1]), 1e-30
            )
            contributions = {
                k: np.array(v["forces"])[atom - 1] for k, v in all_fixed.items()
            }
            contributions.update(
                {"adjustable:" + k: v[atom - 1] for k, v in adjustable.items()}
            )
            atoms.append(
                {
                    "atom": atom,
                    "QM_force": data["QM_force"][i, atom - 1].tolist(),
                    "MM_force": final[atom - 1].tolist(),
                    "residual": residual[i, atom - 1].tolist(),
                    "residual_vector_norm": float(
                        np.linalg.norm(residual[i, atom - 1])
                    ),
                    "contributions": [
                        {
                            "family": k,
                            "vector": v.tolist(),
                            "norm": float(np.linalg.norm(v)),
                            "projection_on_residual_unit": float(v @ unit),
                        }
                        for k, v in contributions.items()
                    ],
                    "sum_component_norms": float(
                        sum(np.linalg.norm(v) for v in contributions.values())
                    ),
                    "norm_of_sum": float(np.linalg.norm(final[atom - 1])),
                }
            )
        record.update(
            id=row["id"],
            adjustable_family_forces={k: v.tolist() for k, v in adjustable.items()},
            limiting_site_analysis=atoms,
        )
        records.append(record)
    sensitivity = {}
    for variant in ("bridge_angle", "SiO_stretch"):
        amended = extend(model, rows, data, variant)
        s = span_sensitivity(p, physical_problem(amended))
        derivatives = [
            {"id": r["id"], **derivative_check(model, r["coordinates"], variant)}
            for r in rows
        ]
        require(
            s["independent"]
            and max(d["force_finite_difference_error"] for d in derivatives) < 1e-6,
            "J68_independent_derivative_required",
        )
        sensitivity[variant] = {
            "sites": references(model, variant),
            "span": s,
            "derivative_checks": derivatives,
        }
    write(
        root / "fixed-family-localization.json",
        signed(
            {
                "rows": records,
                "scope": "All fixed families and full explicit adjustable vector. Magnitudes alone do not identify cause; cancellation, residual direction and projected scalar-potential sensitivities retained.",
                "witness": pin(root / "solves/bound-free-Clarabel.json"),
            }
        ),
    )
    write(root / "baseline-pre-fit-sensitivities.json", signed(sensitivity))
    frozen = {
        "schema": "island.j68-baseline-hypotheses.v1",
        "before_variant_fit": True,
        "activation": "Two-backend original bound-free lower bounds exceed1; fixed baseline and weighted projection localized before variant solves.",
        "localization_sha256": pin(root / "fixed-family-localization.json"),
        "sensitivity_sha256": pin(root / "baseline-pre-fit-sensitivities.json"),
        "variants": [
            {
                "id": v,
                "equation": "E=a sum_s(q_s-q0_s)+b sum_s(q_s-q0_s)^2; F=-sum_s[a+2b(q_s-q0_s)] grad(q_s)",
                "variable": "angle radians" if v == "bridge_angle" else "bond angstrom",
                "units": ["kcal/mol/radian", "kcal/mol/radian^2"]
                if v == "bridge_angle"
                else ["kcal/mol/angstrom", "kcal/mol/angstrom^2"],
                "sites": references(model, v),
                "coefficient_sharing": "One linear and one quadratic coefficient shared across listed source-type environment.",
                "qualification": "Additive polynomial about unchanged source equilibrium; not a propagated equilibrium change; diagnostic only.",
                "original58_bounds": "removed for model-span comparison only",
                "new2_bounds": "none; diagnostic only",
                "selection_reason": "Nonzero independently projected energy/force directions at bridge O5/Si2/Si6 where residuals limit original bound-free model.",
            }
            for v in sensitivity
        ],
        "maximum_variants": 2,
        "combination_authorized": False,
        "solvers": ["clarabel", "scs"],
        "limits": "Identical four limits on protected22 and full26; no fit/validation tuning.",
        "accepted_coefficients": 0,
        "validation_activation": False,
        "new_QM": 0,
        "MD_steps": 0,
    }
    require(
        not (root / "plan/baseline-hypotheses.json").exists(),
        "J68_additive_hypothesis_freeze",
    )
    write(root / "plan/baseline-hypotheses.json", signed(frozen))
    print(
        {
            v: {
                "independent": d["span"]["independent"],
                "outside_old_span_norms": d["span"]["outside_old_span_norms"],
                "capture": d["span"]["target_projected_capture_norm"],
            }
            for v, d in sensitivity.items()
        }
    )
    print("frozen", pin(root / "plan/baseline-hypotheses.json"))
    print("projection", projected["weighted_residual_norm"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    localize(Path.cwd())
