"""Subtract verified fixed contributions from real native QM; never fit a curve alone."""

import argparse
import json
from pathlib import Path

import numpy as np

from island.pcff_completion import pin, require, signed, write
from island.pcff_qm_pilot.hamiltonian import load_model, parts


def residuals(repo, evidence):
    repo, evidence = Path(repo), Path(evidence)
    m = load_model(repo)
    data = json.loads((evidence / "normalized-reference.json").read_text())
    samples = data["samples"]
    records = []
    geometries = {}
    for key, s in samples.items():
        g = s["coordinates"]
        identity = s["job_identity"]
        p = parts(m, g)
        qmf = np.array(s["forces"])
        remainder = qmf - p["fixed_forces"]
        geometries[key] = {
            "coordinates": g,
            "fixed_energy": p["fixed_energy"],
            "fixed_forces": p["fixed_forces"].tolist(),
            "components": p["components"],
            "unknown_energy_jacobian": p["energy_jacobian"].tolist(),
            "unknown_force_jacobian": p["force_jacobian"].tolist(),
        }
        records.append(
            {
                "native_job": key,
                "native_job_identity": identity,
                "QM_energy": s["energy"],
                "fixed_MM_energy": p["fixed_energy"],
                "QM_minus_fixed_energy": s["energy"] - p["fixed_energy"],
                "QM_force_RMS": float(np.sqrt(np.mean(qmf * qmf))),
                "fixed_MM_force_RMS": float(np.sqrt(np.mean(p["fixed_forces"] ** 2))),
                "unfitted_remainder_force_RMS": float(
                    np.sqrt(np.mean(remainder * remainder))
                ),
                "not_a_fitting_error": True,
                "unknowns_physically_zero": False,
            }
        )
    return signed(
        {
            "schema": "island.j60-real-residual-subtraction.v1",
            "hamiltonian": "island_research_class2_sih_disiloxane_j60_v1",
            "relationship": "QM = verified fixed MM + complete explicitly fitted unknown model + model discrepancy",
            "coefficients": None,
            "columns": p["columns"],
            "records": records,
            "geometries": geometries,
            "interpretation": "Raw energy zero is electronic and not an isolated MM cross-term energy. A single training-reference energy offset is required for fitting; unfitted force remainder is not model validation error.",
        }
    )


def family_check(repo, evidence):
    """Reuse all exact common-geometry def2/cc tight-d results; no new QM jobs."""
    from island.pcff_qm_pilot.hamiltonian_analysis import convergence

    evidence = Path(evidence)
    primary = evidence / "normalized-reference.json"
    corrected = evidence / "stages/tight-d-basis/normalized-reference.json"
    a, b = (json.loads(p.read_text()) for p in (primary, corrected))
    extra = {k: v for k, v in a["samples"].items() if k.startswith("def2/")}
    for key, s in extra.items():
        other = b["samples"]["cc5/" + key.split("/")[1]]
        require(s["coordinates"] == other["coordinates"], "family_geometry_mismatch")
        require(
            all(
                s["protocol"][k] == other["protocol"][k]
                for k in ("method", "grid", "thresholds")
            ),
            "family_protocol_conflation",
        )
    combined = {**b, "samples": {**b["samples"], **extra}}
    c = convergence(repo, evidence, combined)
    contrast = next(
        x for x in c["comparisons"] if x["effect"] == "basis_family_sensitivity"
    )
    return signed(
        {
            "schema": "island.j60-independent-family-contrast.v1",
            "input_sha256": {"primary": pin(primary), "tight_d": pin(corrected)},
            "protocol_a": "PBE0/def2-QZVP, including the reused exact J59 dense optimized-center derivative",
            "protocol_b": "PBE0/cc-pV(5+d)Z",
            "fixed_settings": "same geometries, 150x974 grid, SCF thresholds and def2-universal-jkfit; no optimization or method mixing",
            "contrast": contrast,
            "mixed_derivatives": {
                k: c["mixed_derivatives"][k] for k in ("def2", "cc5")
            },
            "new_QM_jobs": 0,
            "interpretation": "Ancillary orbital-family/cardinality/polarization sensitivity. Neither a CBS claim nor a newly selected acceptance threshold.",
        }
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j60"))
    p.add_argument("--family-check", action="store_true")
    a = p.parse_args()
    if a.family_check:
        write(
            a.evidence / "cross-family-reference-check.json",
            family_check(a.repo, a.evidence),
        )
    else:
        write(
            a.evidence / "real-residual-subtraction.json", residuals(a.repo, a.evidence)
        )
