"""Reject J22 receipts that claim pinned-LAMMPS periodic runtime success.

This is a data-control check. It does not execute LAMMPS or OpenMM.
"""

import json
import sys
from pathlib import Path


def validate(path):
    payload = json.loads(Path(path).read_text())
    if payload.get("schema") not in {
        "island_phase_4j22_evidence_v2",
        "island_phase_4j22_correction_receipt_v1",
    }:
        raise ValueError("unsupported J22 receipt schema")
    runtime = payload.get("lammps_periodic_runtime_comparison", {})
    if runtime.get("kspace_package") is False:
        if runtime.get("ewald_pppm_single_point_completed") is True:
            raise ValueError("receipt claims LAMMPS Ewald/PPPM success without KSPACE")
        if runtime.get("status") != "blocked: KSPACE package unavailable":
            raise ValueError("receipt must record blocked KSPACE status")
    correction = payload.get("correction", {})
    if correction and correction.get("pinned_lammps_kspace_package") is False and (
            correction.get("lammps_ewald_pppm_completed") is True
            or correction.get("lammps_openmm_comparison_claimed") is True
    ):
        raise ValueError("correction claims unavailable LAMMPS periodic success")
    if payload.get("production_validated") is True:
        raise ValueError("J22 readiness must remain false")
    if payload.get("simulation_readiness") == "established":
        raise ValueError("J22 readiness must remain not_established")
    return True


if __name__ == "__main__":
    for argument in sys.argv[1:]:
        validate(argument)
        print(f"valid J22 receipt: {argument}")
