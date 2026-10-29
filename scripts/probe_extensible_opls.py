"""Real source comparison for the additive finite external OPLS integration."""

import argparse
import json
import shutil
from copy import deepcopy
from dataclasses import replace
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.dynamics import (
    DynamicsSegmentOptions,
    LangevinOptions,
    create_dynamics_checkpoint,
    initialize_velocities,
    load_dynamics_checkpoint,
    resume_dynamics,
    run_dynamics_segment,
    save_dynamics_checkpoint,
)
from island.forcefields import preparation as native
from island.forcefields.extensible import (
    ExternalOPLSError,
    ExternalOPLSOptions,
    ForceFieldRequest,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.external_pdms import LIMITS
from island.forcefields.opls_external import write_source_manifest
from island.forcefields.opls_external.storage import save_preparation
from island.minimization import minimize_geometry


def file_hash(path):
    return sha256(path.read_bytes()).hexdigest()


def make_real_control_source(xml, license_file, output):
    """Explicit finite control only; no author PDMS model conventions are changed."""
    output.mkdir()
    shutil.copy(xml, output / "source.xml")
    shutil.copy(license_file, output / "LICENSE")
    proof = {
        "file": "LICENSE",
        "sha256": file_hash(output / "LICENSE"),
        "line": 5,
        "raw_text": (output / "LICENSE").read_text().splitlines()[4],
    }
    pin = write_source_manifest(
        output / "source.xml",
        output / "source.json",
        profile_id="island_external_opls_foyer_control_v1",
        provenance={
            "citation": "Foyer-distributed OPLS XML; finite CH software control only",
            "doi": None,
            "charge_method": "unchanged NonbondedForce Atom charge",
            "license_status": "authorized_local_use",
            "license_evidence": proof,
            "origin": "Foyer dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd",
            "conversion": "Unmodified XML; finite NoCutoff software control, no new scientific model",
        },
        model={
            "nonbonded_method": "NoCutoff",
            "periodic": False,
            "improper_policy": "source_defines_no_impropers_for_declared_domain",
            "chemical_state": "neutral_closed_shell",
            "elements": ["C", "H"],
        },
    )
    return output / "source.json", pin


def independent_energy(system, assignment):
    """NumPy equations independent of OpenMM construction, coordinates in A."""
    from itertools import combinations
    from math import acos, cos, pi, sqrt

    xyz = {i: np.asarray(system.coordinates.get(i)) for i in system.topology.sites}
    result = dict.fromkeys(("bond", "angle", "proper_torsion", "lj", "coulomb"), 0.0)
    for r in assignment["bonds"]:
        a, b = (xyz[i] for i in r["sites"])
        p = r["converted"]
        result["bond"] += (
            0.5
            * p["k_kj_mol_angstrom2"]
            * (np.linalg.norm(a - b) - p["length_angstrom"]) ** 2
        )
    for r in assignment["angles"]:
        a, b, c = (xyz[i] for i in r["sites"])
        u, v = a - b, c - b
        p = r["converted"]
        theta = acos(
            float(np.clip(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v), -1, 1))
        )
        result["angle"] += (
            0.5 * p["k_kj_mol_radian2"] * (theta - p["angle_degrees"] * pi / 180) ** 2
        )
    for r in assignment["rb"]:
        a, b, c, d = (xyz[i] for i in r["sites"])
        n1, n2 = np.cross(b - a, c - b), np.cross(c - b, d - c)
        costheta = float(np.dot(n1, n2) / np.linalg.norm(n1) / np.linalg.norm(n2))
        values = r["converted"]["coefficients_kj_mol"]
        if assignment["torsion_form"] == "RBTorsionForce":
            result["proper_torsion"] += sum(
                k * (-costheta) ** i for i, k in enumerate(values)
            )
        else:
            theta = acos(float(np.clip(costheta, -1, 1)))
            result["proper_torsion"] += 0.5 * sum(
                k * (1 + (-1) ** i * cos((i + 1) * theta)) for i, k in enumerate(values)
            )
    pairs = {tuple(r["sites"]): r for r in assignment["pairs"]}
    for i, j in combinations(sorted(xyz), 2):
        a, b = assignment["sites"][i], assignment["sites"][j]
        pair = pairs.get((i, j), {"lj_scale": 1.0, "coulomb_scale": 1.0})
        distance = np.linalg.norm(xyz[i] - xyz[j])
        sigma = sqrt(a["sigma_angstrom"] * b["sigma_angstrom"])
        result["lj"] += (
            pair["lj_scale"]
            * 4
            * sqrt(a["epsilon_kj_mol"] * b["epsilon_kj_mol"])
            * ((sigma / distance) ** 12 - (sigma / distance) ** 6)
        )
        result["coulomb"] += (
            pair["coulomb_scale"]
            * 1389.35457644382
            * a["charge_e"]
            * b["charge_e"]
            / distance
        )
    return result


def independent_rows(xml, assignment, native_assignment):
    """DOM reconstruction independent of production ElementTree/Expat readers."""
    from xml.dom import minidom

    dom = minidom.parse(str(xml))
    count = 0
    for r in assignment["source_row_evidence"]:
        section, key = r["source_record"].split("/")
        tag, index = key[:-1].split("[")
        node = dom.getElementsByTagName(section)[0].getElementsByTagName(tag)[
            int(index)
        ]
        attrs = dict(node.attributes.items())
        if attrs != r["source_attributes"]:
            raise ValueError("Independent DOM row mismatch")
        raw = minidom.parseString(r["raw_text"]).documentElement
        if dict(raw.attributes.items()) != attrs:
            raise ValueError("Raw row evidence mismatch")
        count += 1
    for family in ("bonds", "angles", "rb"):
        rows, refs = assignment[family], native_assignment[family]
        if len(rows) != len(refs):
            raise ValueError("Native interaction count mismatch")
        for a, b in zip(rows, refs, strict=True):
            if (
                tuple(a["sites"]) != tuple(b["sites"])
                or tuple(a["types"]) != tuple(b["types"])
                or a["source_record"] != b["source_record"]
                or a["source_attributes"] != b["source_attributes"]
            ):
                raise ValueError("Native row, type, order or mapping mismatch")
    for i, row in assignment["sites"].items():
        ref = native_assignment["sites"][i]
        if (
            row["type"] != ref["type"]
            or row["source_attributes"] != ref["source_attributes"]
            or row["charge_e"] != ref["charge_e"]
            or row["mass_dalton"] != ref["mass_dalton"]
        ):
            raise ValueError("Native atom mapping mismatch")
    return {
        "raw_rows_reconstructed": count,
        "native_atom_type_order_charge_mass_mapping_equal": True,
        "native_interaction_counts_and_selected_source_rows_equal": True,
    }


def probe(xml, license_file, work):
    path, pin = make_real_control_source(xml, license_file, work / "control_source")
    options = ExternalOPLSOptions(path, pin)
    controls = []
    for name, smiles in (("ethane", "CC"), ("butane", "CCCC")):
        s = from_smiles(smiles, random_seed=2026)
        p = prepare_forcefield(s, ForceFieldRequest("external_opls", options))
        e = create_evaluator(s, p)
        r = e.evaluate()
        n = native.prepare_forcefield(
            s, native.ForceFieldRequest("oplsaa", native.OPLSOptions(xml))
        )
        ref = native.create_evaluator(s, n).evaluate()
        row_checks = independent_rows(
            xml, p.payload["assignment"], n._native.payload["resolved"]
        )
        independent = independent_energy(s, p.payload["assignment"])
        component_errors = {
            k: abs(v - r.energy_components[k]) for k, v in independent.items()
        }
        if max(component_errors.values()) > 1e-8:
            raise ValueError("Independent energy disagreement")
        force_error = max(
            abs(r.forces[i][k] - ref.forces[i][k]) for i in r.forces for k in range(3)
        )
        # Every Cartesian component, two finite-difference spacings.
        fd = {}
        for h in (1e-4, 5e-5):
            errors = []
            for i in s.topology.sites:
                for axis in range(3):
                    plus = deepcopy(s)
                    minus = deepcopy(s)
                    a = s.coordinates.get(i).copy()
                    b = a.copy()
                    a[axis] += h
                    b[axis] -= h
                    plus.coordinates.set(i, a)
                    minus.coordinates.set(i, b)
                    value = -(
                        sum(independent_energy(plus, p.payload["assignment"]).values())
                        - sum(
                            independent_energy(minus, p.payload["assignment"]).values()
                        )
                    ) / (2 * h)
                    errors.append(abs(value - r.forces[i][axis]))
            fd[str(h)] = max(errors)
        if force_error > 1e-8 or max(fd.values()) > 1e-4:
            raise ValueError("Force disagreement")
        minimum = minimize_geometry(s, e)
        if not minimum.converged or not minimum.final_evaluation_verified:
            raise ValueError("Minimization not verified")
        moved = minimum.to_system(s)
        e.validate_system(moved)
        velocities = initialize_velocities(
            moved, temperature_kelvin=10, seed=123
        ).velocities
        opts = LangevinOptions(
            timestep_fs=0.01,
            temperature_kelvin=10,
            friction_per_ps=1,
            thermostat_seed=456,
            steps=4,
            max_evaluations=10,
            max_frames=10,
        )
        whole = run_dynamics_segment(moved, e, velocities, opts)
        first = run_dynamics_segment(moved, e, velocities, replace(opts, steps=2))
        if not whole.completed or not first.completed:
            raise ValueError("NVT control failed")
        checkpoint = create_dynamics_checkpoint(first)
        cp_path = work / f"{name}_checkpoint.json"
        save_dynamics_checkpoint(checkpoint, cp_path)
        restored = load_dynamics_checkpoint(cp_path)
        end = resume_dynamics(
            restored, first.to_system(moved), e, DynamicsSegmentOptions(2, 10, 10)
        )
        same_state = end.final_state == whole.final_state
        same_rng = end.payload["rng"] == whole.payload["rng"]
        if (
            not end.completed
            or not same_state
            or not same_rng
            or restored != checkpoint
        ):
            raise ValueError("Continuation mismatch")
        save_preparation(p, work / f"{name}_preparation.json")
        controls.append(
            {
                "fixture": name,
                "smiles": smiles,
                "atom_count": s.number_of_sites,
                "preparation_identity": p.identity,
                "status": p.payload["status"],
                "energy_kJ_mol": r.potential_energy,
                "energy_components": dict(r.energy_components),
                "source_row_reconstruction": row_checks,
                "native_energy_absolute_difference": abs(
                    r.potential_energy - ref.potential_energy
                ),
                "native_force_max_absolute_difference": force_error,
                "independent_component_errors": component_errors,
                "independent_finite_difference_force_max_errors": fd,
                "minimization": {
                    "converged": minimum.converged,
                    "final_verified": minimum.final_evaluation_verified,
                    "final_energy_kJ_mol": minimum.final_evaluation.potential_energy,
                },
                "finite_NVT_steps": 4,
                "split_continuation_state_exact": same_state,
                "split_continuation_RNG_exact": same_rng,
                "checkpoint_roundtrip_bytes_exact": restored.payload_json
                == checkpoint.payload_json,
                "checkpoint_file_sha256": file_hash(cp_path),
                "split_and_whole_checkpoint_envelopes_identical": False,
                "checkpoint_note": "Execution lineage/counters differ by segment; only physical final state and complete RNG are identical.",
            }
        )
    blocked = []
    fixtures = [
        ("disiloxane_SiH3", from_smiles("[SiH3]O[SiH3]", random_seed=2026)),
        ("HMDS", from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=2026)),
    ]
    for dp in (1, 2, 4):
        fixtures.append(
            (
                f"PDMS_DP{dp}",
                build_linear_polymer(
                    "[*:1]O[Si](C)(C)[*:2]",
                    dp=dp,
                    head_end_group="[*][Si](C)(C)C",
                    tail_end_group="[*]C",
                    add_hydrogens=True,
                    random_seed=2026,
                ),
            )
        )
    for name, s in fixtures:
        for opt, label in (
            (ExternalOPLSOptions(), "no_external_source"),
            (options, "CH_control_source"),
        ):
            p = prepare_forcefield(s, ForceFieldRequest("external_opls", opt))
            try:
                create_evaluator(s, p)
            except ExternalOPLSError as error:
                blocked.append(
                    {
                        "fixture": name,
                        "case": label,
                        "atoms": s.number_of_sites,
                        "graph_valid": p.payload["graph"] is not None,
                        "status": error.status,
                        "diagnostics": error.diagnostics,
                        "execution_profile": p.payload["profile_identity"],
                        "numerical_engine_instantiated": False,
                        "MD_run": False,
                    }
                )
            else:
                raise ValueError("Unexpected PDMS execution")
    return {
        "schema": "island_external_opls_integration_probe_v1",
        "software_integration_passed": True,
        "xml_sha256": file_hash(xml),
        "license_sha256": file_hash(license_file),
        "control_source_manifest_sha256": pin,
        "control_profile_identity": "island_external_opls_foyer_control_v1",
        "controls": controls,
        "PDMS_cases": blocked,
        "PDMS_executable": False,
        "PDMS_acceptance_gate": 1,
        "production_validated": False,
        "simulation_readiness": "not_established",
        "limitations": LIMITS,
        "scope": "Finite CH software controls using unmodified real Foyer source. No new physical parameters or PDMS execution profile.",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xml", required=True, type=Path)
    p.add_argument("--license", required=True, type=Path)
    p.add_argument("--work", required=True, type=Path)
    p.add_argument("--output", required=True, type=Path)
    a = p.parse_args()
    if a.work.exists() or a.output.exists():
        p.error("Work and output must be new paths")
    a.work.mkdir(parents=True)
    report = probe(a.xml, a.license, a.work)
    a.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "software_integration_passed": True,
                "PDMS_executable": False,
                "report_sha256": file_hash(a.output),
            }
        )
    )


if __name__ == "__main__":
    main()
