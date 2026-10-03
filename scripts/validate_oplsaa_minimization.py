"""Declare, then execute bounded OPLS minimization from retained 4G2 inputs."""

import argparse
import time
from dataclasses import asdict
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path
from unittest.mock import patch

import numpy as np
from validate_oplsaa_energy import make_system, reference, reference_value

from island.core.coordinate_provenance import coordinate_hash
from island.evaluation.oplsaa import OPLSSinglePointEvaluator
from island.forcefields.oplsaa import (
    OPLSParameterizationResult,
    load_oplsaa_source,
    parameterize_oplsaa,
)
from island.minimization import MinimizationOptions, minimize_geometry
from island.workflows import storage
from island.workflows.bundle import record, system_data, system_from

CASES = ("butane", "ethanol", "pe3", "peo3")
OPTIONS = MinimizationOptions(
    max_iterations=5000, max_evaluations=10000, force_tolerance=0.1
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def inputs(root, name, source):
    system = system_from(
        storage.decode(storage.read_json(root / f"{name}-system.json"))
    )
    parameters = OPLSParameterizationResult(
        (root / f"{name}-parameters.json").read_text()
    )
    parameters.validate_integrity(system, source)
    return system, parameters


def declaration(root, source):
    prior = storage.read_json(root / "outcomes.json")
    rows = {r["case"]: r for r in prior["cases"]}
    cases = {}
    for name in CASES:
        s, p = inputs(root, name, source)
        if (
            rows[name]["status"] != "passed"
            or p.identity != rows[name]["parameter_identity"]
        ):
            raise ValueError("Retained parameter identity contradicts 4G2 acceptance")
        cases[name] = {
            "files": {
                f"{name}-{suffix}.json": digest(root / f"{name}-{suffix}.json")
                for suffix in ("system", "parameters")
            },
            "parameter_identity": p.identity,
            "initial_coordinate_fingerprint": coordinate_hash(
                {i: s.coordinates.get(i) for i in s.topology.sites}
            ),
        }
    return {
        "schema": "island_opls_minimization_acceptance_v1",
        "source": source.identity,
        "prior_outcomes_sha256": digest(root / "outcomes.json"),
        "cases": cases,
        "options": asdict(OPTIONS),
        "platform": "Reference",
        "reference": "independent pinned Foyer.apply -> ParmEd geometric createSystem from 4G2",
        "energy_atol_kj_mol": 1e-5,
        "energy_rtol": 2e-10,
        "force_atol_kj_mol_angstrom": 1e-5,
        "force_rtol": 2e-10,
        "expected_rejection": {
            "case": "ps3",
            "psmiles": "[*:1]CC(c1ccccc1)[*:2]",
            "dp": 3,
            "template_seed": 2026,
            "assembly_seed": 2026,
            "reason": "Native component neutrality failed",
        },
        "versions": {
            n: version(n) for n in ("numpy", "scipy", "openmm", "foyer", "parmed")
        },
        "outer_retries": 0,
    }


def execute(root, output, xml, plan):
    import openmm as mm

    source = load_oplsaa_source(xml)
    outcomes = []
    for name in CASES:
        row = {"case": name, "status": "failed"}
        start = time.perf_counter()
        try:
            s, p = inputs(root, name, source)
            before = storage.json_bytes(storage.encode(system_data(s)))
            e = OPLSSinglePointEvaluator(s, p, source)
            constructions = []
            constructor = mm.Context

            def context(
                *args, constructor=constructor, constructions=constructions, **kwargs
            ):
                obj = constructor(*args, **kwargs)
                constructions.append(1)
                return obj

            with patch.object(mm, "Context", context), e.open_session() as session:
                minimum = minimize_geometry(s, session, OPTIONS)
            minimum.validate_integrity()
            storage.publish(
                output / f"{name}-minimum.json",
                storage.json_bytes(storage.encode(record(minimum))),
            )
            row.update(
                initial_energy=minimum.initial_energy,
                final_energy=minimum.final_energy,
                initial_fmax=minimum.initial_fmax,
                final_fmax=minimum.final_fmax,
                initial_rms_force=minimum.initial_rms_force,
                final_rms_force=minimum.final_rms_force,
                iterations=minimum.iterations,
                evaluations=minimum.evaluations,
                termination_reason=minimum.termination_reason,
                final_evaluation_verified=minimum.final_evaluation_verified,
                converged=minimum.converged,
                parameter_identity=p.identity,
                initial_coordinate_fingerprint=minimum.initial_evaluation.coordinate_fingerprint,
                final_coordinate_fingerprint=minimum.final_evaluation.coordinate_fingerprint,
                minimization_context_constructions=len(constructions),
                minimization_seconds=time.perf_counter() - start,
            )
            assert len(constructions) == 2, "Expected session and fresh final Context"
            ref, inspection = reference(s, xml)
            ids = sorted(s.topology.sites)
            xyz = np.array([minimum.coordinates[i] for i in ids])
            energy, forces, parts = reference_value(ref, xyz)
            final = minimum.final_evaluation
            actual = {
                k: final.energy_components[k] for k in ("bond", "angle", "rb_proper")
            }
            actual["nonbonded"] = (
                final.energy_components["lj"] + final.energy_components["coulomb"]
            )
            row.update(
                reference_context_constructions=1,
                energy_abs_error=abs(final.potential_energy - energy),
                force_max_abs_error=float(
                    np.max(np.abs(np.array([final.forces[i] for i in ids]) - forces))
                ),
                component_abs_errors={k: abs(actual[k] - parts[k]) for k in parts},
                reference_inspection=inspection,
            )
            assert np.isclose(final.potential_energy, energy, atol=1e-5, rtol=2e-10)
            assert np.allclose(
                [final.forces[i] for i in ids], forces, atol=1e-5, rtol=2e-10
            )
            assert all(
                np.isclose(actual[k], parts[k], atol=1e-5, rtol=2e-10) for k in parts
            )
            assert storage.json_bytes(storage.encode(system_data(s))) == before
            if not (
                minimum.converged
                and minimum.final_evaluation_verified
                and minimum.final_fmax <= 0.1
            ):
                raise ValueError(
                    f"Declared minimization gate unmet: {minimum.termination_reason}"
                )
            copied = minimum.to_system(s)
            copied.validate()
            p.validate_integrity(copied, source)
            storage.publish(
                output / f"{name}-optimized.json",
                storage.json_bytes(storage.encode(system_data(copied))),
            )
            row["status"] = "passed"
        except Exception as error:  # noqa: BLE001 -- bounded acceptance retains every failure
            row["failure"] = f"{type(error).__name__}: {error}"
        row["seconds"] = time.perf_counter() - start
        outcomes.append(row)
        print(
            name,
            row["status"],
            row.get("termination_reason"),
            row.get("final_fmax"),
            row.get("failure", ""),
            flush=True,
        )
        storage.publish(
            output / "outcomes.json",
            storage.json_bytes(
                {"declaration": plan, "cases": outcomes, "complete": False}
            ),
            replace=True,
        )
    try:
        parameterize_oplsaa(make_system("[*:1]CC(c1ccccc1)[*:2]", 3), source)
    except Exception as error:  # noqa: BLE001 -- expected rejection has a specific gate
        status = (
            "expected_rejection"
            if "Native component neutrality failed" in str(error)
            else "failed"
        )
        outcomes.append({"case": "ps3", "status": status, "failure": str(error)})
    else:
        outcomes.append(
            {
                "case": "ps3",
                "status": "failed",
                "failure": "Expected native charge rejection absent",
            }
        )
    # Recheck consumed source bytes after execution; no artifact may change in place.
    assert declaration(root, source) == plan
    complete = all(
        r["status"] == ("expected_rejection" if r["case"] == "ps3" else "passed")
        for r in outcomes
    )
    storage.publish(
        output / "outcomes.json",
        storage.json_bytes(
            {
                "declaration": plan,
                "cases": outcomes,
                "complete": complete,
                "production_validated": False,
                "simulation_readiness": "not_established",
            }
        ),
        replace=True,
    )
    return 0 if complete else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--retained", type=Path, required=True)
    parser.add_argument("--xml", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute-declared", action="store_true")
    args = parser.parse_args()
    source = load_oplsaa_source(args.xml)
    plan = declaration(args.retained, source)
    if not args.execute_declared:
        args.output.mkdir(parents=True, exist_ok=False)
        storage.publish(args.output / "declaration.json", storage.json_bytes(plan))
        return 0
    assert storage.read_json(args.output / "declaration.json") == plan
    assert {p.name for p in args.output.iterdir()} == {"declaration.json"}
    return execute(args.retained, args.output, args.xml, plan)


if __name__ == "__main__":
    raise SystemExit(main())
