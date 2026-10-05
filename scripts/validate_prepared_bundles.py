"""Retained I1 preparation -> bundle -> relocate -> child reconstruction/evaluation."""

import argparse
import os
import shutil
import subprocess
import sys
import time
from contextlib import ExitStack
from pathlib import Path
from unittest.mock import patch

import numpy as np

from island.forcefields import (
    AmberBundleArtifacts,
    PreparedForceFieldSources,
    adopt_forcefield,
    create_evaluator,
    load_prepared_forcefield,
    save_prepared_forcefield,
)
from island.workflows import storage
from island.workflows.bundle import preparation_from, record, system_data, system_from


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    return storage.checksum(Path(path).read_bytes())


def sources(args):
    return PreparedForceFieldSources(opls_xml=args.xml, pcff_frc=args.frc)


def verify_files(root, files):
    for name, checksum in files.items():
        require(
            digest(storage.child(root, name)) == checksum,
            "Retained artifact hash differs: " + name,
        )


def retained(args):
    evidence = storage.read_json(args.evidence)["cases"][args.family]
    verify_files(args.retained, evidence["files"])
    require(
        storage.read_json(args.retained / "outcome.json") == evidence["outcome"],
        "Retained I1 outcome differs",
    )
    require(evidence["outcome"]["passed"], "I1 case did not pass")
    p = storage.read_json(args.retained / "declaration.json")["system"]
    # I1 declaration used ordinary JSON for this field. Restore only its explicitly
    # defined stable-ID coordinate keys; preserve recorded metadata without inference.
    require(
        all(type(k) is str and str(int(k)) == k for k in p["coordinates"]),
        "Noncanonical retained coordinate keys",
    )
    p["coordinates"] = {int(k): v for k, v in p["coordinates"].items()}
    system = system_from(p)
    artifacts = None
    source = None
    if args.family in ("gaff", "gaff2"):
        payload = storage.decode(storage.read_json(args.retained / "preparation.json"))
        checksum = payload["record"]["artifact_sha256"]["result.prmtop"]
        candidates = [
            name
            for name, h in evidence["files"].items()
            if name.endswith("/result.prmtop") and h == checksum
        ]
        require(
            len(candidates) == 1, "Exactly one manifest-listed original prmtop required"
        )
        prmtop = storage.child(args.retained, candidates[0])
        native = preparation_from(system, payload, prmtop)
        artifacts = AmberBundleArtifacts(prmtop)
    elif args.family == "oplsaa":
        from island.forcefields.oplsaa import (
            OPLSParameterizationResult,
            load_oplsaa_source,
        )

        source = load_oplsaa_source(args.xml)
        native = OPLSParameterizationResult(
            (args.retained / "parameters.json").read_text()
        )
    else:
        from island.forcefields.pcff import (
            load_pcff_model,
            load_pcff_parameters,
            load_pcff_source,
        )

        source = load_pcff_source(args.frc)
        assignment = load_pcff_parameters(
            args.retained / "parameters.json", source, system=system
        )
        native = load_pcff_model(
            args.retained / "model.json", assignment, system=system
        )
    prepared = adopt_forcefield(
        system, args.family, native, source=source if args.family == "oplsaa" else None
    )
    require(
        prepared.identity == evidence["outcome"]["wrapper_identity"],
        "Retained facade identity differs",
    )
    require(
        prepared.metadata == evidence["outcome"]["prepared"],
        "Retained facade metadata differs",
    )
    return system, prepared, artifacts, evidence


def compare(result, expected):
    for key in (
        "parameter_fingerprint",
        "model_fingerprint",
        "coordinate_fingerprint",
        "evaluation_fingerprint",
    ):
        require(getattr(result, key) == getattr(expected, key), "Changed " + key)
    require(
        set(result.energy_components) == set(expected.energy_components),
        "Component coverage changed",
    )
    require(set(result.forces) == set(expected.forces), "Force coverage changed")
    energy = abs(result.potential_energy - expected.potential_energy)
    components = max(
        abs(result.energy_components[k] - v)
        for k, v in expected.energy_components.items()
    )
    force = max(
        abs(result.forces[s][i] - expected.forces[s][i])
        for s in expected.forces
        for i in range(3)
    )
    a = [
        result.potential_energy,
        *[result.energy_components[k] for k in sorted(result.energy_components)],
        *[x for s in sorted(result.forces) for x in result.forces[s]],
    ]
    b = [
        expected.potential_energy,
        *[expected.energy_components[k] for k in sorted(expected.energy_components)],
        *[x for s in sorted(expected.forces) for x in expected.forces[s]],
    ]
    require(
        np.allclose(a, b, atol=1e-10, rtol=1e-12), "Numerical bundle comparison failed"
    )
    return {
        "energy_kj_mol": energy,
        "components_kj_mol": components,
        "forces_kj_mol_angstrom": force,
    }


def worker(args):

    from island.evaluation import EvaluationResult

    expected = storage.decode(storage.read_json(args.child / "expected.json"))
    require(os.getpid() != expected["parent_pid"], "Not a separate process")

    def forbidden(*a, **kw):
        raise RuntimeError("Preparation/execution invoked during bundle load")

    # Offline consistency arithmetic is permitted; scientific execution is not.
    with ExitStack() as stack:
        for name in (
            "subprocess.run",
            "subprocess.Popen",
            "openmm.Context",
            "island.forcefields.ambertools.AmberToolsParameterizationEngine.parameterize",
            "island.forcefields.oplsaa.parameterize_oplsaa",
            "island.forcefields.oplsaa.type_atoms",
            "island.forcefields.pcff.type_pcff_atoms",
            "island.forcefields.pcff.assign_automatic_pcff_charges",
            "island.forcefields.pcff.assign_pcff_parameters",
            "island.forcefields.pcff.define_pcff_model",
        ):
            stack.enter_context(patch(name, forbidden))
        loaded = load_prepared_forcefield(
            args.child / "relocated", sources=sources(args)
        )
    system, prepared = loaded.system, loaded.prepared
    require(
        prepared.identity == expected["prepared_identity"], "Facade identity differs"
    )
    require(
        system_data(system) == expected["system"], "Original system/provenance differs"
    )
    evaluator = create_evaluator(system, prepared)
    reference = EvaluationResult(**expected["evaluation"])
    comparisons = [compare(evaluator.evaluate_fresh(), reference)]
    with evaluator.open_session() as session:
        comparisons.extend(
            [
                compare(session.evaluate(), reference),
                compare(session.evaluate_fresh(), reference),
            ]
        )
    row = {
        "passed": True,
        "pid": os.getpid(),
        "parent_pid": expected["parent_pid"],
        "load_execution_guards_passed": True,
        "original_system_equal": True,
        "prepared_identity": prepared.identity,
        "comparisons": comparisons,
        "foyer_imported": "foyer" in sys.modules,
        "rdkit_imported": "rdkit" in sys.modules,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    storage.publish(args.child / "child-outcome.json", storage.json_bytes(row))
    return row


def execute(args):
    system, prepared, artifacts, evidence = retained(args)
    plan = {
        "schema": "island_prepared_bundle_acceptance_v1",
        "family": args.family,
        "input_hashes": evidence["files"],
        "prepared_identity": prepared.identity,
        "native_identity": prepared.metadata["native_identity"],
        "source_identity": prepared.metadata["source"],
        "atol": 1e-10,
        "rtol": 1e-12,
        "units": ["kJ/mol", "kJ/(mol*angstrom)"],
        "execution": "retained records only; no parameterization, QM or new coordinates",
        "charge_provenance": "synthetic provided +/-0.01 e fixture"
        if args.family in ("gaff", "gaff2")
        else "source-native",
        "python": sys.version,
        "outer_retries": 0,
    }
    storage.publish(args.output / "declaration.json", storage.json_bytes(plan))
    start = time.perf_counter()
    before = create_evaluator(system, prepared).evaluate_fresh()
    require(
        np.isclose(
            before.potential_energy,
            evidence["outcome"]["energy_kj_mol"],
            atol=1e-10,
            rtol=1e-12,
        ),
        "I1 retained energy differs",
    )
    storage.publish(
        args.output / "expected.json",
        storage.json_bytes(
            storage.encode(
                {
                    "parent_pid": os.getpid(),
                    "prepared_identity": prepared.identity,
                    "system": system_data(system),
                    "evaluation": record(before),
                }
            )
        ),
    )
    save_prepared_forcefield(
        system,
        prepared,
        args.output / "before-relocation",
        artifacts=artifacts,
        sources=sources(args),
    )
    shutil.move(args.output / "before-relocation", args.output / "relocated")
    command = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--child",
        str(args.output.resolve()),
    ]
    for key in ("xml", "frc"):
        value = getattr(args, key)
        if value is not None:
            command.extend(["--" + key, str(value.resolve())])
    storage.publish(args.output / "command.json", storage.json_bytes(command))
    completed = subprocess.run(
        command,
        cwd=args.output,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        timeout=600,
        check=False,
    )
    storage.publish(args.output / "child.log", completed.stdout.encode())
    require(completed.returncode == 0, "Child failed: " + completed.stdout[-2000:])
    child = storage.read_json(args.output / "child-outcome.json")
    verify_files(args.retained, evidence["files"])
    return {
        "passed": True,
        "family": args.family,
        "seconds": time.perf_counter() - start,
        "child": child,
        "retained_inputs_unchanged": True,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--family", choices=("gaff", "gaff2", "oplsaa", "pcff"))
    p.add_argument("--retained", type=Path)
    p.add_argument(
        "--evidence",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "docs/evidence/phase_4i1.json",
    )
    p.add_argument("--output", type=Path)
    p.add_argument("--xml", type=Path)
    p.add_argument("--frc", type=Path)
    p.add_argument("--child", type=Path)
    args = p.parse_args()
    if args.child:
        worker(args)
        return 0
    if not all((args.family, args.retained, args.output)):
        p.error("--family --retained --output required")
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        row = execute(args)
    except Exception as exc:  # noqa: BLE001 -- retain failures, never substitute inputs
        row = {
            "passed": False,
            "family": args.family,
            "failure": f"{type(exc).__name__}: {exc}",
        }
    storage.publish(args.output / "outcome.json", storage.json_bytes(row))
    storage.publish(
        args.output / "files.json",
        storage.json_bytes(
            {
                str(p.relative_to(args.output)): digest(p)
                for p in args.output.rglob("*")
                if p.is_file()
            }
        ),
    )
    print(row.get("failure", f"{args.family}: passed"))
    return int(not row["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
