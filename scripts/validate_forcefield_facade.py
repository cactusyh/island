"""Bounded fresh-dispatch/adoption fidelity; run each backend in its own environment."""

import argparse
import sys
import time
from copy import deepcopy
from pathlib import Path

import numpy as np

from island.builders import build_linear_polymer
from island.forcefields import (
    AmberToolsOptions,
    ForceFieldRequest,
    OPLSOptions,
    PCFFOptions,
    adopt_forcefield,
    create_evaluator,
    prepare_forcefield,
)
from island.workflows import storage
from island.workflows.bundle import system_data


def execute(args):
    system = build_linear_polymer(
        "[*:1]CC[*:2]",
        dp=3,
        coordinate_method="local_templates",
        template_seed=2026,
        assembly_seed=2026,
    )
    before = deepcopy(system.to_dict())
    if args.family in ("gaff", "gaff2"):
        options = AmberToolsOptions(
            args.family,
            "provided",
            {
                i: (0.01 if n % 2 == 0 else -0.01)
                for n, i in enumerate(sorted(system.topology.sites))
            },
            amberhome=args.amberhome,
            work_root=args.output,
            retain_success_artifacts=True,
            charge_source="Synthetic alternating +/-0.01 e software fixture; not scientifically validated",
        )
    elif args.family == "oplsaa":
        options = OPLSOptions(args.source)
    else:
        options = PCFFOptions(args.source, (0, 0, 1), (0, 0, 1))
    plan = {
        "family": args.family,
        "system": system_data(system),
        "construction_seeds": [2026, 2026],
        "scope": "PE DP3 facade/native fidelity, no dynamics or QM",
        "atol": 1e-10,
        "rtol": 1e-12,
        "platform": "Reference",
        "python": sys.version,
        "provided_charges": "synthetic alternating +/-0.01 e"
        if args.family in ("gaff", "gaff2")
        else None,
        "source_sha256": storage.checksum(args.source.read_bytes())
        if args.source
        else None,
        "pcff_special_pairs": {"lj": [0, 0, 1], "coulomb": [0, 0, 1]}
        if args.family == "pcff"
        else None,
        "amber_timeout_seconds": 600,
        "outer_retries": 0,
    }
    storage.publish(args.output / "declaration.json", storage.json_bytes(plan))
    start = time.perf_counter()
    prepared = prepare_forcefield(system, ForceFieldRequest(args.family, options))
    native = prepared.native_result
    source = prepared.source
    adopted = adopt_forcefield(system, args.family, native, source=source)
    assert adopted.identity == prepared.identity
    if args.family in ("gaff", "gaff2"):
        from island.evaluation import OpenMMSinglePointEvaluator

        direct = OpenMMSinglePointEvaluator(system, native.imported_result)
        storage.publish(
            args.output / "preparation.json",
            storage.json_bytes(
                storage.encode(
                    {
                        "record": dict(native.record),
                        "record_signature": native.record_signature,
                        "import_source": native.imported_result.source,
                        "import_provenance": dict(native.imported_result.provenance),
                    }
                )
            ),
        )
        assert native.record["charge_outcome"]["qm_run"] is False
    elif args.family == "oplsaa":
        from island.evaluation.oplsaa import OPLSSinglePointEvaluator

        direct = OPLSSinglePointEvaluator(system, native, source)
        storage.publish(args.output / "parameters.json", native.json_text.encode())
    else:
        from island.evaluation import PCFFSinglePointEvaluator

        direct = PCFFSinglePointEvaluator(system, native)
        storage.publish(args.output / "model.json", native.json_text.encode())
        storage.publish(
            args.output / "parameters.json", native.assignment.json_text.encode()
        )
    reference = direct.evaluate_fresh()
    errors = []
    for facade in (prepared, adopted):
        evaluator = create_evaluator(system, facade)
        with evaluator.open_session() as session:
            for result in (
                evaluator.evaluate_fresh(),
                session.evaluate(),
                session.evaluate_fresh(),
            ):
                assert result.parameter_fingerprint == reference.parameter_fingerprint
                assert result.model_fingerprint == reference.model_fingerprint
                assert result.evaluation_fingerprint == reference.evaluation_fingerprint
                a = np.array(
                    [
                        result.potential_energy,
                        *result.energy_components.values(),
                        *[v for i in sorted(result.forces) for v in result.forces[i]],
                    ]
                )
                b = np.array(
                    [
                        reference.potential_energy,
                        *reference.energy_components.values(),
                        *[
                            v
                            for i in sorted(reference.forces)
                            for v in reference.forces[i]
                        ],
                    ]
                )
                assert np.allclose(a, b, atol=1e-10, rtol=1e-12)
                errors.append(float(np.max(abs(a - b))))
    assert system.to_dict() == before
    return {
        "passed": True,
        "family": args.family,
        "seconds": time.perf_counter() - start,
        "prepared": prepared.metadata,
        "wrapper_identity": prepared.identity,
        "parameter_fingerprint": reference.parameter_fingerprint,
        "model_fingerprint": reference.model_fingerprint,
        "energy_kj_mol": reference.potential_energy,
        "max_absolute_discrepancy": max(errors),
        "fresh_preparation": True,
        "adoption": True,
        "session_and_fresh": True,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--family", choices=("gaff", "gaff2", "oplsaa", "pcff"), required=True
    )
    p.add_argument("--source", type=Path)
    p.add_argument("--amberhome", type=Path)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        result = execute(args)
    except Exception as exc:  # noqa: BLE001 -- retain failed acceptance, no fallback
        result = {
            "passed": False,
            "family": args.family,
            "failure": f"{type(exc).__name__}: {exc}",
        }
    storage.publish(args.output / "outcome.json", storage.json_bytes(result))
    storage.publish(
        args.output / "files.json",
        storage.json_bytes(
            {
                str(f.relative_to(args.output)): storage.checksum(f.read_bytes())
                for f in args.output.rglob("*")
                if f.is_file()
            }
        ),
    )
    print(result.get("failure", f"{args.family}: passed"))
    return int(not result["passed"])


if __name__ == "__main__":
    raise SystemExit(main())
