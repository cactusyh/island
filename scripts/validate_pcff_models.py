"""Construct model definitions from checksummed H3 records without rewriting them."""

import argparse
import json
import time
from hashlib import sha256
from pathlib import Path

from island.forcefields.pcff import (
    define_pcff_model,
    load_pcff_parameters,
    load_pcff_source,
    save_pcff_model,
    special_pair_policy,
)
from island.workflows.storage import json_bytes, publish


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--assignments", required=True)
    parser.add_argument("--evidence", default="docs/evidence/phase_4h3.json")
    parser.add_argument("--output", required=True)
    parser.add_argument("--lj", required=True, nargs=3, type=float)
    parser.add_argument("--coulomb", required=True, nargs=3, type=float)
    args = parser.parse_args(argv)
    policy = special_pair_policy(lj=args.lj, coulomb=args.coulomb)
    evidence = json.loads(Path(args.evidence).read_text())
    inputs = Path(args.assignments)
    expected = evidence["files"]
    # Validate every original file before consuming any scientific record.
    for name, checksum in expected.items():
        if sha256((inputs / name).read_bytes()).hexdigest() != checksum:
            raise ValueError("Historical input hash mismatch: " + name)
    source = load_pcff_source(args.source)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=False)
    publish(
        out / "declaration.json",
        json_bytes(
            {
                "schema": "island_pcff_model_acceptance_v1",
                "source": source.identity,
                "input_hashes": expected,
                "historical_evidence_sha256": sha256(
                    Path(args.evidence).read_bytes()
                ).hexdigest(),
                "special_pairs": policy,
                "required_cases": [o["name"] for o in evidence["report"]["outcomes"]],
                "gates": [
                    "historical integrity",
                    "model definition complete",
                    "policy-zero provenance",
                    "source hashes unchanged",
                ],
                "no_numerical_or_scientific_attestation": True,
            }
        ),
    )
    outcomes = []
    for case in evidence["report"]["outcomes"]:
        name = case["name"]
        start = time.perf_counter()
        try:
            assignment = load_pcff_parameters(inputs / (name + ".json"), source)
            before = assignment.identity
            result = define_pcff_model(assignment, special_pairs=policy)
            p = result.payload
            if not p["model_definition_complete"]:
                raise ValueError(str(p["diagnostics"]))
            zeros = [t for t in p["terms"] if t["origin"] == "policy_derived_zero"]
            assert len(zeros) == case["counts"][3]
            assert all(
                t["coefficients"] == [0]
                and not t["source_rows"]
                and t["raw_status"] == "missing"
                for t in zeros
            )
            path = out / (name + ".json")
            save_pcff_model(result, path)
            assert (
                sha256((inputs / (name + ".json")).read_bytes()).hexdigest()
                == expected[name + ".json"]
            )
            outcome = {
                "name": name,
                "status": "passed",
                "assignment_identity": before,
                "model_identity": result.identity,
                "model_definition_complete": True,
                "raw_parameter_coverage_complete": p["raw_parameter_coverage_complete"],
                "term_origins": p["term_origins"],
                "terms": len(p["terms"]),
                "special_pairs": len(p["special_pair_inventory"]),
                "bytes": path.stat().st_size,
                "file_sha256": sha256(path.read_bytes()).hexdigest(),
            }
        except Exception as error:  # noqa: BLE001 -- durable per-case outcomes
            outcome = {
                "name": name,
                "status": "failed",
                "error": f"{type(error).__name__}: {error}",
            }
        outcome["seconds"] = time.perf_counter() - start
        outcomes.append(outcome)
        print(json.dumps(outcome), flush=True)
        publish(out / "report.json", json_bytes({"outcomes": outcomes}), replace=True)
    for name, checksum in expected.items():
        assert sha256((inputs / name).read_bytes()).hexdigest() == checksum
    return int(any(o["status"] != "passed" for o in outcomes))


if __name__ == "__main__":
    raise SystemExit(main())
