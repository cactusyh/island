"""Validate retained records and measure cold/cache checks; no force evaluation."""

import argparse
import time
from hashlib import sha256
from pathlib import Path

from island.forcefields import PreparedForceFieldSources, load_prepared_forcefield
from island.forcefields.pcff import (
    clear_pcff_validation_cache,
    load_pcff_model,
    load_pcff_parameters,
    load_pcff_source,
    pcff_validation_cache_info,
)
from island.forcefields.pcff.automatic import graph_system
from island.forcefields.pcff.source import _parsed
from island.workflows import storage


def timed(function):
    start = time.perf_counter()
    result = function()
    return result, time.perf_counter() - start


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--declaration", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    d = storage.read_json(a.declaration)
    for path, expected in d["files"].items():
        if sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Input changed: {path}")
    clear_pcff_validation_cache()
    _parsed.cache_clear()
    loaded, cold = timed(
        lambda: load_prepared_forcefield(
            d["bundle"], sources=PreparedForceFieldSources(pcff_frc=d["source"])
        )
    )
    system, prepared = loaded.system, loaded.prepared
    _, hot = timed(lambda: [prepared.validate_integrity(system) for _ in range(20)])
    short = {
        "cold_seconds": cold,
        "twenty_validation_seconds": hot,
        "native_identity": prepared.native_result.identity,
        "facade_identity": prepared.identity,
        "cache": pcff_validation_cache_info(),
    }
    clear_pcff_validation_cache()
    _parsed.cache_clear()

    def long_load():
        source = load_pcff_source(d["source"])
        assignment = load_pcff_parameters(d["long_parameters"], source)
        system = graph_system(
            assignment.payload["charge_record"]["automatic_typing"]["graph"]
        )
        model = load_pcff_model(d["long_model"], assignment, system=system)
        return system, model

    (long_system, model), cold = timed(long_load)
    _, hot = timed(lambda: [model.validate_integrity(long_system) for _ in range(20)])
    long = {
        "sites": len(long_system.topology.sites),
        "cold_seconds": cold,
        "twenty_validation_seconds": hot,
        "model_identity": model.identity,
        "cache": pcff_validation_cache_info(),
        "role": "historical PE DP50 record integrity, not new strict coverage",
    }
    limits = d["budget_seconds"]
    passed = (
        short["cold_seconds"] <= limits["cold_bundle_load"]
        and short["twenty_validation_seconds"] <= limits["twenty_hot_checks"]
        and long["cold_seconds"] <= limits["long_chain_cold_load_and_validation"]
        and long["twenty_validation_seconds"] <= limits["twenty_long_hot_checks"]
    )
    for path, expected in d["files"].items():
        if sha256(Path(path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Input changed during validation: {path}")
    storage.publish(
        a.output,
        storage.json_bytes(
            {
                "short": short,
                "long": long,
                "budget_passed": passed,
                "input_hashes_unchanged": True,
                "declaration_sha256": sha256(a.declaration.read_bytes()).hexdigest(),
                "cold_definition": "native and source parser caches cleared; Python modules may already be imported",
            }
        ),
    )
    print({"budget_passed": passed, "short": short, "long": long})
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
