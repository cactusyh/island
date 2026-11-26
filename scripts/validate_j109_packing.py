"""Replay the 18 unique J105 geometries with independent per-chain conformers."""

import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from hashlib import sha256
from pathlib import Path


def execute(job):
    from island.workflows.recipe import run_recipe

    row, directory = job
    record = {"name": row["name"], "recipe": row["recipe"], "scope": "geometry_only"}
    try:
        result = run_recipe(
            row["recipe"],
            directory,
            lammps=os.environ["ISLAND_ASSIGNMENT_LAMMPS"],
            packmol=os.environ["ISLAND_PACKMOL"],
            stop_after="packing",
        )
        assert result == {"status": "paused", "stage": "packing"}
        record.update(
            passed=True,
            attempts=json.loads(
                (directory / "packing/packing-attempts.json").read_text()
            ),
        )
    except Exception as exc:  # noqa: BLE001 - complete denominator
        record.update(passed=False, error_type=type(exc).__name__, diagnostic=str(exc))
    return record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("out", type=Path)
    p.add_argument(
        "--panel", type=Path, default=Path("docs/evidence/phase_4j105/panel.json")
    )
    args = p.parse_args()
    rows = json.loads(args.panel.read_text())
    selected = [r for r in rows if r["recipe"]["forcefield"]["name"] == "gaff2"]
    assert len(selected) == 18
    # Geometry inputs must be identical across the old three-field matrix.
    for r in rows:
        representative = next(
            x
            for x in selected
            if x["chemistry"] == r["chemistry"]
            and x["recipe"]["architecture"] == r["recipe"]["architecture"]
        )
        assert {k: v for k, v in r["recipe"].items() if k != "forcefield"} == {
            k: v for k, v in representative["recipe"].items() if k != "forcefield"
        }
    args.out.mkdir(parents=True, exist_ok=False)
    records = []
    with ProcessPoolExecutor(2) as pool:
        futures = [
            pool.submit(execute, (row, args.out / row["name"])) for row in selected
        ]
        for future in as_completed(futures):
            record = future.result()
            records.append(record)
            report = {
                "panel_sha256": sha256(args.panel.read_bytes()).hexdigest(),
                "expected": 18,
                "recorded": len(records),
                "passed": sum(r["passed"] for r in records),
                "scope": "geometry_only",
                "records": records,
            }
            (args.out / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
            print(record["name"], record["passed"], flush=True)
    if not all(r["passed"] for r in records):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
