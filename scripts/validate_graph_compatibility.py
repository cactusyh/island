"""Replay the pinned 25 J107 legacy-only cases with complete outcomes and process bounds."""

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path


def write(path, value):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    os.replace(temporary, path)


def worker(path):
    from island.forcefields.assignment import assign_polymer_parameters
    from island.forcefields.assignment_diagnostics import assignment_failure

    job = json.loads(path.read_text())
    started = time.monotonic()
    result = {"id": job["id"], "field": job["field"], "job": job}
    try:
        a = assign_polymer_parameters(
            job["psmiles"],
            job["field"],
            dp=4,
            field_dir=job["field_dir"],
            seed=2026,
            backend_timeout=job["backend_timeout"],
        )
        result.update(success=True, assignment=asdict(a))
    except Exception as exc:  # noqa: BLE001 - retain every failed input in the denominator
        result.update(success=False, failure=assignment_failure(exc))
    result["seconds"] = time.monotonic() - started
    write(path.with_name(path.name.replace(".input.json", ".result.json")), result)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--pcff-field", type=Path)
    parser.add_argument("--opls-field", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--backend-timeout", type=float, default=30)
    args = parser.parse_args()
    if args.case:
        worker(args.case)
        return
    if not all((args.out, args.pcff_field, args.opls_field)):
        parser.error("out and both field directories are required")
    if args.workers < 1 or not 0 < args.backend_timeout <= args.timeout:
        parser.error("invalid workers or timeout budget")
    source = (
        Path(__file__).resolve().parents[1]
        / "docs/evidence/phase_4j107/repeat-only-cases.json"
    )
    records = json.loads(source.read_text())
    assert len(records) == 25
    from island.forcefields._assignment_cache import source_binding
    from island.forcefields.emc import FIELDS, _field_sources

    directories = {"pcff": args.pcff_field.resolve(), "opls": args.opls_field.resolve()}
    bindings = {}
    for field, directory in directories.items():
        root, binary, _, _ = _field_sources(FIELDS[field], None, directory)
        bindings[field] = source_binding(root, binary, directory)
    manifest = {
        "input_sha256": sha256(source.read_bytes()).hexdigest(),
        "sources": bindings,
        "timeout": args.timeout,
        "backend_timeout": args.backend_timeout,
    }
    args.out.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out / "manifest.json"
    if manifest_path.exists() and json.loads(manifest_path.read_text()) != manifest:
        raise ValueError("Resume source or budget changed")
    write(manifest_path, manifest)

    def execute(record):
        name = record["id"] + "-" + record["field"]
        request, response = (
            args.out / (name + suffix) for suffix in (".input.json", ".result.json")
        )
        job = {k: record[k] for k in ["id", "field", "psmiles"]}
        job.update(
            field_dir=str(directories[job["field"]]),
            backend_timeout=args.backend_timeout,
        )
        write(request, job)
        if response.exists():
            result = json.loads(response.read_text())
            if result["job"] != job:
                raise ValueError("Resume input changed")
            return result
        with (args.out / (name + ".log")).open("w") as log:
            process = subprocess.Popen(
                [sys.executable, __file__, "--case", str(request)],
                stdout=log,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )
            try:
                code = process.wait(timeout=args.timeout)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
                code = "timeout"
        if code or not response.exists():
            result = {
                "id": job["id"],
                "field": job["field"],
                "job": job,
                "success": False,
                "failure": {
                    "exit_code": code,
                    "stage": "worker",
                    "diagnostic": (args.out / (name + ".log")).read_text()[-4000:],
                },
            }
            write(response, result)
        result = json.loads(response.read_text())
        print(name, result["success"], flush=True)
        return result

    with ThreadPoolExecutor(args.workers) as pool:
        results = list(pool.map(execute, records))
    summary = {
        "expected": 25,
        "recorded": len(results),
        "complete": len({(r["id"], r["field"]) for r in results}) == 25,
        "successful": sum(r["success"] for r in results),
        "md": "unmeasured",
        "cases": [{k: r[k] for k in ("id", "field", "success")} for r in results],
    }
    write(args.out / "summary.json", summary)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
