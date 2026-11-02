"""Run and retain fenced software shards/checks, never scientific workflows."""

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j72_c1_fence import assert_installed

REUSE = {
    "ISLAND_J32_TEST_RUN": "/tmp/island-j68-validation/pytest-terminal-1/j32-bundles0",
    "ISLAND_J33_C1_RUN": "/tmp/island-j68-validation/pytest-terminal-1/j33-c10/runs",
    "ISLAND_J34_RUN": "/tmp/island-j68-validation/pytest-terminal-5/j340/runs",
    "ISLAND_J35_RUN": "/tmp/island-j68-validation/pytest-terminal-1/j350/runs",
    "ISLAND_J39_RUN": "/tmp/island-j68-validation/pytest-terminal-3/j39-retained0/run",
    "ISLAND_J40_RUN": "/tmp/island-j68-validation/pytest-terminal-2/j400/run",
    "ISLAND_J40_C2_RUN": "/tmp/island-j68-validation/pytest-terminal-0/j40-c20/run",
}


def run(command, output, env, cwd=None):
    start = time.monotonic()
    with output.open("w") as stream:
        result = subprocess.run(
            command,
            env=env,
            cwd=cwd,
            stdout=stream,
            stderr=subprocess.STDOUT,
            check=False,
        )
    return {
        "command": command,
        "returncode": result.returncode,
        "log": output.name,
        "cwd": str(cwd or Path.cwd()),
        "runtime_seconds": time.monotonic() - start,
    }


def suite(output, workers):
    assert_installed()
    output.mkdir(parents=True, exist_ok=True)
    env = {**os.environ, **REUSE}
    collect = run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-p",
            "scripts.phase_4j72_c1_test_policy",
            "--collect-only",
            "-q",
            "--c1-collection-out",
            str(output / "collection.json"),
        ],
        output / "collection.log",
        env,
    )
    (output / "collection-result.json").write_text(json.dumps(collect, indent=2) + "\n")
    if collect["returncode"]:
        return 1
    nodes = json.loads((output / "collection.json").read_text())["nodes"]
    times = {}
    for case in ET.parse(
        "docs/evidence/phase_4j72/validation/terminal-results.xml"
    ).iter("testcase"):
        node = (
            case.attrib["classname"].replace(".", "/") + ".py::" + case.attrib["name"]
        )
        times[node] = float(case.attrib.get("time", 1))
    groups, loads = [[] for _ in range(workers)], [0.0] * workers
    for node in sorted(nodes, key=lambda n: (-times.get(n, 1), n)):
        i = min(range(workers), key=lambda j: loads[j])
        groups[i].append(node)
        loads[i] += times.get(node, 1)
    (output / "suite-policy.json").write_text(
        json.dumps(
            {
                "workers": workers,
                "threads_per_worker": 1,
                "reused_readonly_fixtures": REUSE,
                "scope_skips_are_not_scientific_passes": True,
                "resources": {
                    n: Path("/sys/fs/cgroup/" + n).read_text()
                    for n in ("memory.max", "memory.current", "memory.stat", "cpu.max")
                },
                "node_shards": groups,
                "estimated_seconds": loads,
            },
            indent=2,
        )
        + "\n"
    )
    processes = []
    for i, group in enumerate(groups):
        log = (output / f"full-{i}.log").open("w")
        command = [
            sys.executable,
            "-m",
            "pytest",
            "-p",
            "scripts.phase_4j72_c1_test_policy",
            "-q",
            *group,
            "--basetemp",
            str(output / f"pytest-{i}"),
            "--junitxml",
            str(output / f"full-{i}.xml"),
        ]
        processes.append(
            (
                subprocess.Popen(
                    command, env=env, stdout=log, stderr=subprocess.STDOUT
                ),
                log,
                command,
            )
        )
    print(
        f"Started {workers} fenced software shards ({len(nodes)} unique nodes)",
        flush=True,
    )
    results = []
    for i, (process, log, command) in enumerate(processes):
        code = process.wait()
        log.close()
        results.append({"shard": i, "returncode": code, "command": command})
        (output / "suite-results.json").write_text(json.dumps(results, indent=2) + "\n")
        print(i, code, flush=True)
    return int(any(r["returncode"] for r in results))


def checks(output):
    assert_installed()
    output.mkdir(parents=True, exist_ok=True)
    parent = Path.cwd().parent / "island-validation"
    previous = json.loads(
        Path(
            "docs/evidence/phase_4j72/validation/historical-and-environment-checks.json"
        ).read_text()
    )
    commands = [(r["name"], r["command"], Path(r["cwd"])) for r in previous]
    commands.append(
        (
            "J72_receipt",
            [sys.executable, "scripts/validate_phase_4j72_receipt.py"],
            Path.cwd(),
        )
    )
    # Check only supported scientific environments; the ambient environment's
    # already disclosed dependency conflicts remain historical advisory facts.
    commands += [
        (
            "ruff",
            [shutil.which("ruff"), "check", "src", "scripts", "tests"],
            Path.cwd(),
        ),
        (
            "compileall",
            [sys.executable, "-m", "compileall", "-q", "src", "scripts", "tests"],
            Path.cwd(),
        ),
        ("diff", ["git", "diff", "--check"], Path.cwd()),
    ]
    results = []
    for name, command, cwd in commands:
        environment = dict(os.environ)
        # The source tree for old receipts stays isolated; do not inject the
        # current src into their historical package path.
        if "historical-base" in str(cwd):
            environment["PYTHONPATH"] = (
                str(parent / "phase4j66-solvers-only") + ":.:src"
            )
        result = {
            "name": name,
            **run(command, output / (name + ".log"), environment, cwd),
        }
        results.append(result)
        (output / "checks.json").write_text(json.dumps(results, indent=2) + "\n")
        print(name, result["returncode"], flush=True)
    return int(any(r["returncode"] for r in results))


def corrections(output, source, workers):
    """Rerun whole affected files, retaining the initial XML without replacement."""
    assert_installed()
    output.mkdir(parents=True, exist_ok=True)
    seen_shards, scheduled, running, results = set(), set(), {}, []
    environment = {**os.environ, **REUSE}
    while len(seen_shards) < 4 or running or scheduled - {r["file"] for r in results}:
        for shard in range(4):
            file = source / f"full-{shard}.xml"
            if shard in seen_shards or not file.is_file():
                continue
            try:
                cases = list(ET.parse(file).iter("testcase"))
            except ET.ParseError:
                continue  # writer has not closed the original XML yet
            seen_shards.add(shard)
            for case in cases:
                if case.find("failure") is not None or case.find("error") is not None:
                    name = case.attrib["classname"].replace(".", "/") + ".py"
                    if name != "tests/test_phase_4j72_c1.py":
                        scheduled.add(name)
        started = {r["file"] for r in results} | set(running)
        for name in sorted(scheduled - started):
            if len(running) >= workers:
                break
            label = Path(name).stem
            log = (output / (label + ".log")).open("w")
            command = [
                sys.executable,
                "-m",
                "pytest",
                "-p",
                "scripts.phase_4j72_c1_test_policy",
                name,
                "-q",
                "--basetemp",
                str(output / (label + "-tmp")),
                "--junitxml",
                str(output / (label + ".xml")),
            ]
            process = subprocess.Popen(
                command, env=environment, stdout=log, stderr=subprocess.STDOUT
            )
            running[name] = process, log, command
        for name, (process, log, command) in list(running.items()):
            code = process.poll()
            if code is not None:
                log.close()
                results.append({"file": name, "returncode": code, "command": command})
                del running[name]
                print(name, code, flush=True)
        (output / "results.json").write_text(
            json.dumps(
                {
                    "source": str(source),
                    "seen_shards": sorted(seen_shards),
                    "affected_files": sorted(scheduled),
                    "results": results,
                    "running": sorted(running),
                    "workers": workers,
                },
                indent=2,
            )
            + "\n"
        )
        if (
            len(seen_shards) == 4
            and not running
            and scheduled == {r["file"] for r in results}
        ):
            break
        time.sleep(1)
    return int(any(r["returncode"] for r in results))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("mode", choices=("suite", "checks", "corrections"))
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--source", type=Path)
    p.add_argument("--workers", type=int, default=4)
    a = p.parse_args()
    if not 1 <= a.workers <= 8:
        p.error("workers must be between one and eight software processes")
    if a.mode == "corrections":
        if a.source is None:
            p.error("corrections requires --source")
        code = corrections(a.output.resolve(), a.source.resolve(), a.workers)
    else:
        code = (
            suite(a.output.resolve(), a.workers)
            if a.mode == "suite"
            else checks(a.output.resolve())
        )
    raise SystemExit(code)
