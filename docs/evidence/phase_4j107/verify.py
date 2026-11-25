"""Verify archived J107 evidence with the Python standard library only."""

import gzip
import hashlib
import importlib.util
import json
import tarfile
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    root = Path(__file__).resolve().parent
    manifest = json.loads((root / "manifest.json").read_text())
    for name, digest in manifest["files_sha256"].items():
        require(
            hashlib.sha256((root / name).read_bytes()).hexdigest() == digest,
            f"Evidence checksum mismatch: {name}",
        )
    for archive, checksums in manifest["archives"].items():
        expected = json.loads((root / checksums).read_text())
        observed = {}
        with tarfile.open(root / archive) as handle:
            for member in handle:
                require(member.isfile(), f"Unexpected archive member: {member.name}")
                require(member.name not in observed, "Duplicate archive member")
                observed[member.name] = hashlib.sha256(
                    handle.extractfile(member).read()
                ).hexdigest()
        require(observed == expected, f"Archive member mismatch: {archive}")

    inputs = json.loads((root / "inputs.json").read_text())
    with gzip.open(root / "results.jsonl.gz", "rt") as handle:
        rows = [json.loads(line) for line in handle]
    repo = root.parents[2]
    script = repo / "scripts/benchmark_final_graphs.py"
    spec = importlib.util.spec_from_file_location("j107_benchmark", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    summary = module.summarize(rows, inputs["jobs"])
    require(
        summary == json.loads((root / "summary.json").read_text()), "Summary mismatch"
    )
    require(summary["complete"] and len(rows) == 3996, "Incomplete final panel")
    require(
        hashlib.sha256(
            (repo / "docs/evidence/reference_assignment.csv.gz").read_bytes()
        ).hexdigest()
        == inputs["panel_sha256"],
        "Original 982-input panel changed",
    )
    key = lambda row: tuple(row[k] for k in ("panel", "id", "field", "mode"))
    jobs = {key(job): job for job in inputs["jobs"]}
    for row in rows:
        require(
            all(row[k] == v for k, v in jobs[key(row)].items()), "Job input changed"
        )
        pinned = inputs["sources"][row["field"]]["field_sha256"]
        if "preflight" in row:
            require(
                row["preflight"]["field_sha256"] == pinned, "Preflight source mismatch"
            )
        if row["assigned"]:
            require(row["field_sha256"] == pinned, "Assigned source mismatch")
        require(row["md"] == "unmeasured", "Coverage must not claim MD")

    reused = [r for r in rows if "reuse" in r]
    require(len(reused) == 1316, "Unexpected native-result reuse count")
    require(all(r["mode"] == "repeat" for r in reused), "Graph results must be rerun")
    with gzip.open(root / "exploratory-incomplete-results.jsonl.gz", "rt") as handle:
        require(sum(1 for _ in handle) == 2702, "Exploratory record count mismatch")
    numerical = json.loads((root / "numerical-and-md.json").read_text())
    require(
        len(numerical) == 11 and all(r["passed"] for r in numerical),
        "Numerical failure",
    )
    print(
        "Verified: 3996 final records, 2702 exploratory records, 11 numerical/MD groups;"
    )
    print("all evidence and archive member checksums match.")


if __name__ == "__main__":
    main()
