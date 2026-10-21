"""Additive J32 receipt, archived real runs, and strict read-only verification."""

import argparse
import json
import subprocess
import tarfile
import tempfile
from pathlib import Path
from types import SimpleNamespace

from island.analysis.network import _body, _digest, _signed
from scripts.phase_4j31_contract import typing_evidence
from scripts.phase_4j32_assignment import checksum
from scripts.phase_4j32_contract import CASES, validate_matrix
from scripts.phase_4j32_real_fixture import require
from scripts.validate_phase_4j31_receipt import junit_counts

BASE = "d212fbaac72583e5a91f1519166c17a65e6e7218"
BRANCH = "codex/phase-4j32-real-multichain-condensed-phase"
FILES = (
    "src/island/analysis/condensed.py",
    "scripts/phase_4j32_real_fixture.py",
    "scripts/phase_4j32_assignment.py",
    "scripts/phase_4j32_contract.py",
    "scripts/validate_phase_4j32.py",
    "scripts/phase_4j32_receipt.py",
    "scripts/validate_phase_4j32_gates.py",
    "tests/test_phase_4j32_condensed.py",
    "tests/test_phase_4j32_bundle.py",
    "docs/phase_4j32.md",
)
SUITES = {"j32", "j31_c1", "j26_j30", "j22_j23", "full_supported"}
CHECKS = {
    "ruff",
    "compileall",
    "pip_check",
    "foyer_pip_check",
    "analysis_pip_check",
    "git_diff_check",
    "historical_receipt",
    "retained_evidence",
}


def extract_archive(archive, destination):
    """Only ordinary relative files/directories; reject links and traversal."""
    with tarfile.open(archive, "r:gz") as stream:
        members = stream.getmembers()
        for m in members:
            require(
                not m.name.startswith("/")
                and ".." not in Path(m.name).parts
                and (m.isfile() or m.isdir()),
                "Unsafe archived member",
            )
        stream.extractall(destination, members=members, filter="data")


def historical_hashes():
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", BASE, "docs/evidence"], text=True
    ).splitlines()
    paths += ["docs/paper/reproducibility_manifest.json", ".gitattributes"]
    hashes = {}
    from hashlib import sha256

    for path in paths:
        data = subprocess.check_output(["git", "show", f"{BASE}:{path}"])
        require(Path(path).read_bytes() == data, "Historical evidence changed: " + path)
        hashes[path] = sha256(data).hexdigest()
    return hashes


def validate_runs(root):
    from island.analysis.condensed import load_bundle, validate_fixture

    matrix = validate_matrix(json.loads((root / "matrix.json").read_text()))
    for name in CASES:
        path = root / name.replace("/", "-")
        row = matrix["cases"][name]["result"]
        require(
            json.loads((path / "result.json").read_text()) == row,
            "Matrix/result mismatch",
        )
        if name != "pcff/network":
            record = load_bundle(
                path / "bundle", expected_identity=row["analysis_identity"]
            )
            require(
                record["bindings"]["final_graph_identity"]
                == row["final_graph_identity"]
                and record["setup"] == row["setup"]
                and record["execution"]["comparison"] == row["comparison"],
                "Case/bundle bindings mismatch",
            )
            for log in ("same-process.log", "fresh-process.log", "analysis.log"):
                require(
                    (path / log).read_text().strip() == row["analysis_identity"],
                    "Missing relocation/reconstruction output",
                )
        else:
            e = json.loads((path / "fixture.json").read_text())
            packed = validate_fixture(
                e, {"system": e["packed_system"], "final_graph": e["final_graph"]}
            )
            replay = typing_evidence(SimpleNamespace(packed=packed), path / "pcff.frc")
            require(
                replay == row["diagnostic"]["typing_evidence"],
                "PCFF typing replay mismatch",
            )
            log = (path / "expected-typing-rejection.log").read_text()
            require(
                "bridge_data" in log
                and "Incomplete automatic typing cannot assign charges" in log,
                "Missing actual PCFF rejection output",
            )
    return matrix


def validate_structure(receipt, root):
    require(receipt["identity"] == _digest(_body(receipt)), "Receipt identity mismatch")
    require(
        receipt["schema"] == "island.phase-4j32-evidence.v1"
        and receipt["base_commit"] == BASE
        and receipt["branch"] == BRANCH,
        "Wrong receipt base/branch/schema",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Readiness promotion",
    )
    require(receipt["remaining_unmet_gates"] == [], "Unmet mandatory gates")
    commit = receipt["validated_implementation_commit"]
    require(len(commit) == 40, "Full implementation SHA required")
    from hashlib import sha256

    require(
        set(receipt["implementation_files"]) == set(FILES), "Implementation inventory"
    )
    for path, expected in receipt["implementation_files"].items():
        original = subprocess.check_output(["git", "show", f"{commit}:{path}"])
        require(
            checksum(path) == expected == sha256(original).hexdigest(),
            "Tested implementation mismatch: " + path,
        )
    require(
        receipt["historical_hashes"] == historical_hashes(), "Historical hash inventory"
    )
    actual = {
        str(p.relative_to(root)): checksum(p)
        for p in sorted((root / "artifacts").rglob("*"))
        if p.is_file()
    }
    require(receipt["artifacts"] == actual, "Artifact hashes or inventory mismatch")
    validation = receipt["validation"]
    require(
        set(validation["suites"]) == SUITES and set(validation["checks"]) == CHECKS,
        "Validation inventory mismatch",
    )
    for name, row in {**validation["suites"], **validation["checks"]}.items():
        require(
            row["exit_code"] == 0
            and row["tested_commit"] == commit
            and row["tested_commit_end"] == commit,
            "Gate did not pass on tested commit: " + name,
        )
        if name in SUITES:
            counts = junit_counts(root / "artifacts" / (name + ".xml"))
            require(
                counts == row["counts"]
                and counts["passed"] > 0
                and counts["errors"] == counts["failures"] == 0,
                "Test count mismatch: " + name,
            )
    run = receipt["acceptance_command"]
    require(
        run["exit_code"] == 0
        and run["tested_commit"] == commit
        and run["tested_commit_end"] == commit,
        "Acceptance executed on another commit",
    )
    validate_matrix(receipt["matrix"])
    return receipt


def validate(path):
    path = Path(path)
    receipt = validate_structure(json.loads(path.read_text()), path.parent)
    with tempfile.TemporaryDirectory(prefix="island-j32-readonly-") as temp:
        root = Path(temp)
        extract_archive(path.parent / "artifacts/acceptance.tar.gz", root)
        require(
            validate_runs(root) == receipt["matrix"], "Receipt/archived matrix mismatch"
        )
    return receipt


def generate(directory, destination):
    """Publish new evidence only; never rewrite historical or existing output."""
    import gzip
    import shutil

    destination.mkdir(parents=True, exist_ok=False)
    artifacts = destination / "artifacts"
    artifacts.mkdir()
    validation = {"suites": {}, "checks": {}}
    for name in sorted(SUITES | CHECKS | {"acceptance"}):
        row = json.loads((directory / (name + "-command.json")).read_text())
        require(row["exit_code"] == 0, "Cannot publish failed gate: " + name)
        for extension in (".log", "-command.json"):
            shutil.copyfile(
                directory / (name + extension), artifacts / (name + extension)
            )
        if name in SUITES:
            shutil.copyfile(directory / (name + ".xml"), artifacts / (name + ".xml"))
            row["counts"] = junit_counts(artifacts / (name + ".xml"))
            validation["suites"][name] = row
        elif name in CHECKS:
            validation["checks"][name] = row
    shutil.copyfile(
        directory / "retained-evidence.json", artifacts / "retained-evidence.json"
    )
    for source in directory.glob("preflight-*"):
        if source.is_file():
            shutil.copyfile(source, artifacts / source.name)
    with tarfile.open(
        artifacts / "software-negative-controls.tar.gz", "w:gz"
    ) as controls:
        for source in sorted((directory / "j32-tmp").rglob("*")):
            if source.is_file() and source.name in {"result.json", "failure.log"}:
                controls.add(
                    source,
                    arcname=str(source.relative_to(directory / "j32-tmp")),
                    recursive=False,
                )
    archive = artifacts / "acceptance.tar.gz"
    run = directory / "acceptance"
    selected = [run / "matrix.json"]
    for case in CASES:
        path = run / case.replace("/", "-")
        selected.extend(p for p in path.iterdir() if p.is_file())
        selected.extend(p for p in (path / "bundle").rglob("*") if p.is_file())
    with (
        archive.open("wb") as raw,
        gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as zipped,
        tarfile.open(fileobj=zipped, mode="w") as stream,
    ):
        for path in sorted(set(selected)):
            info = stream.gettarinfo(str(path), arcname=str(path.relative_to(run)))
            info.uid = info.gid = info.mtime = 0
            info.uname = info.gname = ""
            with path.open("rb") as source:
                stream.addfile(info, source)
    probe = "import sys,json,importlib.metadata as m; out={'python':sys.version,'executable':sys.executable};\nfor n in ('numpy','rdkit','openmm','parmed','pytest','scipy','foyer','MDAnalysis'):\n try: out[n]=m.version(n)\n except m.PackageNotFoundError: out[n]='not installed'\nprint(json.dumps(out))"
    from scripts.validate_phase_4j32_gates import ANALYSIS, FOYER, PYTHON

    environments = {
        name: json.loads(subprocess.check_output([executable, "-c", probe], text=True))
        for name, executable in (
            ("supported", PYTHON),
            ("pinned_foyer", FOYER),
            ("analysis", ANALYSIS),
        )
    }
    receipt = _signed(
        {
            "schema": "island.phase-4j32-evidence.v1",
            "branch": BRANCH,
            "base_commit": BASE,
            "validated_implementation_commit": subprocess.check_output(
                ["git", "rev-parse", "HEAD"], text=True
            ).strip(),
            "implementation_files": {p: checksum(p) for p in FILES},
            "historical_hashes": historical_hashes(),
            "validation": validation,
            "acceptance_command": json.loads(
                (directory / "acceptance-command.json").read_text()
            ),
            "environments": environments,
            "matrix": json.loads((run / "matrix.json").read_text()),
            "artifacts": {
                str(p.relative_to(destination)): checksum(p)
                for p in sorted(artifacts.rglob("*"))
                if p.is_file()
            },
            "remaining_unmet_gates": [],
            "negative_control_scope": "Injected failures and mutated records are software controls, not scientific evidence; see focused JUnit and software-negative-controls.tar.gz",
            "short_run_scope": "Each successful case executes 1000 uninterrupted steps and 500+500 resumed steps; historical integration tests execute their existing short workflows, never equilibrium MD",
            "blocked_or_unexecuted_scope": [
                "PCFF cyclic network: verified expected typing rejection, no dynamics",
                "GAFF: unexecuted; GAFF2 selected for Amber family",
                "GAFF2 disconnected multichain input: outside existing native adapter domain",
                "NPT: unexecuted; no density relaxation claim",
                "Linear system Rg: ambiguous relative images of disconnected components",
                "Equilibrium, Tg, mechanics, diffusion coefficient: not established by 0.5 ps acceptance trajectory",
                "LAMMPS KSPACE/Ewald/PPPM and universal force-field coverage: unvalidated",
            ],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    (destination / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n"
    )
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "receipt",
        nargs="?",
        type=Path,
        default=Path("docs/evidence/phase_4j32/receipt.json"),
    )
    parser.add_argument("--runs", type=Path, help="validate extracted actual runs only")
    parser.add_argument("--generate", type=Path, help="completed validation directory")
    args = parser.parse_args()
    if args.generate:
        generate(args.generate, args.receipt.parent)
    elif args.runs:
        validate_runs(args.runs)
    else:
        validate(args.receipt)
    print("J32 receipt validation: PASS")


if __name__ == "__main__":
    main()
