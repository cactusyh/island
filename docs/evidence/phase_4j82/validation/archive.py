"""Archive exact J82 fenced software checks and immutable validation metadata."""
import json
import lzma
import os
import shutil
import subprocess
from pathlib import Path
from xml.etree import ElementTree as ET

from scripts.phase_4j82_coverage import BASE, BRANCH, IMPLEMENTATION, ROOT, historical, pin, read, scope, write

scope()
repo = Path.cwd()
root = repo / ROOT
v = root / "validation"
v.mkdir(exist_ok=False)
commit = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
write(v / "tested-implementation.json", {"commit": commit, "files": {n: pin(repo / n) for n in IMPLEMENTATION}})
write(v / "repository-state.json", {"base_commit": BASE, "implementation_commit": commit, "branch": BRANCH, "fetched_origin_main": subprocess.check_output(["git", "rev-parse", "origin/main"], text=True).strip(), "main_modified": False, "historical_files_preserved": len(historical(repo)), "unrelated_untracked_files": "preserved; no checkout-copy/source assets used"})
shutil.copyfile("/tmp/j82-terminal-collection.json", v / "collection.json")
(v / "terminal.xml.xz").write_bytes(lzma.compress(Path("/tmp/j82-terminal.xml").read_bytes()))
active = Path(os.environ["ISLAND_C1_POLICY"]).parent
runs = []
for directory in sorted(Path("/tmp").glob("island_j82_*")):
    if not directory.is_dir() or not (directory / "policy.json").is_file() or directory == active:
        continue
    assert (directory / "result.json").is_file(), directory
    dest = v / "fences" / directory.name
    shutil.copytree(directory, dest)
    events = dest / "events.jsonl"
    (dest / "events.jsonl.xz").write_bytes(lzma.compress(events.read_bytes()))
    events.unlink()  # New archive copy only; original fence evidence stays in /tmp.
    runs.append(dest.relative_to(root).as_posix())
write(v / "fenced-runs.json", runs)
logs = v / "logs"
logs.mkdir()
for p in sorted(Path("/tmp").glob("island_j82_*.log")):
    if p.stem != active.name:
        (logs / (p.name + ".xz")).write_bytes(lzma.compress(p.read_bytes()))
checks = []
for name in ("terminal_tests", "ruff", "compileall", "pip_supported", "J81_receipt", "J80_receipt", "J79_receipt", "J78_receipt", "J76_receipt", "C1_receipt", "offline_replay", "relocated_replay", "diff"):
    label = "terminal" if name == "terminal_tests" else name
    revision = "_2" if label in ("terminal", "ruff", "compileall", "offline_replay", "relocated_replay", "diff") else ""
    attempt = "island_j82_final_" + label + revision
    result = read(Path("/tmp") / attempt / "result.json")
    assert result["returncode"] == 0, attempt
    checks.append({"name": name, "returncode": 0, "log": "validation/logs/" + attempt + ".log.xz"})
write(v / "checks.json", checks)
reasons = {
    "island_j82_focused_01": "The synthetic tacticity repeat had indistinguishable terminal substituents; the builder correctly rejected unrealizable CIP states. Corrected the test end caps; did not weaken the builder.",
    "island_j82_lint_01": "Development lint found import formatting, one unused import, nested conditional and missing explicit subprocess check flags. Corrected in new J82 files only.",
}
failures = []
for n in runs:
    r = read(root / n / "result.json")
    if r["returncode"]:
        assert Path(n).name in reasons, n
        failures.append({"run": n, "returncode": r["returncode"], "reason": reasons[Path(n).name], "scientific_execution": False})
write(v / "development-outcomes.json", {"failed_software_attempts": failures, "early_development_stdout": "Initial interactive smoke/focused/lint tool output was not separately redirected; policy, event and result bytes are retained. Every terminal verification has an archived stdout/stderr log.", "scientific_execution": False, "superseded_successful_suite": "island_j82_final_terminal: 978 passed, 58 skipped before final stereo correction; terminal delivery is final_terminal_2"})
cases = list(ET.parse("/tmp/j82-terminal.xml").iter("testcase"))
groups = {}
skips = []
for case in cases:
    group = "J82_focused" if case.attrib["classname"] == "tests.test_phase_4j82" else "regressions"
    counts = groups.setdefault(group, {"passed": 0, "skipped": 0, "failure": 0, "error": 0})
    outcome = next((n for n in ("failure", "error", "skipped") if case.find(n) is not None), "passed")
    counts[outcome] += 1
    if outcome == "skipped":
        skips.append({"test": case.attrib["classname"] + "::" + case.attrib["name"], "reason": case.find("skipped").attrib.get("message", "")})
shutil.copyfile("/tmp/island_j82_relocate.py", v / "relocate.py")
shutil.copyfile("/tmp/island_j82_archive.py", v / "archive.py")
(v / "superseded-terminal.xml.xz").write_bytes(lzma.compress(Path("/tmp/island_j82_terminal_initial.xml").read_bytes()))
write(v / "scope-skips.json", {"explicit_skips": skips, "scientific_checks_passed": 0, "scientific_execution": {"QM": 0, "MD": 0, "converter": 0, "physical_fitting": 0, "energy_force": 0, "reserved_validation": 0}, "synthetic_controls": "Manufactured graph/rule/parameter/audit/receipt controls only; fake native executables are rejected before execution.", "full_scientific_suite": "Not run: no scientific execution authorized."})
write(v / "summary.json", {"tested_implementation_commit": commit, "groups": groups, "unique_terminal_tests": len(cases), "historical_files_preserved": len(historical(repo)), "archived_fenced_attempts": len(runs), "failed_software_attempts": len(failures), "required_checks": len(checks), "native_pcff_unchanged": True, "runtime_registration": False, "scientific_validation": "not_established"})
print(json.dumps(read(v / "summary.json"), sort_keys=True))
