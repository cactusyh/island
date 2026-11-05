"""Relocate offline replay inputs and an explicit tracked code inventory."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from scripts.phase_4j79_literature import ROOT, pin, replay, scope

scope()
repo = Path.cwd()
destination = Path(tempfile.mkdtemp(prefix="island-j79-relocated-"))
tree = subprocess.check_output(
    ["git", "ls-tree", "-r", "--full-tree", "HEAD"], text=True
)
files = [row.split("\t", 1)[1] for row in tree.splitlines()]
inventory = {}
for name in files:
    if not name.startswith(("src/island/", "scripts/")) or Path(name).suffix not in (".py", ".json"):
        continue
    src = repo / name
    target = destination / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, target)
    inventory[name] = pin(src)
shutil.copytree(repo / ROOT, destination / ROOT)
policy = json.loads(Path(os.environ["ISLAND_C1_POLICY"]).read_text())
env = dict(os.environ)
env["PYTHONPATH"] = os.pathsep.join(
    (policy["bootstrap"], str(destination), str(destination / "src"))
)
run = subprocess.run(
    [sys.executable, "-m", "scripts.phase_4j79_literature", "replay"],
    cwd=destination, env=env, capture_output=True, text=True, check=False
)
result = {
    "destination": str(destination),
    "relocated_source_files": inventory,
    "command": run.args,
    "returncode": run.returncode,
    "stdout": run.stdout,
    "stderr": run.stderr,
    "source_PDF_copied": False,
    "network_required": False,
    "original_checkout_required_by_replay": False,
}
if run.returncode == 0:
    assert json.loads(run.stdout) == replay(ROOT)
print(json.dumps(result, indent=2, sort_keys=True))
raise SystemExit(run.returncode)
