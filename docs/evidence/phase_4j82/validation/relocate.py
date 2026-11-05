"""Relocate J82 replay code/evidence; retain the explicitly recorded inherited fence."""
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from scripts.phase_4j82_coverage import ROOT, pin, replay, scope

scope()
repo = Path.cwd()
destination = Path(tempfile.mkdtemp(prefix="island-j82-relocated-"))
tree = subprocess.check_output(["git", "ls-tree", "-r", "--full-tree", "HEAD"], text=True)
inventory = {}
for row in tree.splitlines():
    name = row.split("\t", 1)[1]
    if not name.startswith(("src/island/", "scripts/")) or Path(name).suffix not in (".py", ".json"):
        continue
    src, target = repo / name, destination / name
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, target)
    inventory[name] = pin(src)
shutil.copytree(repo / ROOT, destination / ROOT)
policy = json.loads(Path(os.environ["ISLAND_C1_POLICY"]).read_text())
env = {**os.environ, "PYTHONPATH": os.pathsep.join((policy["bootstrap"], str(destination), str(destination / "src")))}
run = subprocess.run([sys.executable, "-m", "scripts.phase_4j82_coverage", "replay"], cwd=destination, env=env, capture_output=True, text=True, check=False)
if run.returncode == 0:
    assert json.loads(run.stdout) == replay(ROOT)
print(json.dumps({"destination": str(destination), "copied_code_sha256": inventory, "command": run.args, "returncode": run.returncode, "stdout": run.stdout, "stderr": run.stderr, "source_assets_copied": False, "network_required": False, "replay_code_and_evidence_relocated": True, "source_asset_root_required": False, "security_fence_inherited_from_original_checkout": True, "inherited_fence_runtime": policy["runtime"]}, indent=2, sort_keys=True))
raise SystemExit(run.returncode)
