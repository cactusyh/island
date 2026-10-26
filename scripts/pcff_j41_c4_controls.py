"""Execute reserialized C4 negative controls in fresh separate processes."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

from island.exceptions import PCFFError
from island.forcefields.pcff.mixed_sip import load_pcff_sip_distribution
from scripts.phase_4j41_c2 import read, write
from scripts.validate_phase_4j41_c3_receipt import reload_variant
from scripts.validate_phase_4j41_c4_receipt import validate_case

KINDS = (
    "source",
    "companion",
    "mixed_source",
    "graph",
    "order",
    "charge",
    "dependency",
    "coupling",
    "zero",
    "profile",
    "p_label",
)


def worker(folder):
    try:
        load_pcff_sip_distribution(folder)
        v = reload_variant(folder, "pcff.frc")
        validate_case(folder / "methylsilane", v, "mixed_silane")
    except (ValueError, PCFFError) as error:
        print(json.dumps({"rejected": True, "diagnostic": str(error)}))
        return 0
    print(json.dumps({"rejected": False}))
    return 1


def execute(study, output):
    output.mkdir(parents=True, exist_ok=False)
    results = []
    for kind in KINDS:
        folder = output / kind
        folder.mkdir()
        for n in ("pcff.frc", "pcff.rlb", "pcff_templates.dat", "inventory.json"):
            shutil.copyfile(study / "pcff" / n, folder / n)
        case = folder / "methylsilane"
        shutil.copytree(study / "pcff/methylsilane", case)
        p = read(case / "raw-reconstruction.json")
        if kind in ("source", "companion"):
            n = "pcff.frc" if kind == "source" else "pcff_templates.dat"
            with (folder / n).open("ab") as f:
                f.write(b"\n")
        elif kind == "mixed_source":
            shutil.copyfile(
                study / "compass_published/compass_published.frc", folder / "pcff.frc"
            )
        elif kind == "graph":
            p["system_sha256"] = "0" * 64
        elif kind == "order":
            p["atom_order"].reverse()
        elif kind == "charge":
            p["charges"]["partial_charges"][1] = "0"
        elif kind == "dependency":
            p["dependencies"].pop()
        elif kind == "coupling":
            p["missing_couplings"] = []
        elif kind == "zero":
            p["interactions"][0]["resolution"]["values"] = ["0", "0"]
        elif kind == "p_label":
            p["labels"][2] = "p="
        else:
            p["execution_profile"] = "invented_c4_profile"
        write(case / "raw-reconstruction.json", p)
        cmd = [sys.executable, str(Path(__file__).resolve()), "--worker", str(folder)]
        run = subprocess.run(
            cmd, capture_output=True, text=True, check=False, timeout=180
        )
        (folder / "worker.log").write_text(run.stdout + run.stderr)
        result = {
            "mutation": kind,
            "argv": cmd,
            "exit_code": run.returncode,
            "result": json.loads(run.stdout) if run.returncode == 0 else None,
        }
        results.append(result)
        print(kind, run.returncode, flush=True)
    report = {
        "results": results,
        "all_rejected": all(
            r["exit_code"] == 0 and r["result"]["rejected"] for r in results
        ),
    }
    (output / "mutation-audit.json").write_text(json.dumps(report, indent=2) + "\n")
    return 0 if report["all_rejected"] else 1


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--worker", type=Path)
    p.add_argument("--study", type=Path)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    if a.worker:
        raise SystemExit(worker(a.worker))
    if not a.study or not a.output:
        p.error("--study and --output are required")
    raise SystemExit(execute(a.study, a.output))
