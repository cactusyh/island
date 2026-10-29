"""Run explicit local msi2lmp as a failure oracle, never as chemical authority."""

import argparse
import json
import subprocess
from hashlib import sha256
from pathlib import Path

import numpy as np

from island.charge_references.records import unpack
from island.workflows.bundle import system_from
from scripts.validate_pcff_singlepoint import reference_inputs


def run(study, distribution, executable, expected_sha256, output):
    study, distribution, executable, output = map(
        Path, (study, distribution, executable, output)
    )
    executable = executable.resolve()
    assert sha256(executable.read_bytes()).hexdigest() == expected_sha256
    output.mkdir(parents=True, exist_ok=False)
    summary = json.loads((study / "report.json").read_text())
    source = (distribution / "pcff.frc").resolve()
    assert (
        sha256(source.read_bytes()).hexdigest() == summary["source_hashes"]["pcff.frc"]
    )
    cases = []
    for fixture in summary["fixtures"]:
        name = fixture["fixture"]
        folder = output / name
        folder.mkdir()
        p = unpack((study / f"{name}-parameters.json").read_text())
        system = system_from(unpack((study / f"{name}-graph.json").read_text()))
        xyz = np.array(
            [system.coordinates.get(i) for i in sorted(system.topology.sites)]
        )
        reference_inputs(
            system,
            p["typing"]["assignments"],
            p["exact_native_charge_vector"],
            xyz,
            folder,
        )
        argv = [
            str(executable),
            "reference",
            "-class",
            "II",
            "-frc",
            str(source),
            "-p",
            "3",
            "-nocenter",
        ]
        proc = subprocess.run(
            argv, cwd=folder, capture_output=True, text=True, timeout=120, check=False
        )
        (folder / "stdout.log").write_text(proc.stdout)
        (folder / "stderr.log").write_text(proc.stderr)
        first = next(
            (
                s.strip()
                for s in (proc.stdout + proc.stderr).splitlines()
                if "Unable" in s
            ),
            None,
        )
        record = {
            "fixture": name,
            "argv": argv,
            "exit_code": proc.returncode,
            "first_failure": first,
            "status": "parameter_incomplete"
            if proc.returncode
            else "scientifically_unqualified",
            "chemical_authority": False,
            "execution_authorized": False,
        }
        (folder / "command.json").write_text(json.dumps(record, indent=2) + "\n")
        cases.append(record)
    result = {
        "msi2lmp_sha256": expected_sha256,
        "source_sha256": summary["source_hashes"]["pcff.frc"],
        "cases": cases,
        "ignore_used": False,
        "backend_MD_executions": [],
        "profile_identity": None,
        "acceptance_gate": 1,
    }
    (output / "report.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("study", "distribution", "executable", "output"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--executable-sha256", required=True)
    args = parser.parse_args()
    print(
        json.dumps(
            run(
                args.study,
                args.distribution,
                args.executable,
                args.executable_sha256,
                args.output,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
