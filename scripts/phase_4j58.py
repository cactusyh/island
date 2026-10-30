"""Run/replay the bounded J58 Psi4 pilot, with preserved native outputs."""

import argparse
import json
from pathlib import Path

from island.pcff_completion import require, write
from island.pcff_qm_pilot.analysis import diagnostics, identifiability
from island.pcff_qm_pilot.campaign import execute, prepare
from island.pcff_qm_pilot.importer import normalize


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "operation",
        choices=(
            "prepare",
            "execute",
            "normalize",
            "diagnostics",
            "identifiability",
            "recover",
        ),
    )
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--distribution", type=Path)
    p.add_argument("--evidence", type=Path, default=Path("docs/evidence/phase_4j58"))
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--python", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--normalized", type=Path)
    p.add_argument("--forcefield", default="pcff")
    a = p.parse_args()
    require(a.forcefield == "pcff", "wrong_profile_family")
    if a.operation == "prepare":
        print(json.dumps(prepare(a.repo, a.distribution, a.evidence, a.work), indent=2))
    elif a.operation == "execute":
        execute(a.repo, a.evidence, a.work, a.python)
    elif a.operation == "recover":
        from island.pcff_qm_pilot.recovery import recover

        for name in ("angle-angle", "bond-bond"):
            print(json.dumps(recover(a.work / "outputs" / name / "hessian")))
    elif a.operation == "normalize":
        print(
            json.dumps(normalize(a.repo, a.distribution, a.evidence, a.work, a.output))
        )
    elif a.operation == "diagnostics":
        write(a.output, diagnostics(a.evidence, a.work))
    else:
        write(
            a.output, identifiability(a.repo, a.distribution, a.evidence, a.normalized)
        )


if __name__ == "__main__":
    main()
