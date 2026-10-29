"""Build an explicitly modified diagnostic kernel from supplied local LAMMPS.

No download; no change to the original checkout or pinned executable. Replaces
only the three/two erfc polynomial evaluations with libm erfc, in an isolated
copy. Charges, coefficients, Coulomb conversion constant, force expressions and
acceptance bounds are untouched. This executable is NEVER an author-source MD
profile and never grants execution authority to the canonical registry.
"""

import argparse
import json
import shutil
import subprocess
from pathlib import Path

from island.pdms_lammps import pin

PINS = {
    "KSPACE/pair_coul_long.cpp": "bd732a11df7a3691934e2d372723d011e1132b093d1b784cd9b5a1b120ca1661",
    "KSPACE/pair_lj_charmm_coul_long.cpp": "8efd83f9c3569b21598177adbf253c03a62fad0767d2d30147d0fe8f686b42e5",
}
REPLACEMENTS = {
    "KSPACE/pair_lj_charmm_coul_long.cpp": (
        "erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;",
        3,
    ),
    "KSPACE/pair_coul_long.cpp": (
        "erfc = t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * expm2;",
        2,
    ),
}


def patch(root):
    records = {}
    for name, (old, count) in REPLACEMENTS.items():
        path = root / "src" / name
        before = pin(path)
        if name in PINS and before != PINS[name]:
            raise ValueError("Unexpected local LAMMPS source: " + name)
        text = path.read_text()
        if text.count(old) != count:
            raise ValueError("Unexpected kernel expression: " + name)
        replacement = "erfc = std::erfc(grij); // J49 explicit exact-kernel diagnostic"
        path.write_text(text.replace(old, replacement))
        records[name] = {
            "original_sha256": before,
            "modified_sha256": pin(path),
            "original_expression": old,
            "replacement": replacement,
            "occurrences": count,
        }
    return records


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--jobs", type=int, default=8)
    a = p.parse_args(argv)
    a.output.mkdir(parents=True, exist_ok=False)
    source = a.output / "source"
    source.mkdir()
    for name in ("src", "cmake", "lib"):
        shutil.copytree(a.source / name, source / name, symlinks=True)
    records = patch(source)
    commands = [
        [
            "cmake",
            "-S",
            str(source / "cmake"),
            "-B",
            str(a.output / "build"),
            "-D",
            "BUILD_MPI=off",
            "-D",
            "PKG_KSPACE=on",
            "-D",
            "PKG_MOLECULE=on",
            "-D",
            "CMAKE_BUILD_TYPE=Release",
        ],
        ["cmake", "--build", str(a.output / "build"), "-j", str(a.jobs)],
    ]
    for n, command in enumerate(commands):
        with (a.output / f"build-{n}.log").open("w") as log:
            subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True)
    (a.output / "kernel.json").write_text(
        json.dumps(
            {
                "changes": records,
                "commands": commands,
                "executable_sha256": pin(a.output / "build/lmp"),
                "diagnostic_only": True,
                "MD_authority": False,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
