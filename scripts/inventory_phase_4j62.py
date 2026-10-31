"""Inspect actual installed orbital shells without executing QM calculations."""

import importlib.metadata as md
import json
import os
import platform
import sys
from pathlib import Path

import psi4

from island.pcff_completion import digest, pin, signed, write
from island.pcff_qm_pilot.control_engine import molecule, request

root = Path("docs/evidence/phase_4j62").resolve()
# Psi4 writes timer.dat at shutdown; use persistent task scratch, not user files.
scratch = Path.cwd().parent / "island-validation/phase4j62-qm/inventory-scratch"
scratch.mkdir(parents=True, exist_ok=True)
os.chdir(scratch)
plan = json.loads((root / "plan/scientific-plan.json").read_text())
bases = {}
for fixture, spec, names in [
    (
        "SiH",
        {
            **plan["control"]["center"],
            "elements": ["Si", "O", "Si"] + ["H"] * 6,
            "basis": "cc-pv(5+d)z",
            "grid": [150, 974],
            "thresholds": [1e-10, 1e-8],
        },
        ("cc-pv(5+d)z", "cc-pv(q+d)z"),
    ),
    ("HMDS", plan["HMDS"], ("def2-svp", "def2-tzvp")),
]:
    xyz = spec.get("coordinates", spec.get("initial_coordinates"))
    mol = molecule(psi4, request(spec, xyz))
    for name in names:
        b = psi4.core.BasisSet.build(mol, "BASIS", name, quiet=True)
        shells = [
            {
                "center": b.shell_to_center(s),
                "am": b.shell(s).am,
                "pure": b.shell(s).is_pure(),
                "exponents": [b.shell(s).exp(i) for i in range(b.shell(s).nprimitive)],
                "coefficients": [
                    b.shell(s).coef(i) for i in range(b.shell(s).nprimitive)
                ],
            }
            for s in range(b.nshell())
        ]
        path = (
            Path(psi4.core.get_datadir())
            / "basis"
            / (
                {
                    "cc-pv(5+d)z": "cc-pv_5pd_z.gbs",
                    "cc-pv(q+d)z": "cc-pv_qpd_z.gbs",
                }.get(name, name + ".gbs")
            )
        )
        bases[fixture + "/" + name] = {
            "nbf": b.nbf(),
            "nshell": b.nshell(),
            "shells": shells,
            "molecular_shell_identity": digest(shells),
            "source_file": str(path),
            "source_sha256": pin(path),
        }
write(
    root / "engine-environment.json",
    signed(
        {
            "engine": "Psi4",
            "version": psi4.__version__,
            "interpreter": sys.executable,
            "interpreter_sha256": pin(sys.executable),
            "library_sha256": pin(psi4.core.__file__),
            "python": sys.version,
            "platform": platform.platform(),
            "dependencies": {
                n: md.version(n)
                for n in ("numpy", "scipy", "qcelemental", "qcengine", "optking")
            },
            "bases": bases,
            "auxiliary_sha256": pin(
                Path(psi4.core.get_datadir()) / "basis/def2-universal-jkfit.gbs"
            ),
            "threads": 4,
            "memory_GiB": 4,
            "protocols": "plan/scientific-plan.json",
            "other_adapters": "unexecuted",
            "environment_modified": False,
            "basis_composition_only_no_QM_calls": True,
        }
    ),
)
