import importlib.metadata as md
import json
import platform
import sys
from pathlib import Path

import psi4

from island.pcff_completion import digest, pin, signed, write

center = json.loads(Path("docs/evidence/phase_4j59/center-effects.json").read_text())[
    "coordinates"
]
m = psi4.geometry(
    "0 1\n"
    + "\n".join(
        e + " " + " ".join(map(str, x))
        for e, x in zip(["Si", "O", "Si"] + ["H"] * 6, center)
    )
    + "\nunits angstrom\nsymmetry c1\nno_com\nno_reorient"
)
bases = {}
for name in ("def2-qzvp", "def2-qzvpp", "cc-pvqz", "cc-pv5z"):
    b = psi4.core.BasisSet.build(m, "BASIS", name, quiet=True)
    shells = [
        {
            "center": b.shell_to_center(s),
            "angular_momentum": b.shell(s).am,
            "pure": b.shell(s).is_pure(),
            "exponents": [b.shell(s).exp(i) for i in range(b.shell(s).nprimitive)],
            "coefficients": [b.shell(s).coef(i) for i in range(b.shell(s).nprimitive)],
        }
        for s in range(b.nshell())
    ]
    f = Path(psi4.core.get_datadir()) / "basis" / (name + ".gbs")
    bases[name] = {
        "file": str(f),
        "sha256": pin(f),
        "nbf": b.nbf(),
        "nshell": b.nshell(),
        "shells": shells,
        "molecular_shell_identity": digest(shells),
        "per_element_functions": {
            e: sum(
                b.shell(s).nfunction
                for s in range(b.nshell())
                if b.shell_to_center(s) == i
            )
            for i, e in ((0, "Si"), (1, "O"), (3, "H"))
        },
    }
write(
    "docs/evidence/phase_4j60/engine/inventory.json",
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
            "def2_QZVP_QZVPP_molecular_shells_identical": bases["def2-qzvp"][
                "molecular_shell_identity"
            ]
            == bases["def2-qzvpp"]["molecular_shell_identity"],
            "basis_comparison": "cc-pVQZ vs cc-pV5Z is a correlation-consistent cardinal change; no CBS extrapolation or method adequacy claim. def2 QZ labels give identical molecular shells in this installed library.",
            "auxiliary_basis": {
                "file": "def2-universal-jkfit.gbs",
                "sha256": pin(
                    Path(psi4.core.get_datadir()) / "basis/def2-universal-jkfit.gbs"
                ),
            },
            "reference_settings_source": "../plan/reference-plan.json",
            "other_adapters": "unexecuted",
            "environment_modified": False,
        }
    ),
)
