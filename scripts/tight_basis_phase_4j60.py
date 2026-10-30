"""Inspect actual installed tight-d shell additions before the refinement jobs."""

import json
from pathlib import Path

import psi4

from island.pcff_completion import digest, pin, signed, write

p = json.loads(Path("docs/evidence/phase_4j60/plan/reference-plan.json").read_text())
body = "\n".join(
    e + " " + " ".join(map(str, x))
    for e, x in zip(p["elements"], p["geometries"]["center"], strict=True)
)
m = psi4.geometry("0 1\n" + body + "\nunits angstrom\nsymmetry c1\nno_reorient\nno_com")
records = {}
for name, file in (
    ("cc-pvqz", "cc-pvqz.gbs"),
    ("cc-pv(q+d)z", "cc-pv_qpd_z.gbs"),
    ("cc-pv5z", "cc-pv5z.gbs"),
    ("cc-pv(5+d)z", "cc-pv_5pd_z.gbs"),
):
    b = psi4.core.BasisSet.build(m, "BASIS", name, quiet=True)
    path = Path(psi4.core.get_datadir()) / "basis" / file
    shells = [
        {
            "center": b.shell_to_center(s),
            "am": b.shell(s).am,
            "exponents": [b.shell(s).exp(i) for i in range(b.shell(s).nprimitive)],
            "coefficients": [b.shell(s).coef(i) for i in range(b.shell(s).nprimitive)],
        }
        for s in range(b.nshell())
    ]
    records[name] = {
        "basis_file": file,
        "sha256": pin(path),
        "nbf": b.nbf(),
        "nshell": b.nshell(),
        "shells": shells,
        "shell_identity": digest(shells),
        "native_header_citations": "\n".join(path.read_text().splitlines()[:28]),
    }
write(
    "docs/evidence/phase_4j60/engine/tight-d-basis-composition.json",
    signed(
        {
            "records": records,
            "changed_elements": ["Si"],
            "meaning": "One extra pure-d shell per Si adds5 functions, total10 per molecule. This is tight-polarization sensitivity at fixed cardinality; the corrected QZ/5Z pair then tests cardinality separately.",
            "not_CBS_extrapolation": True,
        }
    ),
)
