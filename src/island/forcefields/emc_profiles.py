"""Install explicitly selected typing profiles on a pinned EMC distribution."""

import argparse
import json
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory

from .emc import _field_sources

PCFF = {
    "cff/pcff/pcff.frc": "22974f591c98016c62e6986d8ca831713dbde38a203ed0e5805409b17fb00013",
    "cff/pcff/pcff_templates.dat": "c43b923bab8e68cf534902c79e54de01b502e6e25dd653910ba2aa41f06c45e6",
}


def carboxamide(template):
    # Both the carbonyl carbon and the amide N have explicit native PCFF types.
    rule = """type: n
    ! neutral carboxamide nitrogen (carbonyl substituent is C)
    template: (>N(-C(=O)(-C))(-*)(-*))
    atom_test:1
      formal_charge:0
      aromaticity:NON_AROMATIC
    end_test
end_type

"""
    return template.replace("precedence:", rule + "precedence:").replace(
        "(na(n3m)(n4m)(n_2))", "(na(n)(n3m)(n4m)(n_2))"
    )


def install(profile, directory, *, emc_root=None):
    if profile != "pcff-carboxamide-v1":
        raise ValueError(f"Unknown EMC profile: {profile}")
    field, pins, edits = "pcff", PCFF, {"cff/pcff/pcff_templates.dat": carboxamide}
    _, _, source, hashes = _field_sources(field, emc_root, None)
    if hashes != pins:
        raise ValueError("Installed field differs from the validated source")
    out = Path(directory)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise FileExistsError(out)
    with TemporaryDirectory(dir=out.parent) as tmp:
        stage = Path(tmp) / "fields"
        stage.mkdir()
        result = {}
        for name in hashes:
            path = stage / name
            path.parent.mkdir(parents=True, exist_ok=True)
            raw = (
                edits[name]((source / name).read_text()).encode()
                if name in edits
                else (source / name).read_bytes()
            )
            path.write_bytes(raw)
            result[name] = sha256(raw).hexdigest()
        manifest = {"profile": profile, "base": hashes, "installed": result}
        (stage / "profile.json").write_text(json.dumps(manifest, indent=2) + "\n")
        stage.rename(out)
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("profile")
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    print(json.dumps(install(args.profile, args.directory), indent=2))
