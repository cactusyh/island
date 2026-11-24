"""Install the published PolyJarvis OPLS-AA field into a separate field directory."""

import argparse
import os
import subprocess
from hashlib import sha256
from http.client import HTTPException
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.request import urlopen

from island.exceptions import ForceFieldError

REVISION = "c9fbfcccf617e59a89bda75bc4b166c6426ed4e8"
# SHA-256 of each completed file in the authors' emc_fields/manifest.json.
PATCHED = {
    "opls-aa.define": "83380b23eb8c609f4cc57b82c0e7b964af9df169544bbe3f01a2dfe80a6c9e38",
    "opls-aa.prm": "21d4c647c6cb116f25eabb3f9241e185a069b14797fe0750c2e989cec81caf73",
    "opls-aa.top": "e1e627c72987f1e3081a2cae6ed7120079515bc092c882944cadfa3ce0d26611",
}
STOCK = {
    "opls-aa.define": "834aa4789a92d5ca4f8a768de59959abad1db1e3406a52445d0bc160d83c8422",
    "opls-aa.prm": "19386ed3b4345d834a72f36f7cbae2e1c9fe60378e311e51e76906a48a0590ef",
    "opls-aa.top": "982a10a5e5d8922ed757fe8c0750f71687ba77f86a36f2798cb195f1a67361fe",
}
PATCHES = {
    "opls-aa.define": "ecb64aabc169e54480c8d265024e4c8873317d622c139069f97b5a3a27a2f61b",
    "opls-aa.prm": "5ff5183d8ed944a4ad0bca8733dd633e9d466429d322ec0fd8da343557e9c118",
    "opls-aa.top": "66eeaf1b748cbee7ae588e62eba9e68826a1467d4ef21d341a4fd83804921d17",
}


def install_polyjarvis_opls(destination, *, checkout=None):
    """Build and verify the authors' three field files, then publish together.

    destination is a field directory usable with assign_emc(field_dir=...).
    checkout can point to a local PolyJarvis checkout for offline installation.
    """

    def read(relative):
        if checkout is not None:
            return (Path(checkout) / "emc_fields" / relative).read_bytes()
        raw = f"https://raw.githubusercontent.com/arz-2/PolyJarvis/{REVISION}"
        mirror = f"https://cdn.jsdelivr.net/gh/arz-2/PolyJarvis@{REVISION}"
        try:
            with urlopen(f"{raw}/emc_fields/{relative}", timeout=10) as resource:
                return resource.read()
        except (OSError, HTTPException):
            with urlopen(f"{mirror}/emc_fields/{relative}", timeout=30) as resource:
                return resource.read()

    root = Path(destination).expanduser().resolve()
    target = root / "opls" / "2024"

    def verified():
        return all(
            (target / k).is_file()
            and sha256((target / k).read_bytes()).hexdigest() == v
            for k, v in PATCHED.items()
        )

    if target.exists():
        if verified():
            return root
        raise FileExistsError(f"Different OPLS field already present at {target}")
    target.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(dir=target.parent, prefix=".polyjarvis-") as temp:
        stage = Path(temp) / "2024"
        stage.mkdir()
        for name, pin in PATCHED.items():
            base = "opls-2024-" + name
            stock = read("stock/" + base)
            if sha256(stock).hexdigest() != STOCK[name]:
                raise ForceFieldError(
                    f"PolyJarvis stock field hash mismatch for {name}"
                )
            p = stage / name
            p.write_bytes(stock)
            patch = read(base + ".patch")
            if sha256(patch).hexdigest() != PATCHES[name]:
                raise ForceFieldError(f"PolyJarvis patch hash mismatch for {name}")
            subprocess.run(
                ["patch", "--batch", "--binary", "--fuzz=0", name],
                cwd=stage,
                input=patch,
                check=True,
                capture_output=True,
            )
            if sha256(p.read_bytes()).hexdigest() != pin:
                raise ForceFieldError(
                    f"Patched PolyJarvis field hash mismatch for {name}"
                )
        try:
            os.rename(stage, target)
        except OSError:
            if not verified():
                raise
    return root


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--checkout", type=Path)
    args = parser.parse_args()
    print(install_polyjarvis_opls(args.destination, checkout=args.checkout))
