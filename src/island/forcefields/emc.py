"""EMC assignment of H-capped homopolymers using source-native LAMMPS outputs."""

import os
import platform
import re
import subprocess
from dataclasses import dataclass
from hashlib import sha256
from importlib.util import find_spec
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory

from island.exceptions import ForceFieldError


@dataclass(frozen=True)
class EMCAssignment:
    """Matching LAMMPS data and parameter commands in EMC atom-ID order."""

    family: str
    psmiles: str
    dp: int
    seed: int
    data: str
    parameters: str
    emc_version: str
    field_sha256: dict[str, str]


def _field_input(psmiles):
    from rdkit import Chem

    from island.chemistry import parse_psmiles

    # Re-serialize the parsed graph before inserting it into EMC's input grammar.
    if not isinstance(psmiles, str) or any(c.isspace() for c in psmiles):
        raise ValueError("Expected one repeat-unit SMILES without whitespace")
    unit = parse_psmiles(psmiles)
    mol = unit.mol
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        point.bond_order != 1 for point in unit.attachment_points
    ):
        raise ValueError("EMC requires a connected repeat with single-bond attachments")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        if atom.GetIsotope() or atom.GetNumRadicalElectrons():
            raise ValueError(
                "EMC assignment does not support isotopic or radical sites"
            )
    smiles = Chem.MolToSmiles(
        mol, rootedAtAtom=unit.head.dummy_atom_index, canonical=False
    )
    if "/" in psmiles or "\\" in psmiles:
        raise ValueError("EMC cannot preserve double-bond stereochemistry")
    explicit = Chem.AddHs(mol)
    return (
        smiles,
        explicit.GetNumAtoms() - 2,
        explicit.GetNumBonds() - 2,
        Chem.GetFormalCharge(mol),
    )


def _installation(root):
    if root is None:
        root = os.environ.get("EMC_ROOT")
    if root is None:
        spec = find_spec("pyemc")
        if spec is None or not spec.origin:
            raise ImportError("Install emc-pypi or set EMC_ROOT to an EMC installation")
        root = Path(spec.origin).parent / "emc"
    root = Path(root).expanduser().resolve()
    names = {
        ("Linux", "x86_64"): "emc_linux_x86_64",
        ("Linux", "aarch64"): "emc_linux_aarch64",
        ("Darwin", "x86_64"): "emc_macos",
        ("Darwin", "arm64"): "emc_macos",
        ("Windows", "AMD64"): "emc_win32.exe",
    }
    name = names.get((platform.system(), platform.machine()))
    if not name:
        raise OSError("No bundled EMC executable for this platform")
    return root, root / "bin" / name


def _run(command, work, env):
    result = subprocess.run(
        command,
        check=False,
        cwd=work,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        errors="replace",
    )
    if result.returncode:
        raise ForceFieldError(
            f"EMC {Path(command[0]).name} exited {result.returncode}:\n"
            + result.stdout[-3500:]
        )


def assign_emc(
    psmiles: str,
    family: str,
    *,
    dp: int,
    seed: int = 2026,
    emc_root: str | Path | None = None,
) -> EMCAssignment:
    """Type one H-capped chain and return its full data and coefficient files.

    The provider builds coordinates at 0.1 g/cm³ for assignment. No MD is run.
    family selects the EMC pcff or opls/2024/opls-aa database explicitly.
    """
    if family == "opls":
        family = "oplsaa"
    fields = {"pcff": "pcff", "oplsaa": "opls/2024/opls-aa"}
    if family not in fields:
        raise ValueError("EMC family must be pcff or oplsaa")
    if type(dp) is not int or dp < 1:
        raise ValueError("dp must be a positive integer")
    if type(seed) is not int or not 0 < seed < 2**31:
        raise ValueError("seed must be a positive signed 32-bit integer")
    smiles, atoms_per_repeat, bonds_per_repeat, repeat_charge = _field_input(psmiles)
    root, binary = _installation(emc_root)
    field = fields[family]
    field_files = (
        ("cff/pcff/pcff.frc", "cff/pcff/pcff_templates.dat")
        if family == "pcff"
        else tuple(field + ext for ext in (".define", ".prm", ".top"))
    )
    sources = {
        name: sha256((root / "field" / name).read_bytes()).hexdigest()
        for name in field_files
    }
    script = f"""ITEM OPTIONS
field {field}
number true
density 0.1
ITEM END

ITEM GROUPS
repeat {smiles},1,repeat:2
cap *[H],1,repeat:1,1,repeat:2
ITEM END

ITEM CLUSTERS
chain alternate 1
ITEM END

ITEM POLYMERS
chain
100 repeat,{dp},cap,2
ITEM END
"""
    with TemporaryDirectory(prefix="island-emc-") as directory:
        work = Path(directory)
        (work / "chain.esh").write_text(script)
        env = {**os.environ, "EMC_ROOT": str(root)}
        _run(["perl", root / "scripts/emc.pl", f"-seed={seed}", "chain.esh"], work, env)
        _run([binary, "build.emc"], work, env)
        data = (work / "chain.data").read_text()
        parameters = (work / "chain.params").read_text()

    expected_atoms = dp * atoms_per_repeat + 2
    expected_bonds = dp * bonds_per_repeat + dp + 1
    counts = {
        label: int(n)
        for n, label in re.findall(r"(?m)^\s*(\d+)\s+(atoms|bonds)\s*$", data)
    }
    if counts != {"atoms": expected_atoms, "bonds": expected_bonds}:
        raise ForceFieldError(
            f"EMC output size differs from the H-capped chain: {counts}"
        )
    if not re.search(r"(?m)^pair_coeff\s", parameters) or re.search(
        r"(?i)(?<!\S)[+-]?(?:nan|inf)(?!\S)", parameters + "\n" + data
    ):
        raise ForceFieldError("EMC returned missing or non-finite coefficients")
    atom_rows = (
        data.partition("\nAtoms\n")[2].lstrip().partition("\n\n")[0].splitlines()
    )
    if len(atom_rows) != expected_atoms:
        raise ForceFieldError("EMC returned an incomplete LAMMPS Atoms section")
    net_charge = fsum(float(row.split()[3]) for row in atom_rows)
    # EMC writes partial charges with four decimal places in its full-atom data.
    if abs(net_charge - dp * repeat_charge) > expected_atoms * 0.00005 + 1e-9:
        raise ForceFieldError(
            f"EMC chain charge differs from its formal charge: {net_charge}"
        )
    version = re.match(r"LAMMPS output created by (EMC [^\n]+)", data)
    if not version:
        raise ForceFieldError("Missing EMC version in LAMMPS output")
    return EMCAssignment(
        family, psmiles, dp, seed, data, parameters, version[1], sources
    )
