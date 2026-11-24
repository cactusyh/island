"""EMC assignment with source-native LAMMPS data, coefficients and conventions."""

import os
import platform
import re
import subprocess
from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
from importlib.util import find_spec
from math import fsum, isfinite
from pathlib import Path
from tempfile import TemporaryDirectory

from island.exceptions import ForceFieldError

FIELDS = {
    "pcff": "pcff",
    "pcff_ore": "pcff_ore",
    "opls": "opls/2024/opls-aa",
    "oplsaa": "opls/2024/opls-aa",
    **{
        f"opls/{year}/opls-{mode}": f"opls/{year}/opls-{mode}"
        for year in ("2012", "2017", "2024")
        for mode in ("aa", "ua")
    },
}


@dataclass(frozen=True)
class EMCAssignment:
    """Matching LAMMPS files in EMC atom-ID order. UA chains have methyl end caps."""

    family: str
    psmiles: str
    dp: int
    seed: int
    data: str
    parameters: str
    emc_version: str
    field_sha256: dict[str, str]
    field: str = ""
    emc_smiles: str = ""
    bond_stereo_omitted: bool = False
    nchains: int = 1
    density: float = 0.1
    automatic_ring_depth: bool = False


def _repeat(psmiles, stereo_policy):
    from rdkit import Chem

    from island.chemistry import parse_psmiles

    if not isinstance(psmiles, str) or any(c.isspace() for c in psmiles):
        raise ValueError("Expected one repeat-unit SMILES without whitespace")
    if stereo_policy not in ("strict", "polyjarvis"):
        raise ValueError("stereo_policy must be strict or polyjarvis")
    directional = "/" in psmiles or "\\" in psmiles
    if directional and stereo_policy == "strict":
        raise ValueError("EMC cannot preserve double-bond stereochemistry")
    parsed = psmiles.translate(str.maketrans("", "", "/\\")) if directional else psmiles
    unit = parse_psmiles(parsed)
    mol = unit.mol
    if len(Chem.GetMolFrags(mol)) != 1 or any(
        p.bond_order != 1 for p in unit.attachment_points
    ):
        raise ValueError("EMC requires a connected repeat with single-bond attachments")
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
        if atom.GetNumRadicalElectrons():
            raise ValueError("EMC assignment does not support radical sites")
    return unit, directional


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
    return result.stdout


def _build(
    smiles,
    field,
    root,
    binary,
    field_dir,
    dp,
    seed,
    nchains,
    density,
    *,
    auto_depth=False,
):
    ua = field.endswith("-ua")
    repeat = smiles + (",1,repeat:2, 1,cap:1, 2,cap:1" if ua else ",1,repeat:2")
    cap = "*C" if ua else "*[H],1,repeat:1,1,repeat:2"
    script = f"""ITEM OPTIONS
field {field}
number true
density {density}
ITEM END

ITEM GROUPS
repeat {repeat}
cap {cap}
ITEM END

ITEM CLUSTERS
chain alternate {nchains}
ITEM END

ITEM POLYMERS
chain
100 repeat,{dp},cap,2
ITEM END
"""
    return _build_script(script, field, root, binary, field_dir, seed, auto_depth)


def _build_script(
    script, field, root, binary, field_dir, seed, auto_depth=False, *, psf=False
):
    with TemporaryDirectory(prefix="island-emc-") as directory:
        work = Path(directory)
        (work / "chain.esh").write_text(script)
        env = {**os.environ, "EMC_ROOT": str(root)}
        _run(
            [
                "perl",
                root / "scripts/emc.pl",
                f"-field_location={field_dir}",
                f"-seed={seed}",
                *(["-emc_depth=auto"] if auto_depth else []),
                "chain.esh",
            ],
            work,
            env,
        )
        log = _run([binary, "build.emc"], work, env)
        parameter_file = (
            f"cff/{field}/{field}.frc" if field.startswith("pcff") else field + ".prm"
        )
        imported = re.findall(r"(?m)^Info: importing '([^']+)'", log)
        if (field_dir / parameter_file).resolve() not in {
            Path(path).resolve() for path in imported
        }:
            raise ForceFieldError(
                f"EMC did not confirm loading {field_dir / parameter_file}"
            )
        missing = re.findall(
            r"(?im)^.*(?:increment pair .* not found|no (?:\w+ )*coefficients found for).*$",
            log,
        )
        if missing:
            raise ForceFieldError(
                "EMC reported missing parameters:\n" + "\n".join(missing)
            )
        data, params = (
            (work / "chain.data").read_text(),
            (work / "chain.params").read_text(),
        )
        if psf:
            import gzip

            with gzip.open(work / "chain.psf.gz", "rt") as handle:
                return data, params, handle.read()
        return data, params


def _field_sources(field, emc_root, field_dir):
    root, binary = _installation(emc_root)
    field_dir = Path(field_dir).resolve() if field_dir is not None else root / "field"
    files = (
        (f"cff/{field}/{field}.frc", f"cff/{field}/{field}_templates.dat")
        if field.startswith("pcff")
        else tuple(field + ext for ext in (".define", ".prm", ".top"))
    )
    sources = {
        name: sha256((field_dir / name).read_bytes()).hexdigest() for name in files
    }
    return root, binary, field_dir, sources


def assign_emc(
    psmiles: str,
    family: str,
    *,
    dp: int,
    seed: int = 2026,
    emc_root: str | Path | None = None,
    field_dir: str | Path | None = None,
    stereo_policy: str = "strict",
    nchains: int = 1,
    density: float = 0.1,
) -> EMCAssignment:
    """Assign a specific EMC database. polyjarvis stereo policy permits unspecified E/Z."""
    from rdkit import Chem

    if family not in FIELDS:
        raise ValueError(f"Unknown EMC force field {family!r}")
    if type(dp) is not int or dp < 1 or type(nchains) is not int or nchains < 1:
        raise ValueError("dp and nchains must be positive integers")
    if type(seed) is not int or not 0 < seed < 2**31:
        raise ValueError("seed must be a positive signed 32-bit integer")
    if isinstance(density, bool) or not isfinite(density) or density <= 0:
        raise ValueError("density must be a positive finite number")
    unit, omitted = _repeat(psmiles, stereo_policy)
    mol = unit.mol
    isotopic = any(a.GetIsotope() for a in mol.GetAtoms())
    field = FIELDS[family]
    ua = field.endswith("-ua")
    if ua and isotopic:
        raise ValueError("Isotope resolution requires an all-atom EMC field")
    root, binary, field_dir, sources = _field_sources(field, emc_root, field_dir)
    explicit = Chem.AddHs(mol)
    expected_atoms = (dp * (explicit.GetNumAtoms() - 2) + 2) * nchains
    cycle_rank = mol.GetNumBonds() - mol.GetNumAtoms() + 1
    previous_smiles = None
    auto_depth = False
    for kekule in (False, True):
        # EMC can add excess hydrogens in fused aromatic systems. Retry with
        # the same graph drawn with alternating single and double bonds.
        smiles = Chem.MolToSmiles(
            mol,
            rootedAtAtom=unit.head.dummy_atom_index,
            canonical=False,
            kekuleSmiles=kekule,
        )
        if smiles == previous_smiles:
            continue
        previous_smiles = smiles
        args = smiles, field, root, binary, field_dir, dp, seed, nchains, density
        try:
            data, params = _build(*args, auto_depth=auto_depth)
        except ForceFieldError as exc:
            if auto_depth or "increase scanning depth?" not in str(exc):
                raise
            auto_depth = True
            data, params = _build(*args, auto_depth=True)
        counts = {
            label: int(n)
            for n, label in re.findall(r"(?m)^\s*(\d+)\s+(atoms|bonds)\s*$", data)
        }
        atoms = counts.get("atoms", 0)
        size_ok = (
            (dp * mol.GetNumHeavyAtoms() + 2) * nchains
            <= atoms
            <= expected_atoms + 6 * nchains
            if ua
            else atoms == expected_atoms
        )
        if (
            size_ok
            and counts.get("bonds") == atoms - nchains + dp * nchains * cycle_rank
        ):
            break
    else:
        raise ForceFieldError(
            f"EMC output size differs from the requested chain: {counts}"
        )
    if not re.search(r"(?m)^pair_coeff\s", params) or re.search(
        r"(?i)(?<!\S)[+-]?(?:nan|inf)(?!\S)", params + "\n" + data
    ):
        raise ForceFieldError("EMC returned missing or non-finite coefficients")
    atom_rows = [
        line.split()
        for line in data.partition("\nAtoms\n")[2]
        .lstrip()
        .partition("\n\n")[0]
        .splitlines()
    ]
    per_chain = defaultdict(list)
    for row in atom_rows:
        per_chain[row[1]].append(float(row[3]))
    if (
        len(atom_rows) != atoms
        or len(per_chain) != nchains
        or any(len(charges) != atoms // nchains for charges in per_chain.values())
    ):
        raise ForceFieldError(
            "EMC returned an incomplete chain or LAMMPS Atoms section"
        )
    for charges in per_chain.values():
        net_charge = fsum(charges)
        if (
            abs(net_charge - dp * Chem.GetFormalCharge(mol))
            > len(charges) * 0.00005 + 1e-9
        ):
            raise ForceFieldError(
                f"EMC chain charge differs from its formal charge: {net_charge}"
            )
    if isotopic:
        from ._emc_isotopes import apply_isotope_masses

        data, params = apply_isotope_masses(data, params, psmiles, dp)
    version = re.match(r"LAMMPS output created by (EMC [^\n]+)", data)
    if not version:
        raise ForceFieldError("Missing EMC version in LAMMPS output")
    return EMCAssignment(
        "oplsaa" if family == "opls" else family,
        psmiles,
        dp,
        seed,
        data,
        params,
        version[1],
        sources,
        field,
        smiles,
        omitted,
        nchains,
        density,
        auto_depth,
    )
