"""GAFF assignment on an owned copy of an ISLAND or RDKit molecule."""

from dataclasses import dataclass
from hashlib import sha256
from math import isfinite
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import TYPE_CHECKING

from island.core import MolecularSystem
from island.exceptions import ForceFieldError

if TYPE_CHECKING:
    from rdkit.Chem import Mol

FAMILIES = ("gaff", "gaff2", "gaff2_mod")
QM_CHARGES = ("RESP", "ESP", "Mulliken", "Lowdin")


@dataclass(frozen=True)
class RadonPyAssignment:
    """Native parameterized molecule. site_ids[index] identifies the input atom."""

    mol: "Mol"
    site_ids: tuple[int, ...]
    family: str
    version: str
    database_sha256: str
    charge_method: str | None
    lj14_scale: float
    coulomb14_scale: float


def assign_radonpy(
    system: "MolecularSystem | Mol",
    family: str,
    *,
    charge: str | None = None,
    use_mdl: bool = True,
    retry_mdl: bool = True,
    db_file: str | Path | None = None,
    charge_options: dict | None = None,
) -> RadonPyAssignment:
    """Assign the selected RadonPy parameter library, in native units.

    ISLAND inputs retain stable site IDs. RDKit inputs retain atom order,
    stereo, conformers and preassigned charge properties, with zero-based IDs.
    charge_options are forwarded to RadonPy's charge calculator.
    """
    from rdkit import Chem

    if isinstance(system, Chem.Mol):
        mol = Chem.Mol(system)
        site_ids = tuple(range(mol.GetNumAtoms()))
    else:
        from island.chemistry import system_to_rdkit_graph, to_rdkit
        from island.chemistry._rdkit_stereo import enforce_stored_cip_labels

        if charge in QM_CHARGES or len(system.coordinates):
            graph = to_rdkit(system)
        else:
            graph = system_to_rdkit_graph(system)
            enforce_stored_cip_labels(graph.mol, system, graph.site_id_to_rdkit_index)
        mol = graph.mol
        site_ids = tuple(
            graph.rdkit_index_to_site_id[i] for i in range(mol.GetNumAtoms())
        )
        if not mol.GetNumAtoms() or any(
            not atom.GetAtomicNum() or atom.GetTotalNumHs() for atom in mol.GetAtoms()
        ):
            raise ForceFieldError(
                "RadonPy requires a capped graph with explicit H sites"
            )

    return _assign(
        mol,
        site_ids,
        family,
        charge=charge,
        use_mdl=use_mdl,
        retry_mdl=retry_mdl,
        db_file=db_file,
        charge_options=charge_options,
    )


def _assign(
    mol,
    site_ids,
    family,
    *,
    charge=None,
    use_mdl=True,
    retry_mdl=True,
    db_file=None,
    charge_options=None,
):
    if family not in FAMILIES:
        raise ValueError(f"Unknown RadonPy family {family!r}")
    if charge not in (None, "gasteiger", "zero", *QM_CHARGES):
        raise ValueError(f"Unknown RadonPy charge method {charge!r}")
    if charge is None and charge_options:
        raise ValueError("charge_options requires a charge method")
    try:
        import radonpy
        from radonpy.core import calc
        from radonpy.ff import gaff, gaff2, gaff2_mod
    except ImportError as exc:
        raise ImportError("RadonPy assignment requires radonpy-pypi") from exc

    database = (
        Path(db_file)
        if db_file is not None
        else Path(gaff.__file__).parent / "ff_dat" / f"{family}.json"
    )
    engine = {
        "gaff": gaff.GAFF,
        "gaff2": gaff2.GAFF2,
        "gaff2_mod": gaff2_mod.GAFF2_mod,
    }[family](db_file=str(database))
    if not engine.ff_assign(mol, useMDL=use_mdl, retryMDL=retry_mdl):
        raise ForceFieldError(
            f"RadonPy could not assign {family} to the supplied graph"
        )
    if charge is not None:
        with TemporaryDirectory(prefix="island-radonpy-") as work:
            options = {"work_dir": work, "tmp_dir": work, **(charge_options or {})}
            ok = calc.assign_charges(mol, charge=charge, **options)
        if not ok or any(
            not a.HasProp("AtomicCharge")
            or not isfinite(a.GetDoubleProp("AtomicCharge"))
            for a in mol.GetAtoms()
        ):
            raise ForceFieldError(f"RadonPy could not assign finite {charge} charges")
    return RadonPyAssignment(
        mol,
        site_ids,
        family,
        radonpy.__version__,
        sha256(database.read_bytes()).hexdigest(),
        charge,
        engine.param.lj_c14,
        engine.param.c_c14,
    )
