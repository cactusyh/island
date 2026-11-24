"""GAFF assignment on an owned copy of an ISLAND system using RadonPy."""

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


@dataclass(frozen=True)
class RadonPyAssignment:
    """Native parameterized molecule. site_ids[index] is the ISLAND site ID."""

    mol: "Mol"
    site_ids: tuple[int, ...]
    family: str
    version: str
    database_sha256: str
    charge_method: str | None
    lj14_scale: float
    coulomb14_scale: float


def assign_radonpy(
    system: MolecularSystem,
    family: str,
    *,
    charge: str | None = None,
    use_mdl: bool = True,
    retry_mdl: bool = True,
) -> RadonPyAssignment:
    """Run GAFF or GAFF2 with RadonPy's type rules, parameter tables and charges.

    Native coefficients retain RadonPy units (angstrom, kcal/mol, degrees).
    charge=None requests bonded and LJ assignment without partial charges.
    QM charge methods use the supplied conformer and RadonPy's optional QM stack.
    """
    if family not in ("gaff", "gaff2"):
        raise ValueError("RadonPy family must be gaff or gaff2")
    qm_methods = ("RESP", "ESP", "Mulliken", "Lowdin")
    if charge not in (None, "gasteiger", "zero", *qm_methods):
        raise ValueError(f"Unknown RadonPy charge method {charge!r}")

    try:
        import radonpy
        from radonpy.core import calc
        from radonpy.ff import gaff, gaff2
    except ImportError as exc:
        raise ImportError("RadonPy assignment requires radonpy-pypi") from exc

    from island.chemistry import system_to_rdkit_graph, to_rdkit
    from island.chemistry._rdkit_stereo import enforce_stored_cip_labels

    if charge in qm_methods or len(system.coordinates):
        graph = to_rdkit(system)
    else:
        graph = system_to_rdkit_graph(system)
        enforce_stored_cip_labels(graph.mol, system, graph.site_id_to_rdkit_index)
    mol = graph.mol
    if not mol.GetNumAtoms() or any(
        not atom.GetAtomicNum() or atom.GetTotalNumHs() for atom in mol.GetAtoms()
    ):
        raise ForceFieldError("RadonPy requires a capped graph with explicit H sites")

    database = Path(gaff.__file__).parent / "ff_dat" / f"{family}.json"
    engine_type = gaff.GAFF if family == "gaff" else gaff2.GAFF2
    engine = engine_type(db_file=str(database))
    if not engine.ff_assign(mol, useMDL=use_mdl, retryMDL=retry_mdl):
        raise ForceFieldError(
            f"RadonPy could not assign {family} to the supplied graph"
        )
    if charge is not None:
        with TemporaryDirectory(prefix="island-radonpy-") as work:
            ok = calc.assign_charges(mol, charge=charge, work_dir=work, tmp_dir=work)
        if not ok or any(
            not a.HasProp("AtomicCharge")
            or not isfinite(a.GetDoubleProp("AtomicCharge"))
            for a in mol.GetAtoms()
        ):
            raise ForceFieldError(f"RadonPy could not assign finite {charge} charges")

    return RadonPyAssignment(
        mol,
        tuple(graph.rdkit_index_to_site_id[i] for i in range(mol.GetNumAtoms())),
        family,
        radonpy.__version__,
        sha256(database.read_bytes()).hexdigest(),
        charge,
        engine.param.lj_c14,
        engine.param.c_c14,
    )
