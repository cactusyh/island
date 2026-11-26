"""Checkpointed PCFF DGEBA/PACM cure in a fully periodic orthorhombic cell."""

import fcntl
import json
import os
import shutil
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from hashlib import sha256
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from rdkit import Chem

from .chemistry import AMINE, DGEBA, THIOL, epoxides, monomer, react
from .dynamics import pack, segment
from .network import bond_images, image, select, statistics


@dataclass(frozen=True)
class CureConfig:
    dgeba: int = 8
    pacm: int | None = None
    seed: int = 2026
    temperature: float = 450
    packing_density: float = 0.85
    cutoff: float = 5.0
    min_ring: int = 6
    probability: float = 1
    target: float = 1
    batch_size: int = 2
    pre_steps: int = 4000
    batch_steps: int = 800
    post_steps: int = 4000
    max_idle: int = 5

    petmp: int = 0
    donor_angle_tolerance: float | None = None
    epoxy_smiles: str = DGEBA
    hardener_smiles: str | None = None

    @property
    def template(self):
        return THIOL if self.petmp else AMINE

    @property
    def hardener_count(self):
        return self.petmp or self.pacm

    def __post_init__(self):
        if self.hardener_smiles is None:
            object.__setattr__(self, "hardener_smiles", self.template.hardener)
        for smi in (self.epoxy_smiles, self.hardener_smiles):
            m = Chem.MolFromSmiles(smi)
            if (
                m is None
                or len(Chem.GetMolFrags(m)) != 1
                or any(a.GetAtomicNum() == 0 for a in m.GetAtoms())
            ):
                raise ValueError(
                    "Reaction monomers require one capped, valid SMILES graph"
                )
        if not epoxides(Chem.MolFromSmiles(self.epoxy_smiles)):
            raise ValueError("Epoxy monomer has no supported terminal epoxide")
        hardener = Chem.AddHs(Chem.MolFromSmiles(self.hardener_smiles))
        if not sum(map(len, self.template.donors(hardener).values())):
            raise ValueError("Hardener has no supported donor hydrogen")
        if self.pacm is None:
            object.__setattr__(self, "pacm", 0 if self.petmp else 4)
        if any(type(v) != int or v < 0 for v in (self.pacm, self.petmp)) or bool(
            self.pacm
        ) == bool(self.petmp):
            raise ValueError(
                "Select exactly one positive hardener count, pacm or petmp"
            )
        if self.petmp and self.donor_angle_tolerance is None:
            object.__setattr__(self, "donor_angle_tolerance", 45.0)
        if self.donor_angle_tolerance is not None and (
            not self.petmp or not 0 < self.donor_angle_tolerance <= 180
        ):
            raise ValueError(
                "The PETMP donor angle tolerance must be in (0,180] degrees"
            )
        for v in (
            self.dgeba,
            self.hardener_count,
            self.seed,
            self.min_ring,
            self.batch_size,
            self.pre_steps,
            self.batch_steps,
            self.post_steps,
            self.max_idle,
        ):
            if type(v) != int or v < 1:
                raise ValueError("Counts, steps and seed must be positive integers")
        if (
            self.seed >= 2**31 - 1
            or not all(
                np.isfinite(x) and x > 0
                for x in (
                    self.temperature,
                    self.packing_density,
                    self.cutoff,
                    self.probability,
                    self.target,
                )
            )
            or max(self.target, self.probability) > 1
        ):
            raise ValueError("Invalid cure settings")


def mixture(cfg):
    monomers = [
        monomer(cfg.epoxy_smiles, cfg.seed),
        monomer(cfg.hardener_smiles, cfg.seed),
    ]
    mol = None
    k = 0
    for m, n in zip(monomers, [cfg.dgeba, cfg.hardener_count]):
        for _ in range(n):
            part = Chem.Mol(m)
            for a in part.GetAtoms():
                a.SetIntProp("monomer", k)
            mol = part if mol is None else Chem.CombineMols(mol, part)
            k += 1
    Chem.SanitizeMol(mol)
    return monomers, mol


@contextmanager
def stage(out, index):
    with TemporaryDirectory(dir=out, prefix=".cure-") as tmp:
        p = Path(tmp) / "stage"
        p.mkdir()
        yield p
        p.rename(out / f"stage-{index:04d}")


def _write_stage(work, record, xyz, box, vel):
    np.savez(work / "state.npz", xyz=xyz, box=box, vel=vel)
    record["state_sha256"] = sha256((work / "state.npz").read_bytes()).hexdigest()
    record["artifacts_sha256"] = {
        p.name: sha256(p.read_bytes()).hexdigest()
        for p in work.iterdir()
        if p.is_file()
    }
    (work / "record.json").write_text(
        json.dumps(record, indent=2, allow_nan=False) + "\n"
    )


def bound_cure_inputs(directory, library):
    """Load the actual reaction recipe and its parameter library from a cure."""
    manifest = json.loads((Path(directory) / "config.json").read_text())
    raw = Path(library).read_bytes()
    if sha256(raw).hexdigest() != manifest["library_sha256"]:
        raise ValueError("Parameter library differs from cure checkpoint")
    cfg = CureConfig(**manifest["config"])
    lib = json.loads(raw)
    if lib.get("chemistry", "pacm") != cfg.template.name:
        raise ValueError("Parameter library chemistry disagrees with hardener")
    return cfg, lib


def run_cure(
    directory,
    library,
    config=None,
    *,
    lammps="lmp",
    packmol="packmol",
    stop_after=None,
):
    """Resume at the last published stage. stop_after counts completed MD stages."""
    if stop_after is not None and (type(stop_after) is not int or stop_after < 0):
        raise ValueError("stop_after must be a nonnegative stage count")
    if config is None:
        config = CureConfig()
    out = Path(directory)
    out.mkdir(parents=True, exist_ok=True)
    with (out / ".lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        return _run(out, Path(library), config, lammps, packmol, stop_after)


def _run(out, library, cfg, lammps, packmol, stop_after):
    lammps = shutil.which(lammps)
    if not lammps:
        raise FileNotFoundError("LAMMPS executable missing")
    raw = library.read_bytes()
    lib = json.loads(raw)
    if lib.get("chemistry", "pacm") != cfg.template.name:
        raise ValueError("Parameter library chemistry disagrees with hardener")
    if lib.get("backend") != "radonpy-graph" and (
        cfg.epoxy_smiles != DGEBA or cfg.hardener_smiles != cfg.template.hardener
    ):
        raise ValueError("Custom monomers require the full-graph GAFF backend")
    if lib.get("backend") == "radonpy-graph":
        from .gaff import sources

        version, hashes = sources(lib["field"])
        if version != lib["engine_version"] or hashes != lib["field_sha256"]:
            raise ValueError("GAFF source/version differs from the network library")
        if any(lib[k] != getattr(cfg, k) for k in ("epoxy_smiles", "hardener_smiles")):
            raise ValueError("Reaction monomers differ from the network library")
    settings = asdict(cfg)
    if cfg.epoxy_smiles == DGEBA and cfg.hardener_smiles == cfg.template.hardener:
        settings.pop("epoxy_smiles")
        settings.pop("hardener_smiles")
    if not cfg.petmp:
        settings.pop("petmp")
        settings.pop("donor_angle_tolerance")
    manifest = {
        "protocol": 1,
        "config": settings,
        "library_sha256": sha256(raw).hexdigest(),
    }
    path = out / "config.json"
    if path.exists():
        if json.loads(path.read_text()) != manifest:
            raise ValueError("Cure settings or parameters differ from checkpoint")
    else:
        temp = out / ".config.tmp"
        temp.write_text(json.dumps(manifest, indent=2))
        os.replace(temp, path)
    monomers, mol = mixture(cfg)
    rng = np.random.default_rng(cfg.seed)
    stages = sorted(out.glob("stage-*"))
    records = []
    for i, d in enumerate(stages):
        if d.name != f"stage-{i:04d}":
            raise ValueError("Missing cure stage")
        rec = json.loads((d / "record.json").read_text())
        if sha256((d / "state.npz").read_bytes()).hexdigest() != rec["state_sha256"]:
            raise ValueError("Corrupt cure checkpoint")
        for name, expected in rec.get("artifacts_sha256", {}).items():
            if sha256((d / name).read_bytes()).hexdigest() != expected:
                raise ValueError("Corrupt cure stage artifact: " + name)
        mol = react(
            mol,
            [cfg.template.event(e) for e in rec["events"]],
            cfg.template,
        )
        records.append(rec)
    if stages:
        with np.load(stages[-1] / "state.npz") as saved:
            xyz, box, vel = (saved[k] for k in ("xyz", "box", "vel"))
        rng.bit_generator.state = records[-1]["rng"]
    initial_epoxides = len(epoxides(monomers[0])) * cfg.dgeba
    initial_donor_h = (
        sum(map(len, cfg.template.donors(monomers[1]).values())) * cfg.hardener_count
    )
    target = int(np.ceil(cfg.target * initial_epoxides))
    while not records or records[-1]["phase"] != "post":
        if stop_after is not None and len(records) >= stop_after:
            break
        index = len(records)
        with stage(out, index) as work:
            info = {"events": [], "candidate_count": 0, "rejected": {}}
            reason = None
            idle = 0
            if index == 0:
                exe = shutil.which(packmol)
                if not exe:
                    raise FileNotFoundError("PACKMOL executable missing")
                xyz, box = pack(
                    monomers,
                    [cfg.dgeba, cfg.hardener_count],
                    cfg.packing_density,
                    cfg.seed,
                    exe,
                    work,
                )
                vel = None
                phase, steps = "pre", cfg.pre_steps
            else:
                old = records[-1]
                formed = old["statistics"]["crosslinks"]
                if not epoxides(mol):
                    reason = "epoxide_exhausted"
                elif not any(cfg.template.donors(mol).values()):
                    reason = "thiol_h_exhausted" if cfg.petmp else "amine_h_exhausted"
                elif formed >= target:
                    reason = "target_conversion"
                elif old["idle"] >= cfg.max_idle:
                    probability_rejections = any(
                        r["rejected"].get("probability", 0)
                        for r in records[-cfg.max_idle :]
                    )
                    reason = (
                        "no_accepted_event"
                        if probability_rejections
                        else "no_eligible_contact"
                    )
                if reason:
                    phase, steps = "post", cfg.post_steps
                else:
                    phase, steps = "cure", cfg.batch_steps
                    events, info = select(
                        mol,
                        xyz,
                        box,
                        cfg.cutoff,
                        rng,
                        min(cfg.batch_size, target - formed),
                        cfg.min_ring,
                        cfg.probability,
                        cfg.template,
                        cfg.donor_angle_tolerance,
                    )
                    mol = react(mol, events, cfg.template)
                    idle = 0 if events else old["idle"] + 1
                    for e in events:
                        # Position the transferred H along the old N-H direction at
                        # the acceptor O to avoid an artificial stretched O-H bond.
                        dv = xyz[e.hydrogen] - xyz[e.donor]
                        dv += image(xyz[e.donor], xyz[e.hydrogen], box) * box
                        xyz[e.hydrogen] = (
                            xyz[e.oxygen] + 0.97 * dv / np.linalg.norm(dv)
                        ) % box
            xyz, box, vel, thermo = segment(
                work,
                mol,
                lib,
                xyz,
                box,
                vel,
                steps=steps,
                temp=cfg.temperature,
                seed=int(rng.integers(1, 900000000)),
                lammps=lammps,
                minimize=(not index or bool(info["events"])),
            )
            stat = statistics(
                mol,
                xyz,
                box,
                cfg.dgeba,
                cfg.hardener_count,
                cfg.template,
                initial_epoxides=initial_epoxides,
                initial_donor_h=initial_donor_h,
            )
            if stat["crosslinks"] != initial_epoxides - stat["remaining_epoxides"]:
                raise ValueError("Reaction and functional group counts disagree")
            rec = dict(
                phase=phase,
                index=index,
                **info,
                reason=reason,
                idle=idle,
                rng=rng.bit_generator.state,
                thermo=thermo,
                statistics=stat,
            )
            (work / "bond-images.json").write_text(
                json.dumps(bond_images(mol, xyz, box))
            )
            _write_stage(work, rec, xyz, box, vel)
        records.append(rec)
    return records
