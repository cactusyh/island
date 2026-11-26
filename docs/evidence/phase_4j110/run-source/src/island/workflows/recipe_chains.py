"""Deterministic per-chain recipes and independently embedded mixture packing."""

import shutil
import time
from copy import deepcopy
from dataclasses import replace

import numpy as np

from .recipe_config import positive


def weighted(rng, values, weights, count):
    if (
        not values
        or len(values) != len(weights)
        or any(
            type(p) not in (int, float) or not np.isfinite(p) or p <= 0 for p in weights
        )
    ):
        raise ValueError(
            "Distribution values and positive finite probabilities must align"
        )
    if not np.isclose(sum(weights), 1, atol=1e-10, rtol=0):
        raise ValueError("Distribution probabilities must sum to one")
    return rng.choice(values, size=count, p=weights).tolist()


def expand(recipe):
    """Resolve every DP and sequence before embedding; never modify caller input."""
    if recipe.components is None:
        return [
            (
                "polymer",
                replace(
                    recipe,
                    chains=1,
                    seed=(recipe.seed + 104729 * i - 1) % 900_000_000 + 1,
                ),
            )
            for i in range(recipe.chains)
        ]
    if not recipe.components:
        raise ValueError("Mixture requires at least one component")
    rows, names = [], set()
    rng = np.random.default_rng(recipe.seed)
    for component in recipe.components:
        data = deepcopy(component)
        name = data.pop("name")
        if not isinstance(name, str) or not name or name in names:
            raise ValueError("Component names must be unique nonempty strings")
        names.add(name)
        architecture = data.pop("architecture", "homopolymer")
        if architecture not in ("homopolymer", "block", "branched", "copolymer"):
            raise ValueError("Mixture components must be finite-chain architectures")
        polymer = data.pop("polymer")
        lengths = data.pop("lengths", None)
        distribution = data.pop("length_distribution", None)
        sequence_probabilities = data.pop("sequence_probabilities", None)
        count = data.pop("chains", len(lengths) if lengths is not None else 1)
        positive(count, "component chains", integer=True)
        if data:
            raise ValueError("Unknown component options: " + ", ".join(data))
        if lengths is not None and distribution is not None:
            raise ValueError("Choose explicit lengths or a length distribution")
        if lengths is not None or distribution is not None:
            if any(k in polymer for k in ("dp", "blocks", "sequence")):
                raise ValueError(
                    "Length distribution replaces polymer dp/blocks/sequence"
                )
            if distribution is not None:
                if set(distribution) != {"values", "probabilities"}:
                    raise ValueError(
                        "Length distribution requires values and probabilities"
                    )
                for dp in distribution["values"]:
                    positive(dp, "DP", integer=True)
                lengths = weighted(
                    rng, distribution["values"], distribution["probabilities"], count
                )
            if len(lengths) != count:
                raise ValueError("Explicit lengths must match component chain count")
            for dp in lengths:
                positive(dp, "DP", integer=True)
        else:
            lengths = [None] * count
        for dp in lengths:
            chain = deepcopy(polymer)
            if dp is not None:
                chain["dp"] = dp
            if sequence_probabilities is not None:
                if (
                    architecture != "copolymer"
                    or "blocks" in chain
                    or "sequence" in chain
                ):
                    raise ValueError(
                        "Sequence probabilities require a copolymer with dp"
                    )
                labels = list(chain["repeat_units"])
                if set(sequence_probabilities) != set(labels):
                    raise ValueError("Sequence probabilities must cover every repeat")
                chain["sequence"] = weighted(
                    rng,
                    labels,
                    [sequence_probabilities[k] for k in labels],
                    chain.pop("dp"),
                )
            seed = (recipe.seed + 104729 * len(rows) - 1) % 900_000_000 + 1
            rows.append(
                (
                    name,
                    replace(
                        recipe,
                        architecture=architecture,
                        components=None,
                        polymer=chain,
                        chains=1,
                        seed=seed,
                    ),
                )
            )
    return rows


def pack_chains(recipe, out, packmol):
    from rdkit import Chem

    from island.epoxy.dynamics import pack

    from .recipe import _build_chain, write_json

    recipes = expand(recipe)
    attempts = []
    started = time.monotonic()
    for attempt in range(4):
        trial = out / f"attempt-{attempt:03d}"
        trial.mkdir()
        chains, chain_records = [], []
        for i, (component, cfg) in enumerate(recipes):
            work = trial / f"chain-{i + 1:05d}"
            work.mkdir()
            seed = (cfg.seed + 15485863 * attempt - 1) % 900_000_000 + 1
            before = time.monotonic()
            mol = _build_chain(
                cfg, work, random_coordinates=attempt > 0, conformer_seed=seed
            )
            chains.append(mol)
            xyz = mol.GetConformer().GetPositions()
            masses = np.array([a.GetMass() for a in mol.GetAtoms()])
            center = np.average(xyz, axis=0, weights=masses)
            spec = cfg.polymer_spec()
            chain_records.append(
                {
                    "chain": i + 1,
                    "component": component,
                    "architecture": cfg.architecture,
                    "polymer": cfg.polymer,
                    "dp": spec.degree_of_polymerization,
                    "sequence": list(spec.sequence),
                    "conformer_seed": seed,
                    "atoms": mol.GetNumAtoms(),
                    "build_seconds": time.monotonic() - before,
                    "radius_of_gyration_A": float(
                        np.sqrt(
                            np.average(
                                np.sum((xyz - center) ** 2, axis=1), weights=masses
                            )
                        )
                    ),
                }
            )
        record = {
            "directory": trial.name,
            "conformer_seeds": [r["conformer_seed"] for r in chain_records],
            "random_coordinates": attempt > 0,
            "packing_seed": recipe.seed,
            "density_g_cm3": recipe.density_g_cm3,
        }
        try:
            xyz, box = pack(
                chains,
                [1] * len(chains),
                recipe.density_g_cm3,
                recipe.seed,
                packmol,
                trial,
            )
        except RuntimeError:
            record["status"] = "failed"
            attempts.append(record)
            write_json(out / "packing-attempts.json", attempts)
            text = (trial / "pack.log").read_text()
            if attempt == 3 or not any(
                s in text
                for s in (
                    "unable to put the molecules",
                    "ENDED WITHOUT PERFECT PACKING",
                )
            ):
                raise
        else:
            record["status"] = "completed"
            attempts.append(record)
            write_json(out / "packing-attempts.json", attempts)
            for name in ("pack.in", "pack.log", "packed.xyz"):
                shutil.copyfile(trial / name, out / name)
            # Retain the historical first-chain artifact, plus all per-chain graphs.
            shutil.copyfile(trial / "chain-00001/chain.json", out / "chain.json")
            shutil.copyfile(trial / "monomer0.xyz", out / "monomer0.xyz")
            break
    mol, site_map, offset = None, [], 0
    for i, ((component, _), chain) in enumerate(zip(recipes, chains, strict=True)):
        mol = Chem.Mol(chain) if mol is None else Chem.CombineMols(mol, chain)
        site_map.extend(
            {
                "site_id": offset + j + 1,
                "chain": i + 1,
                "component": component,
                "chain_site_id": j + 1,
            }
            for j in range(chain.GetNumAtoms())
        )
        offset += chain.GetNumAtoms()
    mol.RemoveAllConformers()
    (out / "molecule.json").write_text(Chem.MolsToJSON([mol]))
    np.savez(out / "coordinates.npz", xyz=xyz, box=box)
    write_json(out / "site-map.json", site_map)
    write_json(out / "chains.json", chain_records)
    write_json(
        out / "packing-statistics.json",
        {
            "chains": len(chains),
            "atoms": len(xyz),
            "seconds": time.monotonic() - started,
            "box_A": box.tolist(),
            "independent_conformers": True,
        },
    )
