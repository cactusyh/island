"""Replay published PolyJarvis successes and compare a RadonPy PSMILES panel."""

import argparse
import csv
import json
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from pathlib import Path

from island.forcefields.assignment import assign_polymer_parameters
from island.forcefields.emc import FIELDS
from island.forcefields.polyjarvis_fields import PATCHED


def cases(archives, panel, patched_fields, complete_fields=()):
    repeats = {}
    positives = set()
    for archive in archives:
        for line in archive.read_text().splitlines():
            row = json.loads(line)
            name, smiles = row["monomer_id"], row["smiles"]
            if repeats.setdefault(name, smiles) != smiles:
                raise ValueError(f"Inconsistent source PSMILES for {name}")
            if row["built"] and row["field"] in FIELDS:
                positives.add((name, row["field"]))
    pairs = positives | {(name, field) for field in complete_fields for name in repeats}
    for name, field in sorted(pairs):
        yield name, repeats[name], field, patched_fields, (name, field) in positives
    if panel:
        with panel.open(newline="") as handle:
            for row in csv.DictReader(handle):
                for field in ("gaff", "gaff2"):
                    yield row["monomer_ID"], row["smiles"], field, None, None


def compare_radonpy(smiles, field, row):
    from radonpy.core.poly import polymerize_MolFromSmiles
    from radonpy.ff.gaff import GAFF
    from radonpy.ff.gaff2 import GAFF2
    from rdkit import Chem

    reference = polymerize_MolFromSmiles(smiles, n=3, terminal="[H]")
    engine = GAFF() if field == "gaff" else GAFF2()
    row["reference_assigned"] = reference is not None and engine.ff_assign(
        reference, charge="gasteiger"
    )
    mol = assign_polymer_parameters(
        smiles, field, dp=3, charge="gasteiger", atom_order="radonpy"
    ).mol
    row.update(assigned=True, output_smiles=Chem.MolToSmiles(Chem.RemoveHs(mol)))
    if not row["reference_assigned"]:
        return
    row["reference_smiles"] = Chem.MolToSmiles(Chem.RemoveHs(reference))
    row["atom_parameters_match"] = mol.GetNumAtoms() == reference.GetNumAtoms() and all(
        atom.GetProp("ff_type") == other.GetProp("ff_type")
        and all(
            abs(atom.GetDoubleProp(k) - other.GetDoubleProp(k)) < 1e-10
            for k in ("AtomicCharge", "ff_sigma", "ff_epsilon")
        )
        for atom, other in zip(reference.GetAtoms(), mol.GetAtoms(), strict=True)
    )
    row["bonded_parameters_match"] = (
        mol.GetNumBonds() == reference.GetNumBonds()
        and all(
            (
                bond.GetBeginAtomIdx(),
                bond.GetEndAtomIdx(),
                bond.GetProp("ff_type"),
                bond.GetDoubleProp("ff_k"),
                bond.GetDoubleProp("ff_r0"),
            )
            == (
                other.GetBeginAtomIdx(),
                other.GetEndAtomIdx(),
                other.GetProp("ff_type"),
                other.GetDoubleProp("ff_k"),
                other.GetDoubleProp("ff_r0"),
            )
            for bond, other in zip(reference.GetBonds(), mol.GetBonds(), strict=True)
        )
        and all(
            [t.to_dict() for t in getattr(reference, kind)]
            == [t.to_dict() for t in getattr(mol, kind)]
            for kind in ("angles", "dihedrals", "impropers")
        )
    )


def probe(case):
    name, smiles, field, patched_fields, published = case
    row = {
        "id": name,
        "field": field,
        "psmiles": smiles,
        "assigned": False,
        "published_positive": published,
    }
    try:
        if field in ("gaff", "gaff2"):
            compare_radonpy(smiles, field, row)
        else:
            options = {
                "dp": 4,
                "nchains": 2,
                "density": 0.9,
                "stereo_policy": "polyjarvis",
            }
            if field == "opls/2024/opls-aa":
                options["field_dir"] = patched_fields
            result = assign_polymer_parameters(smiles, field, **options)
            row.update(
                assigned=True,
                output_smiles=result.emc_smiles,
                bond_stereo_omitted=result.bond_stereo_omitted,
            )
    except Exception as exc:  # noqa: BLE001 - keep the provider diagnostic for each input
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--polyjarvis-archives", nargs="+", type=Path, default=[])
    parser.add_argument(
        "--complete-emc-fields",
        nargs="+",
        choices=sorted({*FIELDS.values(), "opls"}),
        default=[],
        help="also test these fields on every archived repeat (opls selects stock 2024 AA)",
    )
    parser.add_argument("--patched-fields", type=Path)
    parser.add_argument("--radonpy-panel", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if not args.polyjarvis_archives and not args.radonpy_panel:
        parser.error("Supply a PolyJarvis archive or a RadonPy panel")
    if args.complete_emc_fields and not args.polyjarvis_archives:
        parser.error("--complete-emc-fields requires --polyjarvis-archives")
    if args.polyjarvis_archives:
        if args.patched_fields is None:
            parser.error("--polyjarvis-archives requires --patched-fields")
        for name, digest in PATCHED.items():
            if (
                sha256(
                    (args.patched_fields / "opls/2024" / name).read_bytes()
                ).hexdigest()
                != digest
            ):
                parser.error(f"Expected PolyJarvis field hash for {name}")
    counts = Counter()
    failed_reference = 0
    columns = (
        "id",
        "field",
        "psmiles",
        "assigned",
        "published_positive",
        "reference_assigned",
        "atom_parameters_match",
        "bonded_parameters_match",
        "output_smiles",
        "reference_smiles",
        "bond_stereo_omitted",
        "error",
    )
    with (
        args.output.open("x", newline="") as handle,
        ProcessPoolExecutor(max_workers=args.workers) as pool,
    ):
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in pool.map(
            probe,
            cases(
                args.polyjarvis_archives,
                args.radonpy_panel,
                args.patched_fields,
                args.complete_emc_fields,
            ),
            chunksize=2,
        ):
            writer.writerow(row)
            handle.flush()
            ok = (
                row["assigned"]
                and not row.get("error")
                and row.get("atom_parameters_match", True)
                and row.get("bonded_parameters_match", True)
            )
            counts[row["field"], bool(ok)] += 1
            if not ok and row["published_positive"] is not False:
                failed_reference += 1
    if not counts:
        parser.error("No applicable force-field cases in the supplied inputs")
    for key, count in sorted(counts.items()):
        print(*key, count)
    return int(failed_reference > 0)


if __name__ == "__main__":
    raise SystemExit(main())
