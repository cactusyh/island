"""Declared real Foyer fidelity acceptance, separate from synthetic unit tests."""

import argparse
import time
from collections import Counter
from copy import deepcopy
from dataclasses import replace
from decimal import Decimal
from math import fsum
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.oplsaa import (
    assign_native_charges,
    load_oplsaa_source,
    save_opls_result,
    type_atoms,
)
from island.forcefields.oplsaa.adapter import check_foyer_installation
from island.forcefields.oplsaa.source import PIN
from island.workflows import storage

MATRIX = [
    ("ethane", "CC", 0),
    ("butane", "CCCC", 0),
    ("ethanol", "CCO", 0),
    ("dimethyl_ether", "COC", 0),
    ("benzene", "c1ccccc1", 0),
    ("pe3", "[*:1]CC[*:2]", 3),
    ("pe50", "[*:1]CC[*:2]", 50),
    ("peo3", "[*:1]CCO[*:2]", 3),
    ("ps3", "[*:1]CC(c1ccccc1)[*:2]", 3),
]


def build(smiles, dp):
    if dp:
        from island.builders import build_linear_polymer

        return build_linear_polymer(smiles, dp=dp, generate_3d=False)
    from rdkit import Chem

    from island.chemistry.rdkit_adapter import from_rdkit

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    # Coordinates are deliberately irrelevant to this experiment. No embedding.
    mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))
    return from_rdkit(mol).system


def upstream(system, xml):
    """Independent Forcefield loader + ParmEd graph, not adapter provider/graph."""
    import parmed
    from foyer import Forcefield

    ff = Forcefield(forcefield_files=str(xml))
    graph = parmed.Structure()
    sites = list(reversed(list(system.topology.sites)))
    mapped = {}
    for sid in sites:
        a = system.topology.sites[sid]
        atom = parmed.Atom(name=a.element, atomic_number=a.atomic_number, mass=a.mass)
        graph.add_atom(atom, "SAME", 1)
        mapped[sid] = atom
    for bond in reversed(list(system.topology.bonds.values())):
        graph.bonds.append(
            parmed.Bond(mapped[bond.site1], mapped[bond.site2], order=bond.order)
        )
    result = ff.run_atomtyping(graph, use_residue_map=False)
    types = {sites[i]: row["atomtype"] for i, row in result.items()}
    charges = {sid: ff.get_parameters("atoms", t)["charge"] for sid, t in types.items()}
    return types, charges


def variant(system):
    from island import Coordinates, MolecularSystem, Topology

    remap = {
        sid: i * 13 + 41 for i, sid in enumerate(reversed(list(system.topology.sites)))
    }
    graph = Topology()
    for sid in reversed(list(system.topology.sites)):
        graph.add_site(replace(deepcopy(system.topology.sites[sid]), id=remap[sid]))
    for b in reversed(list(system.topology.bonds.values())):
        graph.add_bond(
            remap[b.site2], remap[b.site1], order=b.order, aromatic=b.aromatic
        )
    return MolecularSystem(
        graph, Coordinates({sid: (9.0, 7.0, -2.0) for sid in remap.values()})
    ), remap


def execute(root, xml, plan):
    source = load_oplsaa_source(xml)
    rows = []
    for name, smiles, dp in MATRIX:
        row = {"case": name, "status": "failed"}
        started = time.perf_counter()
        try:
            s = build(smiles, dp)
            before = deepcopy(s.to_dict())
            typing = type_atoms(s, source)
            save_opls_result(typing, s, source, root / f"{name}-typing.json")
            types, native = upstream(s, xml)
            assert typing.payload["data"]["assignments"] == types
            row["upstream_typing_agreement"] = True
            if name in ("ethane", "butane", "pe3", "pe50"):
                for sid, atom in s.topology.sites.items():
                    if atom.element == "C":
                        h_count = sum(
                            s.topology.sites[
                                b.site2 if b.site1 == sid else b.site1
                            ].element
                            == "H"
                            for b in s.topology.bonds.values()
                            if sid in (b.site1, b.site2)
                        )
                        assert types[sid] == (
                            "opls_135" if h_count == 3 else "opls_136"
                        )
            if name in ("ethanol", "dimethyl_ether"):
                assert {
                    types[sid]
                    for sid, a in s.topology.sites.items()
                    if a.element == "O"
                } == {"opls_154" if name == "ethanol" else "opls_180"}

            # Decimal extraction is independent of the adapter's charge lookup.
            entries = {
                e.attrib["type"]: Decimal(e.attrib["charge"])
                for e in ET.parse(xml).find("NonbondedForce")
            }
            row["xml_spot_checks"] = {
                t: str(entries[t]) for t in sorted(set(types.values()))
            }
            assert all(native[sid] == float(entries[t]) for sid, t in types.items())
            changed, remap = variant(s)
            changed_t = type_atoms(changed, source)
            assert {
                sid: changed_t.payload["data"]["assignments"][remap[sid]]
                for sid in remap
            } == types
            charges = assign_native_charges(s, typing, source)
            assert {
                sid: r["charge"] for sid, r in charges.payload["charges"].items()
            } == native
            charges.validate_integrity(s, source)
            assert s.to_dict() == before
            save_opls_result(charges, s, source, root / f"{name}-charges.json")
            row.update(
                status="passed",
                sites=s.number_of_sites,
                total_e=fsum(native.values()),
                typing_identity=typing.identity,
                charge_identity=charges.identity,
                type_counts=dict(Counter(types.values())),
                max_charge_difference_e=0.0,
                remapping_and_coordinate_independence=True,
            )
        except Exception as error:  # noqa: BLE001 -- preserve failures, no substitutions
            row["failure"] = f"{type(error).__name__}: {error}"
        row["wall_seconds"] = time.perf_counter() - started
        rows.append(row)
        print(name, row["status"], row.get("failure", ""), flush=True)
        storage.publish(
            root / "outcomes.json",
            storage.json_bytes(
                {
                    "schema": "island_opls_acceptance_v1",
                    "declaration": plan,
                    "cases": rows,
                    "complete": len(rows) == len(MATRIX)
                    and all(r["status"] == "passed" for r in rows),
                    "production_validated": False,
                    "simulation_readiness": "not_established",
                }
            ),
            replace=True,
        )
    return all(r["status"] == "passed" for r in rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--execute-declared", action="store_true")
    args = parser.parse_args()
    _, _, _, versions = check_foyer_installation()
    plan = {
        "schema": "island_opls_declaration_v1",
        "source": load_oplsaa_source(args.xml).identity,
        "pin": PIN,
        "dependencies": versions,
        "matrix": [list(r) for r in MATRIX],
        "checks": [
            "exact types and charges against direct Foyer",
            "alkane CH3 opls_135 / CH2 opls_136; alcohol O opls_154 / ether O opls_180",
            "XML decimal extraction",
            "noncontiguous IDs and insertion order",
            "coordinate independence",
            "nonmutation",
            "component and molecular neutrality",
        ],
        "neutrality_tolerance_e_per_site": 1e-12,
        "outer_retries": 0,
        "embedding": False,
        "QM": False,
        "energy_or_dynamics": False,
    }
    if not args.execute_declared:
        args.output.mkdir(parents=True, exist_ok=False)
        storage.publish(args.output / "declaration.json", storage.json_bytes(plan))
        return 0
    if storage.read_json(args.output / "declaration.json") != plan:
        raise ValueError("Declaration changed")
    if {p.name for p in args.output.iterdir()} != {"declaration.json"}:
        raise ValueError("Refuse existing execution output")
    return 0 if execute(args.output, args.xml, plan) else 1


if __name__ == "__main__":
    raise SystemExit(main())
