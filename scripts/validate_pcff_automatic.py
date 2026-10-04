"""Predeclared real-source graph typing and independent native-charge acceptance."""

import argparse
import importlib.util
import json
import sys
import time
import types
from collections import Counter
from decimal import Decimal
from hashlib import sha256
from pathlib import Path

from island.exceptions import IslandError
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    load_pcff_source,
    save_pcff_automatic_record,
    type_pcff_atoms,
)
from island.forcefields.pcff.automatic import PROFILE
from island.forcefields.pcff.source import PIN
from island.workflows.storage import json_bytes, publish

# Hand-authored expectations, frozen before execution. Patterns are a separate RDKit oracle.
PATTERNS = {
    "c": "[CX4H0]",
    "c1": "[CX4H1]",
    "c2": "[CX4H2]",
    "c3": "[CX4H3]",
    "hc": "[H][CX4]",
    "ho": "[H][OX2H1]",
    "oh": "[OX2H1]",
    "oc": "[OX2H0]",
}
CASES = [
    {"name": "ethane", "smiles": "CC", "counts": {"c3": 2, "hc": 6}},
    {"name": "propane", "smiles": "CCC", "counts": {"c3": 2, "c2": 1, "hc": 8}},
    {"name": "isobutane", "smiles": "CC(C)C", "counts": {"c3": 3, "c1": 1, "hc": 10}},
    {
        "name": "neopentane",
        "smiles": "CC(C)(C)C",
        "counts": {"c3": 4, "c": 1, "hc": 12},
    },
    {
        "name": "ethanol",
        "smiles": "CCO",
        "counts": {"c3": 1, "c2": 1, "oh": 1, "hc": 5, "ho": 1},
    },
    {"name": "dimethyl_ether", "smiles": "COC", "counts": {"c3": 2, "oc": 1, "hc": 6}},
    {
        "name": "PE_DP1",
        "psmiles": "[*:1]CC[*:2]",
        "dp": 1,
        "counts": {"c3": 2, "hc": 6},
    },
    {
        "name": "PE_DP3",
        "psmiles": "[*:1]CC[*:2]",
        "dp": 3,
        "counts": {"c3": 2, "c2": 4, "hc": 14},
    },
    {
        "name": "PEO_DP1",
        "psmiles": "[*:1]CCO[*:2]",
        "dp": 1,
        "counts": {"c3": 1, "c2": 1, "oh": 1, "hc": 5, "ho": 1},
    },
    {
        "name": "PEO_DP3",
        "psmiles": "[*:1]CCO[*:2]",
        "dp": 3,
        "counts": {"c3": 1, "c2": 5, "oh": 1, "oc": 2, "hc": 13, "ho": 1},
    },
    {
        "name": "PE_DP50",
        "psmiles": "[*:1]CC[*:2]",
        "dp": 50,
        "counts": {"c3": 2, "c2": 98, "hc": 202},
    },
]
UNSUPPORTED = {
    "ring": "C1CCCCC1",
    "water": "O",
    "peroxide": "COOC",
    "acetal": "COC(OC)C",
    "unsaturated": "CC=O",
}
ROWS = {
    493: ["1.0", "1", "c", "c", "0.0000", "0.0000"],
    505: ["1.0", "1", "c", "h", "-0.0530", "0.0530"],
    519: ["2.1", "8", "c", "o", "0.1330", "-0.1330"],
    810: ["1.0", "1", "h*", "o", "0.4241", "-0.4241"],
}
EQ = {
    "c": (209, "c"),
    "c1": (212, "c"),
    "c2": (213, "c"),
    "c3": (214, "c"),
    "hc": (244, "h"),
    "ho": (249, "h*"),
    "oh": (298, "o"),
    "oc": (296, "o"),
}
TOL = 1e-12


def build(case):
    if "smiles" in case:
        from island.chemistry import from_smiles

        return from_smiles(case["smiles"], add_hydrogens=True, random_seed=2026)
    from island.builders import build_linear_polymer

    return build_linear_polymer(case["psmiles"], dp=case["dp"], generate_3d=False)


def expected_labels(system):
    from rdkit import Chem

    from island.chemistry.rdkit_graph import system_to_rdkit_graph

    converted = system_to_rdkit_graph(system)
    selected = {}
    for label, pattern in PATTERNS.items():
        for match in converted.mol.GetSubstructMatches(Chem.MolFromSmarts(pattern)):
            sid = converted.rdkit_index_to_site_id[match[0]]
            if sid in selected:
                raise ValueError("Ambiguous independent expectation")
            selected[sid] = label
    if set(selected) != set(system.topology.sites):
        raise ValueError("Incomplete independent expected labels")
    return selected


def reference_charges(raw, system, labels):
    # Fixed audited source rows; no production parsing, selection or summation helper.
    lines = raw.decode().splitlines()
    table = {}
    for line, expected in ROWS.items():
        fields = lines[line - 1].split()
        if fields != expected:
            raise ValueError(f"Changed independent source row {line}")
        table[tuple(fields[2:4])] = [Decimal(x) for x in fields[4:]]
    for label, (line, equiv) in EQ.items():
        fields = lines[line - 1].split()
        if fields[2] != label or fields[4] != equiv:
            raise ValueError("Changed equivalence evidence")
    values = {sid: Decimal(0) for sid in labels}
    for bond in system.topology.bonds.values():
        a, b = bond.site1, bond.site2
        pair = (EQ[labels[a]][1], EQ[labels[b]][1])
        inc = table[pair] if pair in table else table[pair[::-1]][::-1]
        values[a] += inc[0]
        values[b] += inc[1]
    return {k: float(v) for k, v in values.items()}


def external_module(directory):
    """Explicitly execute hash-pinned external GPL code, never bundled or installed."""
    directory = Path(directory)
    for name, expected in PROFILE["external_audit"]["files"].items():
        if sha256((directory / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"External source hash mismatch: {name}")
    for package in ["src", "src.atom_typing", "src.atom_typing.typing"]:
        module = types.ModuleType(package)
        module.__path__ = []
        sys.modules[package] = module
    for name, filename in [
        ("src.atom_typing.typing.typing_functions", "typing_functions.py"),
        ("lunar_pcff_pinned", "PCFF.py"),
    ]:
        spec = importlib.util.spec_from_file_location(name, directory / filename)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return module


def external_types(module, system, prefix):
    # Independent BFS neighborhoods for the declared acyclic supported graph cases.
    adjacency = {sid: set() for sid in system.topology.sites}
    for bond in system.topology.bonds.values():
        adjacency[bond.site1].add(bond.site2)
        adjacency[bond.site2].add(bond.site1)
    count = Counter(s.element for s in system.topology.sites.values())
    formula = "-".join(f"{e}{n}" for e, n in sorted(count.items()))
    atoms = {}
    for sid, site in system.topology.sites.items():
        seen = {sid}
        front = {sid}
        ids = {}
        info = {}
        for depth in range(1, 4):
            front = set().union(*(adjacency[n] for n in front)) - seen
            seen |= front
            ids[depth] = sorted(front)
            info[depth] = [
                [system.topology.sites[n].element, 0, len(adjacency[n])]
                for n in sorted(front)
            ]
        atoms[sid] = types.SimpleNamespace(
            element=site.element,
            nb=len(adjacency[sid]),
            ring=0,
            rings=[],
            molecule=types.SimpleNamespace(formula=formula),
            neighbor_ids=ids,
            neighbor_info=info,
        )
    mm = module.nta(types.SimpleNamespace(atoms=atoms), str(prefix), "PCFF")
    return {
        "assignments": {sid: a.nta_type for sid, a in mm.atoms.items()},
        "per_site_info": {sid: a.nta_info for sid, a in mm.atoms.items()},
        "tally": mm.tally,
        "assumed_table": mm.assumed,
        "supported_types": mm.supported_types,
        "graph_adapter": "independent_acyclic_BFS_v1; ring=0 only for declared acyclic cases",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--lunar-dir",
        help="Explicit local directory with pinned PCFF.py and typing_functions.py",
    )
    args = parser.parse_args(argv)
    root = Path(args.output)
    try:
        root.mkdir()
    except OSError as error:
        print(str(error), file=sys.stderr)
        return 1
    declaration = {
        "schema": "island_pcff_automatic_acceptance_v1",
        "source": PIN,
        "profile": PROFILE,
        "cases": CASES,
        "unsupported": UNSUPPORTED,
        "expected_patterns": PATTERNS,
        "reference_rows": ROWS,
        "reference_equivalences": EQ,
        "charge_tolerance_e": TOL,
        "construction": "SMILES ETKDG seed 2026; PSMILES existing builder 2D, no QM",
        "external_requested": bool(args.lunar_dir),
        "external_expected_differences": "ho -> ho2 only; preserve every upstream label/info/tally",
        "gates": [
            "expected_labels",
            "native_charges",
            "unsupported_rejections",
            "external_comparison_if_requested",
        ],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(root / "declaration.json", json_bytes(declaration))
    outcomes = []
    try:
        source = load_pcff_source(args.source)
        module = external_module(args.lunar_dir) if args.lunar_dir else None
        for case in CASES:
            try:
                system = build(case)
                before = system.to_dict()
                start = time.perf_counter()
                typing = type_pcff_atoms(system, source)
                duration = time.perf_counter() - start
                expected = expected_labels(system)
                if (
                    dict(Counter(expected.values())) != case["counts"]
                    or typing.assignments != expected
                ):
                    raise ValueError("Expected type labels/counts not met")
                charges = assign_automatic_pcff_charges(system, typing)
                reference = reference_charges(source.raw, system, expected)
                actual = charges.charges
                error = max(abs(actual[s] - q) for s, q in reference.items())
                if (
                    error > TOL
                    or abs(charges.payload["native_charge_record"]["total_charge"])
                    > TOL
                ):
                    raise ValueError("Charge gate failed")
                if system.to_dict() != before:
                    raise ValueError("Input mutated")
                comparison = {"status": "unavailable"}
                if module:
                    external = external_types(module, system, root / case["name"])
                    publish(
                        root / f"{case['name']}-external.json", json_bytes(external)
                    )
                    differences = [
                        {
                            "site": sid,
                            "island": expected[sid],
                            "lunar": external["assignments"][sid],
                        }
                        for sid in expected
                        if expected[sid] != external["assignments"][sid]
                    ]
                    comparison = {
                        "status": "executed",
                        "differences": differences,
                        "tally": external["tally"],
                    }
                    if any(
                        d["island"] != "ho" or d["lunar"] != "ho2" for d in differences
                    ):
                        raise ValueError(
                            "Undeclared external disagreement; see retained upstream record"
                        )
                save_pcff_automatic_record(typing, root / f"{case['name']}-typing.json")
                save_pcff_automatic_record(
                    charges, root / f"{case['name']}-charges.json"
                )
                outcomes.append(
                    {
                        "case": case["name"],
                        "status": "passed",
                        "sites": len(expected),
                        "counts": case["counts"],
                        "typing_seconds": duration,
                        "graph_identity": typing.payload["graph_identity"],
                        "typing_identity": typing.identity,
                        "charge_identity": charges.identity,
                        "total_charge_e": charges.payload["native_charge_record"][
                            "total_charge"
                        ],
                        "maximum_charge_error_e": error,
                        "external": comparison,
                    }
                )
            except (
                IslandError,
                ValueError,
                OSError,
                ImportError,
                AttributeError,
                KeyError,
            ) as error:
                outcomes.append(
                    {"case": case["name"], "status": "failed", "error": str(error)}
                )
        for name, smiles in UNSUPPORTED.items():
            try:
                typing = type_pcff_atoms(build({"smiles": smiles}), source)
                save_pcff_automatic_record(typing, root / f"{name}-diagnostic.json")
                outcomes.append(
                    {
                        "case": name,
                        "status": "expected_rejection"
                        if not typing.complete
                        else "failed",
                        "diagnostics": typing.payload["diagnostics"],
                    }
                )
            except (IslandError, ValueError, OSError, ImportError) as error:
                outcomes.append({"case": name, "status": "failed", "error": str(error)})
    except (IslandError, ValueError, OSError, ImportError) as error:
        outcomes.append({"case": "setup", "status": "failed", "error": str(error)})
    passed = len(outcomes) == len(CASES) + len(UNSUPPORTED) and all(
        o["status"] in ("passed", "expected_rejection") for o in outcomes
    )
    report = {
        "schema": "island_pcff_automatic_acceptance_report_v1",
        "passed": passed,
        "outcomes": outcomes,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(root / "report.json", json_bytes(report))
    print(json.dumps(report, indent=2))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
