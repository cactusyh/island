"""Independent byte/line/decimal audit; does not import ISLAND's source parser."""

import argparse
import json
import re
from collections import Counter, defaultdict
from decimal import Decimal
from hashlib import sha256
from pathlib import Path


def audit(directory, report_path, expected_sha256):
    raw = Path(report_path).read_bytes()
    assert sha256(raw).hexdigest() == expected_sha256
    p = json.loads(raw)
    sources = {}
    for name, pin in p["source_hashes"].items():
        raw = (Path(directory) / name).read_bytes()
        assert sha256(raw).hexdigest() == pin
        sources[name] = raw.decode().splitlines()
    citations = list(p["include_order"]) + list(p["active_commands"].values())
    citations += [r["evidence"] for r in p["masses"].values()]
    for row in p["selected_coefficient_rows"]:
        citations.extend([row["data_evidence"], row["selected_settings_evidence"]])
        values = row["selected_settings_evidence"]["raw_text"].split("#", 1)[0].split()
        start = 3 if row["family"] == "nonbonded" else 2
        assert values[start:] == row["original_values"]
        source = row["paper_evidence"]
        paper_values = [
            source["raw_text"].split()[i] for i in source["numeric_columns"]
        ]
        assert list(map(Decimal, paper_values)) == list(map(Decimal, values[start:]))
        v = list(map(Decimal, values[start:]))
        c = row["converted_values"]
        if row["family"] == "nonbonded":
            assert Decimal(c["epsilon_kJ_per_mol"]) == v[0] * Decimal("4.184")
            assert Decimal(c["sigma_nm"]) == v[1] / 10
        elif row["family"] == "bond":
            assert (
                Decimal(c["openmm_k_kJ_per_mol_nm2"])
                == 2 * v[0] * Decimal("4.184") * 100
            )
            assert Decimal(c["r0_nm"]) == v[1] / 10
        elif row["family"] == "angle":
            assert Decimal(c["openmm_k_kJ_per_mol_rad2"]) == 2 * v[0] * Decimal("4.184")
        else:
            assert list(map(Decimal, c["opls_K_kJ_per_mol"])) == [
                x * Decimal("4.184") for x in v
            ]
    for c in citations:
        assert c["source_sha256"] == p["source_hashes"][c["source_file"]]
        assert sources[c["source_file"]][c["line"] - 1] == c["raw_text"]
    init = sources["param/OPLS/opls_system.in.init"]
    commands = [
        line.split("#")[0].split() for line in init if line.split("#")[0].strip()
    ]
    assert ["units", "real"] in commands
    assert ["pair_style", "lj/charmm/coul/long", "9.0", "11.0"] in commands
    assert ["pair_style", "lj/cut/coul/long", "11.0", "11.0"] not in commands
    assert ["pair_modify", "mix", "geometric"] in commands
    assert ["special_bonds", "lj/coul", "0.0", "0.0", "0.5"] in commands
    data = "\n".join(sources["param/OPLS/opls_system.data"])
    atom_block = re.search(
        r"\nAtoms # full\n\n(.*?)\n\nVelocities\n", data, re.DOTALL
    ).group(1)
    atoms = [line.split() for line in atom_block.splitlines() if line.strip()]
    assert len(atoms) == 32448
    q_by_type, q_by_molecule = defaultdict(Counter), defaultdict(Decimal)
    for row in atoms:
        q_by_type[row[2]][row[3]] += 1
        q_by_molecule[row[1]] += Decimal(row[3])
    assert dict(q_by_type) == p["source_graph"]["type_charge_histograms"]
    assert sum(q_by_molecule.values()) == Decimal(
        p["source_graph"]["exact_serialized_total_charge"]
    )
    assert {str(v) for v in q_by_molecule.values()} == {"-4.50E-15"}
    assert len(q_by_molecule) == 64
    assert len(p["selected_coefficient_rows"]) == 17
    assert not p["paper_comparison"]["paper_Qeq_model_reproduced"]
    assert not p["paper_comparison"]["pair_style_matches"]
    assert p["execution_profiles"] == p["backend_executions"] == []
    return {
        "source_files_verified": len(sources),
        "raw_citations_verified": len(citations),
        "coefficient_rows_reconstructed": 17,
        "author_atoms_reconstructed": len(atoms),
        "author_molecules_reconstructed": 64,
        "exact_raw_charge_sum": str(sum(q_by_molecule.values())),
        "policy_and_unit_conversions_verified": True,
        "acceptance_gate": 1,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--sha256", required=True)
    args = parser.parse_args()
    print(json.dumps(audit(args.sources, args.audit, args.sha256)))
