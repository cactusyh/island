"""Independent stdlib-only reconstruction of canonical assignments and raw rows."""

import argparse
import json
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

from scripts.audit_phase_4j47_inputs import audit as raw_audit
from scripts.audit_phase_4j47_inputs import read, require


def identity(value):
    return sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def audit(directory, archive, archive_sha256):
    directory = Path(directory)
    independent = raw_audit(directory, archive, archive_sha256)
    registry = json.loads((directory / "canonical-registry.json").read_text())
    require(
        registry["identity"]
        == identity({k: v for k, v in registry.items() if k != "identity"}),
        "independent_registry_hash",
    )
    selected = next(
        r
        for r in registry["entries"]
        if r["profile_identity"] == "island_external_opls_pdms_xiang2025_dft_local_v1"
    )
    for dp in (0, 2, 4, 6):
        root = directory / f"dp{dp}"
        p = json.loads((root / "canonical-assignment.json").read_text())
        source = read(root / "source-preparation.json")
        system = read(root / "system.json")
        require(
            p["identity"] == identity({k: v for k, v in p.items() if k != "identity"}),
            "independent_canonical_assignment_hash",
        )
        require(
            p["family"] == "opls"
            and p["selected_profile"] == selected["profile_identity"]
            and p["source_kind"] == "external_literature"
            and p["registry_identity"] == registry["identity"]
            and p["source_identity"] == identity(selected["source_hashes"])
            and p["charge_method"] == "author_dft_npa"
            and p["source_evidence"] == source["source"],
            "independent_canonical_source",
        )
        graph = p["final_graph"]
        require(
            p["graph_identity"] == identity(graph), "independent_canonical_graph_hash"
        )
        require(
            graph["atom_order"] == [r["id"] for r in system["sites"]],
            "independent_canonical_order",
        )
        expected_sites = [
            {
                "id": a["id"],
                "element": a["element"],
                "atomic_number": a["atomic_number"],
                "mass": a["mass"],
                "formal_charge": a["formal_charge"],
                "radical_electrons": a["metadata"].get("radical_electrons", 0),
                "isotope": a["metadata"].get("isotope", 0),
                "molecule": a["metadata"].get("chain_id"),
            }
            for a in system["sites"]
        ]
        require(graph["sites"] == expected_sites, "independent_canonical_sites")
        require(
            graph["bonds"]
            == sorted(
                [
                    [b["site1"], b["site2"], b["order"], b["aromatic"]]
                    for b in system["bonds"]
                ]
            ),
            "independent_canonical_bonds",
        )
        require(
            p["typed_graph"]["atom_types"]
            == [
                {
                    k: a[k]
                    for k in ("site_id", "author_type_id", "fragment", "molecule_id")
                }
                for a in source["mapping"]
            ],
            "independent_canonical_typing",
        )
        charge = p["graph_charges"]
        exact = {a["site_id"]: a["charge_e"] for a in source["mapping"]}
        require(
            charge["vector"] == [[i, exact[i]] for i in graph["atom_order"]]
            and charge["unit"] == "elementary_charge"
            and charge["charge_generation_method"] == "author_dft_npa"
            and charge["normalization"] is False
            and charge["redistribution"] is False
            and charge["charge_origin"] == source["charge_policy"],
            "independent_canonical_charge",
        )
        require(
            charge["component_totals"]
            == [
                str(sum((Fraction(exact[i]) for i in c), Fraction()))
                for c in graph["components"]
            ],
            "independent_canonical_components",
        )
        require(
            p["assigned_terms"] == source["interactions"]
            and not p["missing_terms"]
            and p["production_validated"] is False
            and p["simulation_readiness"] == "not_established"
            and p["stages"]["lammps_MD_verified"] is False,
            "independent_canonical_authority",
        )
    result = {
        "schema": "island.j48-independent-inputs.v1",
        "raw_reconstruction": independent,
        "canonical_source_graph_type_charge_and_parameter_mapping_reconstructed": True,
        "registry_identity": registry["identity"],
        "cases": [0, 2, 4, 6],
    }
    return {**result, "identity": identity(result)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--archive", type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(argv)
    a.output.write_text(
        json.dumps(
            audit(a.directory, a.archive, a.archive_sha256), indent=2, sort_keys=True
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
