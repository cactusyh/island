"""Reproduce unified selection using explicit local assets; no MD or discovery."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields.preparation import OPLSOptions, PCFFOptions
from island.forcefields.selection import capability_registry, prepare_polymer_forcefield
from island.forcefields.selection.inventory import digest
from island.forcefields.selection.registry import PDMS


def signed(data):
    return {**data, "identity": digest(data)}


def execute(
    output, *, pcff_source=None, opls_source=None, pdms_archive=None, pdms_sha256=None
):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    cases = {}
    sources = {}
    assets = {"pcff": pcff_source, "opls": opls_source, "pdms": pdms_archive}
    for family, path in assets.items():
        if path is not None:
            sources[family] = {
                "filename": Path(path).name,
                "sha256": sha256(Path(path).read_bytes()).hexdigest(),
            }
    options = {
        "pcff": None
        if pcff_source is None
        else PCFFOptions(pcff_source, lj=(0, 0, 1), coulomb=(0, 0, 1)),
        "opls": None if opls_source is None else OPLSOptions(opls_source),
    }
    for label, smiles in [
        ("pcff-control", "C"),
        ("aliphatic", "CCCC"),
        ("aromatic", "c1ccccc1"),
        ("peo", "COCCOC"),
        ("pdms", "C[Si](C)(C)O[Si](C)(C)C"),
        ("silicon", "C[SiH3]"),
        ("phosphorus", "OP(=O)(O)O"),
        ("multichain", "CC.CCO"),
        ("crosslinked", "C1CC2CCC1C2"),
    ]:
        for family in ("pcff", "opls", "gaff", "gaff2"):
            s = from_smiles(smiles, random_seed=46)
            result = prepare_polymer_forcefield(
                s,
                family,
                charge_method="am1bcc" if family.startswith("gaff") else "native",
                options=options.get(family),
            )
            cases[f"{family}-{label}"] = result.payload
    # Exercise actual PSMILES construction, preserving final end-group identity.
    for repeat, label in [
        ("[*:1]CC[*:2]", "polyethylene"),
        ("[*:1]c1ccc([*:2])cc1", "aromatic-polymer"),
        ("[*:1]COC[*:2]", "peo-polymer"),
    ]:
        s = build_linear_polymer(
            repeat,
            dp=2,
            head_end_group="[*]C",
            tail_end_group="[*]C",
            add_hydrogens=True,
            random_seed=46,
        )
        cases[label] = prepare_polymer_forcefield(
            s, "opls", charge_method="native", options=options["opls"]
        ).payload
    if pdms_archive is not None and pdms_sha256 is not None:
        from island.forcefields.pdms_opls import PDMSOPLSOptions

        for dp in (0, 2, 4):
            if dp:
                s = build_linear_polymer(
                    "[*:1]O[Si](C)(C)[*:2]",
                    dp=dp,
                    head_end_group="[*][Si](C)(C)C",
                    tail_end_group="[*]C",
                    add_hydrogens=True,
                    random_seed=46,
                )
                anchors = (s.metadata["polymer"]["head_site_id"],)
            else:
                s = from_smiles("C[Si](C)(C)O[Si](C)(C)C", random_seed=46)
                anchors = ()
            from island import SimulationBox

            s.box = SimulationBox(40, 40, 40)
            cases[f"external-pdms-dp{dp}"] = prepare_polymer_forcefield(
                s,
                "opls",
                source_profile=PDMS,
                charge_method="author_dft_npa",
                options=PDMSOPLSOptions(
                    pdms_archive, pdms_sha256, e_oxygen_site_ids=anchors
                ),
            ).payload
    report = signed(
        {
            "schema": "island.j46-audit.v1",
            "source_assets": sources,
            "cases": cases,
            "registry": capability_registry(),
            "numerical_md_executed": False,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    (output / "audit.json").write_text(
        json.dumps(report, sort_keys=True, indent=2) + "\n"
    )
    summary = signed(
        {
            "schema": "island.j46-summary.v1",
            "audit_identity": report["identity"],
            "cases": {
                name: {
                    k: p[k]
                    for k in (
                        "status",
                        "profile_identity",
                        "source_kind",
                        "stages",
                        "diagnostics",
                        "blocked_terms",
                        "identity",
                    )
                }
                for name, p in cases.items()
            },
            "source_assets": sources,
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n"
    )
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pcff-source", type=Path)
    parser.add_argument("--opls-source", type=Path)
    parser.add_argument("--pdms-archive", type=Path)
    parser.add_argument("--pdms-sha256")
    args = parser.parse_args(argv)
    report = execute(
        args.output,
        pcff_source=args.pcff_source,
        opls_source=args.opls_source,
        pdms_archive=args.pdms_archive,
        pdms_sha256=args.pdms_sha256,
    )
    print(
        json.dumps(
            {
                "audit_identity": report["identity"],
                "cases": len(report["cases"]),
                "numerical_md_executed": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
