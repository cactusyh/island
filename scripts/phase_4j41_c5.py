"""C5 graph/bundle validation and native-source boundary, without invented rows."""

import argparse
import json
from pathlib import Path

from island.builders import build_linear_polymer
from island.charge_references.records import pack, unpack
from island.chemistry import from_smiles
from island.exceptions import IslandError
from island.forcefields.external_pdms import (
    LIMITS,
    ExternalParameterBundle,
    bundle_template,
    load_external_pdms_bundle,
    requirements,
)
from island.forcefields.pcff.sip_variants import SiPSourceVariant
from island.workflows.bundle import system_data, system_from
from scripts.pcff_j41_c4_reference import reconstruct, trace_audit
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j41_c3 import sources
from scripts.phase_4j41_c4 import charge_vector, convert, mapping

BASE = "cc3794db61f84fa14bedfeca5388945526d43aab"
BRANCH = "codex/phase-4j41-c5-external-pdms-parameter-bundle"
PSMILES = "[*:1]O[Si](C)(C)[*:2]"
CASES = (
    ("siosi_oligomer", "[SiH3]O[SiH3]", None),
    ("hexamethyldisiloxane", "C[Si](C)(C)O[Si](C)(C)C", None),
    ("explicit_h_repeat", PSMILES, 1),
    ("linear_pdms_dp2", PSMILES, 2),
    ("end_group_pdms_dp4", PSMILES, 4),
)


def fixture(case):
    _, text, dp = case
    if dp is None:
        return from_smiles(text, random_seed=2026)
    return build_linear_polymer(
        text,
        dp=dp,
        random_seed=2026,
        add_hydrogens=True,
        head_end_group="[*]C",
        tail_end_group="[*]O",
    )


MUTATIONS = (
    "missing_torsion",
    "missing_cross_term",
    "missing_bond_bond",
    "missing_bond_angle",
    "native_missing_increments",
    "graph",
    "order",
    "atom_id",
    "source_hash",
    "element",
    "component",
    "mixed_source",
    "zero",
    "charge_normalization",
    "charge_vector",
    "units",
    "special_bonds",
    "mixing",
    "authority",
    "evidence",
    "checksum",
    "dependency",
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def mutated_record(original, mutation):
    import copy

    p = copy.deepcopy(original)
    stage = None
    if mutation.startswith("missing_"):
        family = {
            "missing_torsion": "torsion_3",
            "missing_cross_term": "end_bond-torsion_3",
            "missing_bond_bond": "bond-bond",
            "missing_bond_angle": "bond-angle",
        }[mutation]
        key = next(k for k, r in p["coverage"].items() if r["family"] == family)
        p["coverage"].pop(key)
        stage = "coverage_mismatch"
    elif mutation == "native_missing_increments":
        p["charge_origin"] = "external_source_bond_increments"
        stage = "charge_increments_missing"
    elif mutation == "graph":
        p["graph_identity"] = "0" * 64
        stage = "graph_identity_mismatch"
    elif mutation == "order":
        p["atom_order"].reverse()
        stage = "atom_order_mismatch"
    elif mutation == "atom_id":
        p["atoms"][0]["id"] = -99
        stage = "atom_id_mismatch"
    elif mutation in ("source_hash", "mixed_source", "evidence"):
        p["files"]["absent-source.json"] = {
            "sha256": "0" * 64,
            "family": "COMPASS" if mutation == "mixed_source" else p["source_family"],
            "role": "evidence",
            "provenance": "negative control: unavailable asset",
        }
        stage = (
            "mixed_source_family" if mutation == "mixed_source" else "asset_unavailable"
        )
    elif mutation == "element":
        p["atoms"][0]["element"] = "P"
        stage = "atom_element_mismatch"
    elif mutation == "component":
        p["graph"]["molecule_membership"][str(p["atom_order"][0])] = "stale"
        stage = "graph_mismatch"
    elif mutation == "zero":
        key, req = next(iter(p["coverage"].items()))
        p["interactions"][key] = {
            "row": "invented_zero",
            "types": [],
            "dependencies": req["dependencies"],
        }
        p["source_rows"]["invented_zero"] = {
            "family": req["family"],
            "types": [],
            "values": ["0"],
            "units": ["kcal/mol"],
            "provenance": "negative control only",
            "citations": [],
        }
        stage = "coefficient_evidence_missing"
    elif mutation == "charge_normalization":
        p["charge_transformations"] = ["normalize_to_zero"]
        stage = "charge_redistribution_forbidden"
    elif mutation == "charge_vector":
        p["charges"] = [
            {"id": p["atom_order"][0], "value": "0", "evidence_row": "absent"}
        ]
        stage = "charge_vector_evidence_mismatch"
    elif mutation == "units":
        p["units"]["energy"] = "eV"
        stage = "units_mismatch"
    elif mutation == "special_bonds":
        p["nonbonded_policy"]["special_bonds"]["coulomb"] = [0, 0, 0.5]
        stage = "special_bonds_mismatch"
    elif mutation == "mixing":
        p["nonbonded_policy"]["mixing"] = "lorentz_berthelot"
        stage = "lj_mixing_mismatch"
    elif mutation == "authority":
        p["authority"]["execution_authorized"] = True
        stage = "authority_mismatch"
    elif mutation == "dependency":
        k, req = next((k, r) for k, r in p["coverage"].items() if r["dependencies"])
        p["interactions"][k] = {"row": "absent", "types": [], "dependencies": []}
        stage = "equilibrium_dependencies_mismatch"
    elif mutation == "checksum":
        # Alter payload without recomputing the envelope checksum.
        raw = json.loads(pack(p))
        raw["payload"]["profile"] = "forged"
        return json.dumps(raw), "malformed_bundle"
    else:
        raise ValueError("Unknown mutation")
    return pack(p), stage


def negative_controls(system, root):
    p = bundle_template(system)
    result = []
    for mutation in MUTATIONS:
        raw, expected = mutated_record(p, mutation)
        b = ExternalParameterBundle(raw)
        outcome = b.validate(system, assets=root, expected_sha256=b.identity).payload
        stages = [e["stage"] for e in outcome["diagnostics"]]
        if outcome["status"] != "parameter_incomplete" or expected not in stages:
            raise ValueError(f"Mutation escaped its declared gate: {mutation}")
        result.append(
            {"mutation": mutation, "expected_stage": expected, "outcome": outcome}
        )
    return result


def execute(args):
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island.j41-c5-study.v1",
        "base_commit": BASE,
        "branch": BRANCH,
        "public_input": "PSMILES/final graph; no user CAR/MDF",
        "external_bundle_supplied": args.bundle is not None,
        "cases": [],
        "execution_profiles": [],
        "native_pcff_acceptance_gate": 1,
        "software_acceptance_gate": 1,
        "limitations": LIMITS,
    }
    variant = next(v for v in sources(args.distribution) if v.filename == "pcff.frc")
    folder = root / "native-pcff"
    folder.mkdir()
    (folder / variant.filename).write_bytes(variant.raw)
    for n, raw, _ in variant.companions:
        (folder / n).write_bytes(raw)
    write(folder / "source-provenance.json", variant.payload["provenance"])
    if args.oracle_tools:
        import shutil

        shutil.copytree(args.oracle_tools, root / "oracle/tools")
    report["native_source_sha256"] = variant.expected_sha256
    report["companion_hashes"] = variant.payload["companion_hashes"]
    for spec in CASES:
        name = spec[0]
        system = fixture(spec)
        case = root / name
        case.mkdir()
        write(case / "system.json", system_data(system))
        write(case / "requirements.json", requirements(system))
        (case / "incomplete-template.json").write_text(pack(bundle_template(system)))
        bundle_path = (
            args.bundle if args.bundle else case / "external-bundle-not-supplied.json"
        )
        _, validation = load_external_pdms_bundle(
            bundle_path,
            system,
            expected_sha256=args.bundle_sha256,
        )
        write(case / "external-validation.json", validation.payload)
        labels, typing = mapping(system, variant, "siloxane", folder)
        audit = reconstruct(system, variant, labels, typing)
        native_case = folder / name
        native_case.mkdir()
        write(native_case / "raw-reconstruction.json", audit)
        converter = {
            "status": "converter_unavailable",
            "reason": "--oracle-tools not supplied",
        }
        if labels is None:
            converter = {
                "status": "typing_blocked",
                "reason": "native final-graph typing rejected; converter not invoked",
            }
        elif charge_vector(audit) is None:
            converter = {
                "status": "parameter_incomplete",
                "reason": "missing native charge increment; converter not invoked",
            }
        if args.oracle_tools and labels and charge_vector(audit) is not None:
            converter = convert(
                native_case,
                system,
                labels,
                charge_vector(audit),
                variant,
                args.oracle_tools.resolve(),
            )
        write(native_case / "converter-result.json", converter)
        report["cases"].append(
            {
                "name": name,
                "graph_identity": requirements(system)["graph_identity"],
                "atom_count": system.number_of_sites,
                "external_status": validation.payload["status"],
                "native_status": audit["status"],
                "missing_native_couplings": len(audit["missing_couplings"])
                if labels
                else None,
                "converter": converter,
                "execution_authorized": False,
            }
        )
    write(root / "negative-controls.json", negative_controls(fixture(CASES[1]), root))
    report["status"] = (
        "asset_unavailable" if args.bundle is None else "parameter_incomplete"
    )
    write(root / "report.json", report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "native_pcff_acceptance_gate": 1,
                "software_acceptance_gate": 1,
            }
        )
    )
    return 1


def reconstruct_study(root):
    root = Path(root)
    report = read(root / "report.json")
    require(
        report["base_commit"] == BASE and report["limitations"] == LIMITS,
        "Changed base or readiness",
    )
    require(
        report["execution_profiles"] == []
        and report["native_pcff_acceptance_gate"] == 1,
        "Forged authority",
    )
    require(
        not report["external_bundle_supplied"]
        and report["software_acceptance_gate"] == 1
        and report["status"] == "asset_unavailable",
        "Unsupported acceptance result",
    )
    folder = root / "native-pcff"
    variant = SiPSourceVariant(
        "pcff.frc",
        (folder / "pcff.frc").read_bytes(),
        report["native_source_sha256"],
        pack(read(folder / "source-provenance.json")),
        tuple(
            (n, (folder / n).read_bytes(), h)
            for n, h in report["companion_hashes"].items()
        ),
    )
    for spec, item in zip(CASES, report["cases"], strict=True):
        name = spec[0]
        case = root / name
        system = system_from(read(case / "system.json"))
        require(
            read(case / "requirements.json") == requirements(system),
            "Graph/coverage changed",
        )
        require(
            unpack((case / "incomplete-template.json").read_text())
            == bundle_template(system),
            "Template changed",
        )
        _, validation = load_external_pdms_bundle(
            case / "external-bundle-not-supplied.json",
            system,
            expected_sha256=None,
        )
        require(
            read(case / "external-validation.json") == validation.payload,
            "External qualification changed",
        )
        labels, typing = mapping(system, variant, "siloxane", root / "native-pcff")
        expected = reconstruct(system, variant, labels, typing)
        native = root / "native-pcff" / name
        require(
            read(native / "raw-reconstruction.json") == expected,
            "Native rows/charges changed",
        )
        require(
            item["graph_identity"] == requirements(system)["graph_identity"],
            "Stale graph",
        )
        if item["converter"].get("full_converter_calls_checked"):
            from scripts.validate_phase_4j41_c4_receipt import replay_converter

            tools_root = root / "oracle"
            replay_converter(
                tools_root, native, variant, item["converter"]["exit_code"]
            )
            trace = read(native / "trace/command.json")["stderr"]
            require(
                trace_audit(trace, variant)
                == read(native / "full-converter-searches.json"),
                "Converter trace changed",
            )
    require(
        read(root / "negative-controls.json")
        == negative_controls(fixture(CASES[1]), root),
        "Mutations changed",
    )
    return {
        "reconstruction_valid": True,
        "native_pcff_acceptance_gate": 1,
        "cases": len(CASES),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distribution", type=Path)
    parser.add_argument("--oracle-tools", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--bundle-sha256")
    parser.add_argument("--reconstruct", type=Path)
    args = parser.parse_args()
    if args.reconstruct:
        print(json.dumps(reconstruct_study(args.reconstruct)))
        return 0
    if not args.distribution or not args.output:
        parser.error("--distribution and --output are required")
    if args.bundle and not args.bundle_sha256:
        parser.error("--bundle requires an out-of-band --bundle-sha256")
    try:
        return execute(args)
    except (ValueError, TypeError, KeyError, OSError, IslandError) as exc:
        # A failed audit is not a successful blocked study.
        print(
            json.dumps(
                {
                    "status": "parameter_incomplete",
                    "stage": "audit_execution_failure",
                    "reason": str(exc),
                }
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
