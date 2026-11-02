"""Offline J70 replay, retaining original observation/source producers."""

import argparse
import lzma
from hashlib import sha256
from pathlib import Path

from island.pcff_completion import pin, require, signed
from island.pcff_qm_pilot.hmds_precision import read
from scripts.intake_phase_4j70 import ROOT, inputs, parse_native, reports
from scripts.replay_phase_4j68 import close


def static_expected(repo):
    """Reconstruct every declaration and metric from the pinned parent evidence."""
    old = repo / "docs/evidence/phase_4j69/source-to-evaluator-closure.json"
    values = read(old)
    return signed(
        {
            "schema": "island.j70-static-closure.v1",
            "producer": "J69 original source/evaluator numerical evidence; freshly replayed in J70, no invented QM producer",
            "producer_sha256": pin(old),
            "parent_reconstruction": read(
                repo / "docs/evidence/phase_4j69/receipt.json"
            )["numerical_reconstruction"],
            "available_subset_closed": True,
            "complete_source_static_closure": "blocked",
            "supplied_missing_rows": 0,
            "configurations": 26,
            "max_energy_difference_kcal_mol": max(
                v["energy_difference"]
                for r in values["configurations"]
                for v in r["families"].values()
            ),
            "max_force_component_difference_kcal_mol_A": max(
                v["force_component_difference"]
                for r in values["configurations"]
                for v in r["families"].values()
            ),
            "new_QM": 0,
            "MD_steps": 0,
            "full_bundle_comparison": "not_activated_missing_authoritative_rows",
        }
    )


def capability_expected(repo, model, queries):
    """Recheck the public facade and every exact missing interaction, offline."""
    from island.forcefield_policy import prepare_polymer_forcefield
    from island.workflows import storage
    from island.workflows.bundle import system_from

    system = system_from(
        storage.decode(
            storage.read_json(repo / "docs/evidence/phase_4j54/audit/hmds/system.json")[
                "payload"
            ]
        )
    )
    public = prepare_polymer_forcefield(system, forcefield="pcff").payload
    missing = []
    for r in sorted(model["interaction_inventory"], key=lambda r: r["id"]):
        if r["assignment"] is None:
            types = [model["typing"][str(i)] for i in r["sites"]]
            missing.append(
                {
                    **r,
                    "query": r["family"] + ":" + ",".join(types),
                    "source_lookup": queries[r["id"]],
                }
            )
    return signed(
        {
            "schema": "island.j70-canonical-capability.v1",
            "ordinary_api": 'forcefield="pcff"',
            "manual_profile_required": False,
            "manual_charge_method_required": False,
            "selected_profile": None,
            "automatic_assignment": False,
            "registered_bounded_profiles": [],
            "bounded_inspection_API": "island.pcff_pdms.prepare_polymer_forcefield",
            "ordinary_public_API": "island.forcefield_policy.prepare_polymer_forcefield",
            "ordinary_public_status": public["source_assignment"]["status"],
            "ordinary_public_diagnostics": public["diagnostics"],
            "bounded_status": "incomplete",
            "source_status": "parameter_incomplete",
            "missing_rows": missing,
            "ordinary_behavior": "unchanged; no HMDS source-complete domain registered",
            "graph_scope": "exact explicit-H neutral27 HMDS inspection; no DP2 or Si-H calculation",
            "cross_family_fallback": False,
            "QM_fit_authorized": False,
        }
    )


def replay(repo, root, deep=False):
    inputs(repo, root)
    for name, value in reports(repo, root).items():
        close(value, read(root / name), "J70." + name)
    review = read(root / "source-review.json")
    frozen = read(root / "source-inventory.json")
    require(
        [(r["path"], r["sha256"], r["category"]) for r in review["sources"]]
        == [(r["path"], r["sha256"], r["category"]) for r in frozen["sources"]],
        "J70_source_provenance_mutation",
    )
    require(
        review["originals_checked"] == 89
        and review["accepted_rows"] == 0
        and review["network_downloads"] == 0
        and all(r["eligible_native_overlay"] is False for r in review["sources"])
        and review["author_artifact_status"] == "author_artifact_not_supplied",
        "J70_source_declarations",
    )
    from island.pcff_qm_pilot.hmds_hamiltonian import load_model

    model = load_model(repo)
    expected_parsed = {
        (r["path"], r["sha256"], r["category"])
        for r in frozen["sources"]
        if Path(r["path"]).suffix == ".frc" and "pcff" in Path(r["path"]).name.lower()
    }
    require(
        len(review["parsed_PCFF_files"]) == len(expected_parsed)
        and {
            (r["path"], r["sha256"], r["category"]) for r in review["parsed_PCFF_files"]
        }
        == expected_parsed,
        "J70_parsed_source_inventory",
    )
    native_queries = None
    for record in review["parsed_PCFF_files"]:
        raw = lzma.decompress(
            (root / "source-bytes" / (record["sha256"] + ".frc.xz")).read_bytes()
        )
        require(sha256(raw).hexdigest() == record["sha256"], "J70_archived_source_SHA")
        result = parse_native(raw, model)
        close(result, record["result"], "J70.raw_parser")
        if record["sha256"] == model["source_hashes"]["pcff.frc"]:
            native_queries = result["queries"]
    converter = read(root / "converter/report.json")
    require(
        len(converter["attempts"]) == 2
        and converter["schema"] == "island.j70-converter-report.v1"
        and [r["attempt"] for r in converter["attempts"]] == [1, 2]
        and converter["successive_fresh_directory_returncodes"] == [17, 17]
        and converter["missing"] == "endbonddihedral hc c3 sio c3"
        and converter["complete"] is False
        and converter["deterministic_failure"] is True
        and converter["deterministic_complete_output"] == "not_established"
        and converter["full_emission"] == "blocked"
        and converter["new_QM"] == converter["MD_steps"] == 0,
        "J70_converter_status",
    )
    for r in converter["attempts"]:
        close(
            r,
            read(root / f"converter/attempt-{r['attempt']}/attempt.json"),
            "J70.attempt_entry",
        )
        require(
            r["source_sha256"] == model["source_hashes"]["pcff.frc"]
            and r["input_files"]
            == {
                n: pin(repo / "docs/evidence/phase_4j69/converter/input" / n)
                for n in ("system.car", "system.mdf", "system.json", "manifest.json")
            }
            and len(r["argv"]) == 9
            and Path(r["argv"][0]).name == "msi2lmp"
            and r["argv"][1:5] == ["system", "-class", "II", "-frc"]
            and Path(r["argv"][5]).name == "pcff.frc"
            and r["argv"][6:] == ["-print", "3", "-nocenter"],
            "J70_strict_converter_request",
        )
        require(
            r["returncode"] == 17
            and r["schema"] == "island.j70-converter-attempt.v1"
            and r["complete"] is False
            and r["failure"] is None
            and r["ignore_used"] is False
            and r["new_QM"] == r["MD_steps"] == 0,
            "J70_strict_converter_failure",
        )
        for kind in ("stdout", "stderr"):
            raw = lzma.decompress(
                (root / f"converter/attempt-{r['attempt']}/{kind}.log.xz").read_bytes()
            )
            require(
                sha256(raw).hexdigest() == r[kind + "_sha256"], "J70_converter_bytes"
            )
            if kind == "stdout":
                require(
                    b"Unable to find endbonddihedral data for hc c3 sio c3" in raw,
                    "J70_missing_converter_row",
                )
        require(
            r["executable_sha256"]
            == pin(repo / "docs/evidence/phase_4j69/converter/msi2lmp"),
            "J70_converter_producer",
        )
    static = read(root / "static-closure.json")
    close(static_expected(repo), static, "J70.static_report_semantics")
    require(
        static["producer_sha256"]
        == pin(repo / "docs/evidence/phase_4j69/source-to-evaluator-closure.json"),
        "J70_original_static_producer",
    )
    require(
        static["complete_source_static_closure"] == "blocked"
        and static["available_subset_closed"] is True
        and static["new_QM"] == static["MD_steps"] == 0,
        "J70_static_scope",
    )
    if deep:
        from scripts.replay_phase_4j69 import replay as parent

        require(
            parent(repo, repo / "docs/evidence/phase_4j69")[
                "source_available_subset_closed"
            ],
            "J70_parent_static_reconstruction",
        )
    cap = read(root / "canonical-capability.json")
    require(native_queries is not None, "J70_native_queries_unavailable")
    close(
        capability_expected(repo, model, native_queries),
        cap,
        "J70.capability_report_semantics",
    )
    require(
        {r["id"] for r in cap["missing_rows"]}
        == {
            r["id"]
            for r in reports(repo, root)["completeness-matrix.json"]["instances"]
        },
        "J70_canonical_exact_instance_mapping",
    )
    require(
        cap["automatic_assignment"] is False
        and cap["selected_profile"] is None
        and len(cap["missing_rows"]) == 378
        and cap["ordinary_api"] == 'forcefield="pcff"',
        "J70_canonical_fail_closed",
    )
    return {
        "offline": True,
        "instances": 378,
        "observations_retained": 26,
        "source_complete": False,
        "converter_complete": False,
        "parameter_assignment_capability": False,
        "new_QM": 0,
        "MD_steps": 0,
        "production_fits": 0,
        "accepted_coefficients": 0,
        "validation_untouched": True,
    }


def prepare_reports(repo, root):
    from island.forcefield_policy import prepare_polymer_forcefield as ordinary
    from island.pcff_completion import write
    from island.pcff_pdms import PCFFPDMSAssets, prepare_polymer_forcefield
    from island.workflows import storage
    from island.workflows.bundle import system_from
    from scripts.replay_phase_4j69 import replay as parent

    result = parent(repo, repo / "docs/evidence/phase_4j69")
    old = repo / "docs/evidence/phase_4j69/source-to-evaluator-closure.json"
    values = read(old)
    write(
        root / "static-closure.json",
        signed(
            {
                "schema": "island.j70-static-closure.v1",
                "producer": "J69 original source/evaluator numerical evidence; freshly replayed in J70, no invented QM producer",
                "producer_sha256": pin(old),
                "parent_reconstruction": result,
                "available_subset_closed": True,
                "complete_source_static_closure": "blocked",
                "supplied_missing_rows": 0,
                "configurations": 26,
                "max_energy_difference_kcal_mol": max(
                    v["energy_difference"]
                    for r in values["configurations"]
                    for v in r["families"].values()
                ),
                "max_force_component_difference_kcal_mol_A": max(
                    v["force_component_difference"]
                    for r in values["configurations"]
                    for v in r["families"].values()
                ),
                "new_QM": 0,
                "MD_steps": 0,
                "full_bundle_comparison": "not_activated_missing_authoritative_rows",
            }
        ),
    )
    source = read(repo / "docs/evidence/phase_4j69/source-assets.json")["sources"][
        "tools/msi2lmp/frc_files/pcff.frc"
    ]["declared_external_path"]
    system = system_from(
        storage.decode(
            storage.read_json(repo / "docs/evidence/phase_4j54/audit/hmds/system.json")[
                "payload"
            ]
        )
    )
    bounded = prepare_polymer_forcefield(
        system, forcefield="pcff", assets=PCFFPDMSAssets(Path(source).parent)
    ).payload
    public = ordinary(system, forcefield="pcff").payload
    write(
        root / "canonical-capability.json",
        signed(
            {
                "schema": "island.j70-canonical-capability.v1",
                "ordinary_api": 'forcefield="pcff"',
                "manual_profile_required": False,
                "manual_charge_method_required": False,
                "selected_profile": None,
                "automatic_assignment": False,
                "registered_bounded_profiles": [],
                "bounded_inspection_API": "island.pcff_pdms.prepare_polymer_forcefield",
                "ordinary_public_API": "island.forcefield_policy.prepare_polymer_forcefield",
                "ordinary_public_status": public["source_assignment"]["status"],
                "ordinary_public_diagnostics": public["diagnostics"],
                "bounded_status": bounded["parameter_completeness"],
                "source_status": "parameter_incomplete",
                "missing_rows": bounded["missing_terms"],
                "ordinary_behavior": "unchanged; no HMDS source-complete domain registered",
                "graph_scope": "exact explicit-H neutral27 HMDS inspection; no DP2 or Si-H calculation",
                "cross_family_fallback": False,
                "QM_fit_authorized": False,
            }
        ),
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--root", type=Path, default=ROOT)
    p.add_argument("--deep-parent", action="store_true")
    p.add_argument("--prepare-reports", action="store_true")
    a = p.parse_args()
    repo = a.repo.resolve()
    root = a.root.resolve()
    if a.prepare_reports:
        prepare_reports(repo, root)
    print(replay(repo, root, a.deep_parent))
