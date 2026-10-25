"""Validate portable C1 source assessments, blocked gates and historical bytes."""

import argparse
import json
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from island.forcefields.pcff.mixed_sip import (
    load_pcff_mixed_sip_assessment,
    load_pcff_sip_distribution,
    prepare_pcff_mixed_sip,
)
from island.workflows import storage
from island.workflows.bundle import system_from
from scripts.pcff_j5_reference import match, read_source
from scripts.pcff_j41_c1_reference import mixed_inventory
from scripts.pcff_j41_reference import audit_queries
from scripts.validate_phase_4j41_receipt import require, sha, signed
from scripts.validate_phase_4j41_receipt import validate_study as validate_baseline

CASE_NAMES = {
    "silane_control",
    "methylsilane",
    "disiloxane",
    "hexamethyldisiloxane",
    "sila_oligomer",
    "si_n",
    "si_s",
    "si_f",
    "si_cl",
    "si_br",
    "si_i",
    "phosphate_ester",
    "phosphonate",
    "phosphine",
    "phosphine_oxide",
    "pf5",
    "mixed_sip",
}


def read(path):
    return storage.decode(storage.read_json(path))


def validate_study(root):
    root = Path(root)
    d = load_pcff_sip_distribution(root)
    report = read(root / "report.json")
    require(report["schema"] == "island.j41-c1-execution.v1", "Wrong study schema")
    require(report["distribution_hashes"] == d.hashes, "Changed distribution")
    require(
        report["acceptance_gate"] == 1
        and report["production_validated"] is False
        and report["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    require(
        report["ewald_status"] == "ewald_unvalidated"
        and report["pressure_status"] == "pressure_unavailable"
        and report["long_npt_status"] == "blocked",
        "Promoted periodic gates",
    )
    require(
        {r["name"] for r in report["cases"]} == CASE_NAMES
        and len(report["cases"]) == len(CASE_NAMES),
        "Changed fixtures",
    )
    require(
        read(root / "mixed-source-inventory.json") == mixed_inventory(d.source.raw),
        "Changed raw inventory",
    )
    baseline = validate_baseline(root / "baseline")
    require(baseline == report["baseline_validation"], "Changed silane replay")
    baseline_model = read(root / "baseline/report.json")["cases"][0]["model_identity"]
    results = []
    for row in report["cases"]:
        folder = root / row["name"]
        system = system_from(read(folder / "system.json"))
        a = load_pcff_mixed_sip_assessment(folder / "assessment.json", d, system=system)
        p = a.payload
        require(
            row["assessment_identity"] == a.identity
            and row["graph_identity"] == p["graph_identity"]
            and row["source_complete"] == p["complete"]
            and row["phosphorus_environments"] == p["phosphorus_environments"]
            and row["charge_origin"] == p["charge_origin"],
            "Changed source/graph assessment",
        )
        require(
            row["new_dynamics_executed"] is False and row["model_published"] is False,
            "Incomplete tranche executed",
        )
        if p["inspection"]:
            require(
                read(folder / "independent-audit.json")
                == audit_queries(system, a, d.source.raw),
                "Changed independent reconstruction",
            )
        if p["complete"]:
            require(row["name"] == "silane_control", "Unexpected coverage promotion")
            require(
                prepare_pcff_mixed_sip(system, a).identity
                == row["model_identity"]
                == baseline_model,
                "Changed finite silane model",
            )
            require(row["converter"]["exit_code"] == 0, "Failed retained converter")
        else:
            require(
                row["workflow_status"] == "blocked_before_model_publication"
                and not (folder / "model.json").exists(),
                "Published incomplete model",
            )
        results.append(
            {
                "name": row["name"],
                "identity": a.identity,
                "source_complete": p["complete"],
            }
        )
    audit = read(root / "compiled-source-audit.json")
    require(
        audit["raw_rows_checked"] == 4359
        and audit["compiled_queries"] == len(audit["queries"]) == 347
        and all(r["agrees"] for r in audit["queries"]),
        "Incomplete compiled source audit",
    )
    require(
        sha(root / "matcher-leaf.so") == report["matcher_sha256"],
        "Changed compiled matcher",
    )
    raw = read_source(d.source.raw)
    indexed = {r["id"]: r for rows in raw.values() for r in rows}
    for q in audit["queries"]:
        independent = match(
            raw, q["family"], q["namespace"], list(q["query"]), guarded_msi=True
        )
        require(
            independent == q["independent_raw_match"],
            "Changed raw compiled-query reconstruction",
        )
        require(
            independent["status"] == q["native"]["status"], "Changed matcher status"
        )
        if q["compiled_row"] is not None:
            row = indexed[q["compiled_row"]]
            require(
                list(map(str, row["values"])) == q["compiled_source_values"],
                "Changed compiled coefficient values",
            )
        if independent["status"] == "assigned":
            require(
                q["compiled_row"] in independent["selected_ids"],
                "Changed source precedence",
            )
            np.testing.assert_allclose(
                independent["values"],
                q["native"]["normalized_values"],
                atol=1e-10,
                rtol=1e-13,
            )
    return {
        "independent_reconstruction": True,
        "silane_model_identity": baseline_model,
        "cases": results,
        "raw_rows_checked": audit["raw_rows_checked"],
        "compiled_queries": audit["compiled_queries"],
        "acceptance_gate": 1,
    }


def validate(directory, *, repo=None):
    directory = Path(directory)
    repo = Path(repo) if repo else Path(__file__).resolve().parents[1]
    receipt = json.loads((directory / "receipt.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    require(
        receipt == signed(receipt) and receipt["schema"] == "island.j41-c1-receipt.v1",
        "Invalid receipt",
    )
    require(
        manifest == signed(manifest)
        and manifest["schema"] == "island.j41-c1-manifest.v1",
        "Invalid manifest",
    )
    require(
        receipt["manifest_sha256"] == sha(directory / "manifest.json"),
        "Changed manifest",
    )
    require(
        receipt["acceptance_gate"] == 1
        and receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    require(
        set(manifest["files"]) == {"execution.tar.xz", "validation.tar.xz"},
        "Missing evidence",
    )
    for name, checksum in manifest["files"].items():
        require(sha(directory / name) == checksum, "Changed evidence " + name)
    for group in ("historical_files", "implementation_files"):
        require(bool(manifest[group]), "Empty preservation inventory")
        for name, checksum in manifest[group].items():
            require(sha(repo / name) == checksum, "Changed " + group + ": " + name)
    expected = {
        "focused",
        "regressions",
        "full",
        "ruff",
        "compileall",
        "diff_check",
        "pip_supported",
        "pip_foyer",
        "pip_analysis",
        "pip_amber",
        "pip_quantum",
    }
    require(
        set(receipt["verification"]) == expected
        and all(v["exit_code"] == 0 for v in receipt["verification"].values()),
        "Required verification failed",
    )
    require(
        all(
            receipt["verification"][name]["junit"]["passed"] > 0
            and not receipt["verification"][name]["junit"]["failures"]
            and not receipt["verification"][name]["junit"]["errors"]
            for name in ("focused", "regressions", "full")
        ),
        "Invalid JUnit results",
    )
    with tempfile.TemporaryDirectory(prefix="island-j41-c1-") as temp:
        with tarfile.open(directory / "execution.tar.xz") as archive:
            for member in archive.getmembers():
                require(
                    (member.isfile() or member.isdir())
                    and not member.issym()
                    and not member.islnk()
                    and not Path(member.name).is_absolute()
                    and ".." not in Path(member.name).parts,
                    "Unsafe archive",
                )
            archive.extractall(temp, filter="data")
        root = Path(temp)
        bundle = json.loads((root / "bundle-manifest.json").read_text())
        require(bundle == signed(bundle), "Changed portable manifest")
        require(
            set(bundle["files"])
            == {
                str(p.relative_to(root))
                for p in root.rglob("*")
                if p.is_file() and p.name != "bundle-manifest.json"
            },
            "Incomplete raw inventory",
        )
        for name, checksum in bundle["files"].items():
            require(sha(root / name) == checksum, "Changed raw artifact")
        require(
            validate_study(root / "study") == receipt["independent_validation"],
            "Changed independent audit",
        )
    return {
        "receipt_valid": True,
        "historical_files": len(manifest["historical_files"]),
        "acceptance_gate": 1,
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("directory", type=Path)
    p.add_argument("--study", action="store_true")
    a = p.parse_args()
    print(
        json.dumps(
            validate_study(a.directory) if a.study else validate(a.directory), indent=2
        )
    )


if __name__ == "__main__":
    main()
