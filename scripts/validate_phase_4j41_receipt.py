"""Validate J41 portable evidence, raw reconstruction and historical byte hashes."""

import argparse
import hashlib
import json
import tarfile
import tempfile
from pathlib import Path

import numpy as np

from island.forcefields.pcff import (
    load_pcff_sip_assessment,
    load_pcff_source,
    prepare_pcff_sip,
)
from island.workflows import storage
from island.workflows.bundle import system_from
from scripts.pcff_j41_reference import audit_queries, inventory

CASES = {
    "silane": True,
    "methylsilane": False,
    "siloxane": False,
    "phosphate_ester": False,
    "phosphonate": False,
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def signed(record):
    record = {k: v for k, v in record.items() if k != "identity"}
    return {
        **record,
        "identity": hashlib.sha256(
            json.dumps(
                record, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode()
        ).hexdigest(),
    }


def validate_study(root):
    root = Path(root)
    report = storage.decode(storage.read_json(root / "report.json"))
    source = load_pcff_source(root / "pcff.frc")
    require(report["schema"] == "island.j41-execution.v1", "Wrong execution schema")
    require(
        {r["name"]: r["source_complete"] for r in report["cases"]} == CASES
        and len(report["cases"]) == len(CASES),
        "Changed fixture coverage",
    )
    require(
        report["acceptance_gate"] == 1
        and report["production_validated"] is False
        and report["simulation_readiness"] == "not_established",
        "Promoted readiness",
    )
    require(report["source_sha256"] == sha(root / "pcff.frc"), "Changed source")
    for name, checksum in report["distribution_hashes"].items():
        require(sha(root / name) == checksum, "Changed local distribution")
    require(
        storage.decode(storage.read_json(root / "source-inventory.json"))
        == inventory(source.raw),
        "Changed independent inventory",
    )
    results = []
    for row in report["cases"]:
        folder = root / row["name"]
        system = system_from(storage.decode(storage.read_json(folder / "system.json")))
        assessment = load_pcff_sip_assessment(
            folder / "assessment.json", source, system=system
        )
        data = assessment.payload
        require(
            row["assessment_identity"] == assessment.identity
            and row["graph_identity"] == data["graph_identity"]
            and row["diagnostics"] == data["diagnostics"]
            and row["charge_origin"] == data["charge_origin"],
            "Changed final-graph assignment",
        )
        if row["name"] in ("silane", "methylsilane", "siloxane"):
            independent = audit_queries(system, assessment, source.raw)
            require(
                independent
                == storage.decode(storage.read_json(folder / "independent-audit.json")),
                "Changed independent raw reconstruction",
            )
        if not row["source_complete"]:
            require(
                row["dynamics_executed"] is False
                and row["workflow_status"] == "blocked_before_model_publication",
                "Blocked source executed dynamics",
            )
            require(not (folder / "model.json").exists(), "Incomplete model published")
        else:
            model = prepare_pcff_sip(system, assessment)
            require(
                model.identity == row["model_identity"]
                and model.assignment.identity == row["parameter_identity"],
                "Changed model/parameters",
            )
            require(
                (folder / "model.json").read_text() == model.json_text
                and (folder / "parameters.json").read_text()
                == model.assignment.json_text,
                "Changed serialized model",
            )
            r = row["singlepoint"]
            require(
                r["energy_force_atol"] == 1e-5 and r["rtol"] == 2e-10,
                "Widened reference bounds",
            )
            for native, reference in (
                ("native_energies_kj_mol", "reference_energies_kj_mol"),
                ("native_forces_kj_mol_angstrom", "reference_forces_kj_mol_angstrom"),
            ):
                require(
                    np.isfinite(r[native]).all()
                    and np.isfinite(r[reference]).all()
                    and np.allclose(r[native], r[reference], atol=1e-5, rtol=2e-10),
                    "Failed numerical comparison",
                )
            require(
                row["minimization"]["converged"]
                and row["minimization"]["final_fmax"] <= 0.1,
                "Failed minimization",
            )
            require(
                len(row["finite_differences"]) >= 12
                and all(
                    r["passed"] and r["max_force_errors"][-1] <= 1e-5
                    for r in row["finite_differences"]
                ),
                "Failed finite differences",
            )
            require(
                row["dynamics_executed"] is True
                and {r["integrator"] for r in row["dynamics"]} == {"nve", "baoab"},
                "Missing finite dynamics",
            )
            for dynamics in row["dynamics"]:
                moved = root / ("silane-" + dynamics["integrator"])
                require(
                    dynamics["status"] == "passed"
                    and dynamics["rng_equal"]
                    and dynamics["fresh_process_checkpoint_byte_equal"]
                    and dynamics["uninterrupted_state_and_rng_exact"],
                    "Failed continuation",
                )
                require(
                    (moved / "in-process-checkpoint.json").read_bytes()
                    == (moved / "continued-checkpoint.json").read_bytes(),
                    "Checkpoint bytes differ",
                )
                require(
                    sha(moved / "continued-checkpoint.json")
                    == dynamics["final_checkpoint_sha256"],
                    "Changed checkpoint hash",
                )
        results.append(
            {
                "name": row["name"],
                "assessment_identity": assessment.identity,
                "source_complete": data["complete"],
            }
        )
    return {"raw_source_reconstruction": True, "cases": results}


def validate(directory, *, repo=None):
    directory = Path(directory)
    repo = Path(repo) if repo else Path(__file__).resolve().parents[1]
    receipt = json.loads((directory / "receipt.json").read_text())
    manifest = json.loads((directory / "manifest.json").read_text())
    require(
        receipt["schema"] == "island.j41-receipt.v1" and receipt == signed(receipt),
        "Invalid receipt",
    )
    require(
        manifest["schema"] == "island.j41-manifest.v1" and manifest == signed(manifest),
        "Invalid manifest",
    )
    require(
        receipt["production_validated"] is False
        and receipt["simulation_readiness"] == "not_established"
        and receipt["acceptance_gate"] == 1,
        "Promoted readiness",
    )
    require(
        sha(directory / "manifest.json") == receipt["manifest_sha256"],
        "Changed manifest",
    )
    require(
        set(manifest["files"]) == {"execution.tar.xz", "validation.tar.xz"},
        "Incomplete evidence",
    )
    for name, checksum in manifest["files"].items():
        require(sha(directory / name) == checksum, "Changed evidence: " + name)
    for group in ("historical_files", "implementation_files"):
        require(bool(manifest[group]), "Empty preservation inventory")
        for name, checksum in manifest[group].items():
            require(sha(repo / name) == checksum, "Changed " + group + ": " + name)
    with tempfile.TemporaryDirectory(prefix="island-j41-validate-") as temp:
        with tarfile.open(directory / "execution.tar.xz") as archive:
            for member in archive.getmembers():
                require(
                    not member.issym()
                    and not member.islnk()
                    and (member.isfile() or member.isdir())
                    and not Path(member.name).is_absolute()
                    and ".." not in Path(member.name).parts,
                    "Unsafe archive member",
                )
            archive.extractall(temp, filter="data")
        root = Path(temp)
        bundle = json.loads((root / "bundle-manifest.json").read_text())
        require(bundle == signed(bundle), "Changed bundle inventory")
        require(
            set(bundle["files"])
            == {
                str(p.relative_to(root))
                for p in root.rglob("*")
                if p.is_file() and p.name != "bundle-manifest.json"
            },
            "Incomplete portable inventory",
        )
        for name, checksum in bundle["files"].items():
            require(sha(root / name) == checksum, "Changed raw artifact")
        require(
            validate_study(root / "study") == receipt["independent_validation"],
            "Changed raw acceptance",
        )
    require(
        set(receipt["verification"])
        == {
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
        and all(v["exit_code"] == 0 for v in receipt["verification"].values()),
        "Failed verification",
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
