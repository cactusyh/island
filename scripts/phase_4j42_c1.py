"""Reproduce classified external-PDMS boundaries without numerical execution."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack
from island.chemistry import from_smiles
from island.forcefields import (
    ExternalPDMSBlockedError,
    ExternalPDMSOptions,
    ForceFieldRequest,
    create_evaluator,
    prepare_forcefield,
)
from island.forcefields.external_pdms import LIMITS, bundle_template
from island.minimization import minimize_geometry
from island.workflows.bundle import system_data, system_from
from scripts.phase_4j41_c2 import read, write

BASE = "2c48a3eaa69e111279869ff6da12640225f78421"
BRANCH = "codex/phase-4j42-c1-blocked-diagnostics"
CASES = (
    "invalid_valence",
    "missing_sio",
    "charged",
    "radical",
    "phosphorus",
    "asset_unavailable",
    "missing_hash",
    "malformed_bundle",
    "incomplete_bundle",
)


def fixture(case):
    text = {
        "invalid_valence": "[SiH2]O",
        "missing_sio": "[SiH3][SiH3].O",
        "phosphorus": "P",
    }.get(case, "[SiH3]O[SiH3]")
    system = from_smiles(text, random_seed=2026)
    silicon = next(
        (a for a in system.topology.sites.values() if a.element == "Si"), None
    )
    if case == "charged":
        silicon.formal_charge = 1
    elif case == "radical":
        silicon.metadata["radical_electrons"] = 1
    return system


def prepare(system, folder, selection):
    path = None if selection["file"] is None else folder / selection["file"]
    return prepare_forcefield(
        system,
        ForceFieldRequest(
            "external_pdms_class2", ExternalPDMSOptions(path, selection["pin"])
        ),
    )


def boundaries(system, result):
    records = []
    for name in ("create_evaluator", "minimization/evaluator"):
        try:
            evaluator = create_evaluator(system, result)
            if name == "minimization/evaluator":
                minimize_geometry(system, evaluator)
        except ExternalPDMSBlockedError as error:
            records.append(
                {
                    "boundary": name,
                    "exception_type": type(error).__name__,
                    "status": error.status,
                    "result_identity": error.result.identity,
                    "diagnostics": error.diagnostics,
                    "validation_diagnostic": error.validation_diagnostic,
                }
            )
        else:
            raise ValueError("Blocked preparation reached a numerical backend")
    return records


def execute(output):
    root = Path(output)
    root.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island.j42-c1-study.v1",
        "base_commit": BASE,
        "branch": BRANCH,
        "cases": [],
        "software_acceptance_gate": 0,
        "scientific_acceptance_gate": 1,
        "execution_profiles": [],
        "backend_executions": [],
        "limitations": LIMITS,
    }
    for case in CASES:
        folder = root / case
        folder.mkdir()
        system = fixture(case)
        selection = {"file": None, "pin": None}
        if case in ("missing_hash", "malformed_bundle", "incomplete_bundle"):
            raw = (
                pack(bundle_template(system)).encode()
                if case == "incomplete_bundle"
                else b"{"
            )
            (folder / "candidate.json").write_bytes(raw)
            selection = {
                "file": "candidate.json",
                "pin": None if case == "missing_hash" else sha256(raw).hexdigest(),
            }
        result = prepare(system, folder, selection)
        write(folder / "system.json", system_data(system))
        write(folder / "selection.json", selection)
        (folder / "preparation.json").write_text(result.json_text)
        outcome = boundaries(system, result)
        write(folder / "boundaries.json", outcome)
        report["cases"].append(
            {
                "name": case,
                "status": result.payload["status"],
                "preparation_identity": result.identity,
                "diagnostics": result.payload["diagnostics"],
                "boundaries": outcome,
            }
        )
    write(root / "report.json", report)
    return report


def reconstruct(directory):
    root = Path(directory)
    report = read(root / "report.json")
    if [r["name"] for r in report["cases"]] != list(CASES):
        raise ValueError("Case inventory changed")
    for row in report["cases"]:
        folder = root / row["name"]
        system = system_from(read(folder / "system.json"))
        result = prepare(system, folder, read(folder / "selection.json"))
        if (
            result.json_text != (folder / "preparation.json").read_text()
            or result.identity != row["preparation_identity"]
        ):
            raise ValueError("Preparation reconstruction changed")
        outcome = boundaries(system, result)
        if outcome != row["boundaries"] or outcome != read(folder / "boundaries.json"):
            raise ValueError("Classified boundary changed")
    return {
        "reconstructed_cases": len(CASES),
        "classified_boundaries": 2 * len(CASES),
        "software_acceptance_gate": 0,
        "scientific_acceptance_gate": 1,
        "limitations": LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reconstruct", type=Path)
    args = parser.parse_args()
    if bool(args.output) == bool(args.reconstruct):
        parser.error("Specify exactly one of --output or --reconstruct")
    print(
        json.dumps(
            reconstruct(args.reconstruct) if args.reconstruct else execute(args.output)
        )
    )


if __name__ == "__main__":
    main()
