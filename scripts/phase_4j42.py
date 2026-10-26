"""Execute/reconstruct the explicit external-PDMS preparation boundary.

Exit 1 is classified blocked acceptance, not a Python exception. No parameter
path or SHA-256 is discovered. With no supplied asset, no numerical backend runs.
"""

import argparse
import json
from pathlib import Path

from island.charge_references.records import pack
from island.forcefields import (
    ExternalPDMSOptions,
    ForceFieldRequest,
    prepare_forcefield,
)
from island.forcefields.external_pdms import LIMITS, bundle_template
from island.forcefields.external_pdms_adapter import FAMILY
from island.workflows.bundle import system_data, system_from
from scripts.phase_4j41_c2 import read, write
from scripts.phase_4j41_c5 import CASES, fixture

BASE = "3825c838246ab2de75bb0c0a494c457c3153fee9"
BRANCH = "codex/phase-4j42-external-pdms-execution"


def prepare(system, path=None, pin=None):
    return prepare_forcefield(
        system, ForceFieldRequest(FAMILY, ExternalPDMSOptions(path, pin))
    )


def execute(output, *, bundle=None, pin=None, case=None):
    root = Path(output)
    root.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island.j42-study.v1",
        "base_commit": BASE,
        "branch": BRANCH,
        "external_bundle_supplied": bundle is not None,
        "independent_sha256_supplied": pin is not None,
        "execution_profiles": [],
        "backend_executions": [],
        "native_pcff_acceptance_gate": 1,
        "acceptance_gate": 1,
        "limitations": LIMITS,
        "cases": [],
    }
    for spec in CASES:
        if case is not None and case != spec[0]:
            continue
        system = fixture(spec)
        result = prepare(system, bundle, pin)
        folder = root / spec[0]
        folder.mkdir()
        write(folder / "system.json", system_data(system))
        (folder / "preparation.json").write_text(result.json_text)
        # Deliberately empty contract control; never numerical PDMS evidence.
        raw = pack(bundle_template(system))
        path = folder / "incomplete-template.json"
        path.write_text(raw)
        from hashlib import sha256

        control = prepare(system, path, sha256(raw.encode()).hexdigest())
        (folder / "incomplete-preparation.json").write_text(control.json_text)
        p = result.payload
        report["cases"].append(
            {
                "name": spec[0],
                "atoms": system.number_of_sites,
                "graph_identity": p["requirements"]["graph_identity"],
                "status": p["status"],
                "required_interactions": len(p["requirements"]["interactions"]),
                "missing_interactions": len(p["missing_interactions"]),
                "unresolved_dependencies": len(p["unresolved_dependencies"]),
                "preparation_identity": result.identity,
                "incomplete_control_identity": control.identity,
            }
        )
    write(root / "report.json", report)
    return report


def reconstruct(root):
    """Fresh-process reconstruction of absent-asset and incomplete-template runs.

    A supplied bundle is intentionally not rediscovered from an old report.
    Reexecution with an asset requires the caller's explicit path and new pin.
    """
    from hashlib import sha256

    from island.forcefields import ExternalPDMSPreparation

    root = Path(root)
    report = read(root / "report.json")
    if report["external_bundle_supplied"] or report["independent_sha256_supplied"]:
        raise ValueError(
            "Reexecution requires the explicit external path and independent pin"
        )
    for row in report["cases"]:
        folder = root / row["name"]
        system = system_from(read(folder / "system.json"))
        result = prepare(system)
        stored = ExternalPDMSPreparation((folder / "preparation.json").read_text())
        stored.validate_integrity(system)
        if stored != result or stored.identity != row["preparation_identity"]:
            raise ValueError("Stored preparation differs from reconstruction")
        path = folder / "incomplete-template.json"
        control = prepare(system, path, sha256(path.read_bytes()).hexdigest())
        if (
            control.json_text != (folder / "incomplete-preparation.json").read_text()
            or control.identity != row["incomplete_control_identity"]
        ):
            raise ValueError("Incomplete-template reconstruction changed")
    return {
        "reconstructed_cases": len(report["cases"]),
        "complete_execution_profiles": [],
        "acceptance_gate": 1,
        "limitations": LIMITS,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bundle", type=Path)
    parser.add_argument("--bundle-sha256")
    parser.add_argument("--case", choices=[c[0] for c in CASES])
    parser.add_argument("--reconstruct", type=Path)
    args = parser.parse_args()
    if args.reconstruct:
        if args.bundle or args.bundle_sha256 or args.output or args.case:
            parser.error("--reconstruct accepts only its explicit study directory")
        print(json.dumps(reconstruct(args.reconstruct)))
        return 0
    if not args.output:
        parser.error("--output is required")
    print(
        json.dumps(
            execute(
                args.output, bundle=args.bundle, pin=args.bundle_sha256, case=args.case
            )
        )
    )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
