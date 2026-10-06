"""J6 diagnostic adjudication. Success is not complete-model acceptance."""

import argparse
import builtins
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

from pcff_j5_reference import charge_raw, read_source
from validate_pcff_operational_support import write

from island.charge_references.records import unpack
from island.forcefields.pcff import (
    assess_pcff_resolution,
    load_pcff_resolution_assessment,
    load_pcff_source,
    save_pcff_resolution_assessment,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.charges import identity
from island.workflows.bundle import system_from
from island.workflows.storage import decode


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def child(args):
    forbidden = {"openmm", "scipy", "rdkit", "parmed", "foyer"}
    if forbidden.intersection(sys.modules):
        raise RuntimeError(
            "Scientific dependency imported before offline reconstruction"
        )
    original = builtins.__import__

    def blocked(name, *a, **kw):
        if name.split(".")[0] in {"openmm", "scipy", "rdkit", "parmed", "foyer"}:
            raise RuntimeError("Forbidden scientific dependency: " + name)
        return original(name, *a, **kw)

    builtins.__import__ = blocked
    manifest = decode(json.loads((args.child / "manifest.json").read_text()))
    source = load_pcff_source(args.source)
    results = []
    for row in manifest:
        record_path = args.child / row["record"]
        system_path = args.child / row["system"]
        if (
            digest(record_path) != row["record_sha256"]
            or digest(system_path) != row["system_sha256"]
        ):
            raise ValueError("Changed diagnostic inputs")
        system = system_from(decode(json.loads(system_path.read_text())))
        record = load_pcff_resolution_assessment(record_path, source, system=system)
        found = identity(unpack(record.json_text))
        if found != row["identity"]:
            raise ValueError("Changed diagnostic identity")
        results.append(
            {"case": row["case"], "identity": found, "bytes_unchanged": True}
        )
    write(
        args.child.parent / "child-reconstruction.json",
        {
            "results": results,
            "scientific_imports_blocked": True,
            "diagnostic_records_only": True,
        },
    )
    return 0


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--declaration", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--child", type=Path)
    p.add_argument("--require-full-source", action="store_true")
    args = p.parse_args()
    if args.child:
        return child(args)
    if args.declaration is None or args.output is None:
        p.error("--declaration and --output required")
    args.output.mkdir(parents=True, exist_ok=False)
    declaration = json.loads(args.declaration.read_text())
    write(args.output / "declaration.json", declaration)
    if digest(args.source) != declaration["source"]["sha256"]:
        raise ValueError("Wrong source pin")
    source = load_pcff_source(args.source)
    raw = read_source(args.source.read_bytes())
    report = {
        "cases": [],
        "operational_gate": False,
        "full_source_gate": False,
        "diagnostic_gate": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    staging = args.output / "staging"
    staging.mkdir()
    manifest = []
    for case in declaration["cases"]:
        try:
            path = Path(case["retained_input"]["path"])
            if digest(path) != case["retained_input"]["sha256"]:
                raise ValueError("Historical system hash mismatch")
            system = system_from(decode(json.loads(path.read_text())))
            typing = type_pcff_atoms(system, source, profile=declaration["profile"])
            record = assess_pcff_resolution(
                system,
                typing,
                resolution_policy=declaration["resolution_policy"],
                special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
            )
            payload = unpack(record.json_text)
            record_path = staging / (case["name"] + ".assessment.json")
            save_pcff_resolution_assessment(record, record_path)
            system_path = staging / (case["name"] + ".system.json")
            shutil.copyfile(path, system_path)
            reference = charge_raw(raw, system, typing.assignments)
            if reference["complete"] != payload["charges_complete"]:
                raise ValueError("Independent charge completeness mismatch")
            old = json.loads((path.parent / "outcome.json").read_text())
            if (
                typing.identity != old["typing_identity"]
                or payload["native_identities"]["charge_identity"]
                != old["charge_identity"]
            ):
                raise ValueError("Historical typing/charge identity changed")
            if (
                payload["native_identities"]["model_identity"] is not None
                and payload["native_identities"]["model_identity"]
                != old["model_identity"]
            ):
                raise ValueError("Historical model identity changed")
            result = {
                "case": case["name"],
                "identity": identity(payload),
                "native_identities": payload["native_identities"],
                "historical_identities_preserved": True,
                "charges_complete": payload["charges_complete"],
                "native_model_complete": payload["native_model_complete"],
                "counts": payload["counts"],
                "structural_blocker_count": payload["structural_blocker_count"],
                "charge_reference": reference,
                "blockers": [
                    e for e in payload["entries"] if e["structural_model_blocker"]
                ],
            }
            manifest.append(
                {
                    "case": case["name"],
                    "record": record_path.name,
                    "system": system_path.name,
                    "identity": result["identity"],
                    "record_sha256": digest(record_path),
                    "system_sha256": digest(system_path),
                }
            )
        except Exception as error:  # noqa: BLE001 -- retain diagnostic failure
            result = {"case": case["name"], "error": f"{type(error).__name__}: {error}"}
        report["cases"].append(result)
        print(case["name"], result.get("counts", result.get("error")), flush=True)
    write(staging / "manifest.json", manifest)
    relocated = args.output / "relocated"
    staging.rename(relocated)
    cmd = [
        sys.executable,
        str(Path(__file__).resolve()),
        "--source",
        str(args.source.resolve()),
        "--child",
        str(relocated.resolve()),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False, timeout=600)
    (args.output / "child.stdout").write_text(proc.stdout)
    (args.output / "child.stderr").write_text(proc.stderr)
    report["child"] = {
        "command": cmd,
        "returncode": proc.returncode,
        "relocated": True,
        "diagnostic_records_only": True,
    }
    report["diagnostic_gate"] = (
        proc.returncode == 0
        and len(manifest) == len(declaration["cases"])
        and len(manifest) > 0
    )
    write(args.output / "report.json", report)
    return 0 if report["diagnostic_gate"] and not args.require_full_source else 1


if __name__ == "__main__":
    raise SystemExit(main())
