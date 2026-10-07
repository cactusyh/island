"""J16 actual local converter/native boundaries. Default gate requires a new polymer.

No auxiliary patterns are executed. Audit/converter gate success is distinct
from executable or full-source coverage. All declared cases are mandatory.
"""

import argparse
import subprocess
import time
from collections import Counter
from pathlib import Path

import numpy as np
from pcff_distribution_contract import frozen as frozen_contract
from pcff_distribution_reference import labels_for, missing_request
from pcff_j5_reference import charge_raw, read_source
from validate_pcff_singlepoint import data_sections, reference_inputs
from validate_pcff_urethanes import assess, contract

from island.exceptions import PCFFError
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldError,
    prepare_forcefield,
)
from island.forcefields.pcff import (
    inspect_pcff_operational_support,
    load_pcff_source,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.distribution import (
    inspect_pcff_distribution,
    save_pcff_distribution_audit,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY, MSI_POLICY
from island.workflows import storage
from island.workflows.bundle import system_from


def native_case(s, c, source, path):
    t = type_pcff_atoms(s, source, profile="island_pcff_source_graph_v6")
    r = {
        "case": c["name"],
        "typing_complete": t.complete,
        "labels": dict(Counter(t.payload["assignments"].values())),
        "typing_diagnostics": t.payload["diagnostics"],
    }
    inspections = {}
    if t.complete:
        for policy in (MSI_POLICY, COMPATIBILITY_POLICY):
            inspections[policy] = inspect_pcff_operational_support(
                s,
                t,
                resolution_policy=policy,
                special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
            )
        v = inspections[COMPATIBILITY_POLICY]
        r.update(
            charges_complete=v["charges_complete"],
            model_complete=v["model_complete"],
            model_diagnostics=v["model_diagnostics"],
            strict_model_complete=inspections[MSI_POLICY]["model_complete"],
            independent_charge=charge_raw(read_source(source.raw), s, t.assignments),
        )
    try:
        prep = prepare_forcefield(
            s,
            ForceFieldRequest(
                "pcff",
                PCFFOptions(
                    path,
                    (0, 0, 1),
                    (0, 0, 1),
                    typing_profile="island_pcff_source_graph_v6",
                    resolution_policy=COMPATIBILITY_POLICY,
                ),
            ),
        )
    except (PCFFError, PreparedForceFieldError) as e:
        r["public_exception"] = {"type": type(e).__name__, "message": str(e)}
    else:
        r.update(
            public_preparation_complete=True,
            native_identity=prep.native_result.identity,
            prepared_identity=prep.identity,
        )
    return r, t, inspections


def artifacts(root):
    return {
        p.name: storage.checksum(p.read_bytes())
        for p in sorted(root.iterdir())
        if p.is_file()
    }


def check_case(result, c, expected):
    if result.get("execution_error"):
        return ["Unrelated execution error"]
    errors = []
    if c["role"] == "historical_positive":
        if not all(
            result.get(k) is True
            for k in (
                "typing_complete",
                "charges_complete",
                "model_complete",
                "public_preparation_complete",
                "independent_types_and_charges",
            )
        ):
            errors.append("Historical operational control regressed")
    else:
        errors.extend(assess(result, expected[c["name"]]))
    converter = result.get("converter", {})
    if result.get("charges_complete") is True:
        if converter.get("returncode") == 0:
            if not converter.get("charge_vector_preserved") or not converter.get(
                "output_sections"
            ):
                errors.append("Unverified converter output")
        elif type(converter.get("returncode")) is not int or not converter.get(
            "independent_missing"
        ):
            errors.append("Unattributed converter failure")
    elif converter != {
        "state": "blocked_by_native_typing_or_charge",
        "not_executed": True,
    }:
        errors.append("Invalid converter prerequisites")
    return errors


def assess_set(results, declaration, expected):
    names = [c["name"] for c in declaration["cases"]]
    if (
        not names
        or len(set(names)) != len(names)
        or Counter(r.get("case") for r in results) != Counter(names)
    ):
        return ["Missing, duplicate, empty or unexpected mandatory cases"]
    lookup = {c["name"]: c for c in declaration["cases"]}
    return [
        f"{r['case']}: {e}"
        for r in results
        for e in check_case(r, lookup[r["case"]], expected)
    ]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--declaration",
        type=Path,
        default=Path("docs/evidence/phase_4j16_declaration.json"),
    )
    for name in ("frc", "rlb", "templates", "converter", "converter-receipt", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument(
        "--gate",
        choices=("audit", "converter", "polymer", "full-source"),
        default="polymer",
    )
    p.add_argument("--require-full-source", action="store_true")
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=False)
    report = {
        "schema": "island_j16_boundary_acceptance_v1",
        "cases": [],
        "failures": [],
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        d = storage.read_json(a.declaration)
        frozen = frozen_contract("declaration")
        if d != frozen:
            raise ValueError("Changed experiment declaration")
        paths = {"frc_path": a.frc, "rlb_path": a.rlb, "templates_path": a.templates}
        hashes = {
            k: d["distribution"][n]["sha256"]
            for k, n in (
                ("frc", "pcff.frc"),
                ("rlb", "pcff.rlb"),
                ("templates", "pcff_templates.dat"),
            )
        }
        audit = inspect_pcff_distribution(**paths, expected_hashes=hashes)
        save_pcff_distribution_audit(a.output / "distribution.json", audit, **paths)
        build = storage.read_json(a.converter_receipt)
        if (
            build != frozen_contract("converter")
            or build["revision"] != d["local_checkout_revision"]
            or storage.checksum(a.converter.read_bytes()) != build["sha256"]
        ):
            raise ValueError("Converter build identity mismatch")
        storage.publish(a.output / "converter-build.json", storage.json_bytes(build))
        source = load_pcff_source(a.frc)
        raw = read_source(source.raw)
        expected = contract()["cases"]
        for c in d["cases"]:
            work = a.output / c["name"]
            work.mkdir()
            result = {"case": c["name"]}
            try:
                original = Path(c["system_path"])
                if storage.checksum(original.read_bytes()) != c["system_sha256"]:
                    raise ValueError("Retained final-graph checksum mismatch")
                s = system_from(storage.decode(storage.read_json(original)))
                before = s.to_dict()
                result, t, inspections = native_case(s, c, source, a.frc)
                storage.publish(
                    work / "native.json",
                    storage.json_bytes(storage.encode(inspections)),
                )
                if result.get("charges_complete"):
                    labels = labels_for(s, c)
                    q = charge_raw(raw, s, labels)
                    v = inspections[COMPATIBILITY_POLICY]["native_charge_record"][
                        "native_charge_record"
                    ]["partial_charges"]
                    if (
                        labels != t.assignments
                        or not q["complete"]
                        or set(v) != set(labels)
                    ):
                        raise ValueError(
                            "Independent reference labels/coverage mismatch"
                        )
                    np.testing.assert_allclose(
                        [v[i] for i in sorted(v)],
                        [float(q["partial_charges"][i]) for i in sorted(v)],
                        atol=1e-12,
                        rtol=0,
                    )
                    result["independent_types_and_charges"] = True
                    storage.publish(
                        work / "independent-charge.json",
                        storage.json_bytes(storage.encode(q)),
                    )
                    reference_inputs(
                        s,
                        labels,
                        {i: float(q["partial_charges"][i]) for i in labels},
                        np.array([s.coordinates.get(i) for i in sorted(labels)]),
                        work,
                    )
                    command = [
                        str(a.converter.resolve()),
                        "reference",
                        "-class",
                        "II",
                        "-frc",
                        str(a.frc.resolve()),
                        "-p",
                        "3",
                        "-nocenter",
                    ]
                    start = time.perf_counter()
                    proc = subprocess.run(
                        command, cwd=work, capture_output=True, timeout=180, check=False
                    )
                    (work / "stdout.txt").write_bytes(proc.stdout)
                    (work / "stderr.txt").write_bytes(proc.stderr)
                    storage.publish(work / "command.json", storage.json_bytes(command))
                    conversion = {
                        "returncode": proc.returncode,
                        "seconds": time.perf_counter() - start,
                    }
                    if proc.returncode:
                        conversion["independent_missing"] = missing_request(
                            raw, proc.returncode, (proc.stdout + proc.stderr).decode()
                        )
                    else:
                        sections = data_sections(work / "reference.data")
                        atomrows = sorted(
                            sections["Atoms"], key=lambda row: int(row[0])
                        )
                        np.testing.assert_allclose(
                            [float(row[3]) for row in atomrows],
                            [float(q["partial_charges"][i]) for i in sorted(labels)],
                            atol=1e-12,
                            rtol=0,
                        )
                        conversion.update(
                            charge_vector_preserved=True,
                            output_sections={k: len(v) for k, v in sections.items()},
                        )
                    result["converter"] = conversion
                else:
                    result["converter"] = {
                        "state": "blocked_by_native_typing_or_charge",
                        "not_executed": True,
                    }
                if (
                    before != s.to_dict()
                    or storage.checksum(original.read_bytes()) != c["system_sha256"]
                ):
                    raise ValueError("Input mutation")
            except Exception as error:  # noqa: BLE001 -- report and fail, never count arbitrary errors as chemistry
                result["execution_error"] = {
                    "type": type(error).__name__,
                    "message": str(error),
                }
            result["artifacts"] = artifacts(work)
            result["acceptance_errors"] = check_case(result, c, expected)
            report["cases"].append(result)
            print(
                c["name"],
                result.get("converter"),
                result["acceptance_errors"],
                flush=True,
            )
        report["failures"] = assess_set(report["cases"], d, expected)
        # Source bytes and authoritative typing/model records cannot change during a run.
        audit.validate_integrity(**paths)
        report.update(
            distribution_identity=audit.identity,
            source=source.identity,
            declaration_sha256=storage.checksum(a.declaration.read_bytes()),
            converter_sha256=build["sha256"],
        )
    except Exception as error:  # noqa: BLE001 -- preserve failed experiment
        report["failures"].append(type(error).__name__ + ": " + str(error))
    storage.publish(
        a.output / "report.json", storage.json_bytes(storage.encode(report))
    )
    # Reconstruct semantics and inspect real artifacts in one shared current gate.
    # Never promote converter exit status or user-supplied report booleans.
    from check_pcff_distribution import check

    try:
        checked = check(
            a.output, frc=a.frc, rlb=a.rlb, templates=a.templates, converter=a.converter
        )
    except Exception as error:  # noqa: BLE001 -- preserve rejected evidence
        checked = {"error": type(error).__name__ + ": " + str(error)}
    storage.publish(a.output / "checked.json", storage.json_bytes(checked))
    gate = "full-source" if a.require_full_source else a.gate
    key = {
        "audit": "audit_gate",
        "converter": "converter_boundary_gate",
        "polymer": "executable_polymer_gate",
        "full-source": "full_source_complete",
    }[gate]
    return 0 if checked.get(key) is True else 1


if __name__ == "__main__":
    raise SystemExit(main())
