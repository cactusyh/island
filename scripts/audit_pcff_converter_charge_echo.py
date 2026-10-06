"""Guanidinium diagnostic: raw source total .9999 e is supplied, never repaired.

No valid native model is created. msi2lmp reads these explicit CAR charges; it
cannot supply missing increments or establish a different charge convention.
"""

import argparse
from pathlib import Path

import numpy as np
from pcff_j5_reference import charge_raw, read_source
from validate_pcff_singlepoint import data_sections, reference_inputs, run

from island.workflows import storage
from island.workflows.bundle import system_from


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--declaration",
        type=Path,
        default=Path("docs/evidence/phase_4j13_declaration.json"),
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    d = storage.read_json(args.declaration)
    case = next(c for c in d["blocked_cases"] if c["name"] == "guanidinium")
    original = Path(case["system_path"])
    converter = Path(d["reference"]["converter_path"])
    for path, expected in [
        (original, case["system_sha256"]),
        (args.source, d["source"]["sha256"]),
        (converter, d["reference"]["converter_path_sha256"]),
    ]:
        if storage.checksum(path.read_bytes()) != expected:
            raise ValueError("Declared input checksum mismatch: " + str(path))
    system = system_from(storage.decode(storage.read_json(original)))
    labels = {
        i: {"C": "c+", "N": "nr", "H": "h*"}[a.element]
        for i, a in system.topology.sites.items()
    }
    q = charge_raw(read_source(args.source.read_bytes()), system, labels)
    if q["complete"]:
        raise ValueError("Expected retained charge failure; source/convention changed")
    charges = {i: float(v) for i, v in q["partial_charges"].items()}
    ids = sorted(labels)
    args.output.mkdir(exist_ok=False)
    storage.publish(args.output / "declaration.json", storage.json_bytes(d))
    reference_inputs(
        system,
        labels,
        charges,
        np.array([system.coordinates.get(i) for i in ids]),
        args.output,
    )
    result = {
        "independent_charge": q,
        "native_charge_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    try:
        run(
            [
                str(converter.resolve()),
                "reference",
                "-class",
                "II",
                "-frc",
                str(args.source.resolve()),
                "-p",
                "3",
                "-nocenter",
            ],
            args.output,
            "converter",
        )
        values = [
            float(r[3]) for r in data_sections(args.output / "reference.data")["Atoms"]
        ]
        result.update(
            converter_success=True, charge_readback=values, total_readback=sum(values)
        )
        np.testing.assert_allclose(
            values, [charges[i] for i in ids], atol=1e-12, rtol=0
        )
        result["charge_vector_preserved"] = True
    except Exception as e:  # noqa: BLE001 -- retain diagnostic converter failures
        result.update(converter_success=False, error=type(e).__name__ + ": " + str(e))
    storage.publish(
        args.output / "report.json", storage.json_bytes(storage.encode(result))
    )
    return 0 if result.get("charge_vector_preserved") else 1


if __name__ == "__main__":
    raise SystemExit(main())
