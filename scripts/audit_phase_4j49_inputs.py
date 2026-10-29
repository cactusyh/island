"""Stdlib-only source and raw LAMMPS FD reconstruction, outside production code."""

import argparse
import json
from pathlib import Path

from scripts.audit_phase_4j48_inputs import audit as source_audit
from scripts.audit_phase_4j48_inputs import identity, require


def audit(root, archive, archive_sha256):
    root = Path(root)
    sources = source_audit(root, archive, archive_sha256)
    counts = {}
    for dp in (0, 2, 4, 6):
        modes = json.loads((root / f"components-dp{dp}.json").read_text())
        for mode, c in modes.items():
            directory = root / f"diagnostic-dp{dp}-{mode}"
            raw = {}
            for line in (directory / "engine.screen").read_text().splitlines():
                if line.startswith("J49FD ") and "$(" not in line:
                    t = line.split()
                    v = [float(x) * 4.184 for x in t[2:]]
                    raw[int(t[1])] = {
                        "bonded": (v[0] + v[1]) + v[2],
                        "lj": v[3],
                        "coulomb": v[4] + v[5],
                        "total": v[6],
                    }
            r = c["finite_difference"]
            require(
                len(raw) == len(r["raw_records"]) * 2, "Missing raw FD energy sample"
            )
            maxima = {k: 0.0 for k in ("bonded", "lj", "coulomb", "total")}
            for n, row in enumerate(r["raw_records"]):
                for k, v in row["components"].items():
                    plus, minus = raw[2 * n][k], raw[2 * n + 1][k]
                    require(
                        plus == v["positive_energy_kj_mol"]
                        and minus == v["negative_energy_kj_mol"],
                        "Changed raw FD energy",
                    )
                    force = -(plus - minus) / (2 * row["h_angstrom"])
                    require(
                        force == v["central_difference_force"],
                        "Changed FD sign/unit/equation",
                    )
                    error = abs(force - v["analytical_force"])
                    require(error == v["absolute_error"], "Changed raw FD error")
                    require(v["passed"] == (error <= 1e-4), "Changed FD acceptance")
                    maxima[k] = max(maxima[k], error)
            require(
                maxima == r["max_error_by_component"] and r["bound"] == 1e-4,
                "Changed component bounds",
            )
            require(
                r["passed"] == all(v <= 1e-4 for v in maxima.values()),
                "Changed FD gate",
            )
            counts[f"dp{dp}/{mode}"] = len(raw)
    workflows = (
        json.loads((root / "diagnostic-workflows.json").read_text())
        if (root / "diagnostic-workflows.json").exists()
        else {}
    )
    for dp, w in workflows.items():
        directory = root / f"workflow-dp{dp}"
        require(
            (directory / "whole.restart").read_bytes()
            == (directory / "resumed.restart").read_bytes(),
            "Native restart mismatch",
        )
        require(
            (directory / "rng-initial.json").read_bytes()
            == (directory / "rng-continuation.json").read_bytes(),
            "Full RNG state mismatch",
        )
        require(
            (directory / "whole.dump").read_bytes()
            == (directory / "resumed.dump").read_bytes(),
            "Full state dump mismatch",
        )
        require(
            w["checkpoint_bytes_equal"]
            and w["complete_RNG_state_equal"]
            and w["full_state_dump_equal"],
            "Changed workflow acceptance",
        )
    result = {
        "schema": "island.j49-independent-inputs.v1",
        "sources": sources,
        "raw_component_FD_reconstructed": counts,
        "native_workflows_reconstructed": sorted(workflows),
        "outside_production_parser": True,
        "source_author_MD_authority": False,
    }
    return {**result, "identity": identity(result)}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("directory", "archive", "output"):
        p.add_argument("--" + name, type=Path, required=True)
    p.add_argument("--archive-sha256", required=True)
    a = p.parse_args(argv)
    a.output.write_text(
        json.dumps(
            audit(a.directory, a.archive, a.archive_sha256), sort_keys=True, indent=2
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
