"""J4 raw-source diatomic oracle plus durable separate-process reconstruction.

The raw-source selection is independent of native selection. The entire pair
inventory consists of one excluded bonded pair; LAMMPS pair_style zero is
therefore an equivalent reference for these specific declared molecules.
"""

import argparse
import builtins
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from validate_pcff_expanded import independent_charges
from validate_pcff_fallbacks import compare, fmt, lammps, raw_row
from validate_prepared_workflow import hashes

from island.chemistry import from_smiles
from island.forcefields import (
    ForceFieldRequest,
    PCFFOptions,
    PreparedForceFieldSources,
    create_evaluator,
    load_prepared_forcefield,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff.fallbacks import DOMAIN_POLICY
from island.forcefields.pcff.organic_domains import PROFILE_NAME
from island.workflows.storage import checksum, json_bytes, publish


def result_data(r, ids):
    return {
        "energy": r.potential_energy,
        "forces": [list(r.forces[i]) for i in ids],
        "components": dict(r.energy_components),
        "model": r.model_fingerprint,
        "parameter": r.parameter_fingerprint,
        "calculation": r.evaluation_fingerprint,
    }


def child(bundle, source, output):
    importer = builtins.__import__

    def blocked(name, *a, **kw):
        if name.split(".")[0] in ("openmm", "rdkit", "foyer", "scipy", "parmed"):
            raise AssertionError(
                "Scientific dependency during offline reconstruction: " + name
            )
        return importer(name, *a, **kw)

    builtins.__import__ = blocked
    try:
        loaded = load_prepared_forcefield(
            bundle, sources=PreparedForceFieldSources(pcff_frc=source)
        )
        loaded.prepared.validate_integrity(loaded.system)
    finally:
        builtins.__import__ = importer
    evaluator = create_evaluator(loaded.system, loaded.prepared)
    ids = sorted(loaded.system.topology.sites)
    fresh = result_data(evaluator.evaluate_fresh(), ids)
    with evaluator.open_session() as session:
        reused = result_data(session.evaluate(), ids)
        independent = result_data(session.evaluate_fresh(), ids)
    publish(
        output,
        json_bytes(
            {
                "pid": os.getpid(),
                "prepared": loaded.prepared.identity,
                "native": loaded.prepared.native_result.identity,
                "fresh": fresh,
                "session": reused,
                "independent": independent,
                "offline_loading": True,
            }
        ),
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--lammps", type=Path)
    p.add_argument("--child-bundle", type=Path)
    a = p.parse_args()
    if a.child_bundle:
        child(a.child_bundle, a.source, a.output)
        return 0
    if a.lammps is None:
        p.error("--lammps required")
    root = a.output
    root.mkdir(parents=True, exist_ok=False)
    declaration = json.loads(
        Path("docs/evidence/phase_4j4_declaration.json").read_text()
    )
    raw = a.source.read_bytes()
    assert checksum(raw) == declaration["source"]["sha256"]
    publish(
        root / "declaration.json",
        json_bytes(
            {
                "parent": checksum(
                    Path("docs/evidence/phase_4j4_declaration.json").read_bytes()
                ),
                "experiment": declaration["numerics"],
                "lammps_sha256": checksum(a.lammps.read_bytes()),
                "source": declaration["source"],
            }
        ),
    )
    results = []
    for case in declaration["cases"][:4]:
        work = root / case["name"]
        work.mkdir()
        out = {"case": case["name"], "passed": False}
        try:
            s = from_smiles(case["smiles"], random_seed=declaration["seed"])
            ids = sorted(s.topology.sites)
            assert len(ids) == 2 and len(s.topology.bonds) == 1
            # Labels and lookup keys authored from the declared atom/source rows,
            # independently of the production parameter/equivalence selections.
            labels = {
                i: "h"
                if s.topology.sites[i].element == "H"
                else s.topology.sites[i].element.lower()
                for i in ids
            }
            halogen = next(t for t in labels.values() if t != "h")
            # Verify the actual automatic bond column, independently from raw text.
            mapping = {}
            section = None
            for line, text in enumerate(raw.decode().splitlines(), 1):
                w = text.split("!", 1)[0].split()
                if not w:
                    continue
                if w[0].startswith("#"):
                    section = w[0][1:]
                elif (
                    section == "auto_equivalence"
                    and len(w) > 6
                    and w[2] in labels.values()
                ):
                    assert w[2] not in mapping
                    mapping[w[2]] = {"bond": w[5], "line": line, "raw": text}
            assert (
                mapping["h"]["bond"] == "h_"
                and mapping[halogen]["bond"] == halogen + "_"
            )
            key = sorted([halogen + "_", "h_"])
            row = raw_row(raw, "quadratic_bond", key)
            r0, k = row["values"]
            q, charge_rows = independent_charges(raw, s, labels)
            prepared = prepare_forcefield(
                s,
                ForceFieldRequest(
                    "pcff",
                    PCFFOptions(
                        source_path=a.source,
                        typing_profile=PROFILE_NAME,
                        resolution_policy=DOMAIN_POLICY,
                        lj=(0, 0, 1),
                        coulomb=(0, 0, 1),
                    ),
                ),
            )
            evaluator = create_evaluator(s, prepared)
            model = prepared.native_result.payload
            assert (
                len(model["terms"]) == 1
                and model["terms"][0]["family"] == "quadratic_bond"
            )
            term = model["terms"][0]
            assert {r["record_id"] for r in term["source_rows"]} == {
                f"quadratic_bond:cff91_auto:{row['line']}"
            }
            np.testing.assert_allclose(
                term["coefficients"], [r0, k * 4.184], atol=1e-12, rtol=0
            )
            assert {r["site"]: r["charge"] for r in model["nonbonded"]} == q
            xyz = np.array([s.coordinates.get(i) for i in ids])
            out["frames"] = []
            baseline = None
            for n, x in enumerate(
                (
                    xyz,
                    xyz + 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3),
                )
            ):
                r = evaluator.evaluate(dict(zip(ids, x, strict=True)))
                e, f, meta = lammps(
                    a.lammps.resolve(),
                    work / f"frame{n}",
                    x,
                    "bond",
                    "harmonic",
                    fmt([k, r0]),
                )
                check = compare(
                    r.potential_energy, np.array([r.forces[i] for i in ids]), e, f
                )
                np.testing.assert_allclose(
                    r.energy_components["quadratic_bond"], e, atol=1e-5, rtol=2e-10
                )
                assert (
                    abs(sum(r.energy_components.values()) - r.potential_energy) < 1e-10
                )
                fd_errors = []

                def energy(z, k=k, r0=r0):
                    return 4.184 * k * (np.linalg.norm(z[0] - z[1]) - r0) ** 2

                for h in declaration["numerics"]["finite_difference_angstrom"]:
                    fd = np.zeros_like(x)
                    for i in range(2):
                        for j in range(3):
                            hi = x.copy()
                            lo = x.copy()
                            hi[i, j] += h
                            lo[i, j] -= h
                            fd[i, j] = -(energy(hi) - energy(lo)) / (2 * h)
                    fd_errors.append(
                        float(np.max(np.abs(fd - np.array([r.forces[i] for i in ids]))))
                    )
                assert fd_errors[-1] < 1e-5
                out["frames"].append(
                    {
                        **check,
                        **meta,
                        "finite_difference_errors": fd_errors,
                        "components": dict(r.energy_components),
                    }
                )
                if n == 0:
                    baseline = result_data(r, ids)
            out.update(
                raw_bond_row=row,
                automatic_bond_columns=mapping,
                charge_rows=charge_rows,
                pair_inventory={"excluded_bonded_pairs": 1, "remaining": 0},
                parameter_inventory_verified=True,
            )
            save_prepared_forcefield(
                s,
                prepared,
                work / "before",
                sources=PreparedForceFieldSources(pcff_frc=a.source),
            )
            before = hashes(work / "before")
            shutil.move(work / "before", work / "relocated")
            cmd = [
                sys.executable,
                str(Path(__file__).resolve()),
                "--source",
                str(a.source.resolve()),
                "--child-bundle",
                str((work / "relocated").resolve()),
                "--output",
                str((work / "child.json").resolve()),
            ]
            run = subprocess.run(
                cmd, cwd=work, capture_output=True, text=True, timeout=180, check=False
            )
            (work / "child.stdout").write_text(run.stdout)
            (work / "child.stderr").write_text(run.stderr)
            assert run.returncode == 0, run.stderr
            cr = json.loads((work / "child.json").read_text())
            assert cr["pid"] != os.getpid()
            assert (
                cr["prepared"] == prepared.identity
                and cr["native"] == prepared.native_result.identity
            )
            max_delta = 0.0
            for result in (cr["fresh"], cr["session"], cr["independent"]):
                for key in ("model", "parameter", "calculation"):
                    assert result[key] == baseline[key]
                for key in ("energy", "forces"):
                    np.testing.assert_allclose(
                        result[key], baseline[key], atol=1e-10, rtol=1e-12
                    )
                    max_delta = max(
                        max_delta,
                        float(
                            np.max(
                                np.abs(np.array(result[key]) - np.array(baseline[key]))
                            )
                        ),
                    )
                assert result["components"].keys() == baseline["components"].keys()
                np.testing.assert_allclose(
                    [result["components"][k] for k in baseline["components"]],
                    list(baseline["components"].values()),
                    atol=1e-10,
                    rtol=1e-12,
                )
            assert hashes(work / "relocated") == before
            out.update(
                passed=True,
                prepared_identity=prepared.identity,
                model_identity=prepared.native_result.identity,
                bundle_hashes=before,
                child_command=cmd,
                child_pid=cr["pid"],
                parent_pid=os.getpid(),
                reconstruction_max_error=max_delta,
                offline_loading=True,
            )
        except Exception as err:  # noqa: BLE001 -- retain all real failures
            out["error"] = type(err).__name__ + ": " + str(err)
        publish(work / "outcome.json", json_bytes(out))
        results.append(out)
        print(case["name"], out.get("error", "passed"), flush=True)
    source_ok = checksum(a.source.read_bytes()) == checksum(raw)
    report = {
        "cases": results,
        "source_unchanged": source_ok,
        "passed": source_ok and len(results) == 4 and all(c["passed"] for c in results),
        "full_source_complete": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    publish(root / "report.json", json_bytes(report))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
