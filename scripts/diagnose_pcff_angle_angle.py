"""Independent J1 converter-role diagnostic; never replaces original acceptance.

Pinned find_angleangle_data calls find_match, which also full-reverses four
labels. That can move the central atom. This reference resolves center-fixed
source rows independently, retaining both raw converter output and overrides.
"""

import argparse
import shutil
from decimal import Decimal
from pathlib import Path

import numpy as np
from validate_pcff_class2 import independent_rows
from validate_pcff_singlepoint import (
    GROUPS,
    TOL,
    aa_reference_coefficients,
    data_sections,
    digest,
    lammps_input,
    run,
)

from island.charge_references.records import unpack
from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields.pcff import (
    load_pcff_model,
    load_pcff_parameters,
    load_pcff_source,
)
from island.forcefields.pcff.automatic import graph_system
from island.workflows.storage import json_bytes, publish


def source_aa_coefficients(raw, sections, labels):
    table, equiv = independent_rows(raw)
    fixed = {}
    evidence = []
    for improper in sections["Impropers"]:
        a, b, c, d = map(int, improper[2:])
        values = []
        for ids in ((a, b, c, d), (c, b, a, d), (a, b, d, c)):
            types = [labels[i] for i in ids]
            eq = [equiv[t][7] for t in types]
            selected = []
            for query in (types, eq):
                swap = [query[3], query[1], query[2], query[0]]
                selected = [r for r in table["angle-angle"] if r[2] in (query, swap)]
                if selected:
                    break
            if not selected:
                raise ValueError(f"Missing independent AA source row: {types}")
            version = max(r[1] for r in selected)
            selected = [r for r in selected if r[1] == version]
            coeff = {Decimal(r[3][0]) for r in selected}
            if len(coeff) != 1:
                raise ValueError("Ambiguous independent AA rows")
            values.append(float(next(iter(coeff))))
            evidence.append(
                {"ids": ids, "labels": types, "source_lines": [r[0] for r in selected]}
            )
        key = improper[1]
        if key in fixed and fixed[key] != values:
            raise ValueError("Type-specific AA role ambiguity")
        fixed[key] = values
    # H5 correction supplies independent equilibrium angles. Only the three
    # source coefficients are replaced here, in their explicit physical roles.
    commands = []
    for line in aa_reference_coefficients(sections):
        fields = line.split()
        fields[3:6] = [format(v, ".17g") for v in fixed[fields[1]]]
        commands.append(" ".join(fields))
    return commands, evidence


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", required=True)
    p.add_argument("--case", required=True)
    p.add_argument("--lammps", required=True)
    p.add_argument("--output", required=True)
    a = p.parse_args(argv)
    original = Path(a.case).resolve()
    out = Path(a.output).resolve()
    out.mkdir(parents=True, exist_ok=False)
    inputs = [
        "typing.json",
        "charges.json",
        "parameters.json",
        "model.json",
        "reference.data",
        "frame-0/input.lmp",
    ]
    hashes = {n: digest(original / n) for n in inputs}
    declaration = {
        "schema": "island_j1_reference_aa_diagnostic_v1",
        "inputs": hashes,
        "source_sha256": digest(a.source),
        "single_change": "center-fixed AA coefficient selection from raw source; preserve H5 equilibrium correction",
        "tolerances": TOL,
        "timeout_seconds": 180,
        "lammps_sha256": digest(a.lammps),
    }
    publish(out / "declaration.json", json_bytes(declaration))
    source = load_pcff_source(a.source)
    parameter = load_pcff_parameters(original / "parameters.json", source)
    model = load_pcff_model(original / "model.json", parameter)
    auto = unpack(parameter.json_text)["charge_record"]["automatic_typing"]
    system = graph_system(auto["graph"])
    labels = auto["assignments"]
    ids = sorted(labels)
    sections = data_sections(original / "reference.data")
    commands, evidence = source_aa_coefficients(
        source.raw, sections, {k + 1: labels[i] for k, i in enumerate(ids)}
    )
    publish(
        out / "source-aa-override.json",
        json_bytes({"commands": commands, "source_evidence": evidence}),
    )
    # Exact originally declared coordinates, as persisted in the input commands.
    text = (original / "frame-0/input.lmp").read_text()
    xyz = np.array(
        [
            [float(line.split()[j]) for j in (4, 6, 8)]
            for line in text.splitlines()
            if line.startswith("set atom ")
        ]
    )
    evaluator = PCFFSinglePointEvaluator(system, model)
    comparisons = []
    for n, frame in enumerate(
        (xyz, xyz + 0.015 * np.sin(np.arange(xyz.size).reshape(xyz.shape) + 0.3))
    ):
        work = out / f"frame-{n}"
        work.mkdir()
        shutil.copyfile(original / "reference.data", work / "reference.data")
        (work / "input.lmp").write_text(
            lammps_input(frame, [0, 0, 1], [0, 0, 1], commands)
        )
        run([Path(a.lammps).resolve(), "-in", "input.lmp"], work, "lammps")
        line = [
            x
            for x in (work / "lammps.log").read_text().splitlines()
            if x.startswith("ENERGIES ")
        ][-1]
        ref = np.array(list(map(float, line.split()[1:]))) * 4.184
        dump = (
            (work / "forces.dump")
            .read_text()
            .split("ITEM: ATOMS id x y z fx fy fz\n")[1]
        )
        f = (
            np.array(
                [
                    list(map(float, r.split()[4:7]))
                    for r in dump.splitlines()
                    if r.strip()
                ]
            )
            * 4.184
        )
        result = evaluator.evaluate(dict(zip(ids, frame)))
        groups = {**GROUPS, "improper": ["angle-angle", "wilson_out_of_plane"]}
        actual = [result.potential_energy] + [
            sum(result.energy_components[k] for k in g) for g in groups.values()
        ]
        force = np.array([result.forces[i] for i in ids])
        comparisons.append(
            {
                "energy_max_abs": float(np.max(np.abs(actual - ref))),
                "force_max_abs": float(np.max(np.abs(force - f))),
            }
        )
        np.testing.assert_allclose(
            actual, ref, atol=TOL["energy_atol_kj_mol"], rtol=TOL["rtol"]
        )
        np.testing.assert_allclose(
            force, f, atol=TOL["force_atol_kj_mol_angstrom"], rtol=TOL["rtol"]
        )
    if hashes != {n: digest(original / n) for n in inputs}:
        raise ValueError("Historical input mutation")
    publish(
        out / "report.json",
        json_bytes(
            {
                "passed": True,
                "comparisons": comparisons,
                "original_unchanged": True,
                "hashes": {
                    str(p.relative_to(out)): digest(p)
                    for p in out.rglob("*")
                    if p.is_file()
                },
            }
        ),
    )
    print(comparisons)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
