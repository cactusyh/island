"""Executed source audit from retained data. Never invokes QM, MD or fitting."""

import argparse
import json
import lzma
import platform
import resource
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import numpy as np

from island.pcff_completion import inventory, pin, require, signed, write
from island.pcff_qm_pilot.hmds_local_baseline import correction, references
from island.pcff_qm_pilot.hmds_precision import read
from island.pcff_qm_pilot.hmds_source_audit import (
    CONVENTIONS,
    closure,
    projection,
    source_graph,
    source_sum,
)
from scripts.analyze_phase_4j68 import complete_check, reconstruct

ROOT = Path("docs/evidence/phase_4j69")
BASE = "8784955934ef6b750d54b3cfea49d824a1218c76"
REQUIRED = "c957d9f327a025e57052531baed184c110e313aa"
FROZEN = "a72ab77e0e4ebc0ff52a9b7287cfef01569b33aa0394a5dfc8ab2f92c96cf36a"


def inputs(repo, root):
    plan = read(root / "plan.json")
    require(pin(root / "plan.json") == FROZEN, "J69_frozen_plan")
    require(
        all(pin(repo / n) == h for n, h in plan["parent_files"].items()),
        "J69_parent_inputs",
    )
    require(
        pin(root / "source-assets.json") == plan["source_assets_sha256"],
        "J69_source_inventory",
    )
    assets = read(root / "source-assets.json")["sources"]
    raw = {}
    from hashlib import sha256

    for n, a in assets.items():
        if a["status"] == "source_unavailable":
            continue
        p = root / a["archive"]
        require(pin(p) == a["archive_sha256"], "J69_archived_source_mutation")
        raw[n] = lzma.decompress(p.read_bytes())
        require(sha256(raw[n]).hexdigest() == a["sha256"], "J69_source_bytes_mutation")
    return plan, raw


def converter_text(root, name):
    path = root / "converter" / name
    if path.exists():
        return path.read_text()
    raw = lzma.decompress(path.with_suffix(path.suffix + ".xz").read_bytes())
    from hashlib import sha256

    hashes = read(root / "converter/native-output-hashes.json")
    require(
        sha256(raw).hexdigest() == hashes[name], "J69_native_converter_output_bytes"
    )
    return raw.decode()


def localize(model, selected, rows, data, values, extra=None):
    report = []
    for index, row in enumerate(rows):
        families, _ = source_sum(model, selected, row["coordinates"], jets=True)
        fixed = sum(-e.d.reshape(27, 3) for e in families.values())
        adjustable = data["force"][index] @ values[:58]
        adjustable_families = {}
        for name in sorted({n.split(":")[0] for n in data["names"]}):
            columns = [
                j for j, n in enumerate(data["names"]) if n.split(":")[0] == name
            ]
            adjustable_families[name] = (
                np.take(data["force"][index], columns, axis=-1) @ values[columns]
            )
        additions = np.zeros((27, 3))
        if extra:
            _, features = correction(model, row["coordinates"], extra)
            additions = features @ values[58:60]
        residual = fixed + adjustable + additions - data["QM_force"][index]
        atoms = []
        for atom in (5, 2, 6, 3, 10):
            j = atom - 1
            r = residual[j]
            vectors = {k: -e.d.reshape(27, 3)[j] for k, e in families.items()}
            vectors["adjustable58"] = adjustable[j]
            if extra:
                vectors["research_diagnostic"] = additions[j]
            diagnostic = {}
            for variant in ("bridge_angle", "SiO_stretch"):
                _, f = correction(model, row["coordinates"], variant)
                diagnostic[variant] = [projection(f[j, :, k], r) for k in range(2)]
            separated = {
                "fixed:" + k: v for k, v in vectors.items() if k != "adjustable58"
            }
            separated.update(
                {"adjustable:" + k: v[j] for k, v in adjustable_families.items()}
            )
            sum_norms = sum(np.linalg.norm(v) for v in separated.values())
            total = fixed[j] + adjustable[j] + additions[j]
            atoms.append(
                {
                    "atom": atom,
                    "residual": r.tolist(),
                    "residual_norm": float(np.linalg.norm(r)),
                    "fixed_total": fixed[j].tolist(),
                    "total_MM_force": total.tolist(),
                    "QM_force": data["QM_force"][index, j].tolist(),
                    "families": {k: projection(v, r) for k, v in vectors.items()},
                    "fixed_and_adjustable_families": {
                        k: projection(v, r) for k, v in separated.items()
                    },
                    "diagnostic_columns": diagnostic,
                    "sum_family_norms": float(sum_norms),
                    "norm_of_sum": float(np.linalg.norm(total)),
                    "cancellation_difference": float(sum_norms - np.linalg.norm(total)),
                }
            )
        report.append(
            {
                "configuration": row["id"],
                "subset": "protected22" if index < 22 else "new4",
                "atoms": atoms,
            }
        )
    return report


def equilibrium_audit(model, raw, graph):
    lines = raw.decode().splitlines()
    source_rows = read_source_rows(lines)
    selected = {r["id"]: r for r in graph["nodes"]}
    refs = []
    for interaction in ("quartic_bond:2,5", "quartic_bond:5,6", "quartic_angle:2,5,6"):
        node = selected[interaction]
        refs.append(
            {
                "interaction": interaction,
                "source_rows": node["source_rows"],
                "internal_equilibrium": node["internal_values"][0],
                "LAMMPS_equilibrium": node["LAMMPS_coefficients"][0],
                "coefficient_units": "angstrom for bond; radians internally and degrees in LAMMPS input for angle; K2/K3/K4 use radian displacements",
                "dependent_interactions": [
                    r["id"]
                    for r in graph["nodes"]
                    if any(d["interaction"] == interaction for d in r["dependencies"])
                ],
            }
        )
    osi = [
        r
        for r in source_rows
        if r["section"] == "quartic_angle" and r["types"] == ["osi", "sio", "osi"]
    ]
    return signed(
        {
            "schema": "island.j69-conventions.v1",
            "conventions": CONVENTIONS,
            "references": refs,
            "O_Si_O": {
                "source_rows": osi,
                "HMDS_sites": [],
                "status": "not present in exact HMDS graph; each silicon has one oxygen; no numerical O-Si-O coverage inferred",
            },
            "reversal": "Test all source-energy families under full Cartesian rigid rotation and relabelled coordinate/source-site reversal; BA exchange left/right and AA swap outer I/L at fixed center/shared arm",
            "polynomials": "K2*dq^2+K3*dq^3+K4*dq^4, derivative 2*K2*dq+3*K3*dq^2+4*K4*dq^3; no half",
            "equilibrium_propagation": "Every graph dependency edge reads the same audited source equilibrium; no global equilibrium revision made",
            "J68_diagnostics": {
                v: {
                    "equation": "a*sum(q-q0)+b*sum((q-q0)^2)",
                    "force": "-sum((a+2*b*(q-q0))*grad(q))",
                    "references": references(model, v),
                    "origin": "research_derived",
                    "propagates_global_source_equilibrium": False,
                    "source_authority_for_changed_equilibrium": "source_unavailable",
                }
                for v in ("bridge_angle", "SiO_stretch")
            },
            "phase": {
                "source_plus_cos": "pcff.frc:2657",
                "converter_copies_unchanged_phase": "GetParameters.c:342-350",
                "engine_minus_cos": "dihedral_class2.cpp and dihedral_class2.rst",
                "research_transform": "+180 degree yields exact algebraic equality for same Phi; unchanged J60/J63 convention",
                "unresolved": "printed source geometry and absent author artifacts cannot independently establish native PCFF Phi mapping; existing AAT linear-Phi/cos difference retained",
            },
            "missing_terms": graph["missing_interactions"],
        }
    )


def read_source_rows(lines):
    section = None
    out = []
    for i, line in enumerate(lines, 1):
        w = line.split("!", 1)[0].split()
        if not w:
            continue
        if w[0].startswith("#"):
            section = w[0][1:]
        if section in ("quartic_bond", "quartic_angle") and w[0][0].isdigit():
            n = 2 if section == "quartic_bond" else 3
            out.append(
                {
                    "section": section,
                    "row_id": f"{section}:cff91:{i}",
                    "line": i,
                    "types": w[2 : 2 + n],
                    "raw": line,
                    "original_values": w[2 + n :],
                }
            )
    return out


def converter_partial(graph, stdout):
    """Audit only coefficients actually printed before strict conversion failed."""
    headers = {
        "Bond Types and  Parameters": ("quartic_bond", 2),
        "Angle Types and Parameters": ("quartic_angle", 3),
        "BondBond Types and  Parameters": ("bond-bond", 3),
        "BondAngle Types and  Parameters": ("bond-angle", 3),
    }
    family, arity, comparisons = None, None, []
    for line, text in enumerate(stdout.splitlines(), 1):
        t = text.strip()
        if t in headers:
            family, arity = headers[t]
            continue
        if not family:
            continue
        words = t.split()
        if len(words) <= arity:
            continue
        try:
            values = list(map(float, words[arity:]))
        except ValueError:
            family = None
            continue
        labels = words[:arity]
        nodes = [
            n
            for n in graph["nodes"]
            if n["family"] == family and n["site_types"] in (labels, labels[::-1])
        ]
        require(bool(nodes), "J69_converter_printed_site_mapping")
        for node in nodes:
            reverse = node["site_types"] != labels
            expected = list(node["LAMMPS_coefficients"])
            if reverse and family == "bond-angle":
                expected = expected[:2][::-1] + expected[2:][::-1]
            elif reverse and family == "bond-bond":
                expected = expected[:1] + expected[1:][::-1]
            require(len(values) == len(expected), "J69_converter_coefficient_order")
            difference = float(np.max(abs(np.array(values) - expected)))
            require(difference <= 0.0000500001, "J69_converter_partial_mismatch")
            comparisons.append(
                {
                    "interaction": node["id"],
                    "source_rows": node["source_rows"],
                    "stdout_line": line,
                    "printed_types": labels,
                    "reverse": reverse,
                    "printed_values": values,
                    "generated_values": expected,
                    "maximum_difference": difference,
                    "units": "LAMMPS real; degree reference angles, radian angular coefficients",
                    "classification": "within printed four-decimal rounding",
                }
            )
    require(len(comparisons) == 173, "J69_converter_partial_coverage")
    return signed(
        {
            "schema": "island.j69-partial-converter-closure.v1",
            "comparisons": comparisons,
            "complete_export_available": False,
            "returncode": 17,
            "missing_term": "endbonddihedral hc c3 sio c3",
            "engine_execution": False,
            "inference": "only printed bond/angle/BB/BA parameters; absent sections remain unverified, not zero",
        }
    )


def inherited_authority(repo, model):
    """Inventory J60 unit/pair authorities already bound by the frozen model."""
    from hashlib import sha256

    path = repo / "docs/evidence/phase_4j60/hamiltonian-definition.json"
    require(
        pin(path) == model["J60_equation_authority_sha256"],
        "J69_inherited_equation_authority",
    )
    authority = read(path)
    sources = {}
    for name, record in authority["sources"].items():
        relative = (
            "docs/evidence/phase_4j60/sources/"
            + record["file"].replace("/", "__")
            + ".xz"
        )
        archive = repo / relative
        require(
            sha256(lzma.decompress(archive.read_bytes())).hexdigest()
            == record["sha256"],
            "J69_inherited_source_identity",
        )
        sources[name] = {
            **record,
            "retained_archive": relative,
            "archive_sha256": pin(archive),
            "origin": "native_source_inherited_from_frozen_J60_authority",
            "parser": "raw-source archive; declared equation interpretation retained verbatim",
        }
    return signed(
        {
            "schema": "island.j69-inherited-source-authority.v1",
            "parent_sha256": pin(path),
            "sources": sources,
            "conventions": authority["conventions"],
            "new_source_policy": False,
        }
    )


def numerical(repo, root):
    plan, sources = inputs(repo, root)
    _, model, rows, data, _lo, _hi = reconstruct(repo)
    frc = sources["tools/msi2lmp/frc_files/pcff.frc"]
    graph, selected = source_graph(model, frc)
    from island.pcff_qm_pilot.hmds_hamiltonian import accounting

    pairs = [
        r
        for r in accounting(model, rows[0]["coordinates"])["table"]
        if r["family"] == "LJ96_Coulomb"
    ]
    for pair in pairs:
        pair["source_dependencies"] = [
            f"{family}:{i}"
            for i in pair["sites"]
            for family in ("nonbond(9-6)", "charge")
        ]
        pair["origin"] = (
            "derived sixth-power mixing and declared isolated-pair research convention; not independent native pair rows"
        )
    report = read(repo / "docs/evidence/phase_4j54/audit/hmds/report.json")
    # Retain original row provenance verbatim, independently bind raw source lines.
    lines = frc.decode().splitlines()
    for row in report["source_rows"].values():
        require(
            lines[row["location"]["line"] - 1] == row["raw_text"]
            and row["source_sha256"] == report["source_hashes"]["pcff.frc"],
            "J69_raw_row_location",
        )
    inventory_record = signed(
        {
            "schema": "island.j69-row-inventory.v1",
            "rows": report["source_rows"],
            "typing": model["typing"],
            "source_hashes": model["source_hashes"],
            "source_complete": False,
            "missing_interactions": graph["missing_interactions"],
            "authoritative_scope": "hash-pinned LAMMPS-distributed PCFF FRC with cff91 namespace; not original author artifacts or licensed global PCFF authority",
            "parser_files": {
                n: pin(repo / n)
                for n in (
                    "src/island/forcefields/pcff/source.py",
                    "src/island/forcefields/pcff/class2.py",
                    "scripts/pcff_j5_reference.py",
                    "src/island/pcff_qm_pilot/hmds_source_audit.py",
                )
            },
        }
    )
    c = closure(model, selected, rows, plan["finite_difference_step_A"])
    # Explicitly distinguish floating-point closure from historical convention
    # contrasts; no discrepancy is promoted to source authority by its size.
    c = classify_closure(c)
    models = read(repo / "docs/evidence/phase_4j68/executed-comparison.json")["models"]
    localization = {}
    checks = {}
    for key, variant in (
        ("bounded-Clarabel", None),
        ("bound-free-Clarabel", None),
        ("bridge_angle-clarabel", "bridge_angle"),
        ("SiO_stretch-clarabel", "SiO_stretch"),
    ):
        values = np.array(models[key]["parameter_values"])
        localization[key] = localize(model, selected, rows, data, values, variant)
        if variant is None:
            checks[key] = complete_check(model, rows, data, values)
    return {
        "inherited-source-authority.json": inherited_authority(repo, model),
        "pair-source-dependencies.json": signed(
            {
                "schema": "island.j69-pair-dependencies.v1",
                "pairs": pairs,
                "source_pair_exclusion_authority": "source_unavailable; declared research policy unchanged",
                "charges": model["charges"],
                "charge_origin": model["charge_origin"],
                "charge_identity": model["charge_identity"],
                "charge_normalization": False,
            }
        ),
        "converter-partial-closure.json": converter_partial(
            graph, converter_text(root, "stdout.log")
        ),
        "source-row-inventory.json": inventory_record,
        "source-dependency-graph.json": graph,
        "source-to-evaluator-closure.json": c,
        "equilibrium-phase-prefactor-audit.json": equilibrium_audit(model, frc, graph),
        "residual-family-localization.json": signed(
            {
                "models": localization,
                "units": "energy kcal/mol, force kcal/mol/angstrom; feature forces per coefficient unit",
                "interpretation": "projections and cancellation; largest family norm is not causal evidence",
            }
        ),
        "reconstruction.json": signed(
            {
                "dataset_and_producers_sha256": pin(
                    repo / "docs/evidence/phase_4j68/dataset-and-producers.json"
                ),
                "configuration_order": [r["id"] for r in rows],
                "observations": 26,
                "protected": 22,
                "graph_identity": model["graph_identity"],
                "fixed_model_identity": model["identity"],
                "58_column_names": data["names"],
                "one_energy_only_intercept": True,
                "checks": checks,
                "J68_comparison": {
                    k: {n: v[n] for n in ("rho", "lower", "classification")}
                    for k, v in models.items()
                },
                "new_QM": 0,
                "MD_steps": 0,
                "reserved_validation_untouched": True,
            }
        ),
    }


def classify_closure(c):
    for row in c["configurations"]:
        for family in row["families"].values():
            family["classification"] = (
                "roundoff; finite-difference check within frozen derivative tolerance"
            )
    for row in c["phase_controls"]:
        row["classification"] = (
            "convention_difference_under_same_declared_Phi; native geometry authority unresolved; already disclosed by J60/J63"
        )
    return signed({k: v for k, v in c.items() if k != "identity"})


def converter(repo, root):
    """One strict original converter attempt; never ignore missing coefficients."""
    from island.forcefields.pcff.interop import export_msi2lmp_car_mdf
    from island.workflows import storage
    from island.workflows.bundle import system_from

    _, raw = inputs(repo, root)
    _, model, rows, _, _, _ = reconstruct(repo)
    out = root / "converter"
    require(not out.exists(), "J69_no_converter_overwrite")
    out.mkdir()
    start = time.monotonic()
    with tempfile.TemporaryDirectory(prefix="J69-converter-build-") as temp:
        work = Path(temp)
        sources = []
        for n, b in raw.items():
            if n.startswith("tools/msi2lmp/src/"):
                p = work / Path(n).name
                p.write_bytes(b)
                if p.suffix == ".c":
                    sources.append(str(p))
        exe = work / "msi2lmp"
        build_command = ["cc", "-O2", "-o", str(exe), *sorted(sources), "-lm"]
        b = subprocess.run(
            build_command, capture_output=True, text=True, timeout=60, check=False
        )
        write(
            out / "build.json",
            {
                "argv": build_command,
                "returncode": b.returncode,
                "stdout": b.stdout,
                "stderr": b.stderr,
                "compiler": subprocess.check_output(["cc", "--version"], text=True),
            },
        )
        require(b.returncode == 0, "J69_converter_build_failed")
        system = system_from(
            storage.decode(
                storage.read_json(
                    repo / "docs/evidence/phase_4j54/audit/hmds/system.json"
                )["payload"]
            )
        )
        for i, x in enumerate(rows[0]["coordinates"], 1):
            system.coordinates.set(i, x)
        export_msi2lmp_car_mdf(
            system,
            out / "input",
            types={int(k): v for k, v in model["typing"].items()},
            charges={int(k): v for k, v in model["charges"].items()},
            provenance={
                "model_identity": model["identity"],
                "coordinate_id": rows[0]["id"],
            },
        )
        frc = work / "pcff.frc"
        frc.write_bytes(raw["tools/msi2lmp/frc_files/pcff.frc"])
        command = [
            str(exe),
            "system",
            "-class",
            "II",
            "-frc",
            str(frc),
            "-print",
            "3",
            "-nocenter",
        ]
        r = subprocess.run(
            command,
            cwd=out / "input",
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
        )
        (out / "stdout.log").write_text(r.stdout)
        (out / "stderr.log").write_text(r.stderr)
        result = {
            "schema": "island.j69-converter-attempt.v1",
            "argv": command,
            "returncode": r.returncode,
            "executable_sha256": pin(exe),
            "input_files": inventory(out / "input"),
            "source_complete": False,
            "ignore_used": False,
            "new_QM": 0,
            "MD_steps": 0,
            "runtime_seconds": time.monotonic() - start,
            "outcome": "strict conversion failed; no complete coefficient export authority"
            if r.returncode
            else "completed; completeness separately audited",
        }
        write(out / "attempt.json", signed(result))
        shutil.copyfile(exe, out / "msi2lmp")
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--converter", action="store_true")
    a = p.parse_args()
    repo = a.repo.resolve()
    root = repo / ROOT
    if a.converter:
        print(converter(repo, root))
    else:
        start = time.monotonic()
        resource.setrlimit(resource.RLIMIT_AS, (4294967296, 4294967296))
        result = numerical(repo, root)
        for n, v in result.items():
            require(not (root / n).exists(), "J69_no_result_overwrite")
            write(root / n, v)
        write(
            root / "analysis-resources.json",
            {
                "runtime_seconds": time.monotonic() - start,
                "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                "python": platform.python_version(),
                "numpy": np.__version__,
                "workers": 1,
                "new_QM": 0,
                "MD_steps": 0,
                "research_solver_calls": 0,
            },
        )
        print(json.dumps({"status": "completed", "artifacts": list(result)}))
