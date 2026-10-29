"""Explicit exact-erfc numerical control, never source-author MD authority.

Single Nosé–Hoover thermostat (tchain 1), fixed box, 0.25 fs. Velocities are
created once by PCG64; its entire state is preserved and no stochastic integrator
is used. Native restart bytes, including thermostat state, must match exactly.
The author PPPM/default-table profile remains blocked and the canonical run_md
API remains fail-closed. This harness only follows measured diagnostic gates.
"""

import json
import shutil
import subprocess
from pathlib import Path

import numpy as np

from island.electrostatic_audit import bind_inputs
from island.forcefields.pdms_opls import check, export_lammps
from island.pdms_lammps import _execute, pin
from scripts.phase_4j47 import write

FIX = "fix thermostat all nvt temp 300 300 100 tchain 1\n"
CONTROLS = "atom_modify sort 0 0\nneigh_modify every 1 delay 0 check no\n"
TOTAL_STEPS = 40


def run_control(system, prepared, directory, *, engine, static_control, solver):
    """Bounded algorithm experiment after actual static/FD acceptance only."""
    binding = bind_inputs(system, prepared)
    check(
        all(e <= 3e-4 for e in static_control["energy_error_kj_mol"].values())
        and all(
            f <= 3e-4 for f in static_control["force_error_kj_mol_angstrom"].values()
        )
        and static_control["finite_difference"]["passed"],
        "diagnostic_static_gate_failed",
        "scientifically_unqualified",
    )
    engine.validate()
    source = prepared.validate_integrity(system).validate_integrity(system)
    root = Path(directory)
    export_lammps(system, source, root)
    # Re-execute before MD: caller-supplied comparison JSON is never authority.
    from scripts.phase_4j49 import MODES, Reference, analyze_mode

    check(solver == MODES["matched-ewald"], "diagnostic_solver_contract_mutation")
    xyz = np.array(
        [system.coordinates.get(a["site_id"]) for a in source.payload["mapping"]]
    )
    probes = [
        (
            next(
                n for n, a in enumerate(source.payload["mapping"]) if a["element"] == el
            ),
            axis,
        )
        for el in ("H", "O", "Si", "C")
        for axis in range(3)
    ]
    reference = Reference(source.payload, xyz)
    try:
        prefix = (
            (root / "system.in")
            .read_text()
            .replace("read_data system.data", "read_data ../system.data")
            .replace("include coefficients.in", "include ../coefficients.in")
        )
        verified = analyze_mode(
            root / "static-revalidation",
            prefix + solver,
            xyz,
            probes,
            engine,
            reference,
            source,
        )
    finally:
        reference.close()
    check(
        all(e <= 3e-4 for e in verified["energy_error_kj_mol"].values())
        and all(f <= 3e-4 for f in verified["force_error_kj_mol_angstrom"].values())
        and verified["finite_difference"]["passed"],
        "fresh_diagnostic_static_gate_failed",
        "scientifically_unqualified",
    )
    write(root / "static-verification.json", verified)
    rng = np.random.default_rng(2026)
    before = rng.bit_generator.state
    velocity = rng.normal(0, 0.003, (len(source.payload["mapping"]), 3))
    after = rng.bit_generator.state
    # A velocity seed/distribution is an algorithm control, not a charge/model fit.
    probe = np.random.Generator(np.random.PCG64(0))
    probe.bit_generator.state = after
    future_probe = probe.bit_generator.random_raw(16).tolist()
    rng_record = {
        "algorithm": "PCG64",
        "future_random_raw_16_probe": future_probe,
        "numpy_version": np.__version__,
        "before_initial_velocity_generation": before,
        "after_initial_velocity_generation": after,
        "integrator_rng": "none; deterministic single Nose-Hoover",
        "initial_velocities_angstrom_fs_before_COM_removal_and_temperature_rescaling": velocity.tolist(),
    }
    write(root / "rng-initial.json", rng_record)
    initial = (root / "system.in").read_text() + solver + CONTROLS
    initial += "timestep 0.25\nthermo 10\nthermo_style custom step pe ke etotal temp\nthermo_modify format float %.17g\nmin_style cg\nminimize 1e-8 1e-8 100 1000\nreset_timestep 0\n"
    for i, (x, y, z) in enumerate(velocity):
        initial += f"set atom {i + 1} vx {x:.17g} vy {y:.17g} vz {z:.17g}\n"
    initial += "velocity all zero linear\nvelocity all scale 300\nwrite_restart initial.restart\n"
    (root / "initialize.in").write_text(initial)
    check(
        _execute(engine, root, "initialize.in") == 0, "diagnostic_minimization_failed"
    )
    lines = (root / "engine.screen").read_text().splitlines()
    index = next(
        i
        for i, line in enumerate(lines)
        if "Energy initial, next-to-last, final" in line
    )
    initial_energy, _, final_energy = map(float, lines[index + 1].split())
    stop = next(line.strip() for line in lines if "Stopping criterion =" in line)
    check(
        np.isfinite([initial_energy, final_energy]).all()
        and final_energy < initial_energy,
        "minimization_did_not_lower_energy",
        "scientifically_unqualified",
    )
    for name in ("engine.log", "engine.screen", "process.stdout", "process.stderr"):
        (root / name).rename(root / ("initialize-" + name))
    prefix = (
        "units real\natom_style full\nread_restart initial.restart\n"
        + solver
        + CONTROLS
        + FIX
    )
    for name, steps in (("whole", TOTAL_STEPS), ("split", TOTAL_STEPS // 2)):
        script = (
            prefix
            + f"thermo 1\nthermo_style custom step pe ke etotal temp\nthermo_modify format float %.17g\nrun {steps} start 0 stop {TOTAL_STEPS}\nwrite_restart {name}.restart\nwrite_dump all custom {name}.dump id x y z vx vy vz fx fy fz modify sort id format float %.17g\n"
        )
        (root / (name + ".in")).write_text(script)
        check(_execute(engine, root, name + ".in") == 0, "diagnostic_NVT_failed")
        for file in ("engine.log", "engine.screen", "process.stdout", "process.stderr"):
            (root / file).rename(root / (name + "-" + file))
    relocated = root.parent / (root.name + "-relocated")
    shutil.copytree(root, relocated)
    script = (
        "units real\natom_style full\nread_restart split.restart\n"
        + solver
        + CONTROLS
        + FIX
        + f"thermo 1\nthermo_style custom step pe ke etotal temp\nthermo_modify format float %.17g\nrun {TOTAL_STEPS // 2} start 0 stop {TOTAL_STEPS}\nwrite_restart resumed.restart\nwrite_dump all custom resumed.dump id x y z vx vy vz fx fy fz modify sort id format float %.17g\n"
    )
    (relocated / "resume.in").write_text(script)
    # _execute starts a new engine process; use a separate fresh Python process too.
    code = """import sys
import json
from pathlib import Path
import numpy as np
from island.pdms_lammps import _execute,LAMMPSOptions
from scripts.phase_4j47 import write
root=Path(sys.argv[1])
record=json.loads((root/'rng-initial.json').read_text())
assert record['algorithm']=='PCG64' and record['numpy_version']==np.__version__
rng=np.random.Generator(np.random.PCG64(0))
rng.bit_generator.state=record['after_initial_velocity_generation']
result=_execute(LAMMPSOptions(sys.argv[2],sys.argv[3]),root,'resume.in')
record['after_initial_velocity_generation']=rng.bit_generator.state
probe=np.random.Generator(np.random.PCG64(0))
probe.bit_generator.state=rng.bit_generator.state
assert probe.bit_generator.random_raw(16).tolist()==record['future_random_raw_16_probe']
write(root/'rng-continuation.json',record)
raise SystemExit(result)
"""
    import sys

    proc = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(relocated.resolve()),
            str(engine.validate()),
            engine.expected_sha256,
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    check(proc.returncode == 0, "diagnostic_fresh_process_continuation_failed")
    whole, resumed = (
        (root / "whole.restart").read_bytes(),
        (relocated / "resumed.restart").read_bytes(),
    )
    equal = whole == resumed
    mismatch = [i for i, (a, b) in enumerate(zip(whole, resumed)) if a != b]
    if len(whole) != len(resumed):
        mismatch.append(min(len(whole), len(resumed)))
    dump_equal = (root / "whole.dump").read_bytes() == (
        relocated / "resumed.dump"
    ).read_bytes()
    rng_equal = (root / "rng-initial.json").read_bytes() == (
        relocated / "rng-continuation.json"
    ).read_bytes()
    for name in (
        "resume.in",
        "resumed.restart",
        "resumed.dump",
        "engine.screen",
        "engine.log",
        "process.stdout",
        "process.stderr",
        "rng-continuation.json",
    ):
        shutil.copyfile(relocated / name, root / name)
    shutil.rmtree(relocated)
    result = {
        "schema": "island.j49-diagnostic-workflow.v1",
        "input_binding_identity": binding.identity,
        "profile_identity": prepared.payload["selected_profile"],
        "charge_origin": prepared.payload["charge_origin"],
        "source_model_unchanged": True,
        "canonical_author_MD_authority": False,
        "numerical_solver_override": solver,
        "diagnostic_kernel_sha256": engine.expected_sha256,
        "minimization": {
            "executed": True,
            "style": "cg",
            "max_iterations": 100,
            "initial_energy_kj_mol": initial_energy * 4.184,
            "final_energy_kj_mol": final_energy * 4.184,
            "stopping_criterion": stop,
            "converged": "max iterations" not in stop,
            "not_equilibrium": True,
        },
        "NVT": {
            "executed": True,
            "steps": TOTAL_STEPS,
            "split_step": TOTAL_STEPS // 2,
            "dt_fs": 0.25,
            "temperature_kelvin": 300,
            "thermostat": "single Nose-Hoover; tchain 1; damping 100 fs",
        },
        "relocated_separate_process_continuation": True,
        "checkpoint_bytes_equal": equal,
        "native_checkpoint_mismatched_byte_offsets_zero_based": mismatch,
        "whole_checkpoint_sha256": pin(root / "whole.restart"),
        "resumed_checkpoint_sha256": pin(root / "resumed.restart"),
        "complete_RNG_state_equal": rng_equal,
        "RNG_reconstructed_in_fresh_process": True,
        "future_random_stream_probe_equal": True,
        "full_state_dump_equal": dump_equal,
        "bounded_diagnostic_workflow_passed": equal and rng_equal and dump_equal,
        "NPT_executed": False,
        "production_validated": False,
        "simulation_readiness": "not_established",
    }
    write(root / "result.json", result)
    # Re-read the RNG serialization as a check against accidental truncation.
    check(
        json.loads((root / "rng-initial.json").read_text()) == rng_record,
        "RNG_serialization_changed",
    )
    return result
