"""Synthetic PCFF identity/workflow regressions; not scientific parameter data."""

import pytest
from test_pcff_evaluation import bound, parameters, synthetic  # noqa: F401

from island.evaluation import PCFFSinglePointEvaluator
from island.forcefields import (
    PreparedForceFieldSources,
    adopt_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    expanded,
    special_pair_policy,
    type_pcff_atoms,
)
from island.minimization import MinimizationOptions
from island.workflows import PreparedWorkflowConfig, start_prepared_bundle_workflow
from island.workflows import prepared as workflow


@pytest.fixture
def expanded_case(bound, tmp_path, monkeypatch):  # noqa: F811
    system, old, _ = bound
    source = old.assignment.source
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", source.expected_sha256)
    typing = type_pcff_atoms(system, source, profile=expanded.PROFILE_NAME)
    charges = assign_automatic_pcff_charges(system, typing)
    assignment = assign_pcff_parameters(system, typing, charges)
    spec = define_pcff_model(
        assignment, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    assert spec.payload["model_definition_complete"]
    assert not any(t["family"] == "wilson_out_of_plane" for t in spec.payload["terms"])
    prepared = adopt_forcefield(system, "pcff", spec)
    source_path = tmp_path / "synthetic.frc"
    source_path.write_bytes(source.raw)
    sources = PreparedForceFieldSources(pcff_frc=source_path)
    root = save_prepared_forcefield(
        system, prepared, tmp_path / "input", sources=sources
    )
    config = PreparedWorkflowConfig(
        total_steps=4,
        segment_steps=2,
        recording_interval=1,
        max_evaluations_per_segment=4,
        max_frames_per_segment=3,
        minimization=MinimizationOptions(max_iterations=500, max_evaluations=1000),
    )
    return system, prepared, root, tmp_path / "run", sources, config


def test_expanded_workflow_fingerprints_equal_evaluator(expanded_case):
    system, prepared, *_ = expanded_case
    evaluator = PCFFSinglePointEvaluator(system, prepared.native_result)
    result = evaluator.evaluate()
    assert workflow._identities(system, prepared) == (
        result.parameter_fingerprint,
        result.model_fingerprint,
    )


def test_public_expanded_start_after_verified_minimum(expanded_case, monkeypatch):
    _, _, root, run, sources, config = expanded_case
    minimize = workflow.minimize_geometry
    minima = []

    def capture(*args, **kwargs):
        minimum = minimize(*args, **kwargs)
        minima.append(minimum)
        print(
            "minimum:",
            minimum.converged,
            minimum.final_evaluation_verified,
            minimum.iterations,
            minimum.evaluations,
            flush=True,
        )
        return minimum

    monkeypatch.setattr(workflow, "minimize_geometry", capture)
    status = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert minima[0].converged and minima[0].final_evaluation_verified
    assert status["status"] == "paused" and status["accepted_step"] == 2


@pytest.fixture(params=["historical", "expanded", "wilson"])
def model_case(request, bound, monkeypatch):  # noqa: F811
    from dataclasses import replace

    from test_pcff_automatic import explicit

    from island.forcefields.pcff import automatic, model, source

    system, old, _ = bound
    src = old.assignment.source
    if request.param == "historical":
        return system, old
    if request.param == "wilson":
        raw = src.raw.replace(
            b"#equivalence cff91",
            b"1 1 c= 12.011 C 3 synthetic\n#equivalence cff91\n1 1 c= c c c c c",
        )
        sha = source.digest(raw)
        monkeypatch.setitem(source.PIN, "sha256", sha)
        monkeypatch.setitem(automatic.PROFILE, "source_sha256", sha)
        monkeypatch.setitem(model.PROFILE, "frc_sha256", sha)
        src = source.PCFFSource(raw, sha)
        system = explicit([(0, 1)], [2, 2])
        key = (11, 16)
        system.topology.bonds[key] = replace(system.topology.bonds[key], order=2)
        for sid, xyz in zip(
            sorted(system.topology.sites),
            [
                (0, 0, 0),
                (1.4, 0.1, 0.02),
                (-0.6, 0.9, 0.1),
                (-0.6, -0.9, -0.1),
                (2, 1, 0.15),
                (2, -0.9, 0.02),
            ],
        ):
            system.coordinates.set(sid, xyz)
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", src.expected_sha256)
    typing = type_pcff_atoms(system, src, profile=expanded.PROFILE_NAME)
    charges = assign_automatic_pcff_charges(system, typing)
    spec = define_pcff_model(
        assign_pcff_parameters(system, typing, charges),
        special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1)),
    )
    assert spec.payload["model_definition_complete"]
    assert any(t["family"] == "wilson_out_of_plane" for t in spec.payload["terms"]) == (
        request.param == "wilson"
    )
    return system, spec


# Freeze the reviewed evaluation-settings contract independently of the helper.
HISTORICAL_SETTINGS = {
    "implementation": "island_pcff_singlepoint_v1",
    "compatibility_profile": "island_lammps_pcff_acyclic_cho_v1",
    "method": "NoCutoff",
    "periodic": False,
    "constraints": False,
    "switching": False,
    "tail_correction": False,
    "precision": "double",
    "integration_steps": 0,
    "force_sign": "-dE/dR",
    "mixing": "sixth_power",
    "coulomb_constant_nm": 138.935456264,
}


def test_offline_and_evaluator_identities_preserve_reviewed_contract(
    model_case, monkeypatch
):
    import builtins
    from copy import deepcopy

    import openmm

    from island.evaluation.models import fingerprint
    from island.evaluation.pcff_identity import pcff_evaluation_identity
    from island.forcefields.pcff.charges import identity

    system, spec = model_case
    data = spec.payload
    expected_settings = dict(HISTORICAL_SETTINGS)
    if data["schema"] == "island_pcff_source_model_v1":
        expected_settings.update(
            implementation="island_pcff_source_singlepoint_v1",
            compatibility_profile="island_lammps_pcff_source_graph_v1",
        )
    expected_parameter = identity(data)
    expected_model = fingerprint(
        {"specification": expected_parameter, "settings": expected_settings}
    )
    evaluator = PCFFSinglePointEvaluator(system, spec)
    first = evaluator.evaluate()
    with evaluator.open_session() as session:
        assert session.evaluate() == first
        assert session.evaluate_fresh() == first
    assert (first.parameter_fingerprint, first.model_fingerprint) == (
        expected_parameter,
        expected_model,
    )
    assert dict(first.settings) == expected_settings
    import_ = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] == "openmm":
            pytest.fail("offline identity imported OpenMM")
        return import_(name, *args, **kwargs)

    def forbidden(*args, **kwargs):
        pytest.fail("offline identity constructed an OpenMM resource")

    with monkeypatch.context() as mp:
        mp.setattr(builtins, "__import__", blocked)
        mp.setattr(openmm, "Context", forbidden)
        mp.setattr(openmm, "System", forbidden)
        settings, parameter, model = pcff_evaluation_identity(spec, system=system)
        assert (settings, parameter, model) == (
            expected_settings,
            expected_parameter,
            expected_model,
        )
        settings["implementation"] = "caller mutation"
        moved = deepcopy(system)
        moved.coordinates.translate([1, 2, 3])
        assert pcff_evaluation_identity(spec, system=moved)[0] == expected_settings
        prepared = adopt_forcefield(system, "pcff", spec)
        assert workflow._identities(system, prepared) == (parameter, model)


# Synthetic source pins are overridden ONLY in these labelled software-test children.
CHILD_PREFIX = """
from pathlib import Path
import sys
from island.forcefields.pcff import source, automatic, model, expanded
raw=Path(sys.argv[2]).read_bytes(); sha=source.digest(raw)
source.PIN['sha256']=sha
for profile,key in ((automatic.PROFILE,'source_sha256'),(expanded.PROFILE,'source_sha256'),(model.PROFILE,'frc_sha256')):
    profile[key]=sha
from island.forcefields import PreparedForceFieldSources
sources=PreparedForceFieldSources(pcff_frc=Path(sys.argv[2]))
from island.workflows import prepared_workflow_status,read_prepared_workflow_frames,resume_prepared_workflow
"""


def test_relocation_separate_process_and_offline_inspection(expanded_case, monkeypatch):
    import shutil
    import subprocess
    import sys

    from island.workflows import prepared_workflow_status, read_prepared_workflow_frames
    from island.workflows.storage import read_json

    _, prepared, root, run, sources, config = expanded_case
    started = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert started["status"] == "paused" and started["accepted_step"] == 2
    moved = run.with_name("relocated")
    shutil.move(run, moved)
    code = (
        CHILD_PREFIX
        + """
import island.workflows.prepared as workflow

def forbidden(*a, **kw): raise AssertionError('setup repeated on resume')
workflow.initialize_velocities=workflow.minimize_geometry=forbidden
m=resume_prepared_workflow(sys.argv[1],sources=sources)
assert m['status']=='completed' and m['accepted_step']==4
"""
    )
    completed = subprocess.run(
        [sys.executable, "-c", code, str(moved), str(sources.pcff_frc)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    before = {p.name: p.read_bytes() for p in moved.iterdir() if p.is_file()}
    manifest = prepared_workflow_status(moved, sources=sources)
    assert manifest["prepared_identity"] == prepared.identity
    assert [f.step for f in read_prepared_workflow_frames(moved, sources=sources)] == [
        0,
        1,
        2,
        3,
        4,
    ]
    # New interpreter: dependencies genuinely unavailable, not just cached modules.
    code = (
        """
import sys
class Block:
    def find_spec(self,name,*args):
        if name.split('.')[0] in {'openmm','rdkit','parmed','foyer','scipy'}:
            raise AssertionError('offline inspection imported '+name)
sys.meta_path.insert(0,Block())
"""
        + CHILD_PREFIX
        + """
a=prepared_workflow_status(sys.argv[1],sources=sources)
assert a['status']=='completed'
assert len(read_prepared_workflow_frames(sys.argv[1],sources=sources))==5
assert resume_prepared_workflow(sys.argv[1],sources=sources)==a
"""
    )
    child = subprocess.run(
        [sys.executable, "-c", code, str(moved), str(sources.pcff_frc)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stderr
    assert before == {p.name: p.read_bytes() for p in moved.iterdir() if p.is_file()}
    assert read_json(moved / "manifest.json")["payload"]["status"] == "completed"

    # Rechecksummed, internally valid minimum observations from the wrong model
    # must not become authoritative just because they use a historical setting.
    from island.exceptions import WorkflowError
    from island.workflows import storage

    stage = next(s for s in manifest["stages"] if s["stage"] == "minimization")
    path = moved / stage["record"]
    p = storage.decode(read_json(path))
    for key in ("initial_evaluation", "final_evaluation"):
        from island.evaluation.models import fingerprint

        p[key]["settings"] = dict(HISTORICAL_SETTINGS)
        p[key]["model_fingerprint"] = fingerprint(
            {
                "specification": p[key]["parameter_fingerprint"],
                "settings": p[key]["settings"],
            }
        )
    raw = storage.json_bytes(storage.encode(p))
    workflow._minimum_record(p).validate_integrity()
    path.write_bytes(raw)
    manifest["files"][stage["record"]] = {
        "bytes": len(raw),
        "sha256": storage.checksum(raw),
    }
    workflow._manifest(moved, manifest)
    for op in (
        prepared_workflow_status,
        read_prepared_workflow_frames,
        workflow.resume_prepared_workflow,
    ):
        with pytest.raises(WorkflowError, match="Minimum diagnostic model mismatch"):
            op(moved, sources=sources)


def test_wrong_specification_settings_rejected_offline(model_case):
    from island.charge_references.records import pack, unpack
    from island.evaluation.pcff_identity import pcff_evaluation_identity
    from island.exceptions import EvaluationInputError

    system, spec = model_case
    p = unpack(spec.json_text)
    p["compatibility_profile"]["name"] = "unrelated source interpretation"
    with pytest.raises(EvaluationInputError, match="Contradictory"):
        pcff_evaluation_identity(type(spec)(pack(p), spec.assignment), system=system)


@pytest.mark.parametrize(
    "failure", ["minimum_budget", "minimum_final", "backend", "startup", "final"]
)
def test_legitimate_expanded_failures_remain_inspectable(
    expanded_case, monkeypatch, failure
):
    from dataclasses import replace

    from island.exceptions import EvaluationError
    from island.workflows import prepared_workflow_status, read_prepared_workflow_frames

    _, _, root, run, sources, config = expanded_case
    if failure == "minimum_budget":
        config = replace(config, minimization=MinimizationOptions(max_evaluations=1))
    elif failure == "minimum_final":
        from island.evaluation.pcff_session import PCFFEvaluationSession

        def fail(*a, **k):
            raise EvaluationError("synthetic final minimum verification failure")

        monkeypatch.setattr(PCFFEvaluationSession, "evaluate_fresh", fail)
    else:
        execute = workflow.run_dynamics_segment

        def injected(*args, **kwargs):
            session = args[1]
            fresh = session.evaluate_fresh
            calls = []

            def fail(*a, **k):
                raise EvaluationError("synthetic dynamics failure")

            def checked(*a, **k):
                calls.append(1)
                result = fresh(*a, **k)
                if len(calls) == 2:
                    components = dict(result.energy_components)
                    components[next(iter(components))] += 1.0
                    return replace(
                        result,
                        potential_energy=result.potential_energy + 1.0,
                        energy_components=components,
                    )
                return result

            monkeypatch.setattr(
                session,
                "evaluate" if failure == "backend" else "evaluate_fresh",
                checked if failure == "final" else fail,
            )
            return execute(*args, **kwargs)

        monkeypatch.setattr(workflow, "run_dynamics_segment", injected)
    m = start_prepared_bundle_workflow(root, run, config, sources=sources)
    assert m["status"] == "stage_failed" and m["accepted_step"] == 0
    assert m["checkpoint"] is None
    assert prepared_workflow_status(run, sources=sources) == m
    assert read_prepared_workflow_frames(run, sources=sources) == ()
