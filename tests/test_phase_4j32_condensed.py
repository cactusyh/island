"""J32 software negative controls are never scientific acceptance evidence."""

import json
import os
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.condensed import (
    statistics,
    validate_assignment,
    validate_execution,
    validate_fixture,
)
from island.analysis.network import _signed
from island.exceptions import AnalysisError, IslandError
from scripts.phase_4j32_contract import validate_case, validate_matrix
from scripts.phase_4j32_real_fixture import build_fixture, validate_complete_chain


@pytest.fixture(scope="module")
def fixture():
    return build_fixture("linear")


@pytest.fixture(scope="module")
def evidence(fixture):
    return json.loads(json.dumps(fixture.evidence()))


def test_complete_chains_and_packing_replay(evidence):
    packed = validate_fixture(
        evidence,
        {"system": evidence["packed_system"], "final_graph": evidence["final_graph"]},
    )
    assert len(packed.system.topology.sites) == 90
    assert evidence["setup"]["chain_count"] == 3
    assert evidence["setup"]["component_count"] == 3
    assert evidence["setup"]["minimum_interchain_distance_angstrom"] >= 1.8


@pytest.mark.parametrize("mutation", ["incomplete", "disconnected"])
def test_bad_chain_explicitly_rejected(fixture, mutation):
    chain = fixture.chains[0].copy()
    if mutation == "incomplete":
        atom = next(i for i, s in chain.topology.sites.items() if s.element == "H")
        chain.topology.remove_site(atom)
        chain.coordinates._positions.pop(atom)
    else:
        bond = next(iter(chain.topology.bonds.values()))
        chain.topology.remove_bond(bond.site1, bond.site2)
    with pytest.raises((ValueError, IslandError)):
        validate_complete_chain(chain)


@pytest.mark.parametrize(
    "mutation", ["graph", "membership", "mapping", "box", "distance"]
)
def test_stale_packing_inputs_rejected(evidence, mutation):
    changed = deepcopy(evidence)
    if mutation == "graph":
        changed["final_graph"]["identity"] = "0" * 64
    elif mutation == "membership":
        changed["final_graph"]["molecule_membership"]["1"] = "wrong"
    elif mutation == "mapping":
        changed["packing_plan"]["identity"] = "0" * 64
    elif mutation == "box":
        changed["packed_system"]["box"] = None
    else:
        changed["setup"]["minimum_interchain_distance_angstrom"] = 0
    with pytest.raises((ValueError, IslandError)):
        validate_fixture(
            changed,
            {
                "system": evidence["packed_system"],
                "final_graph": evidence["final_graph"],
            },
        )


def test_insufficient_frames_rejected():
    with pytest.raises(AnalysisError, match="Insufficient"):
        statistics([1.0] * 4)


@pytest.fixture(scope="module")
def retained(tmp_path_factory):
    from scripts.phase_4j32_receipt import extract_archive

    candidate = os.environ.get("ISLAND_J32_TEST_RUN")
    if candidate:
        return Path(candidate)
    root = tmp_path_factory.mktemp("j32-retained")
    extract_archive(Path("docs/evidence/phase_4j32/artifacts/acceptance.tar.gz"), root)
    return root


@pytest.fixture(scope="module")
def matrix(retained):
    return json.loads((retained / "matrix.json").read_text())


def test_real_required_matrix(matrix):
    validate_matrix(matrix)


@pytest.mark.parametrize(
    "mutation",
    [
        "all_blocked",
        "success_blocked",
        "missing",
        "pcff_success",
        "pcff_stage",
        "pcff_reason",
        "case_identity",
        "execution",
        "reconstruction",
    ],
)
def test_resigned_matrix_mutations_rejected(matrix, mutation):
    m = deepcopy(matrix)
    cases = m["cases"]
    row = cases["oplsaa/linear"]["result"]
    pcff = cases["pcff/network"]["result"]
    if mutation == "all_blocked":
        for c in cases.values():
            c["result"].update(status="blocked", dynamics_executed=False)
    elif mutation == "success_blocked":
        row.update(status="blocked", dynamics_executed=False)
    elif mutation == "missing":
        del cases["gaff2/network"]
    elif mutation == "pcff_success":
        pcff.update(status="complete", dynamics_executed=True)
    elif mutation == "pcff_stage":
        pcff["diagnostic"]["stage"] = "source_loading"
    elif mutation == "pcff_reason":
        pcff["diagnostic"]["message"] = "Unrelated typing failure"
    elif mutation == "case_identity":
        row["fixture"] = "network"
    elif mutation == "execution":
        row.pop("comparison")
    else:
        row.pop("fresh_process_reconstruction")
    for c in cases.values():
        c["result"] = _signed(c["result"])
    with pytest.raises(ValueError):
        validate_matrix(_signed(m))


@pytest.mark.parametrize("mutation", ["charge", "family", "counts", "verification"])
def test_assignment_mutation_rejected(retained, mutation):
    bundle = retained / "oplsaa-linear/bundle"
    backend = json.loads((bundle / "workflow/manifest.json").read_text())["payload"][
        "backend"
    ]
    assignment = json.loads((bundle / "assignment.json").read_text())
    validate_assignment(assignment, backend)
    if mutation == "charge":
        assignment["atom_types_and_charges"][0]["charge"] += 0.1
    elif mutation == "family":
        assignment["periodic_policy"]["family"] = "GAFF2"
    elif mutation == "counts":
        assignment["angles"] += 1
    else:
        assignment["independent_source_verification"] = None
    with pytest.raises(AnalysisError):
        validate_assignment(assignment, backend)


@pytest.mark.parametrize(
    "injection",
    [
        "oserror",
        "valueerror",
        "islanderror",
        "dependency",
        "executable",
        "source_missing",
        "source_corrupt",
    ],
)
def test_actual_cli_nonzero_preserves_diagnostics(tmp_path, injection):
    output = tmp_path / "diagnostic"
    code = """
import runpy, sys
from pathlib import Path
import scripts.phase_4j32_assignment as assignment
from island.exceptions import IslandError
kind, output = sys.argv[1:]
if kind.startswith('source_'):
    assignment.XML = Path(output).parent / 'bad.xml'
    if kind == 'source_corrupt': assignment.XML.write_text('corrupted software control')
else:
    def fail(*args, **kwargs):
        error = {'oserror': OSError, 'valueerror': ValueError, 'islanderror': IslandError,
                 'dependency': ImportError, 'executable': FileNotFoundError}[kind]
        raise error('injected software negative control: ' + kind)
    assignment.prepare = fail
sys.argv = ['validate_phase_4j32', '--output', output, '--family',
            'oplsaa' if kind.startswith('source_') else 'gaff2', '--kind', 'network']
runpy.run_module('scripts.validate_phase_4j32', run_name='__main__')
"""
    p = subprocess.run(
        [sys.executable, "-c", code, injection, str(output)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode != 0, p.stdout + p.stderr
    result = json.loads((output / "result.json").read_text())
    assert result["acceptance_met"] is False
    assert result["status"] == "failed"
    assert result["diagnostic"]["stage"] == "preparation"
    assert (output / "failure.log").stat().st_size > 0
    with pytest.raises(ValueError):
        validate_case(result["family"] + "/network", result)


def test_execution_claim_bindings(retained):
    from island.periodic_md import load_periodic_md_workflow

    bundle = retained / "oplsaa-linear/bundle"
    workflow = load_periodic_md_workflow(bundle / "workflow").payload
    execution = json.loads((bundle / "execution.json").read_text())
    validate_execution(execution, workflow)
    for key, value in [
        ("short_nvt_duration_ps", 123),
        ("engine", {}),
        ("final_checkpoint_sha256", "0" * 64),
        ("glass_transition", 300),
        ("mechanical_properties", 1),
    ]:
        changed = {**execution, key: value}
        with pytest.raises(AnalysisError):
            validate_execution(changed, workflow)


@pytest.mark.parametrize("outcome", ["unexpected_success", "wrong_reason"])
def test_pcff_cli_rejects_unexpected_outcome(tmp_path, outcome):
    code = """
import runpy,sys
import scripts.phase_4j32_assignment as assignment
from island.exceptions import PCFFError
outcome, output = sys.argv[1:]
def preparation(*args, **kwargs):
    if outcome == 'wrong_reason': raise PCFFError('unrelated software typing control')
    return object()
assignment.prepare = preparation
sys.argv = ['validate_phase_4j32', '--output', output, '--family', 'pcff', '--kind', 'network']
runpy.run_module('scripts.validate_phase_4j32', run_name='__main__')
"""
    target = tmp_path / "negative"
    p = subprocess.run(
        [sys.executable, "-c", code, outcome, str(target)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert p.returncode != 0
    record = json.loads((target / "result.json").read_text())
    assert record["status"] == "failed"
    assert record["acceptance_met"] is False
    assert (target / "failure.log").exists()


def test_complete_network_and_crosslink_plan_replay():
    network = build_fixture("network")
    e = json.loads(json.dumps(network.evidence()))
    validate_fixture(e, {"system": e["packed_system"], "final_graph": e["final_graph"]})
    assert e["setup"]["atom_count"] == 63
    assert e["setup"]["component_count"] == 1
    assert e["setup"]["chain_count"] == 3
    assert e["setup"]["crosslink_count"] == 3
