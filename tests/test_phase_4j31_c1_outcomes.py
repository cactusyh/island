"""Injected software negative controls; none are scientific acceptance evidence.

All injection lives in a fresh test interpreter that monkeypatches callables.
The production CLI has no injection switches, environment hooks or bypass mode.
"""

import json
import os
import subprocess
import sys
from copy import deepcopy
from pathlib import Path

import pytest

from island.analysis.network import _signed
from island.exceptions import AnalysisError
from scripts.phase_4j31_contract import case_errors, require_matrix
from scripts.validate_phase_4j31_receipt import manifest_section, validate_structure

RECEIPT = Path("docs/evidence/phase_4j31.json")
MANIFEST = Path("docs/paper/reproducibility_manifest.json")

# Only this test-controlled child alters functions. Production accepts no hook.
CONTROL = r"""
import json,sys
from pathlib import Path
from types import SimpleNamespace
import scripts.validate_phase_4j31_network as driver
from island.exceptions import IslandError,PCFFError
mode,output=sys.argv[1],Path(sys.argv[2])
family='pcff' if mode.startswith('pcff') or mode.startswith('source') else 'gaff'
if mode in ('oserror','valueerror','islanderror','missing_executable','importerror','permission','pcff_wrong_reason','pcff_wrong_stage','pcff_success','execution_failure'):
    def prepare(*args,**kwargs):
        errors={'oserror':OSError,'valueerror':ValueError,'islanderror':IslandError,'missing_executable':FileNotFoundError,'importerror':ModuleNotFoundError,'permission':PermissionError,'pcff_wrong_reason':PCFFError,'pcff_wrong_stage':PCFFError}
        if mode in errors:
            message=driver.PCFF_MESSAGE if mode=='pcff_wrong_stage' else 'INJECTED SOFTWARE NEGATIVE CONTROL: '+mode
            raise errors[mode](message)
        if mode=='pcff_success':
            return object()
        return SimpleNamespace(identity='a'*64,config=SimpleNamespace(payload={'source':{'family':'GAFF','sha256':'a'*64}}),typed_graph=SimpleNamespace(identity='b'*64),graph_charges=SimpleNamespace(identity='c'*64),parameter_assignment=SimpleNamespace(identity='d'*64))
    driver.prepare_network_backend=prepare
    if mode=='execution_failure':
        driver.independent_assignment_check=lambda backend: {}
        def fail(*args):
            raise ValueError('INJECTED SOFTWARE NEGATIVE CONTROL: numerical/workflow failure')
        driver.execute=fail
elif mode=='missing_dependency':
    original=driver.importlib.util.find_spec
    driver.importlib.util.find_spec=lambda name,*a,**k: None if name=='parmed' else original(name,*a,**k)
elif mode in ('source_missing','source_corrupt'):
    driver.FRC=output.parent/'injected-pcff.frc'
    if mode=='source_corrupt':
        driver.FRC.write_text('INJECTED SOFTWARE NEGATIVE CONTROL: corrupted source')
if mode.startswith('matrix_'):
    # Dispatcher/exit-code tests reuse retained records, never new scientific claims.
    rows=json.loads(Path('docs/evidence/phase_4j31.json').read_text())['force_field_source_completeness']['matrix']
    for kind in ('branched','crosslinked'):
        path=output.parent/('native-pcff-'+kind)
        rows['pcff/'+kind]=driver.worker('pcff',kind,path,None)
    def child(command,**kwargs):
        family=command[command.index('--family')+1]
        kind=command[command.index('--kind')+1]
        path=Path(command[command.index('--output')+1])
        if mode=='matrix_interrupt' and family=='oplsaa':
            raise KeyboardInterrupt('INJECTED SOFTWARE NEGATIVE CONTROL: interrupted matrix')
        if mode=='matrix_missing' and family=='oplsaa':
            raise FileNotFoundError('INJECTED missing Foyer executable')
        row=rows[family+'/'+kind]
        if mode=='matrix_unmet' and family=='gaff':
            row={'family':family,'fixture':kind,'status':'blocked','reason':'INJECTED unrelated error','dynamics_executed':False}
        path.mkdir()
        (path/'result.json').write_text(json.dumps(row))
        return SimpleNamespace(returncode=7 if mode=='matrix_nonzero' and family=='gaff' else 0,stdout='SOFTWARE DISPATCH CONTROL\n',stderr='')
    driver.subprocess.run=child
    sys.argv=['control','--output',str(output)]
else:
    sys.argv=['control','--family',family,'--kind','branched','--output',str(output)]
raise SystemExit(driver.main())
"""


def run_control(tmp_path, mode):
    output = tmp_path / "output"
    process = subprocess.run(
        [sys.executable, "-c", CONTROL, mode, str(output)],
        env={
            **os.environ,
            "PYTHONDONTWRITEBYTECODE": "1",
            "OPENMM_CPU_THREADS": "1",
            "OMP_NUM_THREADS": "1",
        },
        capture_output=True,
        text=True,
        check=False,
    )
    return process, output


@pytest.mark.parametrize(
    "mode,classification",
    (
        ("oserror", "unexpected_preparation_failure"),
        ("valueerror", "unexpected_preparation_failure"),
        ("islanderror", "unexpected_preparation_failure"),
        ("missing_executable", "asset_unavailable"),
        ("importerror", "asset_unavailable"),
        ("permission", "asset_unavailable"),
        ("missing_dependency", "asset_unavailable"),
        ("source_missing", "asset_unavailable"),
        ("source_corrupt", "asset_unavailable"),
        ("pcff_wrong_reason", "unexpected_preparation_failure"),
        ("pcff_wrong_stage", "unexpected_preparation_failure"),
        ("pcff_success", "unexpected_preparation_failure"),
        ("execution_failure", "execution_failure"),
    ),
)
def test_cli_negative_controls_fail_and_preserve_diagnostics(
    tmp_path, mode, classification
):
    process, output = run_control(tmp_path, mode)
    assert process.returncode != 0, process.stdout + process.stderr
    result = json.loads((output / "result.json").read_text())
    assert result["diagnostic"]["classification"] == classification
    assert result["acceptance_met"] is False
    assert result["status"] != "complete"
    assert (output / "failure.log").read_text()
    assert "execution" not in result


@pytest.mark.parametrize(
    "mode,accepted",
    (
        ("matrix_valid", True),
        ("matrix_unmet", False),
        ("matrix_missing", False),
        ("matrix_nonzero", False),
    ),
)
def test_matrix_cli_gates_outcomes_and_keeps_partial_results(tmp_path, mode, accepted):
    process, output = run_control(tmp_path, mode)
    assert (process.returncode == 0) is accepted, process.stdout + process.stderr
    report = json.loads((output / "matrix.json").read_text())
    assert report["acceptance_met"] is accepted
    assert (
        "gaff2/crosslinked" in report["matrix"]
    )  # driver continued after earlier failures
    assert len(list(output.glob("*.log"))) == 8
    if accepted:
        require_matrix(report["matrix"])
    else:
        assert report["case_failures"]
        if mode != "matrix_nonzero":
            assert report["contract_errors"]


def real():
    return json.loads(RECEIPT.read_text()), json.loads(MANIFEST.read_text())


def mutate(receipt, defect):
    matrix = receipt["force_field_source_completeness"]["matrix"]
    if defect in ("all_blocked", "success_blocked"):
        for key, row in matrix.items():
            if defect == "all_blocked" or key == "gaff/branched":
                row.update(
                    status="blocked",
                    dynamics_executed=False,
                    reason="MUTATED SOFTWARE NEGATIVE CONTROL: unrelated failure",
                )
    elif defect == "missing_success":
        matrix.pop("gaff2/crosslinked")
    elif defect == "pcff_success":
        matrix["pcff/branched"] = deepcopy(matrix["gaff/branched"])
        matrix["pcff/branched"]["family"] = "pcff"
    elif defect == "case_identity":
        matrix["gaff/branched"]["fixture"] = "crosslinked"
    elif defect == "pcff_reason":
        matrix["pcff/branched"]["reason"] = "PCFFError: unrelated typing failure"
    elif defect == "pcff_stage":
        matrix["pcff/branched"]["diagnostic"] = {
            "classification": "expected_typing_rejection",
            "stage": "parameters",
        }
    elif defect == "pcff_contradiction":
        matrix["pcff/branched"]["execution"] = {"comparison": {"passed": True}}
    elif defect in ("execution", "reconstruction", "source", "graph"):
        row = matrix["gaff/branched"]
        if defect == "execution":
            row.pop("execution")
        elif defect == "reconstruction":
            row["execution"].pop("separate_process_reconstruction")
        elif defect == "source":
            row.pop("source")
        else:
            row["final_graph_identity"] = "0" * 64
    return _signed({k: v for k, v in receipt.items() if k != "identity"})


DEFECTS = (
    "all_blocked",
    "success_blocked",
    "missing_success",
    "pcff_success",
    "case_identity",
    "pcff_reason",
    "pcff_stage",
    "pcff_contradiction",
    "execution",
    "reconstruction",
    "source",
    "graph",
)


@pytest.mark.parametrize("defect", DEFECTS)
def test_consistently_resigned_real_receipt_mutations_fail(tmp_path, defect):
    receipt, manifest = real()
    receipt = mutate(receipt, defect)
    manifest["phase_4j31"] = manifest_section(receipt)
    with pytest.raises(AnalysisError):
        validate_structure(receipt, manifest)
    rp, mp = tmp_path / "receipt.json", tmp_path / "manifest.json"
    rp.write_text(json.dumps(receipt))
    mp.write_text(json.dumps(manifest))
    child = subprocess.run(
        [sys.executable, "scripts/validate_phase_4j31_receipt.py", str(rp), str(mp)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode != 0
    assert child.stdout + child.stderr


def test_valid_retained_contract_and_all_blocked_never_positive():
    receipt, manifest = real()
    require_matrix(
        receipt["force_field_source_completeness"]["matrix"], legacy_pcff=True
    )
    validate_structure(receipt, manifest)
    row = receipt["force_field_source_completeness"]["matrix"]["gaff/branched"]
    row.update(status="blocked", dynamics_executed=False, reason="unrelated")
    assert case_errors("gaff/branched", row)


def test_real_missing_amber_executables_fail_cli(tmp_path):
    output = tmp_path / "missing-assets"
    child = subprocess.run(
        [
            sys.executable,
            "-m",
            "scripts.validate_phase_4j31_network",
            "--family",
            "gaff",
            "--kind",
            "branched",
            "--output",
            str(output),
            "--amberhome",
            str(tmp_path / "absent-amber"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode != 0
    row = json.loads((output / "result.json").read_text())
    assert row["diagnostic"]["classification"] == "asset_unavailable"
    assert row["diagnostic"]["stage"] == "assets"
    assert "antechamber" in row["reason"]
    assert row["dynamics_executed"] is False


def test_original_full_validator_revalidates_retained_receipt():
    child = subprocess.run(
        [
            sys.executable,
            "scripts/validate_phase_4j31_receipt.py",
            str(RECEIPT),
            str(MANIFEST),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stdout + child.stderr


def test_read_only_relocation_and_fresh_reconstruction(tmp_path):
    child = subprocess.run(
        [
            sys.executable,
            "scripts/validate_phase_4j31_c1_evidence.py",
            "--output",
            str(tmp_path / "retained.json"),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert child.returncode == 0, child.stdout + child.stderr
    report = json.loads((tmp_path / "retained.json").read_text())
    assert len(report["revalidated_cases"]) == 6
    assert all(
        row["observables_match"] and row["fresh_process_reconstruction"] == "PASS"
        for row in report["revalidated_cases"].values()
    )
    assert report["new_dynamics_executed"] is False


def test_interrupted_matrix_preserves_partial_report(tmp_path):
    child, output = run_control(tmp_path, "matrix_interrupt")
    assert child.returncode != 0
    report = json.loads((output / "matrix.json").read_text())
    assert report["interruption"]["classification"] == "interrupted"
    assert (output / "interruption.log").read_text()
    assert set(report["matrix"]) == {"pcff/branched", "pcff/crosslinked"}
    assert report["acceptance_met"] is False
    assert report["contract_errors"]


@pytest.mark.parametrize("family", ("oplsaa", "gaff", "gaff2"))
@pytest.mark.parametrize("kind", ("branched", "crosslinked"))
@pytest.mark.parametrize("mode", ("blocked", "missing"))
def test_every_mandatory_success_is_required(family, kind, mode):
    receipt, _ = real()
    matrix = receipt["force_field_source_completeness"]["matrix"]
    name = family + "/" + kind
    if mode == "missing":
        matrix.pop(name)
    else:
        matrix[name] = {
            "family": family,
            "fixture": kind,
            "status": "blocked",
            "reason": "SOFTWARE NEGATIVE CONTROL: unrelated preparation failure",
            "dynamics_executed": False,
        }
    with pytest.raises(AnalysisError):
        require_matrix(matrix, legacy_pcff=True)
