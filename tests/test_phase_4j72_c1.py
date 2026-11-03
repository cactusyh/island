"""C1 software controls: all executable/module probes are synthetic only."""

import copy
import json
import lzma
import os
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

import pytest

from island.pcff_completion import CompletionError, pin, signed, write
from scripts.phase_4j72_c1_fence import ScopeDenied, assert_installed, load_policy
from scripts.phase_4j72_c1_integrity import (
    BASE,
    BRANCH,
    DECLARATIONS,
    PARENT,
    ROOT,
    SCIENTIFIC_FILES,
    exclusion,
    quarantine,
    read,
    reconstruct,
    reject_quarantine,
    verify_events,
)
from scripts.validate_phase_4j72_c1_receipt import replay, semantics

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def q():
    return quarantine(REPO)


def child(code, *arguments, **kwargs):
    return subprocess.run(
        [sys.executable, "-c", code, *map(str, arguments)],
        capture_output=True,
        text=True,
        check=False,
        **kwargs,
    )


def fake(path, source):
    path.write_text("#!/usr/bin/env python3\n# ISLAND_C1_SYNTHETIC_ONLY\n" + source)
    path.chmod(0o755)
    return path


def test_fence_installed_before_pytest(tmp_path):
    assert_installed()
    policy_path = Path(os.environ["ISLAND_C1_POLICY"])
    policy = read(policy_path)
    (tmp_path / "policy.json").write_bytes(policy_path.read_bytes())
    data = Path(policy["events"]).read_bytes()
    ledger = tmp_path / "events.jsonl"
    ledger.write_bytes(data)
    raw = verify_events([tmp_path])
    (tmp_path / "events.jsonl.xz").write_bytes(lzma.compress(data))
    ledger.unlink()
    assert verify_events([tmp_path]) == raw
    assert raw["C1_new_native_QM_attempts"] == raw["C1_MD_steps"] == 0


def test_synthetic_parent_only_patch_does_not_propagate(tmp_path):
    """Negative control recreates J72's inheritance bug, with fake code only.

    The marker script is explicitly executed by Python. It contains no QM/MD
    library or external executable. It demonstrates the missing Python patch,
    not permission to run a forbidden native name under the corrected policy.
    """
    marker = tmp_path / "negative-control-marker"
    script = fake(
        tmp_path / "synthetic_backend.py",
        f"from pathlib import Path\nPath({str(marker)!r}).write_text('fake-only')\n",
    )
    result = child(
        "import subprocess,sys\n"
        "saved=subprocess.run\n"
        "subprocess.run=lambda *a,**k: (_ for _ in ()).throw(RuntimeError('parent-only'))\n"
        "code='import subprocess,sys; subprocess.run([sys.executable,sys.argv[1]],check=True)'\n"
        "result=saved([sys.executable,'-c',code,sys.argv[1]],check=True)\n",
        script,
    )
    assert result.returncode == 0, result.stderr
    assert marker.read_text() == "fake-only"


@pytest.mark.parametrize(
    "name",
    [
        "sqm",
        "antechamber",
        "sander",
        "pmemd",
        "lammps",
        "lmp_serial",
        "msi2lmp",
        "orca",
    ],
)
def test_forbidden_fake_executable_never_launches(tmp_path, name):
    marker = tmp_path / "ran"
    command = [
        str(fake(tmp_path / name, f"open({str(marker)!r},'w').write('ran')\n")),
        "--original-argument",
    ]
    before = command.copy()
    with pytest.raises(ScopeDenied):
        subprocess.run(command, check=False)
    assert command == before and not marker.exists()
    events = [
        json.loads(line)
        for line in Path(read(os.environ["ISLAND_C1_POLICY"])["events"])
        .read_text()
        .splitlines()
    ]
    record = next(e for e in reversed(events) if e.get("command") == command)
    assert record["executable_identity"]["synthetic_fixture"] is True
    assert record["executable_identity"]["sha256"] == pin(command[0])
    assert record["parent_process"] > 0 and record["timestamp_utc"]
    assert record["environment_policy"]["ISLAND_C1_POLICY"]


def test_child_inherits_policy_even_if_caller_clears_environment(tmp_path):
    marker = tmp_path / "ran"
    exe = fake(tmp_path / "sqm", f"open({str(marker)!r},'w').write('ran')\n")
    code = (
        "import os,subprocess,sys\n"
        "assert os.environ['ISLAND_C1_POLICY']\n"
        "try: subprocess.run([sys.argv[1]])\n"
        "except RuntimeError as e: assert 'J72-C1 scope' in str(e)\n"
        "else: raise AssertionError('launched')\n"
    )
    runtimes = [Path(sys.executable)]
    legacy = REPO.parent / "miniforge3/envs/md_rg/bin/python"
    if legacy.is_file():
        runtimes.append(legacy)
    for runtime in runtimes:
        result = subprocess.run(
            [str(runtime), "-c", code, str(exe)],
            env={"PATH": os.environ["PATH"]},
            capture_output=True,
            text=True,
            check=False,
        )
        assert result.returncode == 0, result.stderr
    assert not marker.exists()
    suite = tmp_path / "test_synthetic_scope.py"
    suite.write_text(
        "import subprocess\n"
        "def test_swallowed_rejection_must_be_scope_skip():\n"
        f" try: subprocess.run([{str(exe)!r}],check=False)\n"
        " except Exception: pass\n"
        " raise AssertionError('a scope denial must bypass application handlers')\n"
    )
    xml = tmp_path / "child-pytest.xml"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "-p",
            "scripts.phase_4j72_c1_test_policy",
            str(suite),
            "-q",
            "--junitxml",
            str(xml),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    cases = list(ET.parse(xml).iter("testcase"))
    assert len(cases) == 1 and cases[0].find("skipped") is not None
    assert not marker.exists()


@pytest.mark.parametrize("api", ["execv", "posix_spawn", "ctypes_execve"])
def test_kernel_fence_blocks_python_hook_bypass(tmp_path, api):
    marker = tmp_path / "ran"
    exe = fake(
        tmp_path / "renamed_backend", f"open({str(marker)!r},'w').write('ran')\n"
    )
    operation = {
        "execv": "os.execv(sys.argv[1], [sys.argv[1], 'exact-argument'])",
        "posix_spawn": "p=os.posix_spawn(sys.argv[1],[sys.argv[1], 'exact-argument'],os.environ); _,s=os.waitpid(p,0); assert os.waitstatus_to_exitcode(s)==-9",
        "ctypes_execve": "import ctypes as C; a=(C.c_char_p*3)(os.fsencode(sys.argv[1]),b'exact-argument',None); e=(C.c_char_p*1)(None); C.CDLL(None).execve(os.fsencode(sys.argv[1]),a,e)",
    }[api]
    result = child("import os,sys\n" + operation, exe)
    assert result.returncode == (0 if api == "posix_spawn" else -9), result.stderr
    assert not marker.exists()
    events = [
        json.loads(line)
        for line in Path(read(os.environ["ISLAND_C1_POLICY"])["events"])
        .read_text()
        .splitlines()
    ]
    record = next(
        e for e in reversed(events) if e.get("command") == [str(exe), "exact-argument"]
    )
    assert record["action"] == "rejected" and record["layer"] == "kernel_exec"
    assert record["executable_identity"]["synthetic_fixture"]


@pytest.mark.parametrize("flag", ["-S", "-I", "-E", "-Es"])
def test_python_startup_bypass_rejected(flag):
    with pytest.raises(ScopeDenied, match="python_startup_bypass"):
        subprocess.run(
            [sys.executable, flag, "-c", "raise AssertionError('not reached')"],
            check=False,
        )


def test_kernel_rejects_python_with_erased_policy(tmp_path):
    marker = tmp_path / "ran"
    result = child(
        "import os,sys; os.execve(sys.executable,[sys.executable,'-c',sys.argv[1]],{})",
        f"open({str(marker)!r},'w').write('ran')",
    )
    assert result.returncode == -9 and not marker.exists()


@pytest.mark.parametrize("name", ["psi4", "pyscf", "lammps", "qcengine"])
def test_find_spec_is_readonly_import_cannot_execute_fake_module(tmp_path, name):
    marker = tmp_path / "import-ran"
    fake(tmp_path / (name + ".py"), f"open({str(marker)!r},'w').write('ran')\n")
    code = (
        "import importlib.util,sys\n"
        "sys.path.insert(0,sys.argv[1])\n"
        f"assert importlib.util.find_spec({name!r}) is not None\n"
        f"try: __import__({name!r})\n"
        "except RuntimeError as e: assert 'scientific_module_execution' in str(e)\n"
        "else: raise AssertionError('module executed')\n"
    )
    result = child(code, tmp_path)
    assert result.returncode == 0, result.stderr
    assert not marker.exists()


def test_synthetic_openmm_constructor_and_integrator_step_blocked(tmp_path):
    marker = tmp_path / "integrator-ran"
    fake(
        tmp_path / "openmm.py",
        (
            "class VerletIntegrator:\n"
            f" def __init__(self,*a): open({str(marker)!r},'w').write('constructor')\n"
            f" def step(self,n): open({str(marker)!r},'w').write('step')\n"
        ),
    )
    result = child(
        "import sys;sys.path.insert(0,sys.argv[1]);import openmm\n"
        "for op in (lambda:openmm.VerletIntegrator(0.001),lambda:openmm.VerletIntegrator.step(None,1)):\n"
        " try: op()\n"
        " except RuntimeError as e: assert 'J72-C1 scope' in str(e)\n"
        " else: raise AssertionError('executed')\n",
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert not marker.exists()


def test_synthetic_integrator_restore_blocked_but_system_xml_allowed(tmp_path):
    marker = tmp_path / "restore-ran"
    fake(
        tmp_path / "openmm.py",
        (
            "def XmlSerializer_deserialize(xml):\n"
            f" open({str(marker)!r},'w').write(xml)\n"
            " return 'fake-system-only'\n"
        ),
    )
    result = child(
        "import sys;sys.path.insert(0,sys.argv[1]);import openmm\n"
        "try: openmm.XmlSerializer_deserialize('<Integrator type=\"fake\"/>')\n"
        "except RuntimeError as e: assert 'MD_integrator_restore' in str(e)\n"
        "else: raise AssertionError('restored')\n",
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert not marker.exists()
    result = child(
        "import sys;sys.path.insert(0,sys.argv[1]);import openmm\n"
        "assert openmm.XmlSerializer_deserialize('<System/>')=='fake-system-only'",
        tmp_path,
    )
    assert result.returncode == 0, result.stderr
    assert marker.read_text() == "<System/>"


def test_bad_policy_stops_before_fixture(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("{}")
    with pytest.raises(RuntimeError, match="hash mismatch"):
        load_policy(path, "0" * 64)
    result = child(
        "import os\nfrom scripts.phase_4j72_c1_fence import load_policy\n"
        "try: load_policy(__import__('sys').argv[1],'0'*64)\n"
        "except BaseException: os._exit(97)\n"
        "raise AssertionError('scientific fixture must not start')",
        path,
    )
    assert result.returncode == 97


def test_quarantine_counts_and_original_byte_identities(q):
    assert len(q["members"]) == 192
    assert [a["case"] for a in q["unintended_attempts"]] == [2, 3, 6, 7]
    assert q["verified_MD_steps"] == 46 and q["possible_MD_steps"] == 48
    assert all(
        a["process_id"] is None and not a["authorized"] and not a["scientific_use"]
        for a in q["unintended_attempts"]
    )
    assert (
        q["original_receipt_sha256"]
        == "66b56babcb7aae17528e314d2cf29557015612a03a96f26838f157da3db14112"
    )


@pytest.mark.parametrize(
    "kind", ["source_path", "bundle_hash", "HMDS_QM", "MD_frames", "derived_source"]
)
def test_quarantine_scientific_reference_rejected_even_resigned(q, kind):
    path = next(n for n in q["members"] if n.endswith("sqm.out"))
    value = {
        "source_path": {"source": "docs/evidence/phase_4j72/scope-deviation/" + path},
        "bundle_hash": {"parameter_bundle_hashes": [q["members"][path]["sha256"]]},
        "HMDS_QM": {"HMDS_QM": [{"producer": path}]},
        "MD_frames": {
            "scientific_trajectory_data": [
                {"source": path.replace("sqm.out", "generation-00000000/record.json")}
            ]
        },
        "derived_source": {"raw_sha256": q["archive_sha256"]},
    }[kind]
    with pytest.raises(CompletionError, match="quarantine"):
        reject_quarantine(signed(value), q)


@pytest.mark.parametrize("index", range(len(SCIENTIFIC_FILES)))
def test_every_scientific_role_excludes_injected_quarantine(q, index):
    docs = {n: read(REPO / n) for n in SCIENTIFIC_FILES}
    name = SCIENTIFIC_FILES[index]
    docs[name]["injected_evidence"] = {"source": q["archive_path"]}
    docs[name] = signed({k: v for k, v in docs[name].items() if k != "identity"})
    with pytest.raises(CompletionError, match="quarantine"):
        exclusion(REPO, q, docs)


def receipt_pair():
    nested = {
        "offline": True,
        "unresolved_instances": 378,
        "duplicate_equivalences": 16,
        "unresolved_unequal_COMPASS_assignments": 2,
        **DECLARATIONS,
    }
    manifest = signed(
        {
            "schema": "island.j72-c1-manifest.v1",
            "base_commit": BASE,
            "implementation_files": {},
            "historical_files_sha256": "1" * 64,
            "artifacts": {},
            "scientific_scope": DECLARATIONS.copy(),
            "quarantine_manifest_sha256": "2" * 64,
        }
    )
    receipt = signed(
        {
            "schema": "island.j72-c1-receipt.v1",
            "base_commit": BASE,
            "required_J72": BASE,
            "branch": BRANCH,
            **DECLARATIONS,
            "manifest_sha256": "3" * 64,
            "quarantine_manifest_sha256": "2" * 64,
            "software_counts": {
                "tests": 1,
                "passed": 1,
                "skipped": 0,
                "failure": 0,
                "error": 0,
            },
            "terminal_unique_tests": 1,
            "historical_files_preserved": 1,
            "fence_verification": {},
            "offline_reconstruction": nested,
        }
    )
    return receipt, manifest


@pytest.mark.parametrize("key", list(DECLARATIONS))
@pytest.mark.parametrize("where", ["receipt", "manifest", "nested"])
def test_resigned_semantic_mutation_never_erases_scope_violation(key, where, tmp_path):
    receipt, manifest = receipt_pair()
    target = (
        manifest["scientific_scope"]
        if where == "manifest"
        else receipt["offline_reconstruction"]
        if where == "nested"
        else receipt
    )
    old = target[key]
    target[key] = (
        not old if type(old) is bool else old + 1 if type(old) is int else "promoted"
    )
    manifest = signed({k: v for k, v in manifest.items() if k != "identity"})
    write(tmp_path / "manifest.json", manifest)
    receipt["manifest_sha256"] = pin(tmp_path / "manifest.json")
    receipt = signed({k: v for k, v in receipt.items() if k != "identity"})
    write(tmp_path / "receipt.json", receipt)
    assert read(tmp_path / "receipt.json")["manifest_sha256"] == pin(
        tmp_path / "manifest.json"
    )
    with pytest.raises(CompletionError):
        semantics(receipt, manifest)


def test_receipt_exact_types_not_bool_as_int():
    receipt, manifest = receipt_pair()
    semantics(receipt, manifest)
    receipt["unintended_external_QM_attempts"] = True
    receipt = signed({k: v for k, v in receipt.items() if k != "identity"})
    with pytest.raises(CompletionError):
        semantics(receipt, manifest)


def test_scientific_result_producers_and_nonactivation_preserved(q):
    report = exclusion(REPO, q)
    assert report["scientific_evidence_contaminated"] is False
    assert report["J70_unresolved_instances_preserved"] == 378
    assert report["equal_duplicate_equivalences_preserved"] == 16
    assert report["unequal_COMPASS_assignments_unresolved"] == 2
    original = read(REPO / PARENT / "receipt.json")
    assert original["new_QM"] == 4 and original["zero_execution_scope_met"] is False
    assert original["reserved_validation_activated"] is False


def test_fresh_process_relocation_offline_replay(tmp_path):
    # Minimal replay repository relocation: no Git tree, external literature
    # or executable backends. The inherited supervisor keeps its separate
    # original policy; scientific replay inputs all come from the relocated tree.
    relocated_repo = tmp_path / "relocated-repository"
    relocated_repo.mkdir()
    for name in ("src", "scripts"):
        shutil.copytree(
            REPO / name,
            relocated_repo / name,
            ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
        )
    parents = read(REPO / ROOT / "frozen-inputs.json")["inputs"]
    extra = (
        PARENT / "scope-deviation/execution-deviation.json",
        PARENT / "scope-deviation/quarantined-j31-regressions.tar.xz",
        PARENT / "execution-scope-report.json",
        PARENT / "validation/start_suite.py",
        PARENT / "validation/suite-results.json",
        PARENT / "validation/start-suite.log",
    )
    for name in [*parents, *map(str, extra)]:
        target = relocated_repo / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((REPO / name).read_bytes())
    relocated = tmp_path / "corrective-evidence"
    relocated.mkdir()
    for name in ("frozen-inputs.json", "historical-files.json"):
        (relocated / name).write_bytes((REPO / ROOT / name).read_bytes())
    for name, value in reconstruct(REPO, relocated).items():
        write(relocated / name, value)
    result = child(
        "import json,sys\nfrom pathlib import Path\n"
        "from scripts.validate_phase_4j72_c1_receipt import replay\n"
        "r=replay(Path(sys.argv[1]),Path(sys.argv[2]));assert r['unintended_external_QM_attempts']==4\n"
        "assert not r['zero_execution_scope_met'];print(json.dumps(r,sort_keys=True))\n",
        relocated_repo,
        relocated,
        cwd=relocated_repo,
        env={
            **os.environ,
            "PYTHONPATH": str(relocated_repo) + ":" + str(relocated_repo / "src"),
        },
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["HMDS_new_QM"] == 0
    altered = read(relocated / "scope-correction.json")
    altered["zero_execution_scope_met"] = True
    write(
        relocated / "scope-correction.json",
        signed({k: v for k, v in altered.items() if k != "identity"}),
    )
    with pytest.raises(CompletionError):
        replay(REPO, relocated)


def test_quarantine_re_signed_manifest_mutation_reconstructs(tmp_path, q):
    value = copy.deepcopy(q)
    value["unintended_attempts"] = []
    value = signed({k: v for k, v in value.items() if k != "identity"})
    for name in ("frozen-inputs.json", "historical-files.json"):
        (tmp_path / name).write_bytes((REPO / ROOT / name).read_bytes())
    for name, report in reconstruct(REPO, tmp_path).items():
        write(tmp_path / name, report)
    write(tmp_path / "quarantine-manifest.json", value)
    with pytest.raises(CompletionError, match="reconstruction:quarantine-manifest"):
        replay(REPO, tmp_path)
