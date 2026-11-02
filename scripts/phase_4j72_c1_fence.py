"""Opt-in C1 software-verification policy; never imported by the public API.

The Linux supervisor mediates exec syscalls independently of these Python
hooks. These hooks preserve the startup policy and fence in-process entry
points. This is a scope fence for the supported test suite, not a sandbox for
hostile Python code (which could modify its own memory or load arbitrary code).
"""

import hashlib
import importlib.abc
import importlib.machinery
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

POLICY_KEY = "ISLAND_C1_POLICY"
HASH_KEY = "ISLAND_C1_POLICY_SHA256"
BLOCKED_MODULES = {
    "psi4",
    "pyscf",
    "xtb",
    "tblite",
    "qcengine",
    "qcelemental",
    "lammps",
    "sander",
}
NATIVE = {
    "sqm",
    "antechamber",
    "sander",
    "pmemd",
    "psi4",
    "lammps",
    "msi2lmp",
    "g16",
    "g09",
    "orca",
    "nwchem",
    "cp2k",
    "xtb",
    "dftb+",
    "qchem",
    "namd",
    "namd2",
    "gmx",
    "grompp",
    "mdrun",
    "tleap",
    "parmchk2",
}
_policy = None
_policy_path = None
_policy_hash = None


class ScopeDenied(RuntimeError):
    """A requested operation was rejected before scientific execution."""


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_policy(path, expected):
    if sha(path) != expected:
        raise RuntimeError("C1 policy hash mismatch")
    policy = json.loads(Path(path).read_text())
    if policy["schema"] != "island.j72-c1-execution-policy.v1":
        raise RuntimeError("C1 policy schema")
    if sha(policy["runtime"]) != policy["runtime_sha256"]:
        raise RuntimeError("C1 startup implementation changed")
    if sha(Path(policy["bootstrap"]) / "sitecustomize.py") != policy["startup_sha256"]:
        raise RuntimeError("C1 startup bootstrap changed")
    return policy


def executable_identity(command, environment, cwd=None):
    """Resolve identity without running the requested binary, including aliases."""
    word = os.fsdecode(command)
    path = Path(cwd or os.getcwd()) / word if "/" in word else None
    if path is None:
        found = shutil.which(word, path=environment.get("PATH", os.defpath))
        path = Path(found) if found else None
    if path is None or not path.is_file():
        return {
            "requested": word,
            "resolved": None,
            "sha256": None,
            "synthetic_fixture": False,
        }
    path = path.resolve()
    stat = path.stat()
    with path.open("rb") as stream:
        synthetic = b"ISLAND_C1_SYNTHETIC_ONLY" in stream.read(256)
    return {
        "requested": word,
        "resolved": str(path),
        "sha256": sha(path),
        "device": stat.st_dev,
        "inode": stat.st_ino,
        "synthetic_fixture": synthetic,
    }


def append_event(policy, event):
    record = {
        "schema": "island.j72-c1-scope-event.v1",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        **event,
    }
    data = (json.dumps(record, sort_keys=True) + "\n").encode()
    fd = os.open(policy["events"], os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
    try:
        # One append, shared by parent and descendants. Failure is fatal.
        if os.write(fd, data) != len(data):
            raise RuntimeError("incomplete C1 event record")
    finally:
        os.close(fd)


def environment_policy(environment):
    return {
        k: environment.get(k)
        for k in (POLICY_KEY, HASH_KEY, "PYTHONPATH", "PYTHONHOME")
    }


def command_reason(argv, environment, policy, identity):
    if not isinstance(argv, (list, tuple)) or not argv:
        return "shell_or_unstructured_command"
    words = [os.fsdecode(w) for w in argv]
    name = Path(words[0]).name.lower()
    if name in NATIVE or name.startswith(("lmp", "pmemd", "sander", "cp2k", "lammps")):
        return "native_scientific_entry_point"
    resolved = identity.get("resolved")
    entry = policy["executables"].get(resolved)
    if not entry or entry["sha256"] != identity.get("sha256"):
        return "executable_not_in_frozen_allowlist"
    if entry["kind"] == "python":
        # Combined switches (e.g. -Es) also bypass PYTHONPATH/sitecustomize.
        for word in words[1:]:
            if word in ("-c", "-m") or not word.startswith("-"):
                break
            if word.startswith("--") or any(c in word[1:] for c in "SEI"):
                return "python_startup_bypass"
        if (
            environment.get(POLICY_KEY) != policy["policy_path"]
            or environment.get(HASH_KEY) != policy["policy_sha256"]
            or environment.get("PYTHONPATH", "").split(os.pathsep)[0]
            != policy["bootstrap"]
            or environment.get("PYTHONHOME")
        ):
            return "python_policy_not_inherited"
        text = " ".join(words)
        if (
            re.search(r"scripts\.(?:validate_)?phase_4j(3[1-9]|40|41)(\b|_)", text)
            and "_receipt" not in text
            and "scripts.phase_4j37_reference" not in text
            and not any(
                w in words
                for w in ("--offline", "--audit-only", "--help", "reconstruct")
            )
        ):
            return "historical_live_workflow_outside_scope"
    elif entry["kind"] == "git":
        if len(words) < 2 or words[1] not in {
            "ls-tree",
            "show",
            "rev-parse",
            "merge-base",
            "diff",
            "status",
            "cat-file",
        }:
            return "git_operation_outside_readonly_inventory"
    return None


def reject(operation, command, reason, identity=None):
    append_event(
        _policy,
        {
            "action": "rejected",
            "layer": "python",
            "operation": operation,
            "command": command,
            "executable_identity": identity,
            "pid": os.getpid(),
            "parent_process": os.getppid(),
            "reason": reason,
            "environment_policy": environment_policy(os.environ),
        },
    )
    raise ScopeDenied("J72-C1 scope: " + reason)


def inherited_environment(environment=None):
    env = dict(os.environ if environment is None else environment)
    env[POLICY_KEY], env[HASH_KEY] = _policy_path, _policy_hash
    paths = env.get("PYTHONPATH", "").split(os.pathsep)
    env["PYTHONPATH"] = os.pathsep.join(
        [_policy["bootstrap"], *[p for p in paths if p != _policy["bootstrap"]]]
    )
    env.pop("PYTHONHOME", None)
    return env


def forbidden(symbol):
    def blocked(*args, **kwargs):
        reject("python_entry_point", [symbol], "QM_MD_or_fitting_entry_point")

    return blocked


def patch_module(module):
    name = module.__name__
    if name.startswith(("openmm", "simtk.openmm")):
        # Restore/clone can otherwise construct native integrators without
        # calling the public constructor. System-only XML stays read-only.
        if hasattr(module, "XmlSerializer_deserialize"):
            deserialize = module.XmlSerializer_deserialize

            def guarded_deserialize(xml, *args, **kwargs):
                if "<Integrator" in str(xml):
                    reject(
                        "integrator_deserialize",
                        [name + ".XmlSerializer_deserialize"],
                        "MD_integrator_restore",
                    )
                return deserialize(xml, *args, **kwargs)

            module.XmlSerializer_deserialize = guarded_deserialize
        if hasattr(module, "XmlSerializer_clone"):
            clone = module.XmlSerializer_clone

            def guarded_clone(obj, *args, **kwargs):
                if any("Integrator" in c.__name__ for c in type(obj).__mro__):
                    reject(
                        "integrator_clone",
                        [name + ".XmlSerializer_clone"],
                        "MD_integrator_clone",
                    )
                return clone(obj, *args, **kwargs)

            module.XmlSerializer_clone = guarded_clone
        for key in list(vars(module)):
            value = getattr(module, key)
            if isinstance(value, type) and "Integrator" in key:
                # Block construction as well as stepping. No real integrators
                # are constructed even when an old test only wants a Context.
                for method in ("__init__", "step", "stepTo"):
                    if hasattr(value, method):
                        setattr(
                            value, method, forbidden(name + "." + key + "." + method)
                        )
            elif (
                "Integrator" in key
                and callable(value)
                and (
                    key.startswith("new_")
                    or key.endswith(
                        ("_step", "_stepTo", "deserializeIntegrator", "cloneIntegrator")
                    )
                )
            ):
                setattr(module, key, forbidden(name + "." + key))
    if name.startswith("island.dynamics"):
        for key in (
            "kick",
            "verlet_drift",
            "baoab_drift",
            "run_nve",
            "run_langevin",
            "run_dynamics_segment",
            "resume_dynamics",
        ):
            if hasattr(module, key):
                setattr(module, key, forbidden(name + "." + key))
    if name == "island.pcff_qm_pilot.hmds_constrained":
        # Small manufactured software controls remain separate. Never solve
        # the real 58-column HMDS fitting problem during corrective replay.
        original = module.solve

        def solve(system, *args, **kwargs):
            if len(system["scale"]) >= 59:
                reject(
                    "scientific_fit",
                    [name + ".solve"],
                    "real_HMDS_fitting_outside_scope",
                )
            return original(system, *args, **kwargs)

        module.solve = solve


class Finder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        deny = fullname.split(".")[0] in BLOCKED_MODULES
        patch = (
            fullname.startswith(("openmm", "simtk.openmm", "island.dynamics"))
            or fullname == "island.pcff_qm_pilot.hmds_constrained"
        )
        if not deny and not patch:
            return None
        spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        if spec is None or spec.loader is None:
            return spec
        original = spec.loader

        class Loader(importlib.abc.Loader):
            def create_module(self, spec):
                if deny:
                    return None
                return (
                    original.create_module(spec)
                    if hasattr(original, "create_module")
                    else None
                )

            def exec_module(self, module):
                if deny:
                    reject("module_import", [fullname], "scientific_module_execution")
                original.exec_module(module)
                patch_module(module)

        spec.loader = Loader()
        return spec


def install():
    global _policy, _policy_path, _policy_hash
    if _policy is not None:
        return
    _policy_path, _policy_hash = os.environ[POLICY_KEY], os.environ[HASH_KEY]
    _policy = load_policy(_policy_path, _policy_hash)
    _policy.update(policy_path=_policy_path, policy_sha256=_policy_hash)
    status = Path("/proc/self/status").read_text()
    if not re.search(r"TracerPid:\s+[1-9]", status) or "Seccomp:\t2" not in status:
        raise RuntimeError("C1 kernel supervisor absent")
    if any(name in sys.modules for name in BLOCKED_MODULES):
        raise RuntimeError("C1 scientific module imported before fence")
    original = subprocess.Popen

    def guarded(argv, *args, **kwargs):
        env = inherited_environment(kwargs.get("env"))
        identity = (
            executable_identity(argv[0], env, kwargs.get("cwd"))
            if isinstance(argv, (list, tuple)) and argv
            else None
        )
        reason = (
            "shell_execution"
            if kwargs.get("shell")
            else command_reason(argv, env, _policy, identity or {})
        )
        if reason:
            reject(
                "subprocess.Popen",
                list(map(os.fsdecode, argv))
                if isinstance(argv, (list, tuple))
                else argv,
                reason,
                identity,
            )
        if kwargs.get("executable") is not None:
            reject(
                "subprocess.Popen", list(map(os.fsdecode, argv)), "alternate_executable"
            )
        kwargs["env"] = env
        return original(argv, *args, **kwargs)

    subprocess.Popen = guarded
    sys.meta_path.insert(0, Finder())
    append_event(
        _policy,
        {
            "action": "installed",
            "layer": "python",
            "pid": os.getpid(),
            "parent_process": os.getppid(),
            "environment_policy": environment_policy(os.environ),
        },
    )


def assert_installed():
    if _policy is None:
        raise RuntimeError("C1 fence must be installed before pytest starts")
