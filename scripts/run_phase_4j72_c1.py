"""Launch software checks only after installing the inherited C1 execution fence.

Example: python -m scripts.run_phase_4j72_c1 --output /tmp/c1 -- -m pytest
  -p scripts.phase_4j72_c1_test_policy tests/test_phase_4j72_c1.py
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

from scripts.phase_4j72_c1_fence import HASH_KEY, POLICY_KEY, load_policy, sha
from scripts.phase_4j72_c1_supervisor import supervise


def configure(output, python):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    bootstrap = output / "bootstrap"
    bootstrap.mkdir()
    runtime = Path(__file__).with_name("phase_4j72_c1_fence.py").resolve()
    repo = runtime.parent.parent
    startup = (
        "import os,sys,importlib.util,datetime\ntry:\n"
        " if not hasattr(datetime,'UTC'): datetime.UTC=datetime.timezone.utc\n"
        f" spec=importlib.util.spec_from_file_location('scripts.phase_4j72_c1_fence',{str(runtime)!r})\n"
        " module=importlib.util.module_from_spec(spec)\n"
        " sys.modules[spec.name]=module\n"
        " spec.loader.exec_module(module)\n"
        " module.install()\n"
        "except BaseException as exc:\n"
        " os.write(2,('C1 fence startup failed: '+repr(exc)+'\\n').encode())\n"
        " os._exit(97)\n"
    )
    (bootstrap / "sitecustomize.py").write_text(startup)
    executables = {}
    roots = [
        repo.parent / "island-validation" / name / "bin/python"
        for name in (
            "phase4g1-env",
            "phase4e2-env",
            "phase4j32-analysis-env",
            "phase4j36-native-env",
            "phase4j36-psi4-env",
            "phase4j66-solver-env",
        )
    ]
    roots += [
        Path(python),
        Path(sys.executable),
        repo.parent / "miniforge3/envs/md_rg/bin/python",
    ]
    for path in roots:
        if path.is_file():
            path = path.resolve()
            executables[str(path)] = {"kind": "python", "sha256": sha(path)}
    for name in ("git", "ruff"):
        found = shutil.which(name)
        if found:
            path = Path(found).resolve()
            executables[str(path)] = {"kind": name, "sha256": sha(path)}
    ruff = Path(python).parent / "ruff"
    if ruff.is_file():
        executables[str(ruff.resolve())] = {"kind": "ruff", "sha256": sha(ruff)}
    policy = {
        "schema": "island.j72-c1-execution-policy.v1",
        "runtime": str(runtime),
        "runtime_sha256": sha(runtime),
        "startup_sha256": sha(bootstrap / "sitecustomize.py"),
        "bootstrap": str(bootstrap),
        "events": str(output / "events.jsonl"),
        "executables": executables,
        "authorized_new_QM": 0,
        "authorized_MD_steps": 0,
        "new_scientific_fitting": 0,
        "native_default": "deny",
        "shell_default": "deny",
        "kernel_gate": "inherited seccomp TRACE execve/execveat + ptrace EXITKILL",
    }
    path = output / "policy.json"
    path.write_text(json.dumps(policy, indent=2, sort_keys=True) + "\n")
    digest = sha(path)
    load_policy(path, digest)
    policy.update(policy_path=str(path), policy_sha256=digest)
    env = {
        **os.environ,
        POLICY_KEY: str(path),
        HASH_KEY: digest,
        "PYTHONPATH": str(bootstrap) + os.pathsep + os.environ.get("PYTHONPATH", ""),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "OPENMM_CPU_THREADS": "1",
    }
    env.pop("PYTHONHOME", None)
    return policy, env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    arguments = args.args[1:] if args.args[:1] == ["--"] else args.args
    if not arguments:
        parser.error("a software command is required")
    # Keep venv path in argv to preserve its site-packages.
    command = [str(Path(args.python).absolute()), *arguments]
    policy, env = configure(args.output, args.python)
    code = supervise(command, env, policy)
    (args.output / "result.json").write_text(
        json.dumps(
            {
                "command": command,
                "returncode": code,
                "policy_sha256": policy["policy_sha256"],
                "events_sha256": sha(policy["events"]),
                "C1_new_native_QM_attempts": 0,
                "C1_MD_steps": 0,
            },
            indent=2,
        )
        + "\n"
    )
    return code if code is not None and code >= 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
