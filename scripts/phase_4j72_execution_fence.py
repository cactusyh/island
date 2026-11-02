"""J72 test-only execution fence, including Python subprocess descendants.

Ordinary package APIs are unchanged. This module is installed only in the
verification environment; guarded scope exclusions are reported as skips.
"""

import importlib.abc
import importlib.machinery
import os
import re
import subprocess
import sys
from pathlib import Path


class ScopeDenied(RuntimeError):
    pass


DANGEROUS = (
    "start_periodic_md_workflow(",
    "resume_periodic_md_workflow(",
    "run_nve(",
    "run_langevin(",
    "run_dynamics_segment(",
    "resume_dynamics(",
    ".step(",
    ".stepTo(",
    "psi4.gradient(",
    "psi4.energy(",
)


def denied(argv, cwd=None):
    if not isinstance(argv, (list, tuple)):
        return True
    if not argv:
        return False
    words = list(map(str, argv))
    exe = Path(words[0]).name.lower()
    if exe in {"msi2lmp", "sqm", "antechamber", "psi4", "sander", "pmemd"}:
        return True
    if "-S" in words:
        return True
    if (exe.startswith("lmp") or exe == "lammps") and "-in" in words:
        p = Path(cwd or ".") / words[words.index("-in") + 1]
        return not p.is_file() or bool(
            re.search(r"^\s*run\s+(?!0\b)", p.read_text(), re.MULTILINE)
        )
    if "-m" in words and words[words.index("-m") + 1] == "scripts.phase_4j37_reference":
        # This pinned entry point only reconstructs retained raw provenance.
        return False
    text = " ".join(words)
    if (
        re.search(r"scripts\.(?:validate_)?phase_4j(3[1-9]|40|41)(\b|_)", text)
        and "_receipt" not in text
        and not any(
            w in words for w in ("--offline", "--audit-only", "--help", "reconstruct")
        )
    ):
        return True
    if "-c" in words and any(s in text for s in DANGEROUS):
        return True
    if "-m" in words:
        module = words[words.index("-m") + 1]
        if module.startswith("scripts.") and "_receipt" not in module:
            p = Path(__file__).parent / (module.rsplit(".", 1)[-1] + ".py")
            if (
                p.is_file()
                and any(s in p.read_text() for s in DANGEROUS)
                and not any(
                    w in words
                    for w in ("--help", "--offline", "--audit-only", "reconstruct")
                )
            ):
                return True
    return False


def fenced_environment(environment=None):
    env = dict(os.environ if environment is None else environment)
    directory = os.environ.get("ISLAND_J72_FENCE_DIRECTORY")
    if directory:
        env["ISLAND_J72_FENCE_DIRECTORY"] = directory
        entries = env.get("PYTHONPATH", "").split(os.pathsep)
        if directory not in entries:
            env["PYTHONPATH"] = directory + os.pathsep + env.get("PYTHONPATH", "")
    return env


def block_steps(module):
    def forbidden(*a, **kw):
        raise ScopeDenied("J72 child fence: MD steps are forbidden")

    for name in dir(module):
        cls = getattr(module, name)
        if isinstance(cls, type) and hasattr(cls, "step") and "Integrator" in name:
            cls.step = forbidden
            if hasattr(cls, "stepTo"):
                cls.stepTo = forbidden


class Fence(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname == "psi4" or fullname.startswith("psi4."):
            # Discovery is read-only. Block the loader before module code runs,
            # while preserving an absent module's ordinary find_spec(None).
            spec = importlib.machinery.PathFinder.find_spec(fullname, path)
            if spec is None:
                return None

            class NoQuantumLoader(importlib.abc.Loader):
                def create_module(self, spec):
                    return None

                def exec_module(self, module):
                    raise ScopeDenied("J72 child fence: QM startup forbidden")

            spec.loader = NoQuantumLoader()
            return spec
        if fullname not in ("openmm", "openmm.openmm"):
            return None
        spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        if spec is None or spec.loader is None:
            return None
        original = spec.loader

        class Loader(importlib.abc.Loader):
            def create_module(self, spec):
                return (
                    original.create_module(spec)
                    if hasattr(original, "create_module")
                    else None
                )

            def exec_module(self, module):
                original.exec_module(module)
                block_steps(module)

        spec.loader = Loader()
        return spec


def install_child_fence():
    if getattr(subprocess, "_island_j72_fenced", False):
        return
    original = subprocess.Popen

    def guarded(argv, *a, **kw):
        if denied(argv, kw.get("cwd")):
            raise ScopeDenied("J72 child fence: native/QM/MD command forbidden")
        kw["env"] = fenced_environment(kw.get("env"))
        return original(argv, *a, **kw)

    subprocess.Popen = guarded
    subprocess._island_j72_fenced = True
    sys.meta_path.insert(0, Fence())
    for name in ("openmm", "openmm.openmm"):
        if name in sys.modules:
            block_steps(sys.modules[name])
