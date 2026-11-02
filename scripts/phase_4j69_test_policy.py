"""Pytest scope policy: retained fixtures and no new QM or dynamics execution.

All collected tests remain in terminal accounting. Tests requiring a forbidden
execution are reported as scope skips, never as passes. Synthetic conic and
static energy/force checks remain enabled.
"""

import importlib
import re
import subprocess
from pathlib import Path

import pytest


def denied_command(argv, cwd=None):
    if not isinstance(argv, (list, tuple)):
        return False
    words = list(map(str, argv))
    if not words:
        return False
    exe = Path(words[0]).name.lower()
    if exe in {"antechamber", "sqm", "sander", "pmemd", "psi4"}:
        return True
    text = " ".join(words)
    if (
        re.search(r"scripts\.(?:validate_)?phase_4j(3[2-9]|40|41)(\b|_)", text)
        and "_receipt" not in text
        and not any(w in words for w in ("--offline", "--audit-only", "--help"))
    ):
        return True
    if "-c" in words and any(
        w in text
        for w in (
            "run_nve(",
            "run_langevin(",
            ".step(",
            "run_dynamics_segment(",
            "resume_dynamics(",
            "psi4.gradient(",
            "psi4.energy(",
        )
    ):
        return True
    if (exe.startswith("lmp") or exe == "lammps") and "-in" in words:
        p = Path(cwd or ".") / words[words.index("-in") + 1]
        if not p.is_file() or re.search(
            r"^\s*run\s+(?!0\b)", p.read_text(), re.MULTILINE
        ):
            return True
    return False


@pytest.fixture(autouse=True)
def no_new_native_or_MD(monkeypatch):
    old = subprocess.Popen

    def popen(argv, *args, **kwargs):
        if denied_command(argv, kwargs.get("cwd")):
            pytest.skip(
                "J69 scope: no new QM/MD; native integration execution not activated"
            )
        return old(argv, *args, **kwargs)

    monkeypatch.setattr(subprocess, "Popen", popen)

    def forbidden(*args, **kwargs):
        pytest.skip(
            "J69 scope: dynamics steps forbidden, including manufactured integrator runs"
        )

    for name in ("_steps", "engine", "langevin", "checkpoints"):
        module = importlib.import_module("island.dynamics." + name)
        for symbol in ("kick", "verlet_drift", "baoab_drift"):
            if hasattr(module, symbol):
                monkeypatch.setattr(module, symbol, forbidden)
    try:
        import openmm
    except ImportError:
        return
    for name in dir(openmm):
        cls = getattr(openmm, name)
        if isinstance(cls, type) and issubclass(cls, openmm.Integrator):
            for symbol in ("step", "stepTo"):
                if hasattr(cls, symbol):
                    monkeypatch.setattr(cls, symbol, forbidden)
