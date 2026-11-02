"""J72 pytest policy: no native/QM/MD descendants; production APIs unchanged."""

import os
import subprocess
from pathlib import Path

import pytest

from scripts import phase_4j72_execution_fence as fence
from scripts.phase_4j69_test_policy import no_new_native_or_MD  # noqa: F401
from scripts.phase_4j72_execution_fence import denied, fenced_environment


def pytest_collection_modifyitems(items):
    for item in items:
        if item.path.name == "test_phase_4j31_execution.py":
            item.add_marker(
                pytest.mark.skip(reason="J72 scope: J31 native/MD subprocess forbidden")
            )


def pytest_configure(config):
    directory = Path(os.environ["ISLAND_J72_FENCE_DIRECTORY"])
    directory.mkdir(parents=True, exist_ok=True)
    module = Path(__file__).with_name("phase_4j72_execution_fence.py").resolve()
    startup = (
        "import os,runpy\ntry:\n runpy.run_path("
        + repr(str(module))
        + ')["install_child_fence"]()\nexcept BaseException:\n os._exit(97)\n'
    )
    (directory / "sitecustomize.py").write_text(startup)
    os.environ.update(fenced_environment())
    # Global parent fence also covers module/session fixtures before autouse.
    fence.install_child_fence()
    fence.ScopeDenied = pytest.skip


@pytest.fixture(autouse=True)
def no_unfenced_descendants(monkeypatch):
    original = subprocess.Popen

    def guarded(argv, *a, **kw):
        if denied(argv, kw.get("cwd")):
            pytest.skip("J72 scope: native/QM/MD subprocess forbidden")
        kw["env"] = fenced_environment(kw.get("env"))
        return original(argv, *a, **kw)

    monkeypatch.setattr(subprocess, "Popen", guarded)
