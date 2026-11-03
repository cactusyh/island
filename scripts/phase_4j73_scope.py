"""J73 evidence-only restrictions layered on the unchanged inherited C1 fence."""

import importlib.abc
import importlib.machinery
import json
import re
import subprocess
import sys
from pathlib import Path

from scripts import phase_4j72_c1_fence as fence

BLOCKED = {
    "openmm",
    "simtk",
    "lammps",
    "psi4",
    "pyscf",
    "sander",
    "xtb",
    "tblite",
    "qcengine",
    "qcelemental",
}
_installed = False


class DiscoveryOnly(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split(".")[0] not in BLOCKED:
            return None
        spec = importlib.machinery.PathFinder.find_spec(fullname, path)
        if spec is None or spec.loader is None:
            return spec

        class Loader(importlib.abc.Loader):
            def create_module(self, spec):
                return None

            def exec_module(self, module):
                fence.reject(
                    "module_import", [fullname], "J73_no_scientific_module_execution"
                )

        spec.loader = Loader()
        return spec


def install():
    global _installed
    fence.assert_installed()
    if _installed:
        return
    if any(n.split(".")[0] in BLOCKED for n in sys.modules):
        raise RuntimeError("J73 scientific module imported before evidence-only gate")
    sys.meta_path.insert(0, DiscoveryOnly())
    original = subprocess.Popen

    def guarded(argv, *args, **kwargs):
        # C1 independently mediates every native exec in all descendants. This
        # additional cooperative Python gate rejects named scientific drivers.
        # It does not try to infer the meaning of arbitrary Python -c text.
        if isinstance(argv, (tuple, list)) and len(argv) > 1:
            word = str(argv[2] if argv[1] == "-m" and len(argv) > 2 else argv[1])
            if re.search(
                r"(?:^|[/._-])(?:run_qm|native_qm|run_md|fit_hmds|hmds_fit|"
                r"parameter_generation|generate_parameters|reserved_validation|msi2lmp)(?:$|[/._-])",
                word,
            ):
                fence.reject(
                    "subprocess.Popen",
                    list(map(str, argv)),
                    "J73_scientific_driver_forbidden",
                )
        return original(argv, *args, **kwargs)

    subprocess.Popen = guarded
    _installed = True


# Loading this explicit pytest plugin happens before test collection. J73 CLI
# entry points also call install; no public ISLAND API imports this module.
install()


def pytest_configure(config):
    install()


def pytest_addoption(parser):
    parser.addoption("--j73-collection-out")


def pytest_collection_finish(session):
    output = session.config.getoption("--j73-collection-out")
    if output:
        Path(output).write_text(
            json.dumps({"nodes": sorted(i.nodeid for i in session.items)}, indent=2)
            + "\n"
        )


def pytest_sessionfinish(session, exitstatus):
    if any(n.split(".")[0] in BLOCKED for n in sys.modules):
        raise RuntimeError("J73 scientific module reached the verification process")
