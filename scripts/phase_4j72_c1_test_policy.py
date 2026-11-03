"""C1 test policy. Execution-dependent checks are scope skips, never passes."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from scripts.phase_4j72_c1_fence import ScopeDenied, assert_installed

_active_node = None


def pytest_configure(config):
    assert_installed()
    module = sys.modules[ScopeDenied.__module__]
    original = module.reject

    def reject(*args, **kwargs):
        try:
            return original(*args, **kwargs)
        except ScopeDenied as error:
            # Some public adapters intentionally catch Exception and return a
            # failed backend record. A scope skip must bypass those handlers,
            # rather than become a misleading pass or assertion failure.
            if _active_node and not _active_node.startswith(
                "tests/test_phase_4j72_c1.py::"
            ):
                pytest.skip("C1 scope: " + str(error))
            raise

    module.reject = reject
    original_popen = subprocess.Popen

    def popen(argv, *args, **kwargs):
        process = original_popen(argv, *args, **kwargs)
        communicate_original = process.communicate

        def communicate(*a, **kw):
            output = communicate_original(*a, **kw)
            if _active_node and not _active_node.startswith(
                "tests/test_phase_4j72_c1.py::"
            ):
                events = [
                    json.loads(line)
                    for line in Path(module._policy["events"]).read_text().splitlines()
                ]
                rejected = [
                    e
                    for e in events
                    if e["action"] == "rejected"
                    and (
                        e.get("pid") == process.pid
                        or e.get("parent_process") == process.pid
                    )
                ]
                if rejected:
                    # A child adapter may swallow ScopeDenied and return a
                    # backend-failed JSON record (even with return code zero).
                    # Retain its complete captured output, then scope-skip the
                    # check instead of presenting it as scientific verification.
                    module.append_event(
                        module._policy,
                        {
                            "action": "child_scope_skip",
                            "layer": "pytest_policy",
                            "command": list(map(os.fsdecode, argv)),
                            "child_pid": process.pid,
                            "returncode": process.returncode,
                            "test_node": _active_node,
                            "stdout": output[0].hex()
                            if isinstance(output[0], bytes)
                            else output[0],
                            "stderr": output[1].hex()
                            if isinstance(output[1], bytes)
                            else output[1],
                            "byte_outputs_encoding": "hex",
                        },
                    )
                    pytest.skip(
                        "C1 scope: descendant scientific operation rejected; command and captured output retained in scope ledger"
                    )
            return output

        process.communicate = communicate
        return process

    subprocess.Popen = popen


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_protocol(item, nextitem):
    global _active_node
    _active_node = item.nodeid
    try:
        yield
    finally:
        _active_node = None


def pytest_addoption(parser):
    parser.addoption("--c1-collection-out")


def pytest_collection_finish(session):
    output = session.config.getoption("--c1-collection-out")
    if output:
        Path(output).write_text(
            json.dumps(
                {"nodes": sorted(item.nodeid for item in session.items)}, indent=2
            )
            + "\n"
        )


def pytest_collection_modifyitems(items):
    for item in items:
        if (
            item.nodeid
            == "tests/test_phase_4j55.py::test_cli_flags_independent[--require-openmm-runtime-0]"
        ):
            item.add_marker(
                pytest.mark.skip(
                    reason="C1 scope: child CLI requires a real OpenMM integrator for runtime qualification"
                )
            )
        if item.path.name in {
            "test_phase_4j31_execution.py",
            "test_ambertools_integration.py",
        }:
            item.add_marker(
                pytest.mark.skip(reason="C1 scope: real QM/MD integration forbidden")
            )
        if item.nodeid.endswith(
            "::test_child_fence_blocks_integrator_and_native_command_before_execution"
        ):
            item.add_marker(
                pytest.mark.skip(
                    reason="C1 scope: historical test constructs a real OpenMM integrator; replaced by synthetic C1 fixture"
                )
            )


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    outcome = yield
    report = outcome.get_result()
    error = call.excinfo.value if call.excinfo is not None else None
    seen, denial = set(), None
    while error is not None and id(error) not in seen:
        seen.add(id(error))
        if isinstance(error, ScopeDenied):
            denial = error
            break
        error = error.__cause__ or error.__context__
    if denial is not None:
        report.outcome = "skipped"
        report.longrepr = (
            str(item.path),
            item.location[1],
            "C1 scope skip: " + str(denial),
        )
