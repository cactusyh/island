"""A checkpoint survives process exit; no evaluator/RNG instance is shared."""

import importlib.util
from pathlib import Path


def test_three_independent_processes(tmp_path):
    path = Path(__file__).resolve().parents[1] / "examples/checkpoint_continuation.py"
    spec = importlib.util.spec_from_file_location("checkpoint_example", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = module.demonstrate(tmp_path)
    assert all(
        r["max_state_difference"] == 0
        and r["split_calls"] == 21
        and r["uninterrupted_calls"] == 17
        for r in report.values()
    )
