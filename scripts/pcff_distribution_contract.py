"""Immutable J16 experiment inputs, independent from mutable outcome receipts."""

from pathlib import Path

from island.workflows import storage

EVIDENCE = Path(__file__).parents[1] / "docs/evidence"
HASHES = {
    "declaration": "5429c30f58529a822c22621f9b741567433335f7c94440ebdec642bcf54cc2f9",
    "converter": "ea0519a80fc2fdce13a898ba20d9ce9344f84f5d586b3b872c6050dc22bc3690",
    "reference_links": "88cada3d489e34ac4e2dcf3b6f3e338f6de5137446b5be8708f5a20550ee00f6",
}


def frozen(name):
    path = EVIDENCE / f"phase_4j16_{name}.json"
    if storage.checksum(path.read_bytes()) != HASHES[name]:
        raise ValueError("Changed J16 experiment contract: " + name)
    return storage.read_json(path)
