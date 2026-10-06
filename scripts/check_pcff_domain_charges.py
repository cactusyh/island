"""Independent source-increment readback for the two declared J3 bundles."""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path("scripts").resolve()))
from validate_pcff_fallbacks import raw_row

from island.forcefields import PreparedForceFieldSources, load_prepared_forcefield
from island.workflows.storage import checksum

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--source", type=Path, required=True)
p.add_argument("--numerical", type=Path, required=True)
p.add_argument("--output", type=Path, required=True)
a = p.parse_args()
root = a.numerical
src = a.source
sources = PreparedForceFieldSources(pcff_frc=src)
rows = []
for name, key in [("hydrogen_sulfide", ["h", "s"]), ("heavy_water", ["h*", "o*"])]:
    p = root / name
    outcome = json.loads((p / "outcome.json").read_text())
    bundle = load_prepared_forcefield(p / "bundle", sources=sources)
    source_row = raw_row(src.read_bytes(), "bond_increments", key)
    assert len(source_row["values"]) == 2
    expected = {
        i: (
            source_row["values"][0] if a.element == "H" else 2 * source_row["values"][1]
        )
        for i, a in bundle.system.topology.sites.items()
    }
    actual = {
        n["site"]: n["charge"]
        for n in bundle.prepared.native_result.payload["nonbonded"]
    }
    assert (
        actual.keys() == expected.keys()
        and max(abs(actual[i] - expected[i]) for i in actual) < 1e-12
    )
    assert (
        bundle.prepared.identity == outcome["prepared_identity"]
        and bundle.prepared.native_result.identity == outcome["model_identity"]
    )
    rows.append(
        {
            "case": name,
            "raw_increment_row": source_row,
            "actual": actual,
            "expected": expected,
            "max_error": max(abs(actual[i] - expected[i]) for i in actual),
            "prepared_identity": bundle.prepared.identity,
            "model_identity": bundle.prepared.native_result.identity,
        }
    )
from island.workflows.storage import publish

publish(
    a.output,
    json.dumps(
        {"cases": rows, "source_sha256": checksum(src.read_bytes())}, indent=2
    ).encode(),
)
print("charges match independent raw rows")
