from pathlib import Path

from island.forcefields.capabilities import benchmark


def test_archived_coverage_and_evidence_levels():
    b = benchmark(Path(__file__).resolve().parents[1] / "docs/evidence")
    assert b["linear"]["pcff"] == {"assigned": 678, "strict": 646, "total": 982}
    assert b["published_positive"] == {"assigned": 611, "total": 611}
    assert b["missing_increments"]["pcff"]["c_1,na"] == 232
    assert len(b["fixed_graph_md"]) == 12 and len(b["periodic_epoxy"]) == 7

    assert b["md_runs_by_architecture"] == {
        "branched_block": 2,
        "branch": 3,
        "block": 7,
    }
