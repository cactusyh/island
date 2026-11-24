import pytest

pytest.importorskip("scipy")
from island.epoxy.convergence import summarize


def test_stalled_seed_and_matched_time_delta():
    def case(seed, alpha):
        record = {
            "dgeba": 4,
            "pacm": 2,
            "seed": seed,
            "elapsed_seconds": 1,
            "statistics": {"epoxy_conversion": alpha, "winding_rank": 1},
        }
        if alpha == 0.75:
            record["volumes"] = {
                k: {
                    t: {"volume_A3": v, "density_g_cm3": 1 / v}
                    for t, v in zip(("12000", "40000"), vol)
                }
                for k, vol in [("before", (100, 100)), ("after", (95, 90))]
            }
        return record

    cases = [case(41, 0.75), case(173, 0.625)]
    failed = case(2026, 0.75)
    failed.update(conditioning_failure={"state": "after", "diagnostic": "LAMMPS exit"})
    del failed["volumes"]["after"]
    cases.append(failed)
    size = summarize(cases)["sizes"][0]
    assert (size["n_total"], size["n_at_target"]) == (3, 2)
    assert size["seeds_with_paired_measurements"] == [41]
    delta = size["windows"]["40000"]["change_from_6ps"]["shrinkage_percent"]
    assert delta["mean"] == pytest.approx(5) and delta["ci95"] is None
    with pytest.raises(ValueError):
        summarize(cases + cases)
