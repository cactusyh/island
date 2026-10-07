"""Apply the retained minimization, relocation and full-RNG continuation gates to J17."""

from validate_pcff_fallback_workflow import main

if __name__ == "__main__":
    raise SystemExit(
        main(
            expected_implementation="island_pcff_typed_graph_singlepoint_v1",
            experiment={
                "schema": "island_j17_typed_graph_workflow_v1",
                "experiment": "External c/h polyethylene, independently reconstructed provided charges",
                "construction": "Frozen PSMILES DP3 final graph; checked J17 bundle, no preparation on resume",
            },
        )
    )
