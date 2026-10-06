"""J3 reuses the established four-step workflow acceptance contract.

The bundle and preceding numerical receipt select the case explicitly. Child
resume uses the unchanged generic acceptance child; no typing or preparation.
"""

from validate_pcff_fallback_workflow import main

if __name__ == "__main__":
    raise SystemExit(
        main(
            experiment={
                "schema": "island_j3_domain_workflow_v1",
                "experiment": "new J3 bounded workflow; independent model verification in separate numerical receipt",
                "construction": "unchanged checked J3 bundle, from declared SMILES and seed 2026",
            }
        )
    )
