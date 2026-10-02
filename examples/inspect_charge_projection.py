"""Inspect an experimental projection offline; no QM or automatic assignment.

python examples/inspect_charge_projection.py /path/to/projections/pe-dp3-seed2026.json
Generate records using: python -m scripts.audit_charge_conservation --help
"""

import argparse
import json

from island.charge_references import load_charge_projection


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("projection")
    args = parser.parse_args()
    payload = load_charge_projection(args.projection).payload
    source = payload["observation"]["data"]
    projection = payload["projection"]
    print(
        json.dumps(
            {
                "case": source["case"],
                "historical_status": source["historical_status"],
                "historical_failure": source["historical_failure"],
                "policy": projection["policy"],
                "raw_stage": source["raw_stage"],
                "raw_total_e": projection["raw_total_e"],
                "projected_total_e": projection["projected_total_e"],
                "uniform_offset_e": projection["uniform_offset_e"],
                "atom_count": projection["atom_count"],
                "observation_identity": payload["observation_identity"],
                "production_validated": False,
                "simulation_readiness": "not_established",
            },
            indent=2,
            allow_nan=False,
        )
    )


if __name__ == "__main__":
    main()
