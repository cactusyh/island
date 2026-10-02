"""Predict experimental PE charges from a frozen model; no QM or dynamics.

python examples/predict_pe_template.py /path/to/conserving.json --dp 20
This graph-only example uses RDKit construction. The prediction/loading APIs
use core/NumPy only when the authoritative system is already available.
"""

import argparse
import json

from island.builders import build_linear_polymer
from island.charge_references import load_pe_record, predict_pe_charges


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model")
    parser.add_argument("--dp", type=int, default=20)
    args = parser.parse_args()
    system = build_linear_polymer("[*:1]CC[*:2]", dp=args.dp, generate_3d=False)
    result = predict_pe_charges(load_pe_record(args.model), system)
    print(json.dumps(result.payload["prediction"], indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
