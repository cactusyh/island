"""Inspect or reuse an explicit fragment template; never rerun its charge backend.

python examples/generic_fragment_charges.py /path/to/template.json --cache-key SHA --dp 4
Use scripts.validate_fragment_charges to declare/run actual bounded calculations.
"""

import argparse
import json

from island.builders import build_linear_polymer
from island.fragment_charges import assign_fragment_charges, load_fragment_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template")
    parser.add_argument("--cache-key", required=True)
    parser.add_argument("--dp", type=int, default=4)
    args = parser.parse_args()
    template = load_fragment_record(args.template, expected_cache_key=args.cache_key)
    p = template.payload
    system = build_linear_polymer(
        p["calculation"]["fragment"]["definition"], dp=args.dp, generate_3d=False
    )
    # No FF claim from a supplied charge vector; AM1-BCC binding remains explicit.
    ff = p["calculation"]["backend"].get("force_field")
    result = assign_fragment_charges(template, system, force_field=ff)
    print(json.dumps(result.payload["data"], indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
