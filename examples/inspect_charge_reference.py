"""Inspect a saved reference without QM, OpenMM, RDKit or ParmEd execution."""

import argparse
import json

from island.charge_references import load_charge_reference
from island.charge_references.audit import summarize

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference")
    args = parser.parse_args()
    reference = load_charge_reference(args.reference)
    summary = summarize(reference)
    print(json.dumps(summary, indent=2, allow_nan=False))
    print("Whole-oligomer AM1-BCC evidence; not a transferable charge assignment.")
