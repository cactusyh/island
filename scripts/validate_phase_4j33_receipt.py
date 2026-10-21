"""Validate additive J33 evidence without rewriting historical receipts."""

import argparse
import json
from hashlib import sha256
from pathlib import Path

from island.exceptions import ValidationError

REQUIRED_SUCCESS = {"oplsaa/linear", "oplsaa/network", "gaff2/network"}
REQUIRED_NEGATIVE = {"pcff/network", "missing_barostat", "invalid_box", "malformed_config"}


def _digest(value):
    body = {k: v for k, v in value.items() if k != "identity"}
    raw = (json.dumps(body, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode()
    return sha256(raw).hexdigest()


def validate(receipt):
    if receipt.get("schema") != "island.j33-receipt.v1" or receipt.get("identity") != _digest(receipt):
        raise ValidationError("J33 receipt identity mismatch")
    if receipt.get("production_validated") is not False or receipt.get("simulation_readiness") != "not_established":
        raise ValidationError("J33 readiness flags were promoted")
    cases = receipt.get("cases")
    if not isinstance(cases, dict) or not REQUIRED_SUCCESS.issubset(cases) or not REQUIRED_NEGATIVE.issubset(cases):
        raise ValidationError("J33 receipt is missing required success or negative-control cases")
    for key in REQUIRED_SUCCESS:
        row = cases[key]
        if row.get("classification") == "asset_unavailable" and receipt.get("acceptance_met") is False:
            continue
        if row.get("status") != "complete" or row.get("npt_executed") is not True or row.get("equilibrium_status") != "not_established":
            raise ValidationError(f"Required NPT case is not complete: {key}")
    pcff = cases["pcff/network"]
    if pcff.get("status") != "blocked" or pcff.get("classification") != "expected_typing_rejection" or pcff.get("npt_executed") is not False:
        raise ValidationError("PCFF negative control is not the retained typing rejection")
    for key in REQUIRED_NEGATIVE - {"pcff/network"}:
        if cases[key].get("status") not in {"blocked", "not_executed", "asset_unavailable"}:
            raise ValidationError(f"Negative control has an unsupported status: {key}")
    complete = all(cases[key].get("status") == "complete" for key in REQUIRED_SUCCESS)
    if receipt.get("acceptance_met") is True and not complete:
        raise ValidationError("J33 acceptance flag contradicts required cases")
    if receipt.get("acceptance_met") is False and complete:
        raise ValidationError("J33 acceptance flag contradicts complete required cases")
    return receipt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("receipt", type=Path)
    args = parser.parse_args()
    try:
        validate(json.loads(args.receipt.read_text()))
    except Exception as error:  # noqa: BLE001 -- CLI reports a validation diagnostic
        print(f"FAIL: {error}")
        return 1
    print("PASS: J33 receipt validated")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
