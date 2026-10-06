"""Independent line accounting and pinned msi2lmp capability audit; fail closed."""

import argparse
from collections import Counter
from decimal import Decimal, InvalidOperation
from pathlib import Path

from island.forcefields.pcff import load_pcff_source
from island.forcefields.pcff.row_coverage import pcff_source_row_ledger
from island.workflows import storage


def independently_index(raw):
    family, namespace = None, None
    rows = {}
    for line, text in enumerate(raw.decode().splitlines(), 1):
        text = text.strip()
        if text.startswith("#"):
            tokens = text.split()
            family = tokens[0][1:]
            namespace = " ".join(tokens[1:])
            continue
        if (
            not text
            or text.startswith(("!", ">", "@"))
            or family in (None, "version", "reference", "end")
        ):
            continue
        tokens = text.split("!", 1)[0].split()
        try:
            version = Decimal(tokens[0])
        except InvalidOperation:
            continue
        if len(tokens) < 2 or not tokens[1].isdigit():
            continue
        if not version.is_finite():
            raise ValueError("Nonfinite independent source version")
        rows[f"{family}:{namespace}:{line}"] = {
            "line": line,
            "family": family,
            "namespace": namespace,
        }
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--require-full-source", action="store_true")
    args = parser.parse_args()
    source = load_pcff_source(args.source)
    ledger = pcff_source_row_ledger(source)
    independent = independently_index(source.raw)
    if set(independent) != {r["id"] for r in ledger["rows"]}:
        raise ValueError("Independent row accounting mismatch")
    ledger["independent_counts"] = dict(
        Counter(r["family"] for r in independent.values())
    )
    ledger["msi2lmp_scope"] = {
        "typing": False,
        "bond_increments": False,
        "auto_equivalence": False,
        "supplied_type_direct_and_ordinary_resolution": True,
        "wildcards": "exact pass then first file-order wildcard; reversal can be orientation dependent; not universal physical precedence",
        "native_runtime_dependency": False,
    }
    storage.publish(Path(args.output), storage.json_bytes(ledger))
    print(
        {
            "rows": len(independent),
            "counts": ledger["counts"],
            "full_source_complete": False,
        }
    )
    return 1 if args.require_full_source else 0


if __name__ == "__main__":
    raise SystemExit(main())
