"""Audit explicitly supplied local files; auxiliary syntax is not runtime authority."""

import argparse
from pathlib import Path

from island.forcefields.pcff.distribution import (
    inspect_pcff_distribution,
    save_pcff_distribution_audit,
)
from island.workflows import storage


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--frc", type=Path, required=True)
    p.add_argument("--rlb", type=Path, required=True)
    p.add_argument("--templates", type=Path, required=True)
    p.add_argument("--declaration", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    d = storage.read_json(a.declaration)
    hashes = {
        k: d["distribution"][n]["sha256"]
        for k, n in (
            ("frc", "pcff.frc"),
            ("rlb", "pcff.rlb"),
            ("templates", "pcff_templates.dat"),
        )
    }
    paths = {"frc_path": a.frc, "rlb_path": a.rlb, "templates_path": a.templates}
    result = inspect_pcff_distribution(**paths, expected_hashes=hashes)
    save_pcff_distribution_audit(a.output, result, **paths)
    data = result.payload
    print(
        {
            "audit_identity": result.identity,
            "frc": data["frc"]["semantic_counts"],
            "numerical_source_unchanged": data["frc"]["identical_to_historical_pin"],
            "templates": data["templates"]["counts"],
            "rlb": data["rlb"],
            "template_labels_without_frc_atom": data["links"][
                "template_labels_without_frc_atom"
            ],
            "frc_labels_without_template": data["links"]["frc_labels_without_template"],
            "runtime_authorization": data["runtime_authorization"],
        }
    )
    return 0  # This command reports inspection only; no executable-coverage gate.


if __name__ == "__main__":
    raise SystemExit(main())
