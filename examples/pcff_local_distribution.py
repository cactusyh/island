"""Explicit local source inspection; this cannot authorize a new PCFF model."""

import argparse

from island.forcefields.pcff import inspect_pcff_distribution


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for role in ("frc", "rlb", "templates"):
        p.add_argument("--" + role, required=True)
        p.add_argument("--" + role + "-sha256", required=True)
    args = p.parse_args()
    record = inspect_pcff_distribution(
        frc_path=args.frc,
        rlb_path=args.rlb,
        templates_path=args.templates,
        expected_hashes={
            r: getattr(args, r + "_sha256") for r in ("frc", "rlb", "templates")
        },
    )
    v = record.payload
    print(
        {
            "identity": record.identity,
            "numerical_source_unchanged": v["frc"]["identical_to_historical_pin"],
            "template_records": v["templates"]["counts"],
            "unsupported_semantics": v["templates"]["unsupported_semantics"],
            "label_differences": {
                k: v["links"][k]
                for k in (
                    "template_labels_without_frc_atom",
                    "frc_labels_without_template",
                )
            },
            "runtime_authorization": v["runtime_authorization"],
        }
    )


if __name__ == "__main__":
    main()
