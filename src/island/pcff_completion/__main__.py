"""Inspect missing Class-II rows or submit a pinned external completion request."""

import argparse
import json
from pathlib import Path

from . import (
    CompletionError,
    CompletionRecord,
    canonical,
    classify,
    derive,
    fit,
    intake,
    load_assessment,
    pin,
    request_template,
    require,
    write,
)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "operation",
        choices=("classify", "template", "derive", "intake", "fit", "validate"),
    )
    parser.add_argument(
        "--bundle", type=Path, required=True, help="J54 exact graph audit bundle"
    )
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--distribution", type=Path, required=True)
    parser.add_argument(
        "--target-id", help="exact missing interaction ID for derivation"
    )
    parser.add_argument(
        "--request", type=Path, help="expert fitting contract outside payload directory"
    )
    parser.add_argument(
        "--request-sha256", help="independently pinned fitting contract"
    )
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--record", type=Path)
    parser.add_argument("--record-sha256")
    parser.add_argument("--json", type=Path, required=True)
    a = parser.parse_args(argv)
    try:
        system, assessment = load_assessment(
            a.bundle, a.manifest_sha256, a.distribution
        )
        if a.operation == "classify":
            value = classify(system, assessment)
        elif a.operation == "template":
            value = request_template(system, assessment, a.target_id)
        elif a.operation == "derive":
            require(a.target_id is not None, "target_id_required")
            value = derive(system, assessment, a.target_id)
        elif a.operation in ("intake", "fit"):
            require(
                a.request is not None and a.data_root is not None,
                "reference_data_unavailable",
            )
            require(pin(a.request) == a.request_sha256, "fitting_request_hash_mutation")
            request = json.loads(a.request.read_text())
            value = (fit if a.operation == "fit" else intake)(
                system, assessment, request, a.data_root
            )
        else:
            require(a.record is not None, "completion_record_required")
            require(pin(a.record) == a.record_sha256, "completion_record_hash_mutation")
            value = CompletionRecord(
                canonical(json.loads(a.record.read_text()))
            ).validate_integrity(system, assessment, data_root=a.data_root)
        # Never overwrite a source, graph bundle, contract, record or data file.
        output = a.json.resolve()
        roots = [a.bundle.resolve(), a.distribution.resolve()]
        if a.data_root is not None:
            roots.append(a.data_root.resolve())
        require(
            not any(output.is_relative_to(r) for r in roots), "output_overwrites_inputs"
        )
        require(
            output not in [p.resolve() for p in (a.request, a.record) if p],
            "output_overwrites_inputs",
        )
        write(output, value)
        print(
            json.dumps(
                {
                    "identity": value.get("identity"),
                    "parameter_completeness": "incomplete",
                    "status": value.get("status", "classified"),
                }
            )
        )
        return 2 if value.get("status") == "blocked" else 0
    except (ValueError, OSError, KeyError, TypeError) as error:
        # Errors are explicit; malformed input cannot produce a fitted parameter.
        print(
            json.dumps(
                {
                    "status": "rejected",
                    "reason": error.reason
                    if isinstance(error, CompletionError)
                    else str(error),
                    "parameter_completeness": "incomplete",
                }
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
