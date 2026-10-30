"""Inspect attestations, prepare scans, intake reference outputs and report closure."""

import argparse
import json
from pathlib import Path

from island.forcefield_policy.claims import ClaimRequirements
from island.pcff_completion import (
    CompletionError,
    canonical,
    load_assessment,
    pin,
    require,
    write,
)
from island.pcff_completion.qm import ADAPTERS, asset_status

from .inputs import generate
from .readiness import closure
from .records import (
    attestations,
    confirmation_inventory,
    digest,
    original_artifact,
    record,
)
from .reference import fit_readiness, intake
from .targets import FIXTURES, inventory, targets


def load_cases(repo, distribution):
    audit = Path(repo) / "docs/evidence/phase_4j54/audit"
    return {
        name: load_assessment(
            audit / name, pin(audit / name / "manifest.json"), distribution
        )
        for name in FIXTURES
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "operation",
        choices=(
            "confirmations",
            "artifacts",
            "targets",
            "generate",
            "validate-reference",
            "readiness",
            "closure",
            "assets",
            "attach-author",
        ),
    )
    p.add_argument("--forcefield", default="pcff")
    p.add_argument("--repo", type=Path, default=Path.cwd())
    p.add_argument("--distribution", type=Path)
    p.add_argument("--fixture", choices=FIXTURES)
    p.add_argument("--query")
    p.add_argument("--output", type=Path)
    p.add_argument("--engine", choices=tuple(ADAPTERS))
    p.add_argument("--executable")
    p.add_argument("--request", type=Path)
    p.add_argument("--request-sha256")
    p.add_argument("--data-root", type=Path)
    p.add_argument("--verification-sha256")
    p.add_argument("--controls", type=Path)
    p.add_argument("--json", type=Path, required=True)
    for claim in (
        "parameter-complete",
        "lammps-runtime",
        "openmm-runtime",
        "cross-engine-equivalence",
        "scientific-validation",
    ):
        p.add_argument("--require-" + claim, action="store_true")
    a = p.parse_args(argv)
    try:
        require(a.forcefield == "pcff", "wrong_profile_family")
        protected = [
            a.repo / "docs/evidence/phase_4j54",
            a.repo / "docs/evidence/phase_4j56",
        ]
        protected += [
            x for x in (a.distribution, a.data_root, a.request, a.controls) if x
        ]
        for dest in (a.json, a.output):
            if dest:
                require(
                    not any(
                        dest.resolve() == x.resolve()
                        or dest.resolve().is_relative_to(x.resolve())
                        for x in protected
                    ),
                    "output_overwrites_inputs",
                )
        if a.operation in ("confirmations", "artifacts"):
            result = confirmation_inventory()
        elif a.operation == "assets":
            result = record(
                {
                    "assets": [
                        asset_status(e, a.executable)
                        for e in ([a.engine] if a.engine else ADAPTERS)
                    ],
                    "QM_executed": False,
                }
            ).payload
        elif a.operation == "attach-author":
            require(a.request and a.data_root, "author_artifact_not_supplied")
            require(
                pin(a.request) == a.request_sha256, "author_submission_hash_mutation"
            )
            att = next(
                (r for r in attestations() if r.payload["query"] == a.query), None
            )
            require(att is not None, "exact_confirmed_query_required")
            submission = json.loads(a.request.read_text())
            result = original_artifact(
                att.payload, submission, a.data_root, expected_sha256=digest(submission)
            ).payload
        else:
            require(a.distribution is not None, "source_assets_unavailable")
            cases = load_cases(a.repo, a.distribution)
            all_targets = inventory(cases)
            selected = [r["target"] for r in all_targets.payload["targets"]]
            if a.fixture:
                s, assessment = cases[a.fixture]
                selected = [r.payload for r in targets(s, assessment, a.fixture)]
            if a.query:
                selected = [r for r in selected if r["query"] == a.query]
                require(len(selected) == 1, "exact_qm_target_required")
            if a.operation == "targets":
                result = record(
                    {"targets": selected, "status": "qm_fit_required"}
                ).payload
            elif a.operation in ("readiness", "closure"):
                names = [a.fixture] if a.fixture else FIXTURES
                result = record(
                    {
                        "fixtures": {n: closure(*cases[n], n).payload for n in names},
                        "fitting_readiness": [
                            fit_readiness(
                                record({k: v for k, v in t.items() if k != "identity"})
                            ).payload
                            for t in selected
                        ],
                        "parameter_completeness": "incomplete",
                        "status": "parameter_incomplete",
                        "cross_engine_equivalence": "not_claimed",
                        "scientific_validation": "not_established",
                        "production_readiness": "not_established",
                        "production_validated": False,
                        "simulation_readiness": "not_established",
                    }
                ).payload
            elif a.operation == "generate":
                require(
                    a.output is not None and not a.output.exists(),
                    "new_output_directory_required",
                )
                bundles = []
                for i, t in enumerate(selected):
                    s, assessment = cases[t["molecular_fixture"]]
                    target = record({k: v for k, v in t.items() if k != "identity"})
                    bundles.append(
                        generate(
                            s,
                            assessment,
                            target,
                            a.output / f"target-{i:02d}",
                            engines=(a.engine,) if a.engine else tuple(ADAPTERS),
                        )
                    )
                result = record(
                    {
                        "bundles": bundles,
                        "QM_executed": False,
                        "parameter_completeness": "incomplete",
                    }
                ).payload
            else:
                require(
                    len(selected) == 1 and a.request and a.data_root,
                    "exact_reference_target_and_data_required",
                )
                require(
                    pin(a.request) == a.request_sha256,
                    "reference_request_hash_mutation",
                )
                t = selected[0]
                target = record({k: v for k, v in t.items() if k != "identity"})
                dataset = intake(
                    *cases[t["molecular_fixture"]],
                    target,
                    json.loads(a.request.read_text()),
                    a.data_root,
                    verification_sha256=a.verification_sha256,
                )
                result = dataset.payload
                if a.controls:
                    readiness = fit_readiness(
                        target,
                        dataset,
                        json.loads(a.controls.read_text()),
                        system=cases[t["molecular_fixture"]][0],
                        assessment=cases[t["molecular_fixture"]][1],
                        root=a.data_root,
                        verification_sha256=a.verification_sha256,
                    )
                    result = record(
                        {
                            "dataset": dataset.payload,
                            **{
                                k: v
                                for k, v in readiness.payload.items()
                                if k != "identity"
                            },
                        }
                    ).payload
        claims = ClaimRequirements(
            parameter_complete=a.require_parameter_complete,
            lammps_runtime=a.require_lammps_runtime,
            openmm_runtime=a.require_openmm_runtime,
            cross_engine_equivalence=a.require_cross_engine_equivalence,
            scientific_validation=a.require_scientific_validation,
        )
        failures = claims.failures(result)
        if any(vars(claims).values()):
            result = record(
                {
                    **{k: v for k, v in result.items() if k != "identity"},
                    "claim_failures": failures,
                }
            ).payload
        write(a.json, result)
        print(canonical({"identity": result["identity"], "claim_failures": failures}))
        return 2 if failures else 0
    except (ValueError, OSError, KeyError, TypeError) as error:
        print(
            canonical(
                {
                    "status": "fit_rejected",
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
