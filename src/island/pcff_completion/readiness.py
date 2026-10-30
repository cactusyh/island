"""Inspect corrected J56 readiness, validate confirmations or intake reference data."""

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

from island.forcefield_policy.claims import ClaimRequirements

from . import (
    CompletionError,
    classify,
    digest,
    load_assessment,
    pin,
    require,
    signed,
    write,
)
from .confirmation import reported_confirmation, validate_primary, validate_reported
from .qm import templates
from .references import reference_intake
from .targets import FIXTURE_GRAPHS, targets


def readiness(system, assessment, fixture):
    c = classify(system, assessment)
    return signed(
        {
            "schema": "island.pcff-fitting-readiness.v1",
            "family": "pcff",
            "candidate_profile": c["candidate_profile"],
            "selected_profile": None,
            "graph_identity": c["graph_identity"],
            "typing_identity": c["typing_identity"],
            "charge_identity": c["charge_identity"],
            "charge_origin": assessment.payload["charge_origin"],
            "source_hashes": c["source_hashes"],
            "author_confirmation": c["author_confirmation"],
            "effective_execution_contract": c["effective_execution_contract"],
            "missing_parameter_blockers": c["missing_parameter_blockers"],
            "targets": [t.payload for t in targets(system, assessment, fixture)],
            "typing_status": "complete",
            "charge_status": "complete",
            "parameter_completeness": "incomplete",
            "parameter_complete": False,
            "lammps_export_status": "blocked",
            "lammps_runtime_status": "blocked",
            "openmm_build_status": "blocked",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "reference_dataset_status": "not_supplied",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "MD_authorized": False,
            "primary_confirmation_artifact_status": "not_supplied",
            "execution_contract_numerical_interpretations": "not_supplied",
            "status": "parameter_incomplete",
        }
    )


@dataclass(frozen=True)
class FittingReadiness:
    """Corrected report; retains the immutable historical source assessment."""

    json_text: str
    assessment: object

    @property
    def payload(self):
        p = json.loads(self.json_text)
        require(
            p.get("identity")
            == digest({k: v for k, v in p.items() if k != "identity"}),
            "readiness_checksum_mutation",
        )
        return p

    def validate_integrity(self, system):
        self.assessment.validate_integrity(system)
        p = self.payload
        fixture = next(
            (k for k, v in FIXTURE_GRAPHS.items() if v == p["graph_identity"]), None
        )
        require(
            p == readiness(system, self.assessment, fixture),
            "readiness_evidence_mutation",
        )
        return p


def prepare_polymer_forcefield(
    system, *, forcefield="pcff", assets=None, audit_si_h=False
):
    """Family-only corrected planning API. No parameter profile is registered."""
    from island.pcff_pdms import prepare_polymer_forcefield as historical_prepare

    from . import canonical

    require(forcefield == "pcff", "wrong_profile_family")
    a = historical_prepare(
        system, forcefield=forcefield, assets=assets, audit_si_h=audit_si_h
    )
    fixture = next(
        (k for k, v in FIXTURE_GRAPHS.items() if v == a.payload["graph_identity"]), None
    )
    require(fixture is not None, "bounded_fixture_required")
    return FittingReadiness(canonical(readiness(system, a, fixture)), a)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "operation",
        choices=("readiness", "confirmations", "reference-intake", "qm-templates"),
    )
    parser.add_argument("--forcefield", default="pcff")
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--manifest-sha256", required=True)
    parser.add_argument("--distribution", type=Path, required=True)
    parser.add_argument(
        "--fixture", required=True, choices=("disiloxane", "hmds", "dp2", "dp4", "dp6")
    )
    parser.add_argument("--query")
    parser.add_argument("--request", type=Path)
    parser.add_argument("--request-sha256")
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--verification-sha256")
    parser.add_argument("--settings", type=Path)
    parser.add_argument("--scan", type=Path)
    parser.add_argument(
        "--engine", choices=("psi4", "pyscf", "orca", "cp2k", "gaussian")
    )
    parser.add_argument("--executable")
    parser.add_argument("--templates-dir", type=Path)
    for name in (
        "parameter-complete",
        "lammps-runtime",
        "openmm-runtime",
        "cross-engine-equivalence",
        "scientific-validation",
    ):
        parser.add_argument("--require-" + name, action="store_true")
    parser.add_argument("--json", type=Path, required=True)
    a = parser.parse_args(argv)
    try:
        require(a.forcefield == "pcff", "wrong_profile_family")
        roots = [a.bundle.resolve(), a.distribution.resolve()]
        if a.data_root:
            roots.append(a.data_root.resolve())
        for output in (a.json, a.templates_dir):
            if output is not None:
                require(
                    not any(output.resolve().is_relative_to(r) for r in roots)
                    and output.resolve()
                    not in [p.resolve() for p in (a.request, a.scan, a.settings) if p],
                    "output_overwrites_inputs",
                )
        s, assessment = load_assessment(a.bundle, a.manifest_sha256, a.distribution)
        if a.operation == "readiness":
            result = readiness(s, assessment, a.fixture)
        elif a.operation == "confirmations":
            if a.request:
                require(
                    pin(a.request) == a.request_sha256,
                    "confirmation_record_hash_mutation",
                )
                record = json.loads(a.request.read_text())
                queries = sorted(
                    {
                        r["query"]
                        for r in classify(s, assessment)["effective_execution_contract"]
                    }
                )
                # Primary records may cover this exact fixture's query subset.
                require(
                    record.get("source_hashes") == assessment.payload["source_hashes"],
                    "confirmation_source_mutation",
                )
                result = validate_primary(record, a.data_root, digest(record), queries)
            else:
                result = validate_reported(reported_confirmation())
        else:
            ts = targets(s, assessment, a.fixture)
            target = next((t for t in ts if t.payload["query"] == a.query), None)
            require(target is not None, "exact_qm_target_required")
            if a.operation == "reference-intake":
                require(a.request and a.data_root, "reference_data_unavailable")
                require(
                    pin(a.request) == a.request_sha256,
                    "reference_request_hash_mutation",
                )
                result = reference_intake(
                    s,
                    assessment,
                    target,
                    json.loads(a.request.read_text()),
                    a.data_root,
                    verification_sha256=a.verification_sha256,
                )
            else:
                require(
                    a.settings and a.scan and a.templates_dir and a.engine,
                    "qm_input_plan_required",
                )
                result = templates(
                    s,
                    assessment,
                    target,
                    json.loads(a.settings.read_text()),
                    json.loads(a.scan.read_text()),
                    a.templates_dir,
                    engine=a.engine,
                    executable=a.executable,
                )
        requirements = ClaimRequirements(
            parameter_complete=a.require_parameter_complete,
            lammps_runtime=a.require_lammps_runtime,
            openmm_runtime=a.require_openmm_runtime,
            cross_engine_equivalence=a.require_cross_engine_equivalence,
            scientific_validation=a.require_scientific_validation,
        )
        failures = requirements.failures(result)
        if any(vars(requirements).values()):
            result = signed(
                {
                    **{k: v for k, v in result.items() if k != "identity"},
                    "claim_failures": failures,
                }
            )
        write(a.json, result)
        print(
            json.dumps(
                {
                    "identity": result["identity"],
                    "status": result["status"],
                    "parameter_completeness": "incomplete",
                }
            )
        )
        return 2 if a.operation == "readiness" or failures else 0
    except (ValueError, OSError, KeyError, TypeError) as e:
        print(
            json.dumps(
                {
                    "status": "rejected",
                    "reason": e.reason if isinstance(e, CompletionError) else str(e),
                    "parameter_completeness": "incomplete",
                }
            )
        )
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
