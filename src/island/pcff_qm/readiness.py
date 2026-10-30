"""J57 family-only planning API; no final PCFF profile or MD authority."""

from dataclasses import dataclass

from island.pcff_completion import _verified, classify, require
from island.pcff_completion.readiness import prepare_polymer_forcefield as prior_prepare
from island.pcff_completion.targets import FIXTURE_GRAPHS

from .records import Record, attestations, confirmation_inventory, record
from .reference import fit_readiness
from .targets import targets


def closure(system, assessment, fixture):
    p = _verified(system, assessment)
    historical = classify(system, assessment)
    confirmations = {r.payload["query"]: r.payload for r in attestations()}
    effective, missing = [], []
    for row in historical["terms"]:
        if row["classification"] == "author_confirmed":
            require(
                row["query"] in confirmations, "attestation_query_coverage_mismatch"
            )
            effective.append(
                {
                    **row,
                    "historical_status": row["status"],
                    "status": "author_confirmed_audited",
                    "classification": "author_confirmed_audited",
                    "confirmation_identity": confirmations[row["query"]]["identity"],
                    "authority": "user_reported_author_confirmation",
                }
            )
        else:
            require(
                row["classification"] == "qm_fit_required", "unexpected_missing_term"
            )
            missing.append(row)
    ts = targets(system, assessment, fixture)
    return record(
        {
            "schema": "island.pcff-parameter-closure.v2",
            "family": "pcff",
            "selected_profile": None,
            "candidate_profile": p["candidate_profile"],
            "graph_identity": p["graph_identity"],
            "typing_identity": p["typing_identity"],
            "charge_identity": p["charge_identity"],
            "charge_origin": p["charge_origin"],
            "source_hashes": p["source_hashes"],
            "source_identity": p["source_record_identity"],
            "author_confirmation": confirmation_inventory(),
            "effective_execution_contract": effective,
            "user_attestation_blocks_conventions": False,
            "original_artifact_absence_blocks_conventions": False,
            "missing_parameter_blockers": missing,
            "required_dataset_targets": [t.payload for t in ts],
            "fitting_readiness": [fit_readiness(t).payload for t in ts],
            "interaction_inventory": p["interaction_inventory"],
            "domain": p["domain"],
            "typing_status": p["typing_status"],
            "charge_status": p["charge_status"],
            "parameter_completeness": "incomplete",
            "parameter_complete": False,
            "lammps_export_status": "blocked",
            "lammps_runtime_status": "blocked",
            "openmm_build_status": "blocked",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "status": "parameter_incomplete",
            "registered_profiles": [],
            "fitted_parameters": [],
            "QM_executed": False,
            "MD_authorized": False,
            "MD_executed": False,
            "completeness_requirements": [
                "28 immutable user-attested conventions bound",
                "19 target terms supplied, fitted and validated",
                "all Class-II dependencies closed",
                "complete typing and charges",
                "complete units and functional forms",
                "complete fit provenance and uncertainty",
                "independent reconstruction passed",
            ],
        }
    )


@dataclass(frozen=True)
class PreparationReport(Record):
    assessment: object

    def validate_integrity(self, system):
        p = self.payload
        fixture = next(
            (k for k, v in FIXTURE_GRAPHS.items() if v == p["graph_identity"]), None
        )
        require(
            p == closure(system, self.assessment, fixture).payload,
            "closure_evidence_mutation",
        )
        return p


def prepare_polymer_forcefield(
    system, *, forcefield="pcff", assets=None, audit_si_h=False
):
    """Automatically inspect the exact final graph; no profile or charge argument.

    Returns an incomplete preparation report, never registers a partial profile.
    """
    require(forcefield == "pcff", "wrong_profile_family")
    previous = prior_prepare(
        system, forcefield=forcefield, assets=assets, audit_si_h=audit_si_h
    )
    a = previous.assessment
    fixture = next(
        (k for k, v in FIXTURE_GRAPHS.items() if v == a.payload["graph_identity"]), None
    )
    require(fixture is not None, "bounded_fixture_required")
    r = closure(system, a, fixture)
    return PreparationReport(r.json_text, a)
