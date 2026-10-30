"""J57 additive identities: user attestation is never an original author file."""

import json
from dataclasses import dataclass
from datetime import date
from pathlib import Path

from island.pcff_completion import (
    canonical,
    digest,
    inventory,
    pin,
    require,
    signed,
)
from island.pcff_completion.confirmation import reported_confirmation

VERSION = "island.pcff-qm-preparation.v1"
BASE = "5fd06ec0f25532ed164ba534b04a52a7d7a0aedb"
PREREQUISITE = "24ad7698f64d62aeb08ed893eb8372bdbbf33819"
BRANCH = "codex/phase-4j57-pcff-pdms-qm-input-generation"
AUDIT_SCOPE = (
    "the 28 conventions and parameter interpretations explicitly confirmed by the user"
)
USER_STATEMENT = (
    "The user has explicitly confirmed that the author has confirmed the 28 "
    "previously blocked conventions and parameter interpretations."
)
ENVIRONMENTS = (
    "c3,sio,c3,hc",
    "c3,sio,osi,sio",
    "osi,sio,c3,hc",
    "osi,sio,osi,sio",
    "sio,osi,sio,c3",
    "sio,osi,sio,hsi",
    "sio,osi,sio,osi",
)
CONFIRMED_FAMILIES = (
    "angle-angle-torsion_1",
    "angle-torsion_3",
    "end_bond-torsion_3",
    "middle_bond-torsion_3",
)
READINESS_STATES = (
    "qm_data_missing",
    "qm_data_available",
    "fit_ready",
    "fit_completed",
    "fit_validated",
    "fit_rejected",
)


def code_identity():
    return digest({p.name: pin(p) for p in sorted(Path(__file__).parent.glob("*.py"))})


@dataclass(frozen=True)
class Record:
    """Immutable serialized value; a checksum alone is not external authority."""

    json_text: str

    @property
    def payload(self):
        p = json.loads(self.json_text)
        require(
            p.get("identity")
            == digest({k: v for k, v in p.items() if k != "identity"}),
            "record_checksum_mutation",
        )
        return p


def record(value):
    return Record(canonical(signed(value)))


def attestations():
    prior = reported_confirmation()
    common = {
        "schema": "island.pcff-user-attested-confirmation.v1",
        "status": "author_confirmed_audited",
        "evidence_kind": "user_attested_author_confirmation",
        "authority": "user_reported_author_confirmation",
        "source_artifact_status": "not_yet_supplied",
        "audit_scope": AUDIT_SCOPE,
        "user_statement": USER_STATEMENT,
        "prior_confirmation_identity": prior["identity"],
        "source_hashes": prior["source_hashes"],
        "audit_result": "user assertion, scope and historical query membership bound; original author evidence not audited",
        "confirms_numerical_fit_values": False,
        "transfer_authority": "exact historical graph domain only",
    }
    return tuple(
        record(
            {
                **common,
                "query": f"{family}:{env}",
                "interaction_family": family,
                "evidence_checksum": digest({**common, "query": f"{family}:{env}"}),
            }
        )
        for family in CONFIRMED_FAMILIES
        for env in ENVIRONMENTS
    )


def confirmation_inventory():
    return signed(
        {
            "schema": "island.pcff-attestation-inventory.v1",
            "confirmations": [r.payload for r in attestations()],
            "original_author_artifacts": [],
            "source_artifact_status": "not_yet_supplied",
            "independently_hashed_author_evidence": [],
            "author_confirmed_audited": 28,
            "author_artifact_fields_not_invented": True,
        }
    )


def validate_attestation(value):
    require(
        any(r.payload == value for r in attestations()), "user_attestation_mutation"
    )
    return value


def original_artifact(attestation, submission, root, *, expected_sha256):
    """Add a separately pinned evidence record; never replace the attestation.

    This checks a user-submitted textual extraction and reconstruction, not the
    authenticity of a sender. Permission is supplied, not inferred from access.
    """
    validate_attestation(attestation)
    require(digest(submission) == expected_sha256, "author_submission_hash_mutation")
    require(
        submission.get("attestation_identity") == attestation["identity"],
        "author_attestation_binding_mutation",
    )
    files = inventory(root)
    require(files == submission.get("files"), "author_artifact_inventory_mutation")
    for key in (
        "filename",
        "date",
        "author_source_identity",
        "location",
        "exact_text",
        "permission_status",
        "independent_reconstruction",
    ):
        require(bool(submission.get(key)), "author_artifact_field_required", field=key)
    date.fromisoformat(submission["date"])
    require(
        submission["permission_status"] in ("licensed", "permission_granted"),
        "author_artifact_permission_required",
    )
    name = submission["filename"]
    require(
        name in files and submission.get("sha256") == files[name],
        "author_artifact_hash_mutation",
    )
    extraction = submission.get("text_extraction_file", name)
    require(extraction in files, "author_artifact_extraction_required")
    require(
        submission["exact_text"] in (Path(root) / extraction).read_text(),
        "author_text_not_in_evidence",
    )
    reconstruction = submission["independent_reconstruction"]
    require(
        isinstance(reconstruction, dict)
        and reconstruction.get("status") == "passed"
        and all(
            reconstruction.get(k)
            for k in ("method", "implementation_version", "result")
        ),
        "independent_reconstruction_required",
    )
    return record(
        {
            "schema": "island.pcff-original-author-artifact.v1",
            "evidence_kind": "original_author_artifact_with_hashed_extraction",
            "parent_attestation_identity": attestation["identity"],
            "submission": submission,
            "submission_checksum": expected_sha256,
            "evidence_checksum": digest(submission),
            "validation_scope": "file integrity and user-submitted reconstruction; sender authenticity not independently established",
        }
    )
