"""Immutable user-reported author confirmation, without invented author text.

The primary author's document/date/values were not supplied in this session.
The user's correction changes planning statuses; it is not a numerical row or
an executable Hamiltonian. A primary artifact may be validated separately.
"""

from datetime import date
from hashlib import sha256
from pathlib import Path

from island.forcefields.pcff.mixed_sip import DISTRIBUTION_PINS

from . import digest, inventory, require, signed

STATES = (
    "author_confirmed",
    "qm_fit_required",
    "derived_external",
    "fit_completed",
    "fit_validated",
    "parameter_incomplete",
)
ARTIFACT = """The author has now confirmed the previously unresolved source conventions and
parameter interpretation. Therefore, the 28 items previously classified as
`requires_author_confirmation` must no longer block the PCFF-PDMS workflow.

Update their status to:

    author_confirmed

The remaining 19 items genuinely require QM fitting because their numerical
values are not supplied by the author and cannot be derived from existing
source rows.
"""


def reported_confirmation():
    """Preserve exactly the supplied assertion and its limited authority."""
    return signed(
        {
            "schema": "island.pcff-author-confirmation.v1",
            "artifact": ARTIFACT,
            "artifact_sha256": sha256(ARTIFACT.encode()).hexdigest(),
            "evidence_checksum": digest(
                {"artifact": ARTIFACT, "kind": "user_reported_author_confirmation"}
            ),
            "source_hashes": DISTRIBUTION_PINS,
            "profile_identity": "island_external_pcff_pdms_class2_v1",
            "artifact_kind": "user_reported_author_confirmation",
            "recorded_date": "2026-10-30",
            "author_confirmation_date": None,
            "primary_author_artifact": None,
            "primary_artifact_status": "not_supplied",
            "scope": "28 previously author-confirmation-classified directed FRC queries only",
            "affected_interaction_families": [
                "angle-angle-torsion_1",
                "angle-torsion_3",
                "end_bond-torsion_3",
                "middle_bond-torsion_3",
            ],
            "status": "author_confirmed",
            "confirmed_interpretations": None,
            "numerical_values": None,
            "execution_authority": False,
            "limitation": "The reported confirmation changes planning status. Exact confirmed "
            "equations, applicability, values and author date require the primary artifact; "
            "none is inferred from family or backend defaults.",
        }
    )


def validate_reported(record):
    require(record == reported_confirmation(), "author_confirmation_mutation")
    return record


def validate_primary(record, root, expected_sha256, queries):
    """Validate an externally pinned author record and every exact interpretation.

    Checks bind evidence, not the identity/truthfulness of an email's sender.
    No package import or validation confers redistribution rights.
    """
    require(digest(record) == expected_sha256, "confirmation_record_hash_mutation")
    require(record.get("family") == "pcff", "wrong_profile_family")
    hashes = inventory(root)
    require(record.get("files") == hashes, "confirmation_inventory_mutation")
    require(record.get("author") and record.get("scope"), "confirmation_scope_required")
    date.fromisoformat(record["date"])
    require(
        record.get("source_hashes") == DISTRIBUTION_PINS, "confirmation_source_mutation"
    )
    entries = record.get("entries", [])
    require(
        len(entries) == len(queries) and {e["query"] for e in entries} == set(queries),
        "confirmation_query_coverage_mismatch",
    )
    for e in entries:
        require(e.get("status") == "author_confirmed", "confirmation_status_mutation")
        require(
            e.get("interpreted_value") is not None
            and e.get("units")
            and e.get("location")
            and e.get("raw_evidence"),
            "confirmed_interpretation_required",
        )
        name = e.get("source_file")
        require(
            name in hashes and e.get("source_sha256") == hashes[name],
            "confirmation_source_hash_mutation",
        )
        require(
            e["raw_evidence"] in (Path(root) / name).read_text(),
            "confirmation_evidence_not_in_artifact",
        )
        require(
            e.get("evidence_checksum")
            == digest({k: v for k, v in e.items() if k != "evidence_checksum"}),
            "confirmation_evidence_mutation",
        )
    return signed(
        {
            "schema": "island.pcff-primary-confirmation-validation.v1",
            "record_sha256": expected_sha256,
            "files": hashes,
            "status": "author_confirmed",
            "entries": entries,
            "execution_authority": False,
        }
    )
