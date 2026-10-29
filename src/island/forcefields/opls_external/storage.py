"""Explicit, data-only preparation replay; a file cannot authorize itself."""

from hashlib import sha256
from pathlib import Path

from island.forcefields.external_pdms import require
from island.forcefields.pdms_literature import pinned_bytes

from .assignment import ExternalOPLSPreparation, prepare_external_opls
from .source import ExternalOPLSError


def save_preparation(preparation, path):
    require(
        type(preparation) is ExternalOPLSPreparation, "external_preparation_required"
    )
    preparation.validate_integrity()
    raw = preparation.json_text.encode()
    with Path(path).open("xb") as stream:
        stream.write(raw)
    return sha256(raw).hexdigest()


def load_preparation(path, expected_sha256, system, options):
    """Replay graph, assets, typing, charges and every row before accepting data."""
    try:
        raw = pinned_bytes(Path(path), expected_sha256)
        replay = prepare_external_opls(system, options)
        require(raw == replay.json_text.encode(), "stale_or_tampered_preparation")
        return replay
    except ExternalOPLSError:
        raise
    except Exception as error:
        raise ExternalOPLSError(
            "parameter_incomplete",
            [{"stage": "preparation_reconstruction_failed", "reason": str(error)}],
        ) from error
