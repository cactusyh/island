"""Opt-in external OPLS sources, distinct from the historical pinned native API."""

from .assignment import (
    FAMILY,
    ExternalOPLSOptions,
    ExternalOPLSPreparation,
    prepare_external_opls,
)
from .graph import extension_graph
from .source import (
    ExternalOPLSError,
    ExternalOPLSSource,
    load_source,
    write_source_manifest,
)

__all__ = [
    "FAMILY",
    "ExternalOPLSError",
    "ExternalOPLSOptions",
    "ExternalOPLSPreparation",
    "ExternalOPLSSource",
    "extension_graph",
    "load_source",
    "prepare_external_opls",
    "write_source_manifest",
]
