"""Offline, audited parameter-bundle adapters; no runtime registration."""

from .adapter import ADAPTERS, BundleIngestion, adapt, proposals
from .intake import capture, ingest

__all__ = ["ADAPTERS", "BundleIngestion", "adapt", "capture", "ingest", "proposals"]
