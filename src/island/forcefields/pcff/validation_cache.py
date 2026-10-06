"""Bounded, process-local reuse of successful immutable-record validation.

No result, live system, Context or external file handle is cached. Mutable caller
systems are checked on every entry. Content or trusted-rule changes invalidate
reuse; failed validations are never cached. Persisted records remain unchanged.
"""

import sys
from collections import OrderedDict
from functools import wraps
from threading import RLock
from types import FunctionType

from island.charge_references.records import unpack

from ._validation_state import callable_state

MAX_ENTRIES = 16
MAX_BYTES = 64 * 1024 * 1024
_ENTRIES = OrderedDict()
_LOCK = RLock()
_STATS = {"hits": 0, "misses": 0, "evictions": 0, "graph_checks": 0}
_MODULES = (
    "source",
    "charges",
    "automatic",
    "expanded",
    "domains",
    "organic_domains",
    "class2",
    "catalog",
    "fallbacks",
    "model",
)


def clear_pcff_validation_cache():
    with _LOCK:
        _ENTRIES.clear()
        _STATS.update(dict.fromkeys(_STATS, 0))


def pcff_validation_cache_info():
    with _LOCK:
        return {
            **_STATS,
            "entries": len(_ENTRIES),
            "bytes": sum(v[2] for v in _ENTRIES.values()),
            "max_entries": MAX_ENTRIES,
            "max_bytes": MAX_BYTES,
        }


def _context():
    constants, functions = [], []
    for name in _MODULES:
        mod = sys.modules.get("island.forcefields.pcff." + name)
        if mod is None:
            constants.append((name, "not_imported"))
            continue
        for key, value in sorted(vars(mod).items()):
            if key.isupper() and type(value) in (
                dict,
                tuple,
                list,
                set,
                str,
                float,
                int,
            ):
                constants.append((name, key, repr(value)))
            if isinstance(value, FunctionType):
                functions.append(callable_state(value))
            elif isinstance(value, type) and value.__module__ == mod.__name__:
                for member in vars(value).values():
                    if isinstance(member, FunctionType):
                        functions.append(callable_state(member))
                    elif isinstance(member, property):
                        functions.append(callable_state(member.fget))
                    elif isinstance(member, (staticmethod, classmethod)):
                        functions.append(callable_state(member.__func__))
    return tuple(constants), tuple(functions)


def _freeze(value):
    # Typed tokens preserve the distinctions enforced by existing pack comparison.
    if type(value) is dict:
        return (
            "map",
            tuple(
                sorted(((_freeze(k), _freeze(v)) for k, v in value.items()), key=repr)
            ),
        )
    if type(value) in (tuple, list):
        return (type(value).__name__, tuple(_freeze(v) for v in value))
    return type(value).__name__, value


def _source_and_texts(record):
    from .source import PCFFSource, require

    texts = [record.json_text]
    owner = record
    if hasattr(owner, "assignment"):
        owner = owner.assignment
        texts.append(owner.json_text)
    source = owner.source
    require(
        type(source) is PCFFSource and type(source.raw) is bytes,
        "Immutable PCFF source required for validation",
    )
    require(all(type(s) is str for s in texts), "Immutable JSON records required")
    return source, tuple(texts)


def _graph(record):
    data = unpack(
        record.assignment.json_text
        if hasattr(record, "assignment")
        else record.json_text
    )
    if "charge_record" in data:
        return data["charge_record"]["automatic_typing"]["graph"]
    if "automatic_typing" in data:
        return data["automatic_typing"]["graph"]
    return data["graph"]


def validated_record(fn):
    """Reuse semantic validation only; always recheck the caller's final graph."""

    @wraps(fn)
    def call(self, system=None):
        from .automatic import chemical_graph
        from .charges import identity
        from .source import require

        source, texts = _source_and_texts(self)
        context = _context()
        key = (type(self), texts, source.raw, source.expected_sha256, context)
        with _LOCK:
            entry = _ENTRIES.get(key)
            if entry is not None:
                _ENTRIES.move_to_end(key)
                _STATS["hits"] += 1
        if entry is None:
            with _LOCK:
                _STATS["misses"] += 1
            source_bytes, source_hash = source.raw, source.expected_sha256
            fn(self, None)
            after_source, after_texts = _source_and_texts(self)
            require(
                after_texts == texts
                and after_source.raw == source_bytes
                and after_source.expected_sha256 == source_hash,
                "Immutable PCFF record changed during validation",
            )
            # Imports performed by cold validation can extend the trusted context.
            final_context = _context()
            prior = {
                (name, key): value
                for item in context[0]
                if len(item) == 3
                for name, key, value in (item,)
            }
            after = {
                (name, key): value
                for item in final_context[0]
                if len(item) == 3
                for name, key, value in (item,)
            }
            require(
                all(after.get(k) == v for k, v in prior.items())
                and all(f in final_context[1] for f in context[1]),
                "Trusted PCFF validation policy changed during validation",
            )
            graph = _freeze(_graph(self))
            content_identity = identity(unpack(self.json_text))
            size = sum(len(t.encode("utf-8")) for t in texts) + len(source.raw)
            entry = graph, content_identity, size
            key = (type(self), texts, source.raw, source.expected_sha256, final_context)
            if size <= MAX_BYTES:
                with _LOCK:
                    _ENTRIES[key] = entry
                    while (
                        len(_ENTRIES) > MAX_ENTRIES
                        or sum(v[2] for v in _ENTRIES.values()) > MAX_BYTES
                    ):
                        _ENTRIES.popitem(last=False)
                        _STATS["evictions"] += 1
        if system is not None:
            actual = _freeze(chemical_graph(system))
            require(actual == entry[0], "PCFF validated record chemical graph mismatch")
            with _LOCK:
                _STATS["graph_checks"] += 1

    return call


def record_identity(record):
    record.validate_integrity()
    source, texts = _source_and_texts(record)
    key = (type(record), texts, source.raw, source.expected_sha256, _context())
    with _LOCK:
        entry = _ENTRIES.get(key)
    if entry is not None:
        return entry[1]
    # Oversized records stay valid but never evade the declared memory ceiling.
    from .charges import identity

    return identity(unpack(record.json_text))


def profile_validation_summary(record):
    """Owned compact evidence for authorization, after full native validation."""
    from copy import deepcopy

    record.validate_integrity()
    source, texts = _source_and_texts(record)
    key = (type(record), texts, source.raw, source.expected_sha256, _context())
    with _LOCK:
        entry = _ENTRIES.get(key)
        summary = entry[3] if entry is not None and len(entry) > 3 else None
    if summary is None:
        model = unpack(record.json_text)
        assignment = unpack(record.assignment.json_text)
        auto = assignment["charge_record"]["automatic_typing"]
        summary = {
            "model": {
                k: model[k]
                for k in (
                    "schema",
                    "source",
                    "compatibility_profile",
                    "model_definition_complete",
                    "diagnostics",
                )
            },
            "assignment": {
                "parameter_coverage_complete": assignment[
                    "parameter_coverage_complete"
                ],
                "charge_record": {
                    "automatic_typing": {
                        "profile": {"name": auto["profile"]["name"]},
                        "assignments": auto["assignments"],
                    }
                },
            },
        }
        summary["model"]["terms"] = [
            {
                "family": t["family"],
                "origin": t["origin"],
                "source_rows": bool(t["source_rows"]),
                "coefficients": t["coefficients"]
                if t["family"] == "wilson_out_of_plane"
                else [],
            }
            for t in model["terms"]
        ]
        if entry is not None:
            with _LOCK:
                _ENTRIES[key] = (*entry[:3], summary)
    return deepcopy(summary)
