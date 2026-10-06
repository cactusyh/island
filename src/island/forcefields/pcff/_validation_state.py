"""In-memory callable state for cache keys, never persisted scientific identity."""

from types import FunctionType


def callable_state(function):
    """Snapshot code, defaults and closure data without invoking user methods.

    Function identity alone misses hot reloads that replace ``__code__`` and
    policy changes inside mutable defaults/closures. Cycles are represented by
    object identity; containers on other branches are still snapshotted fully.
    This is cache invalidation, not authentication against arbitrary Python code.
    """
    return _state(function, set())


def _state(value, active):
    kind = type(value)
    if value is None or kind in (str, bytes, bool, int, float):
        return kind, value
    oid = id(value)
    if oid in active:
        return "cycle", kind, oid
    active.add(oid)
    try:
        if kind is FunctionType:
            cells = []
            for cell in value.__closure__ or ():
                try:
                    contents = cell.cell_contents
                except ValueError:  # Empty cell in an incompletely built closure.
                    cells.append(("empty_cell",))
                else:
                    cells.append(_state(contents, active))
            return (
                value,
                value.__code__,
                _state(value.__defaults__, active),
                _state(value.__kwdefaults__, active),
                tuple(cells),
            )
        if kind is dict:
            return kind, frozenset(
                (_state(k, active), _state(v, active)) for k, v in value.items()
            )
        if kind in (list, tuple):
            return kind, tuple(_state(v, active) for v in value)
        if kind in (set, frozenset):
            return kind, frozenset(_state(v, active) for v in value)
        # Native policies use data containers, not arbitrary mutable objects.
        # Do not call repr, property getters or serialization on external objects.
        return kind, oid
    finally:
        active.remove(oid)
