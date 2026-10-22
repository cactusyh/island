"""Reproducible J33 fixture selection.

The chemistry and packing are the reviewed complete J32 fixtures.  This
module only adds an explicit J33 namespace and derived NPT configuration; it
never copies a trajectory or a thermodynamic measurement.
"""

from dataclasses import dataclass

from island.analysis.npt import NPTConfiguration
from scripts.phase_4j32_real_fixture import build_fixture


@dataclass(frozen=True)
class J33Fixture:
    kind: str
    fixture: object
    npt: NPTConfiguration


def build_j33_fixture(kind: str = "linear", family: str = "OPLS-AA") -> J33Fixture:
    if kind not in {"linear", "network"}:
        raise ValueError("J33 fixture must be linear or network")
    if family not in {"OPLS-AA", "GAFF2", "PCFF"}:
        raise ValueError("J33 fixture family is unsupported")
    fixture = build_fixture(kind)
    backend_identity = fixture.packed.graph.identity
    source_identity = "0" * 64
    config = NPTConfiguration(
        family=family,
        source_identity=source_identity,
        final_graph_identity=backend_identity,
        packing_identity=fixture.packed.plan.identity,
        box_lengths_nm=tuple(float(x) / 10 for x in fixture.packed.system.box.lengths),
    )
    return J33Fixture(kind, fixture, config)


__all__ = ["J33Fixture", "build_j33_fixture"]
