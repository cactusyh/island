"""Synthetic contracts; authentic source/numerical acceptance lives in the J17 runner."""

from copy import deepcopy
from dataclasses import FrozenInstanceError

# Imported pytest fixtures are intentionally named in test signatures.
# ruff: noqa: F811
import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields import (
    ForceFieldRequest,
    ForceFieldRequestError,
    PCFFOptions,
    PreparedBundleError,
    PreparedForceFieldError,
    PreparedForceFieldSources,
    adopt_forcefield,
    create_evaluator,
    load_prepared_forcefield,
    prepare_forcefield,
    save_prepared_forcefield,
)
from island.forcefields.pcff import (
    PCFFGraphCharges,
    PCFFTypedGraph,
    assign_pcff_parameters,
    assign_pcff_source_types,
    assign_typed_pcff_charges,
    bind_pcff_types,
    define_pcff_model,
    expanded,
    load_pcff_graph_record,
    model,
    provide_pcff_charges,
    save_pcff_graph_record,
    source,
    special_pair_policy,
)
from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY, MSI_POLICY
from island.forcefields.pcff.source import PCFFSource, digest


def bind(s, src, labels=None):
    return bind_pcff_types(
        s,
        src,
        labels or {i: a.element.lower() for i, a in s.topology.sites.items()},
        provenance="Synthetic independently supplied labels",
        evidence_references=["software_fixture_not_chemical_acceptance"],
    )


def provide(s, t, values=None, **kw):
    options = {
        "unit": "elementary_charge",
        "provenance": "Synthetic charge software fixture",
        "evidence_references": ["software_fixture"],
        "component_totals": {min(s.topology.sites): 0},
        "total_charge": 0,
        "resolution_policy": COMPATIBILITY_POLICY,
    }
    options.update(kw)
    return provide_pcff_charges(
        s,
        t,
        values if values is not None else dict.fromkeys(s.topology.sites, 0.0),
        **options,
    )


def test_independent_generic_source_family_and_no_perception(synthetic, monkeypatch):
    s = explicit([(0, 1)], [3, 3])
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", synthetic.identity["sha256"])
    t = bind(s, synthetic)
    with pytest.raises(PCFFError, match="contradict"):
        assign_pcff_source_types(s, synthetic, t.assignments, provenance="external")

    def forbidden(*a, **kw):
        raise AssertionError("automatic perception called")

    monkeypatch.setattr(expanded, "recognize", forbidden)
    t.validate_integrity(s)
    q = assign_typed_pcff_charges(s, t, resolution_policy=COMPATIBILITY_POLICY)
    assert q.complete and q.charges[11] == pytest.approx(-0.3)
    assert "automatic_typing" not in q.payload
    assert q.payload["origin"] == "native_increments"
    assert t.payload["validation"]["automatic_perception"] == "not_performed"
    assert t.payload["validation"]["external_chemical_authority"] == "caller_asserted"
    with pytest.raises(FrozenInstanceError):
        t.json_text = "changed"
    t.assignments.clear()
    assert len(t.assignments) == 8


def test_unavailable_automatic_perception_does_not_gate(synthetic, monkeypatch):
    from island import AtomSite, Coordinates, MolecularSystem, Topology
    from island.forcefields.pcff import type_pcff_atoms

    top = Topology()
    for i in (7, 23):
        top.add_site(AtomSite(i, "hydrogen", 1.008, element="H", atomic_number=1))
    top.add_bond(7, 23, order=1)
    s = MolecularSystem(top, Coordinates())
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", synthetic.identity["sha256"])
    assert not type_pcff_atoms(s, synthetic, profile=expanded.PROFILE_NAME).complete
    t = bind(s, synthetic)
    q = provide(s, t)
    assert q.complete  # Explicit h permits H parent; automatic profile has no rule.


@pytest.mark.parametrize(
    "kind",
    [
        "unknown",
        "element",
        "connectivity",
        "hydrogens",
        "missing",
        "extra",
        "bool",
        "nonstring",
    ],
)
def test_invalid_types(synthetic, kind):
    s = explicit([(0, 1)], [3, 3])
    labels = {i: a.element.lower() for i, a in s.topology.sites.items()}
    if kind == "unknown":
        labels[11] = "unknown"
    elif kind == "element":
        labels[11] = "h"
    elif kind == "connectivity":
        labels[100] = "c"
    elif kind == "hydrogens":
        labels[11] = "c2"
    elif kind == "missing":
        labels.pop(11)
    elif kind == "extra":
        labels[999] = "h"
    elif kind == "bool":
        labels[True] = labels.pop(11)
    elif kind == "nonstring":
        labels[11] = 4
    with pytest.raises(PCFFError):
        bind(s, synthetic, labels)


@pytest.mark.parametrize(
    "kind",
    [
        "nan",
        "inf",
        "bool",
        "str",
        "missing",
        "extra",
        "bool_id",
        "total",
        "component",
        "unit",
        "provenance",
        "evidence",
        "formal_contract",
    ],
)
def test_invalid_charges(synthetic, kind):
    s = explicit([(0, 1)], [3, 3])
    t = bind(s, synthetic)
    values = dict.fromkeys(s.topology.sites, 0.0)
    options = {}
    if kind in {"nan", "inf"}:
        values[11] = float(kind)
    elif kind == "bool":
        values[11] = True
    elif kind == "str":
        values[11] = "0"
    elif kind == "missing":
        values.pop(11)
    elif kind == "extra":
        values[999] = 0.0
    elif kind == "bool_id":
        values[True] = values.pop(11)
    elif kind == "total":
        options["total_charge"] = 1
    elif kind == "component":
        values[11] = 0.01
    elif kind == "unit":
        options["unit"] = "coulomb"
    elif kind == "provenance":
        options["provenance"] = " "
    elif kind == "evidence":
        options["evidence_references"] = []
    elif kind == "formal_contract":
        options["component_totals"] = {11: 1}
    with pytest.raises(PCFFError):
        provide(s, t, values, **options)


def test_component_charge_cancellation_cannot_hide_mismatch(synthetic):
    s = explicit([], [4, 4])
    t = bind(s, synthetic)
    values = dict.fromkeys(s.topology.sites, 0.0)
    values[11], values[16] = 1.0, -1.0
    with pytest.raises(PCFFError, match="component charge mismatch"):
        provide(s, t, values, component_totals={11: 0, 16: 0})


def test_exact_vector_and_graph_record_persistence(synthetic, tmp_path):
    s = explicit([(0, 1)], [3, 3])
    t = bind(s, synthetic)
    values = dict.fromkeys(s.topology.sites, 0.0)
    values[11], values[16] = 0.123456789012, -0.123456789012
    q = provide(s, t, values)
    assert q.charges == values
    assert q.payload["charge_data"]["native_increment_availability"] == "not_evaluated"
    for name, record in [("t", t), ("q", q)]:
        save_pcff_graph_record(record, tmp_path / name)
        loaded = load_pcff_graph_record(tmp_path / name, synthetic, system=s)
        assert loaded.identity == record.identity
    values[11] += 1
    assert q.charges[11] == 0.123456789012
    assert q.identity != provide(s, t).identity


@pytest.mark.parametrize("change", ["bond", "formal", "mass", "atom", "isotope"])
def test_graph_mutation_and_coordinates(synthetic, change):
    s = explicit([(0, 1)], [3, 3])
    t = bind(s, synthetic)
    q = provide(s, t)
    for i in s.topology.sites:
        s.coordinates.set(i, [i, 1, 2])
        s.topology.sites[i].metadata["repeat_unit_index"] = 999
    t.validate_integrity(s)
    q.validate_integrity(s)
    if change == "bond":
        s.topology.remove_bond(11, 16)
        s.topology.add_bond(11, 16, order=2)
    elif change == "formal":
        s.topology.sites[11].formal_charge = 1
    elif change == "mass":
        s.topology.sites[11].mass += 1
    elif change == "atom":
        s.topology.remove_site(100)
    else:
        s.topology.sites[11].metadata["isotope"] = 13
    with pytest.raises(PCFFError):
        t.validate_integrity(s)
    with pytest.raises(PCFFError):
        q.validate_integrity(s)


def test_rechecksummed_records_and_source_mismatch(synthetic):
    s = explicit([(0, 1)], [3, 3])
    t = bind(s, synthetic)
    q = provide(s, t)
    p = t.payload
    p["entries"][11]["environment"]["degree"] = 999
    with pytest.raises(PCFFError, match="Contradictory"):
        PCFFTypedGraph(pack(p), synthetic).validate_integrity()
    p = q.payload
    p["charge_data"]["partial_charges"][11] = 0.123
    with pytest.raises(PCFFError, match="Contradictory"):
        PCFFGraphCharges(pack(p), synthetic).validate_integrity()
    p = q.payload
    p["typed_graph"]["assignments"][11] = "c3"
    with pytest.raises(PCFFError):
        PCFFGraphCharges(pack(p), synthetic).validate_integrity()
    other = PCFFSource(synthetic.raw + b"\n", digest(synthetic.raw + b"\n"))
    with pytest.raises(PCFFError):
        PCFFTypedGraph(t.json_text, other).validate_integrity()


def test_provided_parameters_without_native_increments(parameters, monkeypatch):
    s = explicit([(0, 1)], [3, 3])
    raw = parameters.source.raw.replace(b"1.0 1 c h -0.1 0.1\n", b"")
    sha = digest(raw)
    monkeypatch.setitem(source.PIN, "sha256", sha)
    monkeypatch.setitem(model.PROFILE, "frc_sha256", sha)
    src = PCFFSource(raw, sha)
    t = bind(s, src)
    native = assign_typed_pcff_charges(s, t, resolution_policy=COMPATIBILITY_POLICY)
    assert not native.complete
    q = provide(s, t)
    a = assign_pcff_parameters(s, t, q, resolution_policy=COMPATIBILITY_POLICY)
    m = define_pcff_model(
        a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    assert m.payload["model_definition_complete"]
    assert "native_charge_record" not in a.payload["charge_record"]
    assert all(
        t["origin"] != "source_row" or t["source_rows"] for t in m.payload["terms"]
    )
    with pytest.raises(PCFFError, match="Complete typing and charges"):
        assign_pcff_parameters(s, t, native, resolution_policy=COMPATIBILITY_POLICY)
    # Valid supplied charges never manufacture missing source interactions.
    strict = provide(s, t, resolution_policy=MSI_POLICY)
    a = assign_pcff_parameters(s, t, strict, resolution_policy=MSI_POLICY)
    m = define_pcff_model(
        a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    assert not m.payload["model_definition_complete"]
    with pytest.raises(PreparedForceFieldError, match="Incomplete PCFF model"):
        adopt_forcefield(s, "pcff", m)


def test_public_preparation_bundle_session_and_tampering(parameters, tmp_path):
    s = explicit([(0, 1)], [3, 3])
    # Nonsingular software coordinates; no chemistry evidence claimed.
    import numpy as np

    rng = np.random.default_rng(78)
    for i in s.topology.sites:
        s.coordinates.set(i, rng.normal(size=3))
    src = parameters.source
    t = bind(s, src)
    q = provide(s, t)
    path = tmp_path / "source.frc"
    path.write_bytes(src.raw)
    options = PCFFOptions(
        path,
        (0, 0, 1),
        (0, 0, 1),
        resolution_policy=COMPATIBILITY_POLICY,
        typed_graph=t,
        graph_charges=q,
    )
    prepared = prepare_forcefield(s, ForceFieldRequest("pcff", options))
    assert prepared.metadata["charge_method"] == "external_types_provided_v1"
    sources = PreparedForceFieldSources(pcff_frc=path)
    save_prepared_forcefield(s, prepared, tmp_path / "bundle", sources=sources)
    relocated = tmp_path / "relocated"
    (tmp_path / "bundle").rename(relocated)
    path.rename(tmp_path / "moved.frc")
    sources = PreparedForceFieldSources(pcff_frc=tmp_path / "moved.frc")
    loaded = load_prepared_forcefield(relocated, sources=sources)
    assert loaded.prepared.identity == prepared.identity
    evaluator = create_evaluator(loaded.system, loaded.prepared)
    with evaluator.open_session() as session:
        assert session.evaluate().potential_energy == pytest.approx(
            session.evaluate_fresh().potential_energy
        )
    import json

    from island.workflows.storage import checksum, json_bytes

    data = unpack((relocated / "assignment.json").read_text())
    data["charge_record"]["charge_data"]["partial_charges"][11] = 0.1
    raw = pack(data).encode()
    (relocated / "assignment.json").write_bytes(raw)
    env = json.loads((relocated / "manifest.json").read_text())
    env["payload"]["files"]["assignment"]["sha256"] = checksum(raw)
    env["sha256"] = checksum(json_bytes(env["payload"]))
    (relocated / "manifest.json").write_bytes(json_bytes(env))
    with pytest.raises(PreparedBundleError, match="Contradictory"):
        load_prepared_forcefield(relocated, sources=sources)
    bad = deepcopy(s)
    bad.topology.sites[11].formal_charge = 1
    with pytest.raises(PreparedForceFieldError):
        prepared.validate_integrity(bad)
    with pytest.raises(ForceFieldRequestError):
        PCFFOptions(path, (0, 0, 1), (0, 0, 1), typed_graph=t)


@pytest.mark.parametrize(
    "element,label,number,order", [("O", "o=", 8, 2), ("N", "nt", 7, 3)]
)
def test_historical_wrong_chemistry(
    synthetic, monkeypatch, element, label, number, order
):
    from island import AtomSite, Coordinates, MolecularSystem, Topology
    from island.forcefields.pcff.expanded import MASSES

    raw = synthetic.raw.replace(
        b"#equivalence cff91",
        (
            f"1 1 {label} {MASSES[element]} {element} 1 synthetic\n#equivalence cff91"
        ).encode(),
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    src = PCFFSource(raw, digest(raw))
    top = Topology()
    for i in (1, 2):
        top.add_site(
            AtomSite(i, element, MASSES[element], element=element, atomic_number=number)
        )
    top.add_bond(1, 2, order=order)
    s = MolecularSystem(top, Coordinates())
    with pytest.raises(PCFFError, match="constraint"):
        bind(s, src, {1: label, 2: label})


@pytest.mark.parametrize("mutation", ["formal", "double", "aromatic"])
def test_contradictions_on_fresh_external_binding(synthetic, mutation):
    s = explicit([(0, 1)], [3, 3])
    if mutation == "formal":
        s.topology.sites[11].formal_charge = 1
    else:
        s.topology.remove_bond(11, 16)
        s.topology.add_bond(
            11,
            16,
            order=2 if mutation == "double" else 1.5,
            aromatic=mutation == "aromatic",
        )
    with pytest.raises(PCFFError):
        bind(s, synthetic)


@pytest.mark.parametrize("change", ["types", "charges"])
def test_external_input_changes_invalidate_checkpoint(parameters, change):
    import numpy as np

    from island.dynamics import (
        DynamicsOptions,
        DynamicsSegmentOptions,
        create_dynamics_checkpoint,
        resume_dynamics,
        run_dynamics_segment,
    )
    from island.evaluation import PCFFSinglePointEvaluator
    from island.exceptions import DynamicsCheckpointCompatibilityError

    s = explicit([(0, 1)], [3, 3])
    xyz = np.array(
        [
            [0, 0, 0],
            [1.4, 0.02, 0],
            [-0.3, 1, 0.1],
            [-0.3, -0.6, 0.9],
            [-0.3, -0.6, -0.9],
            [1.7, 1, 0.1],
            [1.7, -0.6, 0.9],
            [1.7, -0.6, -0.9],
        ]
    )
    for i, pos in zip(sorted(s.topology.sites), xyz, strict=True):
        s.coordinates.set(i, pos)
    t = bind(s, parameters.source)
    q = provide(s, t)

    def evaluator(typing, charges):
        a = assign_pcff_parameters(
            s, typing, charges, resolution_policy=COMPATIBILITY_POLICY
        )
        spec = define_pcff_model(
            a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
        )
        return PCFFSinglePointEvaluator(s, spec)

    e = evaluator(t, q)
    with e.open_session() as session:
        segment = run_dynamics_segment(
            s,
            session,
            dict.fromkeys(s.topology.sites, (0.0, 0.0, 0.0)),
            DynamicsOptions(
                timestep_fs=0.001, steps=1, max_evaluations=3, max_frames=2
            ),
        )
    checkpoint = create_dynamics_checkpoint(segment)
    if change == "types":
        labels = {
            i: "c3" if a.element == "C" else "hc" for i, a in s.topology.sites.items()
        }
        t = bind(s, parameters.source, labels)
        q = provide(s, t)
    else:
        values = q.charges
        values[11] = 0.01
        values[16] = -0.01
        q = provide(s, t, values)
    changed = evaluator(t, q)
    assert e.parameter_fingerprint != changed.parameter_fingerprint
    with (
        changed.open_session() as session,
        pytest.raises(DynamicsCheckpointCompatibilityError),
    ):
        resume_dynamics(checkpoint, s, session, DynamicsSegmentOptions(1, 3, 2))
