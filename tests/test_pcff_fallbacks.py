"""Explicitly synthetic software tests; real-source evidence is kept separately."""

from copy import deepcopy
from dataclasses import replace

import numpy as np
import pytest
from test_pcff_automatic import explicit, synthetic  # noqa: F401
from test_pcff_model import parameters  # noqa: F401

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.forcefields.pcff import (
    assign_automatic_pcff_charges,
    assign_pcff_parameters,
    define_pcff_model,
    load_pcff_parameters,
    save_pcff_parameters,
    special_pair_policy,
    type_pcff_atoms,
)
from island.forcefields.pcff.fallbacks import POLICY, lookup
from island.forcefields.pcff.source import PCFFSource, digest


def row(family, types, values, version="1", rid="x", namespace="cff91_auto"):
    return {
        "id": rid,
        "section": family,
        "namespace": namespace,
        "types": types,
        "normalized_values": values,
        "version": version,
    }


def test_lookup_versions_specificity_conflicts_and_zero():
    f = "quadratic_angle"
    rows = {
        str(i): r
        for i, r in enumerate(
            [
                row(f, ["*", "c_", "*"], [1.0, 4.0], rid="generic"),
                row(f, ["h_", "c_", "*1"], [1.2, 0.0], rid="old"),
                row(f, ["h_", "c_", "*1"], [1.3, 0.0], version="2", rid="new"),
            ]
        )
    }
    q = lookup(f, ["h_", "c_", "o_"], "cff91_auto", rows, [])
    assert q["status"] == "assigned" and q["normalized_values"] == [1.3, 0.0]
    assert len(q["candidates"]) >= 3
    rows["conflict"] = row(f, ["*2", "c_", "o_"], [1.4, 6.0], rid="conflict")
    assert (
        lookup(f, ["h_", "c_", "o_"], "cff91_auto", rows, [])["status"] == "ambiguous"
    )
    assert lookup(f, ["h_", "n_", "o_"], "cff91_auto", rows, [])["status"] == "missing"
    assert (
        lookup(
            f, ["h_", "c_", "o_"], "cff91_auto", dict(reversed(list(rows.items()))), []
        )["status"]
        == "ambiguous"
    )


def test_directional_blocks_and_fixed_centers():
    rows = {"ba": row("bond-angle", ["a", "b", "c"], [2.0, 7.0], namespace="cff91")}
    result = lookup("bond-angle", ["c", "b", "a"], "cff91", rows, [])
    assert result["normalized_values"] == [7.0, 2.0]
    assert result["selected"][0]["permutation"] == [2, 1, 0]
    rows = {"aa": row("angle-angle", ["a", "b", "c", "d"], [3.0], namespace="cff91")}
    assert (
        lookup("angle-angle", ["d", "b", "c", "a"], "cff91", rows, [])["status"]
        == "assigned"
    )
    assert (
        lookup("angle-angle", ["d", "c", "b", "a"], "cff91", rows, [])["status"]
        == "missing"
    )
    rows = {"oop": row("wilson_out_of_plane", ["*", "center", "*", "*"], [3.0, 0.0])}
    assert (
        lookup(
            "wilson_out_of_plane", ["a", "center", "b", "c"], "cff91_auto", rows, []
        )["status"]
        == "assigned"
    )
    assert (
        lookup(
            "wilson_out_of_plane", ["center", "a", "b", "c"], "cff91_auto", rows, []
        )["status"]
        == "missing"
    )


@pytest.fixture
def fallback_source(parameters, monkeypatch):  # noqa: F811
    from island.forcefields.pcff import expanded, model, source

    raw = parameters.source.raw.decode()
    start = raw.index("#quartic_bond")
    end = raw.index("#", start + 1)
    raw = (
        raw[:start]
        + """#quadratic_bond cff91_auto
1.0 1 c c 1.53 300
1.0 1 c h 1.1 200
1.0 1 h h 0.75 250
"""
        + raw[end:]
    )
    raw = raw.replace(
        "#auto_equivalence cff91_auto",
        "#auto_equivalence cff91_auto\n9 1 c3 c c c c c c c c c\n9 1 hc h h h h h h h h h",
    )
    raw = raw.encode()
    for mapping, key in [
        (source.PIN, "sha256"),
        (expanded.PROFILE, "source_sha256"),
        (model.PROFILE, "frc_sha256"),
    ]:
        monkeypatch.setitem(mapping, key, digest(raw))
    return PCFFSource(raw, digest(raw))


def make(fallback_source):
    m = explicit([(0, 1)], [3, 3])
    for i, sid in enumerate(sorted(m.topology.sites)):
        m.coordinates.set(sid, [i * 0.3, (i % 3) * 0.7, (i % 2) * 0.5])
    t = type_pcff_atoms(m, fallback_source, profile="island_pcff_source_graph_v1")
    q = assign_automatic_pcff_charges(m, t, resolution_policy=POLICY)
    a = assign_pcff_parameters(m, t, q, resolution_policy=POLICY)
    spec = define_pcff_model(
        a, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    return m, t, q, a, spec


def test_complete_fallback_dependencies_persistence_and_tampering(
    fallback_source, tmp_path
):
    m, t, q, a, spec = make(fallback_source)
    assert spec.payload["model_definition_complete"]
    assert any(x["family"] == "quadratic_bond" for x in a.payload["assignments"])
    for term in a.payload["assignments"]:
        for dep in term["dependencies"]:
            if dep["assignment_id"].startswith("quartic_bond"):
                assert all(
                    "quadratic_bond:cff91_auto:" in r for r in dep["source_rows"]
                )
    with pytest.raises(PCFFError, match="policy mismatch"):
        assign_pcff_parameters(m, t, q)
    old = assign_pcff_parameters(m, t, assign_automatic_pcff_charges(m, t))
    assert not define_pcff_model(
        old, special_pairs=spec.payload["special_pairs"]
    ).payload["model_definition_complete"]
    save_pcff_parameters(a, tmp_path / "assignment.json")
    assert (
        load_pcff_parameters(
            tmp_path / "assignment.json", fallback_source, system=m
        ).identity
        == a.identity
    )
    for field, value in [
        ("resolution_policy", "invented"),
        ("normalized_values", [1, 999]),
    ]:
        data = unpack(a.json_text)
        if field == "resolution_policy":
            data[field] = value
        else:
            next(x for x in data["assignments"] if x["family"] == "quadratic_bond")[
                field
            ] = value
        with pytest.raises(PCFFError):
            type(a)(pack(data), fallback_source).validate_integrity(m)
    changed = deepcopy(m)
    changed.topology.sites[11].mass += 1
    with pytest.raises(PCFFError):
        a.validate_integrity(changed)
    # An absent coupling never becomes zero merely because its base is quadratic.
    raw = fallback_source.raw.replace(
        b"#bond-angle cff91", b"#uninterpreted_bond_angle cff91"
    )
    from island.forcefields.pcff import expanded, model, source

    # Use a local patch context: no historical identities are changed.
    with pytest.MonkeyPatch.context() as patch:
        for mapping, key in [
            (source.PIN, "sha256"),
            (expanded.PROFILE, "source_sha256"),
            (model.PROFILE, "frc_sha256"),
        ]:
            patch.setitem(mapping, key, digest(raw))
        *_, incomplete = make(PCFFSource(raw, digest(raw)))
        assert not incomplete.payload["model_definition_complete"]
        assert any("bond-angle" in d["id"] for d in incomplete.payload["diagnostics"])


@pytest.mark.parametrize(
    "family,coeff,xyz",
    [
        ("quadratic_bond", (1.3, 27.0), [[0, 0, 0], [1.5, 0.2, 0.1]]),
        ("quadratic_angle", (1.7, 12.0), [[0.1, 1.3, 0.2], [0, 0, 0], [1.5, 0.2, 0.1]]),
        (
            "torsion_1",
            (2.3, 3.0, 0.7),
            [[0.1, 1.3, 0.2], [0, 0, 0], [1.5, 0.2, 0.1], [2, 1.1, 1.4]],
        ),
    ],
)
def test_lower_order_compilation_fresh_fd_and_reversal(family, coeff, xyz):
    pytest.importorskip("openmm")
    from test_pcff_evaluation import compiled

    from island.forcefields.pcff.terms import (
        FallbackClass2Term,
        finite_difference_forces,
    )

    coords = dict(enumerate(np.array(xyz)))
    ids = list(coords)
    term = FallbackClass2Term(family, tuple(ids), coeff)
    data = {
        "schema": "island_pcff_source_model_v2",
        "terms": [
            {
                "family": family,
                "sites": ids,
                "coefficients": list(coeff),
                "equilibria": [],
            }
        ],
        "nonbonded": [{"site": i, "rmin": 1, "epsilon": 0, "charge": 0} for i in ids],
        "special_pair_inventory": [],
    }
    energy, forces = compiled(data, xyz)
    assert energy == pytest.approx(term.energy(coords), abs=1e-12)
    errors = []
    for h in (1e-3, 1e-4, 1e-5):
        fd = finite_difference_forces(term, coords, displacement=h)
        errors.append(np.max(np.abs(forces - np.array(list(fd.values())))))
    assert errors[-1] < 1e-7
    if family != "quadratic_bond":
        assert errors[-1] < errors[0]
    data["terms"][0]["sites"] = ids[::-1]
    reverse, rf = compiled(data, xyz)
    assert reverse == pytest.approx(energy, abs=1e-12)
    np.testing.assert_allclose(rf, forces, atol=1e-11)


def test_request_policy_and_native_identity(fallback_source):
    from island.evaluation.pcff_identity import pcff_evaluation_identity
    from island.forcefields import (
        ForceFieldRequestError,
        PCFFOptions,
        adopt_forcefield,
        create_evaluator,
    )

    for value in ("unknown", True, {}):
        with pytest.raises(ForceFieldRequestError):
            PCFFOptions("x", (0, 0, 1), (0, 0, 1), resolution_policy=value)
    m, _t, _q, _a, spec = make(fallback_source)
    settings, param, model = pcff_evaluation_identity(spec, system=m)
    assert settings["implementation"] == "island_pcff_fallback_singlepoint_v1"
    pytest.importorskip("openmm")
    p = adopt_forcefield(m, "pcff", spec)
    e = create_evaluator(m, p)
    assert (e.parameter_fingerprint, e.model_fingerprint) == (param, model)
    # Caller data ownership.
    payload = spec.payload
    payload["terms"].clear()
    assert spec.payload["terms"]
    assert (
        replace(
            PCFFOptions(
                "x",
                (0, 0, 1),
                (0, 0, 1),
                typing_profile="island_pcff_source_graph_v2",
                resolution_policy=POLICY,
            )
        ).resolution_policy
        == POLICY
    )


def test_automatic_charge_column_orientation_and_missing(synthetic, monkeypatch):  # noqa: F811
    from island.forcefields.pcff import expanded, source

    # Remove ordinary c/h rows, supply distinct charge-only automatic labels.
    raw = synthetic.raw.decode()
    lines = []
    section = ""
    for line in raw.splitlines():
        if line.startswith("#"):
            section = line.split()[0]
        fields = line.split()
        if section == "#bond_increments" and fields and fields[0][0].isdigit():
            continue
        lines.append(line)
    raw = (
        "\n".join(lines)
        .replace(
            "#auto_equivalence cff91_auto",
            "#auto_equivalence cff91_auto\n9 1 c3 c cq c c c c c c c\n9 1 hc h hq h h h h h h h",
        )
        .replace(
            "#end",
            "#bond_increments cff91_auto\n9 1 cq cq 0 0\n9 1 hq cq 0.1 -0.1\n#end",
        )
        .encode()
    )
    monkeypatch.setitem(source.PIN, "sha256", digest(raw))
    monkeypatch.setitem(expanded.PROFILE, "source_sha256", digest(raw))
    src = PCFFSource(raw, digest(raw))
    m = explicit([(0, 1)], [3, 3])
    typing = type_pcff_atoms(m, src, profile="island_pcff_source_graph_v1")
    assert not assign_automatic_pcff_charges(m, typing).complete
    q = assign_automatic_pcff_charges(m, typing, resolution_policy=POLICY)
    assert q.complete
    assert q.charges[11] == pytest.approx(-0.3)
    assert q.charges[100] == pytest.approx(0.1)
    assert all(
        x["path"] == "auto_equivalence.bond_increment"
        for x in q.payload["native_charge_record"]["contributions"]
    )
    data = unpack(q.json_text)
    data["resolution_policy"] = "renormalize"
    with pytest.raises(PCFFError):
        type(q)(pack(data), src).validate_integrity()


def test_offline_bundle_fallback_identity_and_coordinate_replacement(
    fallback_source, tmp_path, monkeypatch
):
    import builtins

    from island.evaluation.pcff_identity import pcff_evaluation_identity
    from island.forcefields import (
        PreparedForceFieldSources,
        adopt_forcefield,
        load_prepared_forcefield,
        save_prepared_forcefield,
    )

    m, *_, spec = make(fallback_source)
    frc = tmp_path / "synthetic.frc"
    frc.write_bytes(fallback_source.raw)
    sources = PreparedForceFieldSources(pcff_frc=frc)
    original = adopt_forcefield(m, "pcff", spec)
    save_prepared_forcefield(m, original, tmp_path / "bundle", sources=sources)
    importer = builtins.__import__

    def blocked(name, *args, **kwargs):
        if name.split(".")[0] in ("openmm", "scipy", "rdkit", "foyer", "parmed"):
            raise AssertionError("offline scientific dependency: " + name)
        return importer(name, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(builtins, "__import__", blocked)
        loaded = load_prepared_forcefield(tmp_path / "bundle", sources=sources)
        assert loaded.prepared.identity == original.identity
        a = pcff_evaluation_identity(spec, system=m)
        moved = deepcopy(m)
        moved.coordinates.translate([1, 2, 3])
        assert (
            pcff_evaluation_identity(loaded.prepared.native_result, system=moved) == a
        )
    loaded.system.metadata["caller"] = "new"
    assert "caller" not in loaded.system.metadata
