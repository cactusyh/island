"""Periodic identity, semantic binding, native-family and durability controls."""

import builtins
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from island.exceptions import ValidationError
from island.graph.final import pack
from island.periodic import (
    PeriodicForceFieldConfig,
    PeriodicParameterizedSystem,
    _signed,
    build_openmm_periodic_system,
    export_lammps_data,
    export_lammps_input,
    load_periodic_backend_bundle,
    prepare_periodic_forcefield,
    save_periodic_backend_bundle,
)

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from periodic_j22_fixtures import RETAINED, retained_fixture


@pytest.fixture(scope="module", params=("pcff", "oplsaa", "gaff", "gaff2"))
def native(request):
    if not (RETAINED / request.param / "relocated").exists():
        pytest.skip("Retained real-source J4I2 bundle unavailable")
    return retained_fixture(request.param)


@pytest.fixture(scope="module")
def pcff():
    return retained_fixture("pcff")[0]


def config_kwargs():
    return {
        "family": "PCFF",
        "source": {"family": "PCFF", "sha256": "a" * 64},
        "box_lengths": (30, 30, 30),
        "nonbonded_cutoff": 10,
        "lj_mixing_rule": "sixthpower_9_6",
        "one_four_scaling": (1, 1),
        "electrostatics_method": "pme",
        "pme_tolerance": 1e-8,
        "provenance": "software contract test",
    }


@pytest.mark.parametrize(
    "field,value",
    [
        ("family", "unknown"),
        ("source", {"family": "GAFF", "sha256": "a" * 64}),
        ("source", {"sha256": "z" * 64}),
        ("box_lengths", (30, float("nan"), 30)),
        ("box_lengths", (30, 0, 30)),
        ("periodic", (True, False, True)),
        ("periodic", (1, 1, 1)),
        ("nonbonded_cutoff", 15),
        ("nonbonded_cutoff", True),
        ("lj_mixing_rule", "arithmetic"),
        ("one_four_scaling", (1, True)),
        ("electrostatics_method", "cutoff"),
        ("pme_tolerance", float("inf")),
        ("switching_policy", "automatic"),
        ("units", "metal"),
        ("neighbor_list_policy", "unknown"),
        ("excluded_pair_policy", "guess"),
    ],
)
def test_config_rejection(field, value):
    kw = config_kwargs()
    kw[field] = value
    with pytest.raises(ValidationError):
        PeriodicForceFieldConfig(**kw)


def test_configuration_owned_and_identity():
    kw = config_kwargs()
    c = PeriodicForceFieldConfig(**kw)
    identity = c.identity
    kw["source"]["sha256"] = "b" * 64
    assert c.identity == identity
    q = c.payload
    q["nonbonded_cutoff"] = 9
    with pytest.raises(ValidationError):
        PeriodicForceFieldConfig.from_json(pack(q))
    assert PeriodicForceFieldConfig.from_json(c.json_text).identity == identity


def test_real_family_binding_and_export(native):
    b, _ = native
    b.validate_integrity()
    assert b.final_graph.periodic
    assert b.parameter_assignment.payload["diagnostics"][
        "executable_potential_complete"
    ]
    assert b.payload["production_validated"] is False
    assert b.payload["simulation_readiness"] == "not_established"
    assert sum(b.graph_charges.charges.values()) == pytest.approx(
        b.graph_charges.payload["total_charge"], abs=1e-12
    )
    data = export_lammps_data(b)
    script = export_lammps_input(b)
    assert data == export_lammps_data(b)
    assert "boundary p p p" in script and "kspace_style pppm" in script
    assert "table 0" in script and "special_bonds lj 0 0" in script
    atoms = data.split("Atoms # full\n\n")[1].split("\n\n")[0].splitlines()
    assert [int(r.split()[0]) for r in atoms] == sorted(b.system.topology.sites)
    assert [float(r.split()[3]) for r in atoms] == [
        b.graph_charges.charges[i] for i in sorted(b.system.topology.sites)
    ]
    assert all(int(r.split()[1]) == 1 for r in atoms)


def test_openmm_periodic(native):
    pytest.importorskip("openmm")
    b, _ = native
    m = build_openmm_periodic_system(b)
    assert m.getNumParticles() == len(b.graph_charges.charges)
    assert m.usesPeriodicBoundaryConditions()
    from openmm import NonbondedForce

    nb = next(f for f in m.getForces() if isinstance(f, NonbondedForce))
    assert nb.getNonbondedMethod() == NonbondedForce.PME
    assert nb.getNumExceptions() == len(b.payload["numerical_data"]["pairs"])


def test_optional_dependency_error(pcff, monkeypatch):
    original = builtins.__import__

    def forbidden(name, *args, **kw):
        if name == "openmm":
            raise ImportError("blocked dependency")
        return original(name, *args, **kw)

    monkeypatch.setattr(builtins, "__import__", forbidden)
    with pytest.raises(ValidationError, match="optional dependency openmm"):
        build_openmm_periodic_system(pcff)


@pytest.mark.parametrize(
    "key",
    (
        "typed_graph",
        "graph_charges",
        "parameter_assignment",
        "config",
        "numerical_data",
    ),
)
def test_resigned_nested_mutation(pcff, key):
    p = pcff.payload
    if key == "typed_graph":
        p[key]["atom_types"]["1"] = "wrong"
    elif key == "graph_charges":
        p[key]["charges"]["1"] += 0.1
    elif key == "parameter_assignment":
        p[key]["source"]["sha256"] = "b" * 64
    elif key == "config":
        p[key]["pme_tolerance"] = 0.1
    else:
        p[key]["pcff_terms"].pop()
    if key != "numerical_data":
        p[key] = json.loads(
            _signed({k: v for k, v in p[key].items() if k != "identity"})
        )["payload"]
    raw = _signed({k: v for k, v in p.items() if k != "identity"})
    with pytest.raises(ValidationError, match="contradicts"):
        PeriodicParameterizedSystem(raw, pcff._prepared).validate_integrity()


def test_stale_public_records(pcff):
    p = pcff.typed_graph.payload
    p["atom_types"]["1"] = "wrong"
    from island.graph import UnifiedTypedGraph

    t = UnifiedTypedGraph(_signed({k: v for k, v in p.items() if k != "identity"}))
    with pytest.raises(ValidationError, match="typed_graph"):
        prepare_periodic_forcefield(
            pcff.system, pcff.final_graph, pcff._prepared, pcff.config, typed_graph=t
        )


def test_graph_mutation(pcff):
    system = pcff.system
    system.box.lx += 1
    with pytest.raises(ValidationError):
        prepare_periodic_forcefield(
            system, pcff.final_graph, pcff._prepared, pcff.config
        )
    p = pcff.payload
    p["system"]["sites"][0]["formal_charge"] = 1
    with pytest.raises(ValidationError):
        PeriodicParameterizedSystem(
            _signed({k: v for k, v in p.items() if k != "identity"}), pcff._prepared
        ).validate_integrity()


def test_relocated_native_bundle(native, tmp_path):
    b, sources = native
    from island.forcefields import AmberBundleArtifacts

    family = b._prepared._family
    artifacts = (
        AmberBundleArtifacts(RETAINED / family / "relocated/result.prmtop")
        if family in ("gaff", "gaff2")
        else None
    )
    save_periodic_backend_bundle(
        b, tmp_path / "original", sources=sources, artifacts=artifacts
    )
    shutil.move(tmp_path / "original", tmp_path / "moved")
    loaded = load_periodic_backend_bundle(
        tmp_path / "moved", sources=sources, expected_identity=b.identity
    )
    assert loaded.identity == b.identity
    assert loaded.final_graph.payload == b.final_graph.payload
    with pytest.raises(ValidationError, match="Unexpected"):
        load_periodic_backend_bundle(
            tmp_path / "moved", sources=sources, expected_identity="b" * 64
        )
    code = """
import builtins,sys
old=builtins.__import__
def guarded(name, globals=None, locals=None, fromlist=(), level=0):
    if level == 0 and name.split('.')[0] in {'rdkit','openmm','scipy','foyer'}:raise ImportError('blocked scientific dependency '+name)
    return old(name,globals,locals,fromlist,level)
builtins.__import__=guarded
from island.periodic import load_periodic_backend_bundle
from island.forcefields import PreparedForceFieldSources
sources=PreparedForceFieldSources(**__import__('json').loads(sys.argv[2]))
b=load_periodic_backend_bundle(sys.argv[1],sources=sources,expected_identity=sys.argv[3])
assert not any(k.split('.')[0] in {'rdkit','openmm','scipy','foyer'} for k in sys.modules)
print(b.identity)
"""
    source_args = (
        {}
        if sources is None
        else {k: str(v) for k, v in sources.__dict__.items() if v is not None}
    )
    env = {**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src")}
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            code,
            str(tmp_path / "moved"),
            json.dumps(source_args),
            b.identity,
        ],
        capture_output=True,
        text=True,
        env=env,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == b.identity
    with (tmp_path / "moved/system.data").open("a") as stream:
        stream.write("tampered")
    with pytest.raises(ValidationError, match="Tampered periodic export"):
        load_periodic_backend_bundle(tmp_path / "moved", sources=sources)


@pytest.mark.parametrize("mutation", ("scale", "source", "family"))
def test_configuration_native_policy_binding(pcff, mutation):
    p = pcff.config.payload
    if mutation == "scale":
        p["one_four_scaling"] = [0.5, 0.5]
    elif mutation == "source":
        p["source"]["sha256"] = "b" * 64
    else:
        p["family"] = "GAFF"
        p["source"]["family"] = "GAFF"
        p["lj_mixing_rule"] = "lorentz_berthelot_12_6"
    config = PeriodicForceFieldConfig(
        **{k: v for k, v in p.items() if k not in {"schema", "identity"}}
    )
    with pytest.raises(ValidationError, match="scaling|source|family"):
        prepare_periodic_forcefield(
            pcff.system, pcff.final_graph, pcff._prepared, config
        )


def test_export_preserves_existing_file(pcff, tmp_path):
    target = tmp_path / "system.data"
    target.write_text("prior valid file")
    with pytest.raises(ValidationError, match="exists"):
        export_lammps_data(pcff, target)
    assert target.read_text() == "prior valid file"


def test_missing_cross_term_rejects_before_lammps(pcff):
    from island._periodic_lammps import _pcff

    data = pcff.payload["numerical_data"]
    data["pcff_terms"] = [
        row for row in data["pcff_terms"] if row["family"] != "bond-bond_1_3"
    ]
    with pytest.raises(ValidationError, match="Missing PCFF cross term"):
        _pcff(data)


@pytest.mark.parametrize("origin", ("native_increments", "provided"))
def test_j17_external_origins(origin):
    from periodic_j22_fixtures import FRC

    from island.core import SimulationBox
    from island.forcefields import adopt_forcefield
    from island.forcefields.pcff import (
        assign_pcff_parameters,
        assign_typed_pcff_charges,
        bind_pcff_types,
        define_pcff_model,
        load_pcff_source,
        provide_pcff_charges,
        special_pair_policy,
    )
    from island.forcefields.pcff.fallbacks import COMPATIBILITY_POLICY
    from island.graph import build_psmiles_graph, final_graph

    s, _g, _ = build_psmiles_graph("[*:1]CC[*:2]", dp=3, random_seed=2026)
    source = load_pcff_source(FRC)
    typing = bind_pcff_types(
        s,
        source,
        {i: a.element.lower() for i, a in s.topology.sites.items()},
        provenance="J22 regression of J17 explicit c/h binding",
        evidence_references=["docs/phase_4j17.md"],
    )
    if origin == "native_increments":
        charges = assign_typed_pcff_charges(
            s, typing, resolution_policy=COMPATIBILITY_POLICY
        )
    else:
        charges = provide_pcff_charges(
            s,
            typing,
            dict.fromkeys(s.topology.sites, 0.0),
            unit="elementary_charge",
            provenance="Software charge-origin control, not scientific acceptance",
            evidence_references=["software fixture"],
            component_totals={min(s.topology.sites): 0.0},
            total_charge=0.0,
            resolution_policy=COMPATIBILITY_POLICY,
        )
    assignment = assign_pcff_parameters(
        s, typing, charges, resolution_policy=COMPATIBILITY_POLICY
    )
    model = define_pcff_model(
        assignment, special_pairs=special_pair_policy(lj=(0, 0, 1), coulomb=(0, 0, 1))
    )
    prepared = adopt_forcefield(s, "pcff", model)
    s.box = SimulationBox(30, 30, 30)
    config = PeriodicForceFieldConfig(**{**config_kwargs(), "source": source.identity})
    result = prepare_periodic_forcefield(s, final_graph(s), prepared, config)
    assert result.graph_charges.payload["charge_origin"] == origin
    assert result.graph_charges.charges == charges.charges
    assert result.typed_graph.payload["automatic_perception"] == "not_performed"
