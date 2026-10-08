"""Additive periodic backends from validated, family-owned native assignments.

No typing, charge generation, parameter fitting, packing or integration occurs
here. The explicit chemical view differs from the supplied system only by box.
"""

import json
import os
import shutil
import tempfile
from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
from math import isfinite
from pathlib import Path

from island.exceptions import ValidationError
from island.graph import FinalChemicalGraph
from island.graph.final import _digest, _system_data, _system_from_payload, pack, unpack

CONFIG_SCHEMA = "island_periodic_forcefield_config_v1"
BACKEND_SCHEMA = "island_periodic_parameterized_system_v1"
EXPORT_SCHEMA = "island_periodic_export_v1"
BUNDLE_SCHEMA = "island_periodic_backend_bundle_v1"
FAMILIES = {"PCFF", "OPLS-AA", "GAFF", "GAFF2"}
FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}


def _require(condition, message):
    if not condition:
        raise ValidationError(message)


def _finite(value, label, positive=False):
    _require(
        type(value) in (int, float) and isfinite(value),
        f"{label} must be finite and nonboolean",
    )
    if positive:
        _require(value > 0, f"{label} must be positive")


def _signed(payload):
    payload = json.loads(json.dumps(payload, allow_nan=False))
    return pack({**payload, "identity": _digest(payload)})


def _payload(raw, schema):
    p = unpack(raw)
    _require(p.get("schema") == schema, "Unsupported periodic schema")
    _require(
        p.get("identity") == _digest({k: v for k, v in p.items() if k != "identity"}),
        "Periodic record identity mismatch",
    )
    return p


@dataclass(frozen=True, init=False)
class PeriodicForceFieldConfig:
    """Owned JSON configuration. Distances are angstroms; public energies kJ/mol.

    ``one_four_scaling`` is (LJ, Coulomb). It must match the native policy.
    PME maps to LAMMPS PPPM and OpenMM PME; Ewald maps to Ewald in both.
    """

    json_text: str

    def __init__(
        self,
        *,
        family,
        source,
        box_lengths,
        nonbonded_cutoff,
        lj_mixing_rule,
        one_four_scaling,
        electrostatics_method,
        pme_tolerance,
        provenance,
        evidence=(),
        periodic=(True, True, True),
        neighbor_list_policy="bin_skin_2_angstrom_every_1_check",
        excluded_pair_policy="native_shortest_path",
        switching_policy="none",
        units="angstrom_kj_mol_elementary_charge",
        deterministic_export=True,
    ):
        p = {
            "schema": CONFIG_SCHEMA,
            "family": family,
            "source": deepcopy(source),
            "box_lengths": list(box_lengths),
            "periodic": list(periodic),
            "nonbonded_cutoff": nonbonded_cutoff,
            "lj_mixing_rule": lj_mixing_rule,
            "one_four_scaling": list(one_four_scaling),
            "electrostatics_method": electrostatics_method,
            "pme_tolerance": pme_tolerance,
            "provenance": provenance,
            "evidence": list(evidence),
            "neighbor_list_policy": neighbor_list_policy,
            "excluded_pair_policy": excluded_pair_policy,
            "switching_policy": switching_policy,
            "units": units,
            "deterministic_export": deterministic_export,
        }
        self._validate(p)
        object.__setattr__(self, "json_text", _signed(p))

    @staticmethod
    def _validate(p):
        from island.graph.unified import _source

        _require(p["family"] in FAMILIES, "Unsupported periodic force-field family")
        _source(p["source"], p["family"])
        _require(
            len(p["box_lengths"]) == 3, "Orthorhombic three-dimensional box required"
        )
        for value in p["box_lengths"]:
            _finite(value, "Box length", True)
        _require(
            p["periodic"] == [True, True, True]
            and all(type(v) is bool for v in p["periodic"]),
            "Three periodic boundaries required",
        )
        _finite(p["nonbonded_cutoff"], "Nonbonded cutoff", True)
        _require(
            2 * p["nonbonded_cutoff"] < min(p["box_lengths"]),
            "Cutoff must be less than half the shortest box length",
        )
        _finite(p["pme_tolerance"], "PME/Ewald tolerance", True)
        _require(p["pme_tolerance"] < 1, "PME/Ewald tolerance must be below one")
        _require(
            p["electrostatics_method"] in {"pme", "ewald"},
            "Explicit PME or Ewald electrostatics required",
        )
        expected = {
            "PCFF": "sixthpower_9_6",
            "OPLS-AA": "geometric_12_6",
            "GAFF": "lorentz_berthelot_12_6",
            "GAFF2": "lorentz_berthelot_12_6",
        }
        _require(
            p["lj_mixing_rule"] == expected[p["family"]],
            "LJ mixing rule contradicts family conventions",
        )
        _require(
            len(p["one_four_scaling"]) == 2,
            "Explicit LJ and Coulomb 1-4 scales required",
        )
        for value in p["one_four_scaling"]:
            _finite(value, "1-4 scale")
            _require(0 <= value <= 1, "Invalid 1-4 scale")
        for key, value in {
            "neighbor_list_policy": "bin_skin_2_angstrom_every_1_check",
            "excluded_pair_policy": "native_shortest_path",
            "switching_policy": "none",
            "units": "angstrom_kj_mol_elementary_charge",
        }.items():
            _require(p[key] == value, f"Unsupported periodic {key}")
        _require(p["deterministic_export"] is True, "Deterministic export required")
        _require(
            type(p["provenance"]) is str and p["provenance"].strip(),
            "Periodic provenance required",
        )
        _require(
            type(p["evidence"]) is list
            and all(type(x) is str and x.strip() for x in p["evidence"]),
            "Malformed periodic evidence",
        )

    @property
    def payload(self):
        p = _payload(self.json_text, CONFIG_SCHEMA)
        self._validate(p)
        return p

    @property
    def identity(self):
        return self.payload["identity"]

    @classmethod
    def from_json(cls, raw):
        p = _payload(raw, CONFIG_SCHEMA)
        p.pop("identity")
        p.pop("schema")
        result = cls(**p)
        _require(
            unpack(result.json_text) == unpack(raw),
            "Unexpected periodic configuration fields",
        )
        return result


def periodic_chemical_view(system):
    """Explicit box-free view for validating existing native chemical records.

    Atom IDs, metadata, topology, charges and coordinates are copied unchanged.
    This does not generate an assignment or alter the authoritative final graph.
    """
    result = system.copy()
    result.box = None
    return result


@dataclass(frozen=True)
class PeriodicParameterizedSystem:
    """Owned records plus a validated native witness, never a mutable OpenMM object."""

    json_text: str
    _prepared: object

    def __post_init__(self):
        object.__setattr__(self, "_prepared", deepcopy(self._prepared))

    @property
    def payload(self):
        return _payload(self.json_text, BACKEND_SCHEMA)

    @property
    def identity(self):
        self.validate_integrity()
        return self.payload["identity"]

    @property
    def system(self):
        payload = self.payload["system"]
        result = _system_from_payload(payload)
        if payload.get("impropers"):
            from island.core import Improper

            result.topology.impropers = [
                Improper(**row) for row in payload["impropers"]
            ]
        result.topology.rebuild_derived_interactions()
        return result

    @property
    def final_graph(self):
        return FinalChemicalGraph(pack(self.payload["final_graph"]))

    @property
    def config(self):
        return PeriodicForceFieldConfig.from_json(pack(self.payload["config"]))

    @property
    def typed_graph(self):
        from island.graph import UnifiedTypedGraph

        return UnifiedTypedGraph(pack(self.payload["typed_graph"]))

    @property
    def graph_charges(self):
        from island.graph import UnifiedGraphCharges

        return UnifiedGraphCharges(pack(self.payload["graph_charges"]))

    @property
    def parameter_assignment(self):
        from island.graph import UnifiedParameterAssignment

        return UnifiedParameterAssignment(pack(self.payload["parameter_assignment"]))

    def validate_integrity(self):
        p = self.payload
        system, graph, config = self.system, self.final_graph, self.config
        graph.validate_integrity(system)
        expected = _definition(system, graph, self._prepared, config)
        _require(
            p == unpack(_signed(expected)),
            "Periodic backend contradicts native assignment or nested identities",
        )


def _definition(system, graph, prepared, config):
    from island._periodic_native import binding
    from island.forcefields import PreparedForceField

    _require(type(graph) is FinalChemicalGraph, "FinalChemicalGraph required")
    _require(
        type(config) is PeriodicForceFieldConfig, "PeriodicForceFieldConfig required"
    )
    _require(
        type(prepared) is PreparedForceField,
        "Validated native PreparedForceField required",
    )
    graph.validate_integrity(system)
    p = config.payload
    _require(
        graph.periodic and system.box is not None,
        "Periodic final graph with a box required",
    )
    _require(
        list(system.box.lengths) == p["box_lengths"]
        and list(system.box.periodic) == p["periodic"],
        "Periodic graph/configuration box mismatch",
    )
    prepared.validate_integrity(periodic_chemical_view(system))
    typed, charges, assignment, data = binding(system, graph, prepared, p)
    return dict(
        schema=BACKEND_SCHEMA,
        system=_system_data(system),
        final_graph=graph.payload,
        config=p,
        typed_graph=typed.payload,
        graph_charges=charges.payload,
        parameter_assignment=assignment.payload,
        numerical_data=data,
        native_preparation_identity=prepared.identity,
        final_graph_identity=graph.identity,
        typed_graph_identity=typed.identity,
        charge_identity=charges.identity,
        assignment_identity=assignment.identity,
        config_identity=config.identity,
        **FLAGS,
    )


def prepare_periodic_forcefield(
    system,
    final_graph,
    prepared,
    config,
    *,
    typed_graph=None,
    graph_charges=None,
    parameter_assignment=None,
):
    """Bind an existing native assignment to this exact periodic graph.

    Supplied unified records are optional identity assertions and must equal the
    native-derived binding. Native validation cannot be bypassed by re-signing JSON.
    """
    p = _definition(system, final_graph, prepared, config)
    for record, key in (
        (typed_graph, "typed_graph"),
        (graph_charges, "graph_charges"),
        (parameter_assignment, "parameter_assignment"),
    ):
        if record is not None:
            _require(record.payload == p[key], f"Stale or contradictory periodic {key}")
    return PeriodicParameterizedSystem(_signed(p), prepared)


def build_openmm_periodic_system(backend):
    """Construct a fresh OpenMM System; optional dependency is loaded only here."""
    from island._periodic_openmm import build

    backend.validate_integrity()
    try:
        import openmm as mm
    except ImportError as exc:
        raise ValidationError(
            "Periodic OpenMM construction requires optional dependency openmm"
        ) from exc
    return build(backend.payload["numerical_data"], backend.config.payload, mm)


def _export_texts(backend):
    from island._periodic_lammps import render

    backend.validate_integrity()
    return render(backend)


def _publish_file(path, raw):
    target = Path(path)
    _require(
        not target.exists() and not target.is_symlink(),
        f"Destination already exists: {target}",
    )
    fd, name = tempfile.mkstemp(prefix=".island-periodic-", dir=target.parent)
    try:
        with os.fdopen(fd, "wb") as stream:
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(name, target)  # Exclusive, atomic publication; never overwrite.
    finally:
        os.unlink(name)


def export_lammps_data(backend, path=None):
    text = _export_texts(backend)[0]
    if path is not None:
        _publish_file(path, text.encode())
    return text


def export_lammps_input(backend, path=None):
    """Input reads ``system.data`` in its working directory and performs run 0."""
    text = _export_texts(backend)[1]
    if path is not None:
        _publish_file(path, text.encode())
    return text


def periodic_export_record(backend):
    data, script = _export_texts(backend)
    return unpack(
        _signed(
            {
                "schema": EXPORT_SCHEMA,
                "backend_identity": backend.identity,
                "final_graph_identity": backend.final_graph.identity,
                "typed_graph_identity": backend.typed_graph.identity,
                "charge_identity": backend.graph_charges.identity,
                "parameter_assignment_identity": backend.parameter_assignment.identity,
                "source": backend.config.payload["source"],
                "configuration_identity": backend.config.identity,
                "files": {
                    "system.data": sha256(data.encode()).hexdigest(),
                    "in.periodic": sha256(script.encode()).hexdigest(),
                },
            }
        )
    )


def save_periodic_backend_bundle(backend, path, *, sources=None, artifacts=None):
    """Stage, reconstruct and exclusively publish a relocatable native bundle."""
    from island.forcefields import save_prepared_forcefield

    backend.validate_integrity()
    target = Path(path)
    _require(
        not target.exists() and not target.is_symlink(),
        "Periodic bundle destination exists",
    )
    stage = Path(tempfile.mkdtemp(prefix=".island-periodic-", dir=target.parent))
    try:
        save_prepared_forcefield(
            backend._prepared._system,
            backend._prepared,
            stage / "native",
            sources=sources,
            artifacts=artifacts,
        )
        record = periodic_export_record(backend)
        _publish_file(stage / "backend.json", backend.json_text.encode())
        _publish_file(stage / "system.data", export_lammps_data(backend).encode())
        _publish_file(stage / "in.periodic", export_lammps_input(backend).encode())
        _publish_file(
            stage / "manifest.json",
            _signed(
                dict(
                    schema=BUNDLE_SCHEMA,
                    backend_identity=backend.identity,
                    export=record,
                    openmm_reconstruction="island_periodic_openmm_v1",
                    **FLAGS,
                )
            ).encode(),
        )
        load_periodic_backend_bundle(
            stage, sources=sources, expected_identity=backend.identity
        )
        target.mkdir()  # Reserve the name exclusively; no preexisting directory replaced.
        try:
            os.replace(stage, target)
        except BaseException:
            target.rmdir()
            raise
    finally:
        if stage.exists():
            shutil.rmtree(stage)


def load_periodic_backend_bundle(path, *, sources=None, expected_identity=None):
    """Revalidate retained native source rows; never regenerate chemistry or packing."""
    from island.forcefields import load_prepared_forcefield

    root = Path(path)
    manifest = _payload((root / "manifest.json").read_text(), BUNDLE_SCHEMA)
    native = load_prepared_forcefield(root / "native", sources=sources)
    result = PeriodicParameterizedSystem(
        (root / "backend.json").read_text(), native.prepared
    )
    result.validate_integrity()
    _require(
        manifest["backend_identity"] == result.identity,
        "Stale periodic bundle identity",
    )
    if expected_identity is not None:
        _require(
            result.identity == expected_identity, "Unexpected periodic backend identity"
        )
    export = periodic_export_record(result)
    _require(manifest["export"] == export, "Tampered periodic export record")
    for name, digest in export["files"].items():
        _require(
            sha256((root / name).read_bytes()).hexdigest() == digest,
            f"Tampered periodic export: {name}",
        )
    _require(
        manifest
        == unpack(
            _signed(
                dict(
                    schema=BUNDLE_SCHEMA,
                    backend_identity=result.identity,
                    export=export,
                    openmm_reconstruction="island_periodic_openmm_v1",
                    **FLAGS,
                )
            )
        ),
        "Tampered periodic manifest",
    )
    return result
