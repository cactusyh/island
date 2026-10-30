"""J54 external Class-II PDMS candidate audit; no executable profile is registered.

This additive facade leaves historical PCFF/OPLS and J55 identities unchanged.
Only a single pinned PCFF distribution is queried. Paper and COMPASS rows are
separate comparisons, never fallback parameters. Numerical equivalence is not a
parameter-completeness prerequisite.
"""

import json
import os
from collections import Counter
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

from island.charge_references.records import pack
from island.forcefield_policy.claims import ClaimRequirements, assess
from island.forcefields.pcff.mixed_sip import (
    DISTRIBUTION_PINS,
    load_pcff_sip_distribution,
)
from island.forcefields.pdms_source_completion import assess_pdms_source_completion
from island.forcefields.selection.inventory import canonical, digest
from island.workflows.bundle import system_data

SCHEMA = "island.external-pcff-pdms-class2-assessment.v1"
CANDIDATE = "island_external_pcff_pdms_class2_v1"
REVISION = "e891a3e10973c1a729e391a0aefaa02fd70f8c0f"
REGISTERED_PROFILES = ()
_ASSETS = ContextVar("island_pcff_pdms_sources", default=None)


class Class2AuditError(ValueError):
    def __init__(self, status, reason, **details):
        self.status = status
        self.diagnostics = [{"reason": reason, **details}]
        super().__init__(status + ": " + canonical(self.diagnostics))


def require(ok, status, reason, **details):
    if not ok:
        raise Class2AuditError(status, reason, **details)


def code_identity():
    return digest(
        {
            p.name: sha256(p.read_bytes()).hexdigest()
            for p in sorted(Path(__file__).parent.glob("*.py"))
        }
    )


@dataclass(frozen=True)
class PCFFPDMSAssets:
    """Administrator-installed original distribution; never a parameter overlay."""

    directory: object = None

    def load(self):
        require(
            self.directory is not None,
            "asset_unavailable",
            "pcff_distribution_not_supplied",
        )
        root = Path(self.directory)
        for name, expected in DISTRIBUTION_PINS.items():
            path = root / name
            require(
                path.is_file(),
                "asset_unavailable",
                "incomplete_source_archive",
                file=name,
            )
            require(
                sha256(path.read_bytes()).hexdigest() == expected,
                "source_incompatible",
                "source_hash_mutation",
                file=name,
            )
        return load_pcff_sip_distribution(root)


@contextmanager
def use_pcff_pdms_assets(assets):
    require(
        type(assets) is PCFFPDMSAssets,
        "asset_unavailable",
        "invalid_asset_configuration",
    )
    token = _ASSETS.set(assets)
    try:
        yield
    finally:
        _ASSETS.reset(token)


def configured_assets():
    return _ASSETS.get() or PCFFPDMSAssets(os.environ.get("ISLAND_PCFF_PDMS_SOURCES"))


def domain_diagnostics(system):
    """Narrow graph-domain check; no molecule names, coordinates or OPLS types."""
    system.validate()
    adj = {i: [] for i in system.topology.sites}
    bad = []
    if system.box is not None:
        bad.append(
            {
                "reason": "finite_graph_audit_only; periodic conventions are not established"
            }
        )
    for b in system.topology.bonds.values():
        adj[b.site1].append(b.site2)
        adj[b.site2].append(b.site1)
        if b.order != 1 or b.aromatic:
            bad.append(
                {"sites": [b.site1, b.site2], "reason": "single_bond_domain_only"}
            )
    elements = {i: a.element for i, a in system.topology.sites.items()}
    for i, a in system.topology.sites.items():
        neighbors = sorted(elements[j] for j in adj[i])
        good = (
            (
                a.element == "Si"
                and neighbors in (["C", "C", "C", "O"], ["C", "C", "O", "O"])
            )
            or (a.element == "O" and neighbors == ["Si", "Si"])
            or (a.element == "C" and neighbors == ["H", "H", "H", "Si"])
            or (a.element == "H" and neighbors == ["C"])
        )
        if (
            not good
            or a.formal_charge
            or a.metadata.get("radical_electrons", 0)
            or a.metadata.get("isotope", 0)
        ):
            bad.append(
                {
                    "site": i,
                    "element": a.element,
                    "neighbors": neighbors,
                    "reason": "unsupported_atom_environment_or_terminal_group",
                }
            )
    # Finite, linear, TMS-capped siloxane skeletons only. No cyclic/network claim.
    heavy = [i for i, e in elements.items() if e in ("Si", "O")]
    edges = sum(elements[j] in ("Si", "O") for i in heavy for j in adj[i]) // 2
    if (
        not heavy
        or edges != len(heavy) - 1
        or len(system.topology.connected_components()) != 1
    ):
        bad.append({"reason": "finite_single_linear_tms_capped_domain_required"})
    return sorted(bad, key=canonical)


def _source_rows(p, distribution):
    """Wrap exact source records with versioned conversions and internal citations."""
    found = {}
    refs = {
        s["namespace"]: s
        for s in distribution.source.inventory["sections"]
        if s["name"] == "reference"
    }

    def add(record, kind, fields, values, units, converted=None):
        ref = refs[record["reference"]]
        output_units = (
            units
            if converted is None
            else [
                "radian"
                if u == "degree"
                else u.replace("E", "kJ/mol")
                .replace("L", "angstrom")
                .replace("A", "radian")
                for u in units
            ]
        )
        if converted is not None and len(converted) == 2 * len(output_units):
            output_units *= 2
        value = {
            "id": record["id"],
            "family": kind,
            "functional_form": kind,
            "source_file": "pcff.frc",
            "source_sha256": DISTRIBUTION_PINS["pcff.frc"],
            "citation": {
                "repository": "https://github.com/lammps/lammps",
                "revision": REVISION,
                "path": "tools/msi2lmp/frc_files/pcff.frc",
                "reference": record["reference"],
                "raw_reference": ref["lines"],
            },
            "doi": None,
            "doi_status": "no DOI bound to this distribution row; not replaced with a paper DOI",
            "location": {"line": record["line"], "namespace": record["namespace"]},
            "raw_text": record["raw"],
            "version": record["version"],
            "labels": record.get("types", [record.get("data", {}).get("type")]),
            "fields": fields,
            "original_values": values,
            "original_units": units,
            "unit_symbols": {"E": "kcal/mol", "L": "angstrom", "A": "radian"},
            "converted_values": values if converted is None else converted,
            "converted_units": output_units,
            "symmetric_expansion": "duplicate explicit symmetric coefficient block only"
            if converted is not None and len(converted) != len(values)
            else None,
            "conversion_formula": "identity for raw cells; E[kJ/mol]=4.184*E[kcal/mol]; phase/equilibrium radians=degree*(pi/180); angular displacements are radians",
            "conversion_implementation": "island.forcefields.pcff.class2:v1",
            "conversion_code_sha256": sha256(
                Path(
                    __import__(
                        "island.forcefields.pcff.class2", fromlist=["x"]
                    ).__file__
                ).read_bytes()
            ).hexdigest(),
            "authority": "single_distribution_candidate; not executable; not Sun1997 or COMPASS overlay",
        }
        found[record["id"]] = {**value, "evidence_checksum": digest(value)}

    from island.forcefields.pcff.class2 import inspect_pcff_class2

    catalog = inspect_pcff_class2(distribution.source)["records"]
    for rid, row in p["source_rows"].items():
        raw = catalog[rid]
        add(
            raw,
            row["family"],
            row["fields"],
            row["original_values"],
            row["units"],
            [str(x) for x in raw["normalized_values"]],
        )
    for entry in p["typing"]["entries"].values():
        r = entry["source_atom_record"]["record"]
        w = r["raw"].split()
        add(
            r,
            "atom_typing_mass",
            ["mass", "element", "connections"],
            w[3:6],
            ["dalton", "element_symbol", "count"],
        )
    for decision in p["queries"].values():
        for proof in decision.get("equivalence_evidence", []):
            r = proof["record"]
            add(
                r,
                "equivalence",
                ["nonbond", "bond", "angle", "torsion", "out_of_plane"],
                r["raw"].split()[3:],
                ["atom_type"] * 5,
            )
    for contribution in p["native_charge_record"]["contributions"]:
        for proof in contribution["equivalence_records"]:
            r = proof["record"]
            add(
                r,
                "equivalence",
                ["nonbond", "bond", "angle", "torsion", "out_of_plane"],
                r["raw"].split()[3:],
                ["atom_type"] * 5,
            )
        for match in contribution["matches"]:
            r = match["selection"]["record"]
            r = {**r, "types": r["data"]["types"]}
            add(
                r, "bond_increments", ["dq_i", "dq_j"], r["raw"].split()[4:], ["e", "e"]
            )
    return {k: found[k] for k in sorted(found)}


def _conventions(distribution):
    """Raw equations are source facts; missing physical policies stay missing."""
    sections = distribution.source.inventory["sections"]
    evidence = {}
    for s in sections:
        if s["name"] not in ("reference", "version") and s["namespace"] == "cff91":
            lines = [r for r in s["lines"] if r["raw"].startswith((">", "@"))]
            if lines:
                evidence[s["name"]] = {
                    "source_file": "pcff.frc",
                    "source_sha256": DISTRIBUTION_PINS["pcff.frc"],
                    "lines": lines,
                }
    gaps = {
        "exclusions": "source-specific 1-2/1-3 exclusions absent in this FRC",
        "coulomb_14": "source-specific 1-4 Coulomb rule absent in this FRC",
        "lj_14": "source-specific 1-4 LJ rule absent in this FRC",
        "cutoff": "no declared finite/periodic cutoff or switching policy",
        "long_range": "no declared Coulomb solver, dielectric or periodic image policy",
        "torsion_geometry": "source Phi geometry/phase mapping to executable backend not established",
        "angle_angle_torsion_form": "linear Phi source comment; no authorized cos(Phi) replacement",
        "torsion_torsion_applicability": "source declares torsion-torsion_1 equation but supplies no rows or absence rule",
        "improper_applicability": "structural degree-three Wilson enumeration is a converter rule; domain-specific source omission authority not established",
        "source_authorization": "distribution identity is pinned; no new license/redistribution grant inferred for an external profile",
    }
    conventions = {
        name: {"value": None, "evidence": [], "required": True, "reason": reason}
        for name, reason in gaps.items()
    }
    # Preserve exact printed equations, including the unconverted torsion form.
    for name, value in {
        "chemical_domain": "audit-only linear explicit-H TMS-capped methyl-PDMS; disiloxane is a separate Si-H control",
        "charge_origin": "one pinned FRC bond increments; zero atomic base; no redistribution/normalization",
        "units": {
            "source_energy": "kcal/mol",
            "length": "angstrom",
            "equilibrium_angle": "degree",
            "angular_displacement": "radian",
            "charge": "e",
            "authority": "existing pinned FRC converter interpretation; not independent paper authority",
        },
        "functional_forms": evidence,
        "torsion_convention": evidence.get("torsion_3"),
        "lj_mixing": evidence.get("nonbond(9-6)"),
        "provenance": {"revision": REVISION, "overlay": False},
        "source_hashes": DISTRIBUTION_PINS,
        "evidence": evidence,
    }.items():
        conventions[name] = {
            "value": value,
            "evidence": [
                {
                    "source_sha256": DISTRIBUTION_PINS["pcff.frc"],
                    "authority": "pinned distribution and labelled converter interpretation",
                }
            ],
        }
    return conventions


def _assess(system, assets, *, audit_si_h=False):
    distribution = assets.load()
    diagnostics = domain_diagnostics(system)
    # Si-H is audited only as an explicitly separate source-covered diagnostic.
    if audit_si_h:
        require(
            Counter(a.element for a in system.topology.sites.values())
            == {"Si": 2, "O": 1, "H": 6}
            and all(
                d.get("neighbors") == ["H", "H", "H", "O"]
                if d.get("element") == "Si"
                else d.get("neighbors") == ["Si"]
                for d in diagnostics
            ),
            "typing_ambiguous",
            "silane_control_only",
        )
    elif diagnostics:
        raise Class2AuditError(
            "typing_ambiguous", "outside_declared_domain", environments=diagnostics
        )
    p = assess_pdms_source_completion(system, distribution).payload
    source_rows = _source_rows(p, distribution)
    required, assigned = [], []
    for i, entry in p["typing"]["entries"].items():
        rid = entry["source_atom_record"]["record"]["id"]
        for kind in ("atom_type", "mass", "charge"):
            key = f"{kind}:{i}"
            row_ids = [rid]
            if kind == "charge":
                row_ids = sorted(
                    {
                        m["selection"]["record"]["id"]
                        for c in p["native_charge_record"]["contributions"]
                        if i in c["sites"]
                        for m in c["matches"]
                    }
                )
            required.append(
                {
                    "id": key,
                    "family": kind,
                    "sites": [i],
                    "dependencies": [] if kind == "atom_type" else [f"atom_type:{i}"],
                }
            )
            assigned.append(
                {
                    "id": key,
                    "status": "assigned",
                    "row_ids": row_ids,
                    "evidence": [source_rows[r]["evidence_checksum"] for r in row_ids],
                    "value": p["exact_native_charge_vector"][i]
                    if kind == "charge"
                    else entry["type"]
                    if kind == "atom_type"
                    else entry["source_atom_record"]["record"]["data"]["mass"],
                    "unit": {
                        "charge": "e",
                        "atom_type": "source_label",
                        "mass": "dalton",
                    }[kind],
                    "charge_record_identity": digest(p["native_charge_record"])
                    if kind == "charge"
                    else None,
                }
            )
    for rid, row in sorted(p["interactions"].items()):
        required.append(
            {
                "id": rid,
                "family": row["family"],
                "sites": row["sites"],
                "dependencies": row["dependencies"],
            }
        )
        decision = p["queries"][row["query"]]
        if row["status"] == "assigned":
            assigned.append(
                {
                    "id": rid,
                    "status": "assigned",
                    "query": row["query"],
                    "row_ids": sorted({x["record_id"] for x in decision["selected"]}),
                    "values": decision["normalized_values"],
                    "units": "kJ/mol,angstrom,radian; source functional form retained",
                    "evidence": [
                        source_rows[x["record_id"]]["evidence_checksum"]
                        for x in decision["selected"]
                    ],
                }
            )
    conventions = _conventions(distribution)
    report = assess(
        required,
        assigned,
        conventions,
        typing_complete=True,
        charge_complete=p["native_charge_record"]["complete"],
        source_valid=True,
    )
    typed = p["typing"]["assignments"]
    missing = [
        {
            **row,
            "query": p["interactions"][row["id"]]["query"],
            "source_lookup": p["queries"][p["interactions"][row["id"]]["query"]],
        }
        for row in report["missing_terms"]
    ]
    report.update(
        schema=SCHEMA,
        family="pcff",
        source_kind="external_source_candidate",
        candidate_profile=CANDIDATE,
        selected_profile=None,
        registered_profiles=[],
        graph_identity=p["requirements"]["graph_identity"],
        atom_order=p["requirements"]["atom_order"],
        execution_input_identity=digest(system_data(system)),
        code_identity=code_identity(),
        source_hashes=DISTRIBUTION_PINS,
        typing=typed,
        typing_identity=p["typing_identity"],
        source_record_identity=p["identity"],
        source_rows=source_rows,
        source_to_interaction={r["id"]: r["row_ids"] for r in assigned},
        missing_terms=missing,
        conventions=conventions,
        coverage=p["coverage"],
        charge_origin="single_pinned_pcff_frc_bond_increments_candidate",
        charge_record=p["native_charge_record"],
        charge_vector=p["exact_native_charge_vector"],
        charge_identity=digest(p["native_charge_record"]),
        normalization=False,
        redistribution=False,
        improper_structural_inventory=p["inapplicable"],
        domain={
            "kind": "linear_tms_capped_explicit_H_methyl_pdms",
            "claim": "candidate audit only",
            "SiH_control": audit_si_h,
            "diagnostics": diagnostics,
            "excluded": [
                "Si-OH",
                "Si-O-C",
                "cyclic",
                "crosslinked",
                "phosphorus",
                "open_repeat_attachments",
            ],
            "C-C": "not used by these methyl-only graphs; no C-C coverage extrapolation",
        },
        lammps_export_status="blocked",
        lammps_runtime_status="blocked",
        openmm_build_status="blocked",
        cross_engine_equivalence="not_claimed",
        scientific_validation="not_established",
        production_readiness="not_established",
        production_validated=False,
        simulation_readiness="not_established",
        MD_authorized=False,
        backend_block_reason="missing source parameters/couplings and conventions; no numerical equivalence prerequisite",
        numerical_checks={
            k: "not_run"
            for k in (
                "static",
                "energy",
                "force",
                "finite_difference",
                "minimization",
                "NVT",
                "NPT",
            )
        },
    )
    return json.loads(canonical({**report, "identity": digest(report)}))


@dataclass(frozen=True)
class Class2Assessment:
    json_text: str
    assets: PCFFPDMSAssets
    audit_si_h: bool = False

    def __post_init__(self):
        p = self.payload
        require(
            p["schema"] == SCHEMA
            and p["family"] == "pcff"
            and p["candidate_profile"] == CANDIDATE
            and p["selected_profile"] is None
            and p["registered_profiles"] == []
            and p["parameter_completeness"] == "incomplete"
            and not p["parameter_complete"]
            and not p["MD_authorized"]
            and p["source_hashes"] == DISTRIBUTION_PINS
            and p["code_identity"] == code_identity()
            and p["identity"] == digest({k: v for k, v in p.items() if k != "identity"})
            and p["production_validated"] is False
            and p["scientific_validation"]
            == p["production_readiness"]
            == p["simulation_readiness"]
            == "not_established",
            "source_incompatible",
            "candidate_identity_or_authority_mutation",
        )

    @property
    def payload(self):
        return json.loads(self.json_text)

    def validate_integrity(self, system):
        expected = _assess(system, self.assets, audit_si_h=self.audit_si_h)
        require(
            self.payload == expected,
            "source_incompatible",
            "graph_charge_rows_units_or_authority_mutation",
        )
        return expected

    def require(self, claims=None):
        # Claims are independent; this record cannot authorize MD or register itself.
        self.assets.load()
        claims = (
            ClaimRequirements(parameter_complete=True) if claims is None else claims
        )
        claims.enforce(self.payload)
        return self


def prepare_polymer_forcefield(
    system, forcefield="pcff", *, assets=None, profile=None, audit_si_h=False
):
    """Family-only PCFF candidate inspection; source installation is administrative.

    Incomplete candidates return exact diagnostics and no selected profile. No
    native PCFF fallback or cross-family overlay is attempted for this domain.
    """
    require(
        forcefield == "pcff",
        "source_incompatible",
        "wrong_profile_family",
        family=forcefield,
    )
    require(
        profile in (None, CANDIDATE),
        "source_incompatible",
        "wrong_profile_family_or_identity",
        profile=profile,
    )
    assets = configured_assets() if assets is None else assets
    report = _assess(system, assets, audit_si_h=audit_si_h)
    return Class2Assessment(canonical(report), assets, audit_si_h)


def prepare_polymer(
    *, psmiles, forcefield="pcff", dp, head_end_group, tail_end_group, **kwargs
):
    from island.builders import build_linear_polymer

    system = build_linear_polymer(
        psmiles,
        dp=dp,
        head_end_group=head_end_group,
        tail_end_group=tail_end_group,
        add_hydrogens=True,
        random_seed=2026,
    )
    return system, prepare_polymer_forcefield(system, forcefield=forcefield, **kwargs)


def export_lammps(system, prepared, directory):
    prepared.validate_integrity(system)
    prepared.require()
    raise Class2AuditError(
        "asset_unavailable", "no_qualified_external_class2_export_adapter"
    )


def build_openmm(system, prepared):
    prepared.validate_integrity(system)
    prepared.require()
    raise Class2AuditError(
        "asset_unavailable", "no_qualified_external_class2_openmm_adapter"
    )


def save(system, prepared, directory):
    """Save a replayable data audit, never an executable MD bundle."""
    prepared.validate_integrity(system)
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=False)
    (root / "system.json").write_text(pack(system_data(system)))
    (root / "report.json").write_text(prepared.json_text)
    manifest = {
        "schema": SCHEMA,
        "code_identity": code_identity(),
        "source_hashes": DISTRIBUTION_PINS,
        "audit_si_h": prepared.audit_si_h,
        "report_identity": prepared.payload["identity"],
        "files": {
            n: sha256((root / n).read_bytes()).hexdigest()
            for n in ("system.json", "report.json")
        },
    }
    (root / "manifest.json").write_text(
        canonical({**manifest, "identity": digest(manifest)})
    )
    return sha256((root / "manifest.json").read_bytes()).hexdigest()


def reconstruct(directory, expected_sha256, *, assets=None):
    from island.charge_references.records import unpack
    from island.workflows.bundle import system_from

    root = Path(directory)
    require(
        sha256((root / "manifest.json").read_bytes()).hexdigest() == expected_sha256,
        "source_incompatible",
        "bundle_manifest_mutation",
    )
    manifest = json.loads((root / "manifest.json").read_text())
    require(
        manifest["identity"]
        == digest({k: v for k, v in manifest.items() if k != "identity"}),
        "source_incompatible",
        "manifest_identity",
    )
    require(
        manifest["code_identity"] == code_identity()
        and manifest["source_hashes"] == DISTRIBUTION_PINS,
        "source_incompatible",
        "source_or_code_identity",
    )
    require(
        set(manifest["files"]) == {"system.json", "report.json"},
        "source_incompatible",
        "bundle_membership",
    )
    for name, expected in manifest["files"].items():
        require(
            sha256((root / name).read_bytes()).hexdigest() == expected,
            "source_incompatible",
            "bundle_file_mutation",
            file=name,
        )
    system = system_from(unpack((root / "system.json").read_text()))
    result = prepare_polymer_forcefield(
        system, forcefield="pcff", assets=assets, audit_si_h=manifest["audit_si_h"]
    )
    require(
        result.json_text == (root / "report.json").read_text()
        and result.payload["identity"] == manifest["report_identity"],
        "source_incompatible",
        "reconstructed_parameters_changed",
    )
    return system, result
