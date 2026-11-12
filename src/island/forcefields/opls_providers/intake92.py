"""J92 source-bound native-real OPLS intake. No engine registration or exports."""

from decimal import Decimal
from pathlib import Path

from island.exceptions import IslandError
from island.forcefields.assignment_core.schema import identity, is_hash, require
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.bundle_adapters.schema import keys, loads, text

from .intake91 import outside_git
from .lammps92 import expand_pairs, parse_data, parse_settings, row_metadata
from .polymer import finite_graph
from .registry import _match_key, _result, external_file, fence, obligations, pinned
from .schema import FLAGS, CanonicalOPLS, number

SOURCE_DIR = Path(
    "/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-opls-bundles/pe5"
)
BACKENDS = ("lammps", "openmm")
FILES = (
    "Lig.lammps.lmp",
    "Lig.lammps.in",
    "request.json",
    "provider_metadata.json",
    "independent_audit.json",
)
PARTIAL_FILES = ("LIG.lammps.lmp", "LIG.lammps.in")
READY = {**FLAGS, "automatic_assignment": False}
AUDIT_FIELDS = (
    "schema",
    "status",
    "independent",
    "reviewer",
    "authority",
    "source_binding",
    "graph_identity",
    "settings_identity",
    "expansion_hash",
    "claim_hashes",
    "improper_policy_identity",
)


def source_binding(manifest):
    # Avoid the impossible manifest/audit self-hash cycle. The manifest is pinned
    # separately; the audit binds every OTHER source, request and provider file.
    return identity(
        sorted(
            (f["path"], f["sha256"], f["size"])
            for f in manifest["files"]
            if f["path"] != "independent_audit.json"
        )
    )


def _bundle_manifest(m):
    keys(m, ("schema", "bundle_id", "family", "files"), "manifest_schema_required")
    require(
        m["schema"] == "island.j92.lammps-bundle.v1" and m["family"] == "OPLS-AA",
        "OPLS_manifest_required",
    )
    text(m["bundle_id"], "bundle_id_required")
    require(isinstance(m["files"], list), "complete_file_manifest_required")
    for f in m["files"]:
        keys(f, ("path", "sha256", "size"), "file_manifest_schema")
        require(
            is_hash(f["sha256"]) and type(f["size"]) is int and f["size"] > 0,
            "file_hash_and_size_required",
        )
        require(f["path"] in FILES, "unlisted_source_role")
    require(
        len(m["files"]) == len(FILES) and {f["path"] for f in m["files"]} == set(FILES),
        "complete_file_manifest_required",
    )


def _request(q):
    keys(
        q,
        (
            "schema",
            "target",
            "family",
            "settings",
            "backends",
            "atom_count",
            "graph_identity",
            "atom_mapping",
        ),
        "complete_graph_request_required",
    )
    require(
        q["schema"] == "island.j92.request.v1"
        and q["target"] == "finite_PE_five_mer"
        and q["family"] == "OPLS-AA",
        "exact_PE5_request_required",
    )
    require(
        isinstance(q["backends"], list)
        and q["backends"]
        and all(b in BACKENDS for b in q["backends"])
        and len(q["backends"]) == len(set(q["backends"])),
        "only_OpenMM_LAMMPS_backends",
    )
    require(
        isinstance(q["settings"], dict)
        and q["settings"].get("intended_export_engine") in q["backends"],
        "backend_scope_mismatch",
    )
    require(
        type(q["settings"].get("repeat_count")) is int
        and q["settings"]["repeat_count"] == 5,
        "exact_repeat_count_five_required",
    )
    graph = finite_graph(q["settings"])
    from rdkit import Chem

    repeat = Chem.MolFromSmiles(q["settings"]["psmiles"])
    require(
        repeat is not None and Chem.MolToSmiles(repeat) == "*CC*",
        "exact_PE_repeat_unit_required",
    )
    require(q["graph_identity"] == graph.identity, "exact_graph_identity_mismatch")
    require(
        type(q["atom_count"]) is int and q["atom_count"] == len(graph.atom_order),
        "requested_atom_count_mismatch",
    )
    mol = molecule(graph)
    require(
        all(a.GetSymbol() in ("C", "H") for a in mol.GetAtoms())
        and all(b.GetBondTypeAsDouble() == 1 for b in mol.GetBonds()),
        "PE_single_bond_CH_scope_only",
    )
    require(
        isinstance(q["atom_mapping"], list)
        and len(q["atom_mapping"]) == len(graph.atom_order),
        "explicit_atom_mapping_required",
    )
    require(
        all(
            isinstance(p, list)
            and len(p) == 2
            and all(type(v) is int and v > 0 for v in p)
            for p in q["atom_mapping"]
        ),
        "invalid_atom_mapping",
    )
    mapping = dict(q["atom_mapping"])
    require(
        len(mapping) == len(graph.atom_order)
        and set(mapping.values()) == set(graph.atom_order),
        "atom_mapping_not_bijective",
    )
    return graph, mapping


def _provider(p):
    keys(
        p,
        (
            "schema",
            "provider",
            "version",
            "family",
            "source_authority",
            "license",
            "provenance",
            "atom_typing_origin",
            "charge_origin",
            "synthetic",
            "improper_policy",
        ),
        "complete_provider_metadata_required",
    )
    require(
        p["schema"] == "island.j92.provider.v1"
        and p["provider"] == "LigParGen"
        and p["family"] == "OPLS-AA",
        "LigParGen_OPLS_only",
    )
    for k in (
        "version",
        "source_authority",
        "provenance",
        "atom_typing_origin",
        "charge_origin",
    ):
        text(p[k], "provider_metadata_missing:" + k)
        require(
            p[k].lower() not in ("unknown", "missing", "unresolved"),
            "provider_metadata_unresolved:" + k,
        )
    keys(
        p["license"],
        ("identifier", "authorization"),
        "license_and_authorization_required",
    )
    text(p["license"]["identifier"], "license_identity_missing")
    require(
        p["license"]["authorization"] == "authorized_local_use", "authorization_missing"
    )
    require(type(p["synthetic"]) is bool, "synthetic_status_required")
    policy = p["improper_policy"]
    keys(policy, ("center_position", "omissions"), "explicit_improper_policy_required")
    require(
        type(policy["center_position"]) is int
        and policy["center_position"] in range(4),
        "improper_center_order_required",
    )
    require(
        isinstance(policy["omissions"], list), "explicit_improper_omissions_required"
    )
    for row in policy["omissions"]:
        keys(row, ("sites", "reason"), "improper_omission_schema")
        require(
            isinstance(row["sites"], list)
            and len(row["sites"]) == 4
            and all(type(i) is int and i > 0 for i in row["sites"]),
            "invalid_improper_omission_sites",
        )
        text(row["reason"], "improper_omission_evidence_required")


def _coverage(graph, mapping, parsed, provider, expansion):
    needed, _ = obligations(graph)
    mol = molecule(graph)
    rows = parsed["rows"]
    atoms = {r["id"]: r for r in rows["Atoms"]}
    require(set(atoms) == set(mapping), "source_atom_mapping_inventory_mismatch")
    require(
        len({a["tokens"][1] for a in atoms.values()}) == 1,
        "multiple_source_molecules_rejected",
    )
    masses = {r["id"]: number(r["tokens"][1]) for r in rows["Masses"]}
    for i, atom in atoms.items():
        expected = {"C": Decimal("12.011"), "H": Decimal("1.008")}[
            mol.GetAtomWithIdx(mapping[i] - 1).GetSymbol()
        ]
        require(
            masses[int(atom["tokens"][2])] == expected,
            "source_element_mass_mapping_mismatch",
        )
    charge = sum((number(a["tokens"][3]) for a in atoms.values()), Decimal(0))
    require(
        charge == sum(a.GetFormalCharge() for a in mol.GetAtoms()),
        "charge_conservation_mismatch",
    )
    present = {name: {} for name in needed}
    for i, atom in atoms.items():
        for name in ("atom_typing", "charge", "nonbonded"):
            present[name][_match_key(name, [mapping[i]])] = atom["source"]
    for family, section in (
        ("bond", "Bonds"),
        ("angle", "Angles"),
        ("torsion", "Dihedrals"),
    ):
        for row in rows.get(section, []):
            sites = [mapping[int(i)] for i in row["tokens"][2:]]
            key = _match_key(family, sites)
            require(
                key not in present[family], "duplicate_interaction_instance:" + family
            )
            require(
                key in {_match_key(family, s) for s in needed[family]},
                "source_interaction_not_in_requested_graph:" + family,
            )
            present[family][key] = row["source"]
    require(
        set(present["bond"]) == {_match_key("bond", s) for s in needed["bond"]},
        "source_exact_bond_graph_mismatch",
    )
    policy = provider["improper_policy"]

    def improper_key(sites):
        require(
            len(set(sites)) == 4 and set(sites) <= set(mapping),
            "invalid_improper_mapping",
        )
        sites = [mapping[i] for i in sites]
        c = sites[policy["center_position"]]
        peripheral = tuple(sorted(v for v in sites if v != c))
        return ("improper", c, peripheral)

    allowed = {_match_key("improper", s) for s in needed["improper"]}
    for row in rows.get("Impropers", []):
        key = improper_key([int(i) for i in row["tokens"][2:]])
        require(
            key in allowed and key not in present["improper"],
            "improper_order_or_duplicate_mismatch",
        )
        present["improper"][key] = row["source"]
    for i, row in enumerate(policy["omissions"]):
        key = improper_key(row["sites"])
        require(
            key in allowed and key not in present["improper"],
            "improper_omission_conflict",
        )
        present["improper"][key] = {
            "path": "provider_metadata.json",
            "pointer": f"/improper_policy/omissions/{i}",
            "status": "explicit_not_applicable",
            "reason": row["reason"],
        }
    for p in expansion["pair14"]:
        sites = [mapping[i] for i in p["sites"]]
        present["interaction_14"][_match_key("interaction_14", sorted(sites))] = {
            "expansion_hash": expansion["expansion_hash"]
        }
    present["topology"][_match_key("topology", [1])] = {
        "graph_identity": graph.identity
    }
    report, missing, instances = {}, [], []
    for family, sites in needed.items():
        for s in sites:
            key = _match_key(family, s)
            hit = key in present[family]
            instances.append(
                {
                    "family": family,
                    "sites": list(s),
                    "status": "evidenced" if hit else "missing",
                    "source": present[family].get(key),
                }
            )
            if not hit:
                missing.append(
                    {
                        "family": family,
                        "sites": list(s),
                        "reason": "no_evidenced_source_term",
                    }
                )
        report[family] = {
            "required": len(sites),
            "covered": sum(_match_key(family, s) in present[family] for s in sites),
        }
    # Independently compare expanded graph-distance pairs to inherited inventory.
    _, excluded = obligations(graph)
    mapped_exclusions = {
        tuple(sorted(mapping[i] for i in p["sites"])) for p in expansion["exclusions"]
    }
    require(mapped_exclusions == set(excluded), "exclusion_expansion_not_complete")
    require(
        report["interaction_14"]["required"] == report["interaction_14"]["covered"],
        "pair14_expansion_not_complete",
    )
    report["exclusion"] = {"required": len(excluded), "covered": len(mapped_exclusions)}
    instances.extend(
        {
            "family": "exclusion",
            "sites": list(p),
            "status": "evidenced",
            "source": {"expansion_hash": expansion["expansion_hash"]},
        }
        for p in sorted(mapped_exclusions)
    )
    return report, missing, instances, str(charge)


def _partial_two_file_intake(root, result, data_sha256, settings_sha256):
    """Inspect two pinned provider files without qualifying an unknown graph."""
    require(is_hash(data_sha256), "pinned_LAMMPS_data_sha256_required")
    require(is_hash(settings_sha256), "pinned_LAMMPS_settings_sha256_required")
    paths = {name: external_file(root, name) for name in PARTIAL_FILES}
    raw = {
        PARTIAL_FILES[0]: pinned(paths[PARTIAL_FILES[0]], data_sha256),
        PARTIAL_FILES[1]: pinned(paths[PARTIAL_FILES[1]], settings_sha256),
    }
    specs = {
        name: {"path": name, "sha256": sha, "size": len(raw[name])}
        for name, sha in zip(PARTIAL_FILES, (data_sha256, settings_sha256), strict=True)
    }
    parsed = parse_data(raw[PARTIAL_FILES[0]], specs[PARTIAL_FILES[0]])
    settings = parse_settings(
        raw[PARTIAL_FILES[1]], specs[PARTIAL_FILES[1]], PARTIAL_FILES[0]
    )
    rows = parsed["rows"]
    atoms = rows["Atoms"]
    masses = sorted({row["tokens"][1] for row in rows["Masses"]}, key=number)
    charges = [number(row["tokens"][3]) for row in atoms]
    result.update(
        decision="parameter_bundle_incomplete",
        provider=parsed["provider"],
        provider_version=parsed["provider_version"],
        source_files_inspected=2,
        source_hashes={
            PARTIAL_FILES[0]: data_sha256,
            PARTIAL_FILES[1]: settings_sha256,
        },
        source_files=[specs[name] for name in PARTIAL_FILES],
        source_rows_inspected=sum(len(value) for value in rows.values())
        + len(settings["commands"])
        + len(settings["inert_commands"]),
        conventions=settings,
        parameter_summary={
            "counts": parsed["counts"],
            "atom_count": len(atoms),
            "atom_type_count": len({int(row["tokens"][2]) for row in atoms}),
            "mass_values": masses,
            "charge_count": len(charges),
            "charge_sum": str(sum(charges, Decimal(0))),
            "charge_min": str(min(charges)),
            "charge_max": str(max(charges)),
            "sections": {name: len(value) for name, value in sorted(rows.items())},
        },
        missing_terms=[
            {
                "family": family,
                "sites": None,
                "reason": "exact_graph_obligations_unavailable",
            }
            for family in (
                "atom_typing",
                "bond",
                "angle",
                "torsion",
                "improper",
                "charge",
                "nonbonded",
                "interaction_14",
                "exclusion",
                "topology",
            )
        ],
        blockers=[
            "manifest_json_missing",
            "request_json_missing",
            "provider_metadata_json_missing",
            "input_structure_missing",
            "independent_audit_json_missing",
            "terminal_caps_missing",
            "input_PSMILES_missing",
            "exact_graph_identity_missing",
            "tacticity_missing",
            "special_bonds_expansion_blocked_without_exact_graph",
        ]
        + ["convention_blocked:" + item for item in settings["blockers"]],
    )
    return _result(result)


def intake(
    source_dir=SOURCE_DIR,
    *,
    manifest_sha256=None,
    data_sha256=None,
    settings_sha256=None,
):
    fence()
    r = {
        "schema": "island.j92.intake-result.v1",
        "decision": "source_bundle_required",
        "provider": None,
        "provider_version": None,
        "requested_graph": None,
        "requested_settings": None,
        "qualified_graph": None,
        "source_files_inspected": 0,
        "source_rows_inspected": 0,
        "source_hashes": {},
        "source_files": [],
        "manifest_sha256": manifest_sha256,
        "data_sha256": data_sha256,
        "settings_sha256": settings_sha256,
        "backends": list(BACKENDS),
        "synthetic": False,
        "coverage": {},
        "missing_terms": [],
        "interaction_instances": [],
        "blockers": [],
        "conventions": {},
        "pair14_count": 0,
        "exclusion_count": 0,
        "pair_expansion": None,
        "audit_status": "not_inspected",
        "scientific_graphs_qualified": 0,
        **READY,
    }
    try:
        root = outside_git(source_dir)
        if not root.is_dir():
            r["blockers"] = ["explicit_external_source_directory_missing"]
            return _result(r)
        if not is_hash(manifest_sha256) or not (root / "manifest.json").is_file():
            if all((root / name).is_file() for name in PARTIAL_FILES):
                return _partial_two_file_intake(root, r, data_sha256, settings_sha256)
            r["decision"] = "parameter_bundle_incomplete"
            r["blockers"] = ["pinned_source_manifest_required"]
            return _result(r)
        r["decision"] = "parameter_bundle_incomplete"
        raw_manifest = pinned(root / "manifest.json", manifest_sha256)
        manifest = loads(raw_manifest)
        _bundle_manifest(manifest)
        actual = set()
        for p in root.rglob("*"):
            require(not p.is_symlink(), "bundle_symlink_rejected")
            if not p.is_dir():
                require(p.is_file(), "regular_source_files_required")
                actual.add(p.relative_to(root).as_posix())
        require(actual == {*FILES, "manifest.json"}, "missing_or_unlisted_bundle_file")
        raw = {}
        r["source_hashes"] = {
            "manifest.json": manifest_sha256,
            **{f["path"]: f["sha256"] for f in manifest["files"]},
        }
        for f in manifest["files"]:
            r["source_files_inspected"] += 1
            raw[f["path"]] = pinned(external_file(root, f["path"]), f["sha256"])
            require(len(raw[f["path"]]) == f["size"], "source_size_mutation")
            r["source_files"].append(f)
        q = loads(raw["request.json"])
        provider = loads(raw["provider_metadata.json"])
        audit = loads(raw["independent_audit.json"])
        graph, mapping = _request(q)
        _provider(provider)
        keys(audit, AUDIT_FIELDS, "independent_audit_metadata_required")
        require(
            audit["schema"] == "island.j92.audit.v1"
            and audit["status"] in ("pending", "approved"),
            "audit_status_required",
        )
        r.update(
            provider=provider["provider"],
            provider_version=provider["version"],
            requested_graph=graph.identity,
            requested_settings=q["settings"],
            synthetic=provider["synthetic"],
            backends=q["backends"],
            source_binding=source_binding(manifest),
        )
        needed, excluded = obligations(graph)
        required = {**needed, "exclusion": excluded}
        r["coverage"] = {
            name: {"required": len(sites), "covered": 0, "status": "not_verified"}
            for name, sites in required.items()
        }
        r["missing_terms"] = [
            {
                "family": name,
                "sites": list(s),
                "reason": "source_or_convention_not_verified",
            }
            for name, sites in required.items()
            for s in sites
        ]
        specs = {f["path"]: f for f in manifest["files"]}
        parsed = parse_data(raw["Lig.lammps.lmp"], specs["Lig.lammps.lmp"])
        settings = parse_settings(
            raw["Lig.lammps.in"], specs["Lig.lammps.in"], "Lig.lammps.lmp"
        )
        require(
            parsed["provider_version"] == provider["version"],
            "provider_version_mismatch",
        )
        r["source_rows_inspected"] = (
            sum(len(rr) for rr in parsed["rows"].values())
            + len(settings["commands"])
            + len(settings["inert_commands"])
        )
        r["conventions"] = settings
        require(
            not settings["blockers"],
            "convention_blocked:" + ",".join(settings["blockers"]),
        )
        atoms = [row["id"] for row in parsed["rows"]["Atoms"]]
        bonds = [
            tuple(int(i) for i in row["tokens"][2:])
            for row in parsed["rows"].get("Bonds", [])
        ]
        expansion = expand_pairs(atoms, bonds, graph.identity, settings)
        coverage, missing, instances, charge = _coverage(
            graph, mapping, parsed, provider, expansion
        )
        for item in instances:
            if (
                item["source"]
                and item["source"].get("path") == "provider_metadata.json"
            ):
                item["source"]["sha256"] = specs["provider_metadata.json"]["sha256"]
        r.update(
            coverage=coverage,
            missing_terms=missing,
            interaction_instances=instances,
            pair_expansion=expansion,
            pair14_count=len(expansion["pair14"]),
            exclusion_count=len(expansion["exclusions"]),
            partial_charge_sum=charge,
            record_claims=row_metadata(parsed),
        )
        audited = False
        if audit["status"] == "approved":
            require(audit["independent"] is True, "independent_audit_required")
            for k in ("reviewer", "authority"):
                text(audit[k], "audit_authority_required")
                require(
                    audit[k].lower() not in ("unknown", "missing", "none"),
                    "audit_authority_required",
                )
            expected = {
                "source_binding": source_binding(manifest),
                "graph_identity": graph.identity,
                "settings_identity": identity(settings),
                "expansion_hash": expansion["expansion_hash"],
                "claim_hashes": sorted(c["claim_hash"] for c in r["record_claims"]),
                "improper_policy_identity": identity(provider["improper_policy"]),
            }
            for k, v in expected.items():
                require(audit[k] == v, "independent_audit_mismatch:" + k)
            audited = True
        else:
            require(
                audit["independent"] is False
                and all(
                    audit[k] is None
                    for k in AUDIT_FIELDS
                    if k not in ("schema", "status", "independent")
                ),
                "pending_audit_must_not_claim_approval",
            )
        r["audit_status"] = "independent_audit_verified" if audited else "pending"
        r["blockers"] = (
            ["missing_required_terms"]
            if missing
            else ([] if audited else ["independent_audit_pending"])
        )
        r["decision"] = (
            "parameter_bundle_incomplete"
            if missing
            else "audited_exact_graph"
            if audited
            else "provisional_parameterization"
        )
        if audited and not missing and not provider["synthetic"]:
            r["qualified_graph"] = graph.identity
            r["scientific_graphs_qualified"] = 1
        canonical = CanonicalOPLS.create(
            {
                "schema": "island.j92.native-real-opls.v1",
                "graph_identity": graph.identity,
                "atom_mapping": q["atom_mapping"],
                "source_hashes": r["source_hashes"],
                "data": parsed,
                "settings": settings,
                "pair_expansion": expansion,
                "improper_policy": provider["improper_policy"],
                **READY,
            }
        )
        r["canonical_identity"] = canonical.identity
        return _result(r, canonical)
    except (
        ValueError,
        TypeError,
        KeyError,
        OSError,
        UnicodeError,
        IslandError,
    ) as error:
        r["decision"] = "parameter_bundle_incomplete"
        r["blockers"] = [
            str(error)
            if type(error) is ValueError
            else type(error).__name__ + "_invalid_input"
        ]
        return _result(r)
