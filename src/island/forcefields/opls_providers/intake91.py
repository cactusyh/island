"""J91 explicit external intake. Reuse J90 parsers; never discover providers."""

from pathlib import Path

from island.exceptions import IslandError
from island.forcefields.assignment_core.schema import data, identity, is_hash, require
from island.forcefields.assignment_core.topology import molecule
from island.forcefields.bundle_adapters.schema import keys, loads, text

from .export import render
from .parsers import claim_metadata, pointer
from .polymer import finite_graph
from .registry import (
    OPLSProviderRegistry,
    _audit,
    _coverage,
    _manifest,
    _match_key,
    _result,
    _validate,
    external_file,
    fence,
    obligations,
    pinned,
)
from .schema import FLAGS, FORMS, CanonicalOPLS, row_key

CONTROLS = ("request", "manifest", "audit")
DEFINITIONS = {
    "sigma": "zero_of_12_6_pair_potential",
    "epsilon": "positive_well_depth_minimum_minus_epsilon",
    "minimum": "r_min=2^(1/6)*sigma",
}
CONVERSION = {
    "version": "island.j91.lammps-real.v1",
    "source_units": {"energy": "kJ/mol", "length": "nm"},
    "target_units": {"energy": "kcal/mol", "length": "Angstrom"},
    "formula": {
        "energy": "E_target=E_source/4.184",
        "length": "r_target=10*r_source",
        "bond_k": "K_target=k_source/(2*4.184*100)",
        "angle_improper_k": "K_target=k_source/(2*4.184)",
        "periodic_torsion": "K_target=k_source/4.184; n,phase,order unchanged",
    },
}
CONVERSION = {**CONVERSION, "conversion_hash": identity(CONVERSION)}
REQUEST_FIELDS = (
    "schema",
    "provider",
    "settings",
    "graph_identity",
    "source_authority",
    "self_contained_finite_molecule",
    "expected_conventions",
    "expected_forms",
    "nonbonded_definitions",
)
AUDIT_FIELDS = (
    "schema",
    "request_sha256",
    "source_audit",
    "pair_claim_hashes",
    "exclusion_claim_hashes",
    "export_conversions",
)


def outside_git(path):
    path = Path(path).absolute()
    require(
        not any(p.is_symlink() for p in (path, *path.parents)),
        "source_symlink_or_path_substitution",
    )
    require(
        not any((p / ".git").exists() for p in (path, *path.parents)),
        "licensed_files_must_remain_outside_Git_worktree",
    )
    return path


def input_identity(source_dir=None, **pins):
    return {
        "directory_supplied": source_dir is not None,
        "pins": {n: pins.get(n + "_sha256") for n in CONTROLS},
    }


def pair_claims(raw, file, pairs, exclusions):
    """Exact JSON pointers or physical lines, including explicit exclusions."""
    if file["format"] == "canonical_json":
        return tuple(
            [
                {"sites": list(p), "source": pointer(file, f"/{name}/{i}")}
                for i, p in enumerate(values)
            ]
            for name, values in (("pairs14", pairs), ("exclusions", exclusions))
        )
    require(file["format"] == "gromacs_itp", "self_contained_source_required")
    result = {"pairs": [], "exclusions": []}
    section = None
    for line, original in enumerate(raw.decode().splitlines(), 1):
        content = original.split(";", 1)[0].strip()
        if content.startswith("["):
            section = content.strip("[] ")
        elif content and section in result:
            fields = content.split()
            sites = [int(n) for n in (fields[:2] if section == "pairs" else fields)]
            for b in sites[1:]:
                result[section].append(
                    {
                        "sites": [sites[0], b],
                        "source": pointer(file, f"line:{line};section:{section}"),
                    }
                )
    return result["pairs"], result["exclusions"]


def _audit_metadata(a, request_pin, manifest_pin, graph, convention):
    keys(a, AUDIT_FIELDS, "complete_independent_audit_required")
    require(
        a["schema"] == "island.j91.audit.v1" and a["request_sha256"] == request_pin,
        "audit_request_pin_mismatch",
    )
    s = a["source_audit"]
    keys(
        s,
        (
            "schema",
            "reviewer",
            "authority_basis",
            "independent",
            "manifest_sha256",
            "graph_identity",
            "convention_identity",
            "claim_hashes",
            "atom_claim_hashes",
        ),
        "complete_independent_source_audit_required",
    )
    require(
        s["schema"] == "island.opls-provider-audit.v1" and s["independent"] is True,
        "independent_audit_required",
    )
    for k in ("reviewer", "authority_basis"):
        text(s[k], "audit_authority_required")
        require(
            s[k].lower() not in ("unknown", "none", "missing"),
            "audit_authority_required",
        )
    require(
        s["manifest_sha256"] == manifest_pin
        and s["graph_identity"] == graph
        and s["convention_identity"] == convention,
        "audit_scope_mismatch",
    )
    for obj, names in (
        (s, ("claim_hashes", "atom_claim_hashes")),
        (a, ("pair_claim_hashes", "exclusion_claim_hashes")),
    ):
        for name in names:
            require(
                isinstance(obj[name], list) and all(is_hash(h) for h in obj[name]),
                "audit_claim_hashes_required",
            )
    require(
        isinstance(a["export_conversions"], list)
        and (a["export_conversions"] == [] or a["export_conversions"] == [CONVERSION]),
        "unsupported_or_mutated_conversion_audit",
    )


def _instances(needed, canonical, atoms, records, pair_records, exclusion_records):
    present = {_match_key(r["family"], r["sites"]): r for r in records}
    by_atom = {a["site"]: a for a in atoms}
    results = []
    for family, interactions in needed.items():
        for sites in interactions:
            if family == "atom_typing":
                matches = [by_atom[sites[0]]] if sites[0] in by_atom else []
            elif family in ("interaction_14", "exclusion"):
                rows = pair_records if family == "interaction_14" else exclusion_records
                matches = [r for r in rows if sorted(r["sites"]) == sorted(sites)]
            elif family == "topology":
                matches = (
                    [{"source": canonical["graph_source"]}]
                    if canonical["topology_complete"]
                    else []
                )
            else:
                row = present.get(_match_key(family, sites))
                matches = [row] if row is not None else []
            results.append(
                {
                    "family": family,
                    "sites": list(sites),
                    "status": "source_record_present" if matches else "missing",
                    "sources": [r["source"] for r in matches],
                    "claim_hashes": [identity(r) for r in matches],
                }
            )
    return results


def intake(
    source_dir=None, *, request_sha256=None, manifest_sha256=None, audit_sha256=None
):
    """No parameter bytes read until *all* control metadata and audit pins pass."""
    fence()
    pins = dict(
        zip(CONTROLS, (request_sha256, manifest_sha256, audit_sha256), strict=True)
    )
    r = {
        "schema": "island.j91.intake-result.v1",
        "decision": "source_bundle_required",
        "provider": None,
        "provider_version": None,
        "source_version": None,
        "requested_graph": None,
        "qualified_graph": None,
        "settings": None,
        "source_files_inspected": 0,
        "source_rows_inspected": 0,
        "source_files": [],
        "source_hashes": [],
        "control_hashes": pins,
        "scientific_graphs_covered": 0,
        "synthetic": False,
        "blocking_reasons": [],
        "convention_blockers": [],
        "missing_terms": [],
        "coverage": {},
        "interaction_instances": [],
        "record_claims": [],
        "atom_claims": [],
        "pair_claims": [],
        "exclusion_claims": [],
        "manifest_sha256": manifest_sha256,
        "audit_sha256": audit_sha256,
        "audit_status": "not_inspected",
        "scope": "one_exact_finite_capped_graph_only",
        "metadata_only_replay_requires_external_originals": True,
        "export_conversions": [],
        **FLAGS,
    }
    if source_dir is None or not all(is_hash(p) for p in pins.values()):
        r["blocking_reasons"] = [
            "explicit_external_directory_and_three_control_pins_required"
        ]
        return _result(r)
    metadata_complete = False
    try:
        root = outside_git(source_dir)
        require(root.is_dir(), "source_directory_required")
        paths = {n: external_file(root, n + ".json") for n in CONTROLS}
        require(
            all(p.is_file() for p in paths.values()), "complete_control_files_required"
        )
        q, m, a = [loads(pinned(paths[n], pins[n])) for n in CONTROLS]
        keys(q, REQUEST_FIELDS, "complete_request_metadata_required")
        require(q["schema"] == "island.j91.request.v1", "request_schema_required")
        require(
            q["provider"] in ("polypargen", "ligpargen"),
            "J91_single_polypargen_or_ligpargen_only",
        )
        text(q["source_authority"], "source_authority_required")
        require(
            q["self_contained_finite_molecule"] is True,
            "self_contained_finite_molecule_required",
        )
        graph = finite_graph(q["settings"])
        require(
            type(q["settings"]["repeat_count"]) is int
            and q["settings"]["repeat_count"] > 0,
            "positive_integer_repeat_count_required",
        )
        require(
            {atom.GetSymbol() for atom in molecule(graph).GetAtoms()}
            <= {"C", "H", "O"},
            "J91_initial_CHO_scope_only_no_PDMS",
        )
        require(
            q["graph_identity"] == graph.identity,
            "request_exact_graph_identity_mismatch",
        )
        blockers = _manifest(m, q["provider"], graph)
        require(
            m["license"]["authorization"] == "authorized_local_use",
            "source_authorization_required",
        )
        keys(
            q["expected_conventions"],
            tuple(m["conventions"]),
            "complete_requested_conventions_required",
        )
        for k, v in q["expected_conventions"].items():
            text(v, "requested_convention_missing:" + k)
            if v != m["conventions"][k]:
                blockers.append("requested_convention_mismatch:" + k)
        keys(
            q["nonbonded_definitions"],
            tuple(DEFINITIONS),
            "sigma_epsilon_definitions_required",
        )
        if q["nonbonded_definitions"] != DEFINITIONS:
            blockers.append("sigma_epsilon_minimum_convention_mismatch")
        keys(
            q["expected_forms"],
            ("bond", "angle", "torsion", "improper", "charge", "nonbonded"),
            "complete_structured_expected_forms_required",
        )
        for family, form in q["expected_forms"].items():
            require(
                isinstance(form, dict)
                and form in [v for v in FORMS.values() if v["family"] == family],
                "unsupported_expected_form:" + family,
            )
        _audit_metadata(
            a,
            request_sha256,
            manifest_sha256,
            graph.identity,
            identity(m["conventions"]),
        )
        # Complete metadata validated before directory census or source reads.
        metadata_complete = True
        r.update(
            provider=q["provider"],
            provider_version=m["provider_version"],
            source_version=m["source_version"],
            requested_graph=graph.identity,
            settings=q["settings"],
            synthetic=m["synthetic"],
            license=m["license"],
            provenance=m["provenance"],
            source_authority=q["source_authority"],
            audit_reference={
                "sha256": audit_sha256,
                "reviewer": a["source_audit"]["reviewer"],
                "authority_basis": a["source_audit"]["authority_basis"],
            },
            export_conversions=a["export_conversions"],
        )
        needed, excluded = obligations(graph)
        r["required_interactions"] = {
            k: [list(s) for s in v] for k, v in needed.items()
        } | {"exclusion": [list(p) for p in excluded]}
        expected = {n + ".json" for n in CONTROLS}
        listed = {f["path"] for f in m["files"]}
        require(not expected & listed, "control_source_role_collision")
        actual = set()
        for p in sorted(root.rglob("*")):
            require(not p.is_symlink(), "source_symlink_or_path_substitution")
            if not p.is_dir():
                require(p.is_file(), "regular_source_files_required")
                actual.add(p.relative_to(root).as_posix())
        require(actual == expected | listed, "unlisted_or_missing_bundle_files")
        r["convention_blockers"] = sorted(set(blockers))
        if blockers:
            r["decision"] = "parameter_bundle_incomplete"
            r["blocking_reasons"] = r["convention_blockers"]
            return _blocked(r)
        files = []
        # All bytes pinned first, including opaque attachments. No parsing earlier.
        for f in m["files"]:
            r["source_files_inspected"] += 1
            raw = pinned(external_file(root, f["path"]), f["sha256"])
            require(len(raw) == f["size"], "source_size_mutation")
            files.append((f, raw))
            r["source_files"].append(f)
        r["source_hashes"] = sorted(
            {*pins.values(), *(f["sha256"] for f in m["files"])}
        )
        for f, raw in files:
            if f["role"] == "graph":
                require(
                    f["format"] == "chemical_graph_json" and loads(raw) == data(graph),
                    "source_graph_identity_order_stereo_caps_mismatch",
                )
        sources = [(f, raw) for f, raw in files if f["role"] == "parameters"]
        if not sources:
            r["decision"] = "provenance_only"
            r["blocking_reasons"] = ["no_parameter_records"]
            return _blocked(r)
        f, raw = sources[0]
        parser = OPLSProviderRegistry()._providers[q["provider"]]
        atoms, records, pairs, exclusions = parser.parse(
            raw, f, graph, m["conventions"]
        )
        pairs_meta, exclusions_meta = pair_claims(raw, f, pairs, exclusions)
        for omission_file, omission_raw in files:
            if omission_file["role"] != "omissions":
                continue
            require(
                omission_file["format"] == "improper_omissions_json",
                "unsupported_omission_format",
            )
            omissions = loads(omission_raw)
            require(isinstance(omissions, list), "omission_records_array_required")
            for i, row in enumerate(omissions):
                keys(
                    row,
                    ("sites", "reason", "location"),
                    "explicit_improper_omission_required",
                )
                text(row["location"], "exact_source_location_required")
                records.append(
                    {
                        "id": row_key("improper", row["sites"]),
                        "family": "improper",
                        "sites": row["sites"],
                        "form": None,
                        "values": [],
                        "applicability": "not_applicable",
                        "reason": row["reason"],
                        "source": pointer(omission_file, f"/{i}; " + row["location"]),
                    }
                )
        r["source_rows_inspected"] = (
            len(atoms) + len(records) + len(pairs_meta) + len(exclusions_meta)
        )
        c, diagnostics = _validate(
            graph, atoms, records, pairs, exclusions, needed, excluded, m
        )
        for row in records:
            if (
                row["form"] is not None
                and FORMS[row["form"]] != q["expected_forms"][row["family"]]
            ):
                blockers.append("source_requested_form_mismatch:" + row["family"])
        r["convention_blockers"] = sorted(set(blockers))
        c.update(
            schema="island.canonical-opls.v1",
            provider=q["provider"],
            provider_version=m["provider_version"],
            source_version=m["source_version"],
            profile_id=m["profile_id"],
            family="OPLS-AA",
            graph_identity=graph.identity,
            graph_source=pointer(
                next(f for f in m["files"] if f["role"] == "graph"), "/"
            ),
            atom_order=list(graph.atom_order),
            settings_identity=identity(q["settings"]),
            conventions=m["conventions"],
            source_hashes=r["source_hashes"],
            manifest_sha256=manifest_sha256,
            license=m["license"],
            provenance=m["provenance"],
            synthetic=m["synthetic"],
            **FLAGS,
        )
        # Replace J90's derived exclusion list with only the explicit source list.
        c["exclusions"] = sorted([sorted(p) for p in exclusions])
        r["record_claims"] = [claim_metadata(row) for row in c["records"]]
        r["atom_claims"] = [
            {
                "site": atom["site"],
                "source": atom["source"],
                "claim_sha256": identity(atom),
            }
            for atom in c["atoms"]
        ]
        r["pair_claims"], r["exclusion_claims"] = pairs_meta, exclusions_meta
        _audit(a["source_audit"], r, graph.identity, identity(m["conventions"]))
        for name, rows in (
            ("pair_claim_hashes", pairs_meta),
            ("exclusion_claim_hashes", exclusions_meta),
        ):
            require(
                sorted(a[name]) == sorted(identity(row) for row in rows),
                "independent_audit_pair_or_exclusion_claim_mismatch",
            )
        r["audit_status"] = "independent_exact_graph_source_and_request_pins_verified"
        r["coverage"], r["missing_terms"] = _coverage(needed, c)
        missing_exclusions = [
            list(p) for p in excluded if list(p) not in c["exclusions"]
        ]
        r["coverage"]["exclusion"] = {
            "required": len(excluded),
            "covered": len(excluded) - len(missing_exclusions),
            "status": "incomplete" if missing_exclusions else "complete",
        }
        r["missing_terms"] += [
            {
                "family": "exclusion",
                "sites": p,
                "reason": "explicit_source_exclusion_missing",
            }
            for p in missing_exclusions
        ]
        r["interaction_instances"] = _instances(
            {**needed, "exclusion": excluded},
            c,
            atoms,
            records,
            pairs_meta,
            exclusions_meta,
        )
        if r["missing_terms"]:
            blockers.append("explicit_missing_terms_no_substitution")
        blockers.extend(diagnostics)
        r["blocking_reasons"] = sorted(set(blockers))
        r["decision"] = (
            "parameter_bundle_incomplete" if blockers else "audited_exact_graph"
        )
        if not blockers and not m["synthetic"]:
            r["qualified_graph"] = graph.identity
            r["scientific_graphs_covered"] = 1
        canonical = CanonicalOPLS.create(c)
        r["canonical_identity"] = canonical.identity
        return _result(r, canonical)
    except (
        ValueError,
        TypeError,
        KeyError,
        OSError,
        UnicodeError,
        AttributeError,
        IndexError,
        IslandError,
    ):
        # Errors must not echo arbitrary proprietary source literals into receipts.
        import sys

        error = sys.exception()
        reason = (
            str(error)
            if type(error) is ValueError
            else type(error).__name__ + "_invalid_input"
        )
        r["blocking_reasons"] = ["intake_blocked:" + reason]
        r["decision"] = (
            "parameter_bundle_incomplete"
            if metadata_complete
            else "source_bundle_required"
        )
        return _blocked(r)


def _blocked(r):
    required = r.get("required_interactions", {})
    r["coverage"] = {
        k: {"required": len(v), "covered": 0, "status": "blocked"}
        for k, v in required.items()
    }
    r["missing_terms"] = [
        {"family": k, "sites": s, "reason": "source_or_convention_not_established"}
        for k, values in required.items()
        for s in values
    ]
    r["interaction_instances"] = [
        {
            "family": t["family"],
            "sites": t["sites"],
            "status": "unverified",
            "sources": [],
            "claim_hashes": [],
        }
        for t in r["missing_terms"]
    ]
    return _result(r)


def export(source_dir, *, destination, engine, **pins):
    fence()
    dest = outside_git(destination)
    require(not dest.exists(), "export_destination_must_be_new")
    require(
        not dest.is_relative_to(Path(source_dir).absolute()),
        "export_must_not_modify_source_bundle",
    )
    result = intake(source_dir, **pins)  # Rehash external originals before each export.
    r = result.report
    require(
        result.canonical is not None
        and r["decision"] in ("audited_exact_graph", "provisional_parameterization"),
        "incomplete_result_export_forbidden",
    )
    require(
        engine == r["settings"]["intended_export_engine"],
        "intended_export_engine_mismatch",
    )
    if engine == "lammps":
        require(
            r["export_conversions"] == [CONVERSION],
            "independent_conversion_audit_required",
        )
        require(
            all(
                row["family"] != "torsion" or row["form"] == "torsion_periodic"
                for row in result.canonical.payload["records"]
            ),
            "torsion_conversion_not_authorized",
        )
    raw = render(result.canonical.payload, r, engine).encode()
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("xb") as stream:
        stream.write(raw)
    import hashlib

    return {
        "sha256": hashlib.sha256(raw).hexdigest(),
        "source_hashes": r["source_hashes"],
        "graph_identity": r["requested_graph"],
        "decision": r["decision"],
        **FLAGS,
    }
