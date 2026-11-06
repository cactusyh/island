"""One explicitly selected OPLS provider per exact finite graph; no runtime hook."""

from dataclasses import dataclass, field
from decimal import Decimal
from itertools import combinations
from pathlib import Path

from island.forcefields.assignment_core.schema import data, identity, is_hash, require
from island.forcefields.assignment_core.topology import as_topology, molecule
from island.forcefields.bundle_adapters.intake import _read
from island.forcefields.bundle_adapters.schema import keys, loads, relative_path, text
from island.forcefields.parameters.inventory import derive_interaction_inventory

from .parsers import canonical_json, claim_metadata, gromacs, pointer
from .polymer import finite_graph
from .schema import (
    DIMENSIONS,
    FLAGS,
    PROVIDERS,
    CanonicalOPLS,
    convention_blockers,
    number,
    row_key,
    sealed,
    validate_row,
)

MANIFEST_FIELDS = (
    "schema",
    "provider",
    "provider_version",
    "source_version",
    "profile_id",
    "family",
    "graph_identity",
    "files",
    "license",
    "provenance",
    "conventions",
    "synthetic",
)


def fence():
    # Same fence as J81-J89. Install via run_phase_4j72_c1 before external IO.
    from scripts.phase_4j76_source_resolution import scope

    scope()


def external_file(root, name):
    parts = relative_path(name).parts
    candidate = root
    for part in parts:
        candidate /= part
        require(not candidate.is_symlink(), "source_symlink_or_path_substitution")
    require(candidate.resolve().is_relative_to(root), "source_root_escape")
    return candidate


def pinned(path, checksum):
    path = Path(path).absolute()
    require(
        not any(p.is_symlink() for p in (path, *path.parents)),
        "source_symlink_or_path_substitution",
    )
    return _read(path, checksum)


def obligations(graph):
    """Reuse the existing graph-only interaction inventory, not element heuristics."""
    topology = as_topology(graph)
    inv = derive_interaction_inventory(topology)
    needed = {name: [] for name in DIMENSIONS}
    for site in graph.atom_order:
        for name in ("atom_typing", "charge", "nonbonded"):
            needed[name].append((site,))
    needed["bond"] = list(inv.bonds)
    needed["angle"] = list(inv.angles)
    needed["torsion"] = list(inv.proper_torsions)
    needed["improper"] = [
        (a, center, b, c)
        for center in graph.atom_order
        for a, b, c in combinations(sorted(topology.neighbors(center)), 3)
    ]
    pairs, excluded = [], []
    for a in graph.atom_order:
        frontier, visited = {a}, {a}
        for distance in range(1, 4):
            frontier = {v for u in frontier for v in topology.neighbors(u)} - visited
            visited |= frontier
            target = excluded if distance <= 2 else pairs
            target.extend((a, b) for b in sorted(frontier) if a < b)
    needed["interaction_14"] = sorted(pairs)
    needed["topology"] = [(1,)]
    return needed, sorted(excluded)


def _pairs(values, atom_order):
    require(isinstance(values, list), "pair_array_required")
    output = []
    for pair in values:
        require(
            isinstance(pair, list)
            and len(pair) == 2
            and all(type(i) is int and i in atom_order for i in pair)
            and pair[0] != pair[1],
            "invalid_pair_sites",
        )
        output.append(tuple(sorted(pair)))
    require(len(output) == len(set(output)), "duplicate_pair_record")
    return sorted(output)


def _manifest(m, provider, graph):
    keys(m, MANIFEST_FIELDS, "provider_manifest_schema")
    require(
        m["schema"] == "island.opls-provider-bundle.v1", "unsupported_provider_manifest"
    )
    require(m["provider"] == provider, "conflicting_provider_identity_no_merge")
    require(m["family"] == "OPLS-AA", "OPLS_AA_only_no_family_fallback")
    require(
        m["graph_identity"] == graph.identity,
        "topology_changed_provider_reassignment_required",
    )
    for k in ("provider_version", "source_version", "profile_id"):
        text(m[k], "provider_metadata_missing:" + k)
    require(type(m["synthetic"]) is bool, "synthetic_flag_required")
    keys(
        m["license"],
        ("identifier", "authorization", "redistribution"),
        "license_authorization_required",
    )
    for v in m["license"].values():
        text(v, "license_metadata_missing")
    keys(
        m["provenance"],
        ("citation", "acquisition", "source_uri"),
        "provider_provenance_required",
    )
    for v in m["provenance"].values():
        text(v, "provider_provenance_missing")
    require(isinstance(m["files"], list) and m["files"], "file_manifest_required")
    for file in m["files"]:
        keys(file, ("path", "sha256", "size", "role", "format"), "provider_file_schema")
        relative_path(file["path"])
        require(
            is_hash(file["sha256"]) and type(file["size"]) is int and file["size"] > 0,
            "file_hash_and_size_required",
        )
        require(
            file["role"] in ("graph", "parameters", "omissions", "opaque"),
            "unsupported_file_role",
        )
    require(
        len({f["path"] for f in m["files"]}) == len(m["files"]), "duplicate_source_path"
    )
    require(
        sum(f["role"] == "graph" for f in m["files"]) == 1,
        "exactly_one_pinned_graph_required",
    )
    require(
        sum(f["role"] == "parameters" for f in m["files"]) <= 1,
        "multiple_parameter_sources_no_silent_merge",
    )
    return convention_blockers(m["conventions"])


@dataclass(frozen=True)
class ProviderResult:
    report_json: str
    canonical: CanonicalOPLS | None = field(default=None, repr=False)

    @property
    def report(self):
        return loads(self.report_json)


def _result(value, canonical=None):
    import json

    return ProviderResult(json.dumps(sealed(value), sort_keys=True), canonical)


class OPLSProvider:
    """Offline interface: subclasses never contact a server or launch a provider."""

    provider_id = ""
    formats = ()

    def parse(self, raw, file, graph, conventions):
        require(file["format"] in self.formats, "unsupported_provider_format")
        if file["format"] == "gromacs_itp":
            return gromacs(raw, file, graph.identity, conventions)
        return canonical_json(
            raw, file, graph.identity, database=file["format"] == "oplsaa_db_snapshot"
        )


class PolyParGenProvider(OPLSProvider):
    provider_id = "polypargen"
    formats = ("gromacs_itp", "canonical_json")


class LigParGenProvider(OPLSProvider):
    provider_id = "ligpargen"
    formats = ("gromacs_itp", "canonical_json")


class OPLSAADBProvider(OPLSProvider):
    provider_id = "oplsaa_db"
    formats = ("oplsaa_db_snapshot",)


class AuditedBundleProvider(OPLSProvider):
    provider_id = "audited_bundle"
    formats = ("canonical_json",)


class OPLSProviderRegistry:
    """Parser interface registry, deliberately separate from runtime registration."""

    def __init__(self):
        self._providers = {
            p.provider_id: p
            for p in (
                PolyParGenProvider(),
                LigParGenProvider(),
                OPLSAADBProvider(),
                AuditedBundleProvider(),
            )
        }

    @property
    def providers(self):
        return tuple(sorted(self._providers))

    def assign(
        self,
        settings,
        *,
        provider,
        manifest=None,
        manifest_sha256=None,
        audit=None,
        audit_sha256=None,
    ):
        fence()
        require(
            isinstance(provider, str) and provider in PROVIDERS,
            "one_explicit_OPLS_provider_required_no_merge",
        )
        graph = finite_graph(settings)
        needed, excluded = obligations(graph)
        report = {
            "schema": "island.opls-provider-result.v1",
            "provider": provider,
            "settings_identity": identity(settings),
            "graph_identity": graph.identity,
            "scope": "this_finite_capped_oligomer_only",
            "decision": "source_bundle_required",
            "source_files_inspected": 0,
            "source_rows_inspected": 0,
            "scientific_graphs_covered": 0,
            "synthetic": False,
            "source_hashes": [],
            "source_files": [],
            "record_claims": [],
            "manifest_sha256": manifest_sha256,
            "audit_sha256": audit_sha256,
            "blocking_reasons": [],
            "missing_terms": [],
            "coverage": {},
            "audit_status": "not_supplied",
            **FLAGS,
        }
        # No unpinned source is opened, including when just one pin is absent.
        if manifest is None or not is_hash(manifest_sha256):
            report["blocking_reasons"].append(
                "pinned_user_supplied_source_bundle_required"
            )
            return self._blocked(report, needed)
        try:
            m = loads(pinned(manifest, manifest_sha256))
            blockers = _manifest(m, provider, graph)
            report["synthetic"] = m["synthetic"]
            report["provider_version"] = m["provider_version"]
            report["source_version"] = m["source_version"]
            report["license"] = m["license"]
            report["provenance"] = m["provenance"]
            if m["license"]["authorization"] != "authorized_local_use":
                report["blocking_reasons"].append(
                    "explicit_source_authorization_required"
                )
                return self._blocked(report, needed)
            root = Path(manifest).absolute().parent
            files = []
            # Rehash *all* source bytes, including opaque attachments, before parsing.
            for file in m["files"]:
                raw = pinned(external_file(root, file["path"]), file["sha256"])
                require(len(raw) == file["size"], "source_size_mutation")
                files.append((file, raw))
                report["source_files_inspected"] += 1
                report["source_files"].append(file)
            report["source_hashes"] = sorted(
                {manifest_sha256, *(f["sha256"] for f in m["files"])}
            )
            for file, raw in files:
                if file["role"] == "graph":
                    require(
                        file["format"] == "chemical_graph_json"
                        and loads(raw) == data(graph),
                        "provider_graph_bytes_atom_order_stereo_or_caps_mismatch",
                    )
            records_files = [(f, raw) for f, raw in files if f["role"] == "parameters"]
            if not records_files:
                report["decision"] = "provenance_only"
                report["blocking_reasons"] = [
                    "no_parseable_parameter_records",
                    *blockers,
                ]
                return self._blocked(report, needed)
            file, raw = records_files[0]
            atoms, records, pairs, extra_exclusions = self._providers[provider].parse(
                raw, file, graph, m["conventions"]
            )
            for file, raw in files:
                if file["role"] != "omissions":
                    continue
                require(
                    file["format"] == "improper_omissions_json",
                    "unsupported_omission_format",
                )
                omissions = loads(raw)
                require(isinstance(omissions, list), "omission_records_array_required")
                for i, row in enumerate(omissions):
                    keys(
                        row,
                        ("sites", "reason", "location"),
                        "explicit_improper_omission_required",
                    )
                    records.append(
                        {
                            "id": row_key("improper", row["sites"]),
                            "family": "improper",
                            "sites": row["sites"],
                            "form": None,
                            "values": [],
                            "applicability": "not_applicable",
                            "reason": row["reason"],
                            "source": pointer(file, f"/{i}; " + row["location"]),
                        }
                    )
            report["source_rows_inspected"] = len(atoms) + len(records)
            canonical, diagnostics = _validate(
                graph, atoms, records, pairs, extra_exclusions, needed, excluded, m
            )
            blockers.extend(diagnostics)
            report["record_claims"] = [claim_metadata(r) for r in canonical["records"]]
            report["atom_claims"] = [
                {"site": a["site"], "source": a["source"], "claim_sha256": identity(a)}
                for a in canonical["atoms"]
            ]
            canonical.update(
                {
                    "schema": "island.canonical-opls.v1",
                    "provider": provider,
                    "provider_version": m["provider_version"],
                    "source_version": m["source_version"],
                    "profile_id": m["profile_id"],
                    "family": "OPLS-AA",
                    "graph_identity": graph.identity,
                    "atom_order": list(graph.atom_order),
                    "settings_identity": identity(settings),
                    "conventions": m["conventions"],
                    "source_hashes": report["source_hashes"],
                    "manifest_sha256": manifest_sha256,
                    "license": m["license"],
                    "provenance": m["provenance"],
                    "synthetic": m["synthetic"],
                    **FLAGS,
                }
            )
            audited = False
            if audit is not None or audit_sha256 is not None:
                a = loads(pinned(audit, audit_sha256))
                _audit(a, report, graph.identity, identity(m["conventions"]))
                audited = True
                report["source_hashes"] = sorted(
                    {*report["source_hashes"], audit_sha256}
                )
                canonical["source_hashes"] = report["source_hashes"]
                report["audit_status"] = "independent_exact_graph_audit_pin_verified"
            if provider == "audited_bundle" and not audited:
                report["blocking_reasons"] = [
                    "independent_audit_required_for_audited_provider"
                ]
                return self._blocked(report, needed)
            report["coverage"], report["missing_terms"] = _coverage(needed, canonical)
            if report["missing_terms"]:
                blockers.append("explicit_missing_terms_no_substitution")
            report["blocking_reasons"] = sorted(set(blockers))
            report["decision"] = (
                "parameter_bundle_incomplete"
                if blockers
                else "audited_exact_graph"
                if audited
                else "provisional_parameterization"
            )
            report["scientific_graphs_covered"] = int(
                not blockers and not m["synthetic"]
            )
            payload = CanonicalOPLS.create(canonical)
            report["canonical_identity"] = payload.identity
            return _result(report, payload)
        except (ValueError, TypeError, KeyError, OSError, UnicodeError) as error:
            report["decision"] = (
                "source_bundle_required"
                if str(error).startswith("source_file_unavailable")
                else "parameter_bundle_incomplete"
            )
            report["blocking_reasons"] = ["source_validation_blocked:" + str(error)]
            return self._blocked(report, needed)

    @staticmethod
    def _blocked(report, needed):
        report["coverage"] = {
            name: {"required": len(sites), "covered": 0, "status": "blocked"}
            for name, sites in needed.items()
        }
        report["missing_terms"] = [
            {
                "family": name,
                "sites": list(s),
                "reason": "source_or_convention_not_established",
            }
            for name, sites in needed.items()
            for s in sites
        ]
        return _result(report)


def _validate(
    graph, atoms, records, pairs, extra_exclusions, needed, excluded, manifest
):
    require(
        isinstance(atoms, list) and isinstance(records, list),
        "canonical_arrays_required",
    )
    mol = molecule(graph)
    files = {f["path"]: f for f in manifest["files"]}
    diagnostics = []
    for atom in atoms:
        keys(
            atom,
            ("site", "atom_type", "element", "mass", "source"),
            "canonical_atom_schema",
        )
        site = atom["site"]
        require(type(site) is int and site in graph.atom_order, "invalid_atom_mapping")
        text(atom["atom_type"], "atom_type_missing")
        require("*" not in atom["atom_type"], "wildcard_typing_forbidden")
        expected = mol.GetAtomWithIdx(site - 1)
        require(
            atom["element"] in (None, expected.GetSymbol()),
            "source_element_mapping_mismatch",
        )
        require(number(atom["mass"]) > 0, "invalid_mass")
        # Loose mass identity check detects swapped heavy/H maps, never types atoms.
        require(
            abs(number(atom["mass"]) - Decimal(str(expected.GetMass())))
            < Decimal("0.15"),
            "source_mass_atom_mapping_mismatch",
        )
        atom["element"] = expected.GetSymbol()
        keys(
            atom["source"],
            ("path", "sha256", "pointer"),
            "atom_source_pointer_required",
        )
        text(atom["source"]["pointer"], "atom_source_pointer_required")
    require(
        len({a["site"] for a in atoms}) == len(atoms),
        "ambiguous_or_duplicate_atom_assignment",
    )
    for row in records:
        validate_row(row, graph.atom_order)
        require(
            row["source"]["path"] in files
            and row["source"]["sha256"] == files[row["source"]["path"]]["sha256"],
            "unbound_source_row_hash",
        )
        # Validate arbitrary improper order by center and neighbor set; preserve
        # the original peripheral order for its signed functional form.
        family, sites = row["family"], tuple(row["sites"])
        valid = {row_key(family, s) for s in needed[family]}
        if family == "improper":
            valid = {(s[1], tuple(sorted((s[0], s[2], s[3])))) for s in needed[family]}
            require(
                (sites[1], tuple(sorted((sites[0], sites[2], sites[3])))) in valid,
                "improper_topology_or_ordering_mismatch",
            )
        else:
            require(row_key(family, sites) in valid, "interaction_not_in_exact_graph")
    for atom in atoms:
        require(
            atom["source"]["path"] in files
            and atom["source"]["sha256"] == files[atom["source"]["path"]]["sha256"],
            "unbound_atom_source_hash",
        )
    require(
        len({r["id"] for r in records}) == len(records), "duplicate_record_identity"
    )
    keys_seen = [_match_key(r["family"], r["sites"]) for r in records]
    require(
        len(keys_seen) == len(set(keys_seen)),
        "conflicting_or_duplicate_parameter_candidates",
    )
    pairs = _pairs(pairs, graph.atom_order)
    extras = _pairs(extra_exclusions, graph.atom_order)
    require(set(extras) <= set(excluded), "extra_exclusions_change_exact_graph_policy")
    require(set(pairs) <= set(needed["interaction_14"]), "invalid_14_pair_topology")
    charges = [number(r["values"][0][0]) for r in records if r["family"] == "charge"]
    if len(charges) == len(graph.atom_order):
        # No adjusted/renormalized charges. Decimal source sum must be exact.
        total = sum(charges, Decimal(0))
        formal = sum(a.GetFormalCharge() for a in mol.GetAtoms())
        if total != formal:
            diagnostics.append("charge_conservation_failed")
    source_bonds = {tuple(sorted(r["sites"])) for r in records if r["family"] == "bond"}
    topology_complete = source_bonds == {tuple(sorted(s)) for s in needed["bond"]} and {
        a["site"] for a in atoms
    } == set(graph.atom_order)
    return {
        "atoms": sorted(atoms, key=lambda a: a["site"]),
        "records": sorted(records, key=lambda r: r["id"]),
        "pairs14": [list(p) for p in pairs],
        "exclusions": [list(p) for p in excluded],
        "topology_complete": topology_complete,
    }, diagnostics


def _match_key(family, sites):
    if family == "improper":
        return (family, sites[1], tuple(sorted((sites[0], sites[2], sites[3]))))
    return row_key(family, sites)


def _coverage(needed, canonical):
    present = {_match_key(r["family"], r["sites"]) for r in canonical["records"]}
    missing, coverage = [], {}
    for family, sites in needed.items():
        count = 0
        for s in sites:
            if family == "atom_typing":
                found = s[0] in {a["site"] for a in canonical["atoms"]}
            elif family == "interaction_14":
                found = list(s) in canonical["pairs14"]
            elif family == "topology":
                found = canonical["topology_complete"]
            else:
                found = _match_key(family, s) in present
            count += bool(found)
            if not found:
                missing.append(
                    {
                        "family": family,
                        "sites": list(s),
                        "reason": "exact_source_term_missing_no_fallback",
                    }
                )
        coverage[family] = {
            "required": len(sites),
            "covered": count,
            "status": "complete" if count == len(sites) else "incomplete",
        }
    return coverage, missing


def _audit(audit, report, graph, conventions):
    keys(
        audit,
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
        "independent_audit_schema",
    )
    require(
        audit["schema"] == "island.opls-provider-audit.v1"
        and audit["independent"] is True,
        "independent_audit_required",
    )
    for name in ("reviewer", "authority_basis"):
        text(audit[name], "trusted_audit_identity_required")
        require(
            audit[name].lower() not in ("unknown", "none", "missing"),
            "trusted_audit_identity_required",
        )
    require(
        audit["manifest_sha256"] == report["manifest_sha256"]
        and audit["graph_identity"] == graph
        and audit["convention_identity"] == conventions,
        "audit_scope_manifest_or_convention_mismatch",
    )
    for key, rows in (
        ("claim_hashes", report["record_claims"]),
        ("atom_claim_hashes", report["atom_claims"]),
    ):
        require(
            isinstance(audit[key], list)
            and sorted(audit[key]) == sorted(r["claim_sha256"] for r in rows),
            "independent_audit_claim_coverage_mismatch",
        )
