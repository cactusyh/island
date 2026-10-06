"""Explicit local FRC comparison and source-bound queries, never a mixed model."""

from collections import Counter
from dataclasses import dataclass
from decimal import Decimal
from itertools import combinations
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.exceptions import PCFFError
from island.workflows.storage import publish

from .catalog import AUTO_ROLES, EXTRA, catalog_from_inventory
from .charges import bond_selection, identity
from .class2 import FAMILIES, reverse_values
from .fallbacks import DOMAIN_POLICY, lookup, orientations, resolve, validate_policy
from .source import (
    FLAGS,
    PIN,
    PCFFSource,
    boundary,
    digest,
    parse_frc,
    records,
    require,
    select,
)

VARIANT_SCHEMA = "island_pcff_source_variant_inventory_v1"
COMPARISON_SCHEMA = "island_pcff_source_comparison_v1"
RESOLUTION_SCHEMA = "island_pcff_variant_resolution_v1"
CANDIDATE_PROFILE = "island_pcff_frc_candidate_v1"
POLICY = "island_pcff_explicit_source_queries_v1"
PROVENANCE_FIELDS = {
    "repository",
    "revision",
    "source_path",
    "family",
    "license_status",
    "citation",
}
STATES = {
    "source_row_present",
    "source_row_absent",
    "source_variant_only",
    "wildcard_ambiguous",
    "equilibrium_dependency_missing",
    "policy_derived_zero",
    "native_charge_incomplete",
    "model_incomplete",
}


def sha_valid(value):
    return (
        type(value) is str
        and len(value) == 64
        and all(c in "0123456789abcdef" for c in value)
    )


def decimal_token(token):
    try:
        value = Decimal(token)
    except Exception:  # noqa: BLE001 -- chemical tokens are deliberately not numbers
        return token
    require(value.is_finite(), "Nonfinite FRC token")
    return str(value.normalize()) if value else "0"


def variant_inventory(raw, expected_sha256, provenance):
    require(
        sha_valid(expected_sha256) and digest(raw) == expected_sha256,
        "Variant source hash mismatch",
    )
    require(
        type(provenance) is dict and set(provenance) == PROVENANCE_FIELDS,
        "Complete explicit source provenance required",
    )
    require(
        all(type(v) is str and v.strip() for v in provenance.values()),
        "Invalid source provenance",
    )
    require(
        provenance["family"] in ("PCFF", "PCFF-IFF", "unreviewed PCFF candidate"),
        "Explicit PCFF-family source classification required",
    )
    require(
        expected_sha256 != PIN["sha256"] or provenance["family"] == "PCFF",
        "Pinned PCFF bytes contradict declared source family",
    )
    inv = parse_frc(raw)
    profile = PIN["profile"] if expected_sha256 == PIN["sha256"] else CANDIDATE_PROFILE
    key = {"sha256": expected_sha256, "profile": profile}
    catalog = catalog_from_inventory(inv, key)
    numerical = {r["id"]: r for r in catalog["records"]}
    rows, sections = [], []
    for section in inv["sections"]:
        name, namespace = section["name"], section["namespace"]
        semantic = {r["line"]: r for r in section["records"]}
        ids = []
        for line in section["lines"]:
            text = line["raw"].strip()
            if not text or text.startswith(("!", ">", "@")):
                continue
            fields = text.split("!", 1)[0].split()
            if name in ("version", "reference", "end"):
                continue
            if not fields or not fields[0][0].isdigit():
                continue
            require(
                len(fields) >= 2 and fields[1].isdigit(),
                f"Invalid FRC record at {line['line']}",
            )
            version = Decimal(fields[0])
            require(version.is_finite() and version >= 0, "Invalid variant row version")
            rid = f"{name}:{namespace}:{line['line']}"
            parsed = semantic.get(line["line"])
            number = numerical.get(rid)
            arity = (
                FAMILIES[name][0]
                if name in FAMILIES
                else EXTRA[name][0]
                if name in EXTRA
                else (
                    2
                    if name == "bond_increments"
                    else 1
                    if name in ("atom_types", "equivalence", "auto_equivalence")
                    else 0
                )
            )
            tokens = fields[2 : 2 + arity] + [
                decimal_token(t) for t in fields[2 + arity :]
            ]
            content = {
                "section": name,
                "namespace": namespace,
                "version": decimal_token(fields[0]),
                "reference": fields[1],
                "tokens": tokens,
            }
            row = {
                "id": rid,
                "line": line["line"],
                "content": content,
                "row_identity": identity(content),
                "types": fields[2 : 2 + arity],
                "interpretation": "semantic"
                if parsed
                else "numerical"
                if number
                else "raw_only",
            }
            if parsed:
                row["semantic"] = parsed["data"]
            if number:
                row["numerical"] = {k: v for k, v in number.items() if k != "raw"}
            rows.append(row)
            ids.append(rid)
        sections.append(
            {
                "name": name,
                "namespace": namespace,
                "line": section["line"],
                "record_ids": ids,
                "record_count": len(ids),
                "interpretation": "supported_record_form"
                if name in FAMILIES or name in EXTRA
                else section["interpretation"],
                "raw_lines": section["lines"]
                if name not in FAMILIES and name not in EXTRA and not semantic
                else [],
            }
        )
    counts = dict(Counter(r["content"]["section"] for r in rows))
    return {
        "schema": VARIANT_SCHEMA,
        "selection": key,
        "provenance": provenance,
        "declarations": inv["declarations"],
        "sections": sections,
        "rows": rows,
        "counts": counts,
        "atom_labels": sorted(
            {
                r["semantic"]["type"]
                for r in rows
                if r["content"]["section"] == "atom_types"
            }
        ),
        "native_assignment_authorized": expected_sha256 == PIN["sha256"],
        "semantics": "Inventory and queries only; candidate profiles do not authorize typing, charge assignment or model assembly",
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFSourceVariant:
    """Owned source bytes and typed provenance; runtime never downloads a file."""

    raw: bytes
    expected_sha256: str
    provenance_json: str

    @boundary
    def validate_integrity(self):
        require(type(self.raw) is bytes, "Immutable source bytes required")
        variant_inventory(self.raw, self.expected_sha256, unpack(self.provenance_json))

    @property
    def payload(self):
        self.validate_integrity()
        return variant_inventory(
            self.raw, self.expected_sha256, unpack(self.provenance_json)
        )

    @property
    def selection(self):
        return self.payload["selection"]

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def load_pcff_source_variant(path, *, expected_sha256, provenance):
    result = PCFFSourceVariant(
        Path(path).read_bytes(), expected_sha256, pack(provenance)
    )
    result.validate_integrity()
    return result


def selections(variants):
    require(
        type(variants) in (list, tuple) and bool(variants),
        "Nonempty explicit variants required",
    )
    require(
        all(type(v) is PCFFSourceVariant for v in variants),
        "Expected PCFF source variants",
    )
    payloads = [v.payload for v in variants]
    require(
        len({identity(p) for p in payloads}) == len(payloads),
        "Duplicate variant identity",
    )
    return payloads


def compare_data(payloads):
    payloads = sorted(payloads, key=identity)
    differences = []
    for a, b in combinations(payloads, 2):
        ar, br = (
            {r["row_identity"]: r for r in a["rows"]},
            {r["row_identity"]: r for r in b["rows"]},
        )
        differences.append(
            {
                "left": identity(a),
                "right": identity(b),
                "left_selection": a["selection"],
                "right_selection": b["selection"],
                "same_source_bytes": a["selection"]["sha256"]
                == b["selection"]["sha256"],
                "left_only": [ar[k] for k in sorted(ar.keys() - br.keys())],
                "right_only": [br[k] for k in sorted(br.keys() - ar.keys())],
                "common_content_count": len(ar.keys() & br.keys()),
                "duplicate_candidates_preserved_in_inventories": True,
            }
        )
    return {
        "schema": COMPARISON_SCHEMA,
        "variants": sorted(payloads, key=identity),
        "differences": sorted(differences, key=lambda d: (d["left"], d["right"])),
        "policy": "No mixed rows or source authority inferred from a row difference",
        "full_source_complete": False,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFSourceComparison:
    json_text: str
    variants: tuple

    @boundary
    def validate_integrity(self):
        require(type(self.json_text) is str, "Data-only comparison required")
        p = unpack(self.json_text)
        require(
            p["schema"] == COMPARISON_SCHEMA, "Unsupported source-comparison schema"
        )
        require(
            pack(p) == pack(compare_data(selections(self.variants))),
            "Contradictory PCFF source comparison",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def compare_pcff_sources(variants):
    variants = tuple(variants)
    result = PCFFSourceComparison(pack(compare_data(selections(variants))), variants)
    result.validate_integrity()
    return result


@boundary
def save_pcff_source_comparison(result, path):
    require(type(result) is PCFFSourceComparison, "Expected source comparison")
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_source_comparison(path, *, variants):
    result = PCFFSourceComparison(Path(path).read_text(), tuple(variants))
    result.validate_integrity()
    return result


@dataclass(frozen=True)
class PCFFVariantSelection:
    sha256: str
    profile: str

    def payload(self):
        require(
            sha_valid(self.sha256) and type(self.profile) is str and bool(self.profile),
            "Explicit source hash/profile required",
        )
        return {"sha256": self.sha256, "profile": self.profile}


@dataclass(frozen=True)
class PCFFVariantPolicy:
    primary: PCFFVariantSelection
    fallbacks: tuple = ()
    resolution_policy: str = DOMAIN_POLICY

    def payload(self):
        require(
            type(self.primary) is PCFFVariantSelection
            and type(self.fallbacks) is tuple,
            "Typed source selection and explicit fallback tuple required",
        )
        require(
            all(type(f) is PCFFVariantSelection for f in self.fallbacks),
            "Typed fallback selections required",
        )
        items = [self.primary.payload()] + [f.payload() for f in self.fallbacks]
        require(
            len({pack(x) for x in items}) == len(items), "Repeated source selection"
        )
        validate_policy(self.resolution_policy)
        return {
            "name": POLICY,
            "primary": items[0],
            "fallbacks": items[1:],
            "resolution_policy": self.resolution_policy,
            "fallback_semantics": "Explicit source-bound diagnostic query on missing rows only; never an assembled mixed-source model",
        }


def chosen_variant(variants, selection):
    found = [v for v in variants if v.selection == selection]
    require(bool(found), "Explicit source hash/profile unavailable")
    # Equal-byte mirrors may have different acquisition provenance; retain every origin.
    found.sort(key=lambda v: v.identity)
    return found[0], [v.identity for v in found]


@boundary
def select_native_pcff_source(variants, selection):
    require(
        type(selection) is PCFFVariantSelection,
        "Explicit typed source selection required",
    )
    selections(tuple(variants))
    source, _ = chosen_variant(variants, selection.payload())
    native = PCFFSource(source.raw, source.expected_sha256)
    native.require_assignment()
    return native


def increment_query(inventory, labels):
    inc = records(inventory, "bond_increments")
    paths = [("direct", list(labels), [])]
    for kind, role in (("equivalence", "bond"), ("auto_equivalence", "bond_increment")):
        eq = records(inventory, kind)
        evidence = [select([r for r in eq if r["data"]["type"] == t]) for t in labels]
        if all(evidence):
            paths.append(
                (
                    kind + "." + role,
                    [e["record"]["data"]["families"][role] for e in evidence],
                    evidence,
                )
            )
    searches = []
    for path, types, evidence in paths:
        try:
            result = bond_selection(types, inc)
        except PCFFError as error:
            return {
                "status": "ambiguous",
                "reason": str(error),
                "searches": searches,
                "candidate_ids": [
                    r["id"] for r in inc if r["data"]["types"] in (types, types[::-1])
                ],
            }
        searches.append(
            {
                "path": path,
                "types": types,
                "equivalence_evidence": evidence,
                "candidate_ids": [
                    r["id"] for r in inc if r["data"]["types"] in (types, types[::-1])
                ],
            }
        )
        if result:
            return {
                "status": "assigned",
                "path": path,
                "resolved_types": types,
                "selected": result,
                "normalized_values": result["increments"],
                "searches": searches,
            }
    return {"status": "missing", "searches": searches}


def query_context(variant):
    data = variant.payload
    return {
        "variant": data,
        "identity": identity(data),
        "inventory": parse_frc(variant.raw),
        "catalog": {
            r["numerical"]["id"]: r["numerical"]
            for r in data["rows"]
            if "numerical" in r
        },
    }


def query_one(variant, family, namespace, labels, resolution_policy):
    return query_compiled(
        query_context(variant), family, namespace, labels, resolution_policy
    )


def query_compiled(context, family, namespace, labels, resolution_policy):
    """Internal owned source data validated at the public boundary."""
    inv = context["inventory"]
    data = context["variant"]
    if family == "bond_increments":
        require(
            namespace == "cff91_auto" and len(labels) == 2,
            "Unsupported increment namespace/arity",
        )
        forward = increment_query(inv, labels)
        backward = increment_query(inv, labels[::-1])
        reversal = forward["status"] == backward["status"] and (
            forward["status"] != "assigned"
            or forward["normalized_values"] == backward["normalized_values"][::-1]
        )
    else:
        indexed = context["catalog"]
        require(
            namespace in ("cff91", "cff91_auto"),
            "Unsupported lookup namespace; raw inventory still preserved",
        )
        require(
            family in FAMILIES if namespace == "cff91" else family in AUTO_ROLES,
            "Unsupported query family",
        )
        arity = FAMILIES[family][0] if namespace == "cff91" else len(AUTO_ROLES[family])
        require(len(labels) == arity, "Wrong source-query arity")
        permutation = orientations(family, len(labels))[-1]
        reversed_labels = [labels[k] for k in permutation]

        def run(types):
            return (
                resolve(
                    family, types, indexed, inv, policy=resolution_policy, trace=True
                )
                if namespace == "cff91"
                else lookup(
                    family,
                    types,
                    namespace,
                    indexed,
                    records(inv, "auto_equivalence"),
                    policy=resolution_policy,
                    trace=True,
                )
            )

        forward, backward = run(labels), run(reversed_labels)
        inverse_values = backward.get("normalized_values")
        if inverse_values is not None and permutation == tuple(
            reversed(range(len(labels)))
        ):
            inverse_values = reverse_values(family, inverse_values)
        reversal = forward["status"] == backward["status"] and (
            forward["status"] != "assigned"
            or forward["normalized_values"] == inverse_values
        )
    ambiguous = (
        forward["status"] == "ambiguous"
        or backward["status"] == "ambiguous"
        or not reversal
    )
    state = (
        "wildcard_ambiguous"
        if ambiguous
        else "source_row_present"
        if forward["status"] == "assigned"
        else "source_row_absent"
    )
    nonzero_wilson = (
        family == "wilson_out_of_plane"
        and forward["status"] == "assigned"
        and forward["normalized_values"][1] != 0
    )
    if nonzero_wilson:
        state = "model_incomplete"
    if family == "bond_increments" and ambiguous:
        state = "native_charge_incomplete"
    return {
        "source": data["selection"],
        "variant_identity": context["identity"],
        "classification": state,
        "forward": forward,
        "reverse": backward,
        "physical_reversal_invariant": forward["status"] == "assigned"
        and reversal
        and not ambiguous
        and not nonzero_wilson,
        "row_present": forward["status"] == "assigned",
        "interpretation_diagnostic": "Nonzero Wilson equilibrium signed-permutation semantics unresolved"
        if nonzero_wilson
        else None,
        "coefficient_convention": "Values transformed into original physical roles before exact comparison",
        "native_assignment_authorized": data["selection"]["sha256"] == PIN["sha256"],
        "query_only": True,
    }


def resolution_data(variants, policy_data, request):
    require(
        type(request) is dict
        and set(request) == {"family", "namespace", "types", "sites"},
        "Exact source query required",
    )
    labels = request["types"]
    require(
        type(labels) is list and labels and all(type(t) is str and t for t in labels),
        "Exact type labels required",
    )
    require(
        request["sites"] is None
        or (
            type(request["sites"]) is list
            and len(request["sites"]) == len(labels)
            and all(type(s) is int for s in request["sites"])
            and len(set(request["sites"])) == len(labels)
        ),
        "Exact ordered stable-site coverage required",
    )
    require(
        type(request["family"]) is str and type(request["namespace"]) is str,
        "Explicit family/namespace required",
    )
    for selection in [policy_data["primary"]] + policy_data["fallbacks"]:
        chosen_variant(variants, selection)
    decisions = []
    selected = None
    for selection in [policy_data["primary"]] + policy_data["fallbacks"]:
        variant, mirrors = chosen_variant(variants, selection)
        q = query_one(
            variant,
            request["family"],
            request["namespace"],
            labels,
            policy_data["resolution_policy"],
        )
        q["mirror_identities"] = mirrors
        decisions.append(q)
        if q["classification"] == "source_row_present":
            selected = selection
            break
        if q["classification"] != "source_row_absent":
            break
    outcome = decisions[-1]["classification"]
    if selected and selected != policy_data["primary"]:
        outcome = "source_variant_only"
    return {
        "schema": RESOLUTION_SCHEMA,
        "policy": policy_data,
        "request": request,
        "available_variant_identities": sorted(v.identity for v in variants),
        "attempts": decisions,
        "selected_source": selected,
        "classification": outcome,
        "mixed_source_model": False,
        "model_assembly_authorized": False,
        **FLAGS,
    }


@dataclass(frozen=True)
class PCFFVariantResolution:
    json_text: str
    variants: tuple

    @boundary
    def validate_integrity(self):
        p = unpack(self.json_text)
        require(
            p["schema"] == RESOLUTION_SCHEMA, "Unsupported variant-resolution schema"
        )
        policy = p["policy"]
        expected_policy = PCFFVariantPolicy(
            PCFFVariantSelection(**policy["primary"]),
            tuple(PCFFVariantSelection(**f) for f in policy["fallbacks"]),
            policy["resolution_policy"],
        ).payload()
        require(pack(policy) == pack(expected_policy), "Contradictory source policy")
        selections(self.variants)
        require(
            pack(p) == pack(resolution_data(self.variants, policy, p["request"])),
            "Contradictory variant resolution",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def resolve_pcff_variant_record(
    variants, *, policy, family, types, namespace="cff91", sites=None
):
    require(type(policy) is PCFFVariantPolicy, "Explicit PCFF variant policy required")
    require(type(types) in (list, tuple), "Exact type sequence required")
    require(
        sites is None or type(sites) in (list, tuple), "Exact site sequence required"
    )
    variants = tuple(variants)
    selections(variants)
    request = {
        "family": family,
        "namespace": namespace,
        "types": list(types),
        "sites": None if sites is None else list(sites),
    }
    result = PCFFVariantResolution(
        pack(resolution_data(variants, policy.payload(), request)), variants
    )
    result.validate_integrity()
    return result


@boundary
def save_pcff_variant_resolution(result, path):
    require(type(result) is PCFFVariantResolution, "Expected variant resolution")
    result.validate_integrity()
    publish(Path(path), result.json_text.encode())


@boundary
def load_pcff_variant_resolution(path, *, variants):
    result = PCFFVariantResolution(Path(path).read_text(), tuple(variants))
    result.validate_integrity()
    return result
