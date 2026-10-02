"""Data-only fragment identities, cap bookkeeping and external assignments."""

from dataclasses import dataclass
from functools import wraps
from math import fsum, isfinite
from pathlib import Path

from island.charge_references.records import pack, unpack
from island.exceptions import FragmentChargeError
from island.workflows import storage
from island.workflows.bundle import system_data, system_from

FLAGS = {"production_validated": False, "simulation_readiness": "not_established"}
CAP_POLICY = "two_hydrogen_caps_v1"
H_POLICY = "verified_parent_group_equal_v1"
TRANSFER = "removed_cap_to_bonded_parent_v1"


def need(ok, message):
    if not ok:
        raise FragmentChargeError(message)


def boundary(fn):
    @wraps(fn)
    def call(*args, **kwargs):
        try:
            return fn(*args, **kwargs)
        except FragmentChargeError:
            raise
        except Exception as error:
            raise FragmentChargeError(f"{fn.__name__}: {error}") from error

    return call


def identity(value):
    return storage.checksum(storage.json_bytes(storage.encode(value)))


def number(x):
    return type(x) in (int, float) and isfinite(x)


def atom_identity(atom):
    need(
        atom.formal_charge == 0
        and not atom.metadata.get("isotope")
        and not atom.metadata.get("radical_electrons", 0)
        and not atom.metadata.get("cip_label")
        and atom.metadata.get("chiral_tag", "CHI_UNSPECIFIED") == "CHI_UNSPECIFIED",
        "Charged, isotopic, radical or assigned stereochemical sites unsupported",
    )
    need(number(atom.mass) and atom.mass > 0, "Invalid mass")
    return (
        atom.element,
        atom.atomic_number,
        atom.formal_charge,
        atom.mass,
        atom.metadata.get("aromatic", False),
    )


def edges(system):
    return {
        tuple(sorted((b.site1, b.site2))): (b.order, b.aromatic)
        for b in system.topology.bonds.values()
    }


def fragment_data(p):
    need(
        set(p)
        == {
            "schema",
            "system",
            "definition",
            "seed",
            "construction_version",
            "mapping",
            "data",
        },
        "Invalid capped fragment fields",
    )
    need(
        type(p["definition"]) is str
        and p["definition"]
        and type(p["seed"]) is int
        and 0 <= p["seed"] < 2**31,
        "Invalid definition/seed",
    )
    need(
        type(p["construction_version"]) is str and p["construction_version"],
        "Missing construction version",
    )
    system = system_from(p["system"])
    system.coordinates.validate(system.topology)
    need(
        len(system.topology.connected_components()) == 1 and system.number_of_sites > 0,
        "Require connected fragment",
    )
    need(
        system.metadata.get("fragment_psmiles") == p["definition"]
        and system.metadata.get("fragment_seed") == p["seed"],
        "Fragment construction provenance mismatch",
    )
    for site, atom in system.topology.sites.items():
        atom_identity(atom)
        xyz = system.coordinates.get(site)
        need(
            len(xyz) == 3 and all(isfinite(float(v)) for v in xyz),
            "Invalid fragment coordinates",
        )
        bonds = [
            b for b in system.topology.bonds.values() if site in (b.site1, b.site2)
        ]
        valence = fsum(b.order for b in bonds)
        allowed = {
            "H": {1},
            "C": {4},
            "N": {3},
            "O": {2},
            "F": {1},
            "Cl": {1},
            "Br": {1},
            "S": {2, 4, 6},
            "P": {3, 5},
        }
        need(atom.element in allowed, "Unsupported element in fragment contract")
        if atom.metadata.get("aromatic"):
            need(
                (atom.element == "C" and valence in {4, 4.5})
                or (atom.element == "N" and valence in {3, 4, 4.5})
                or (atom.element == "S" and valence == 3),
                "Invalid aromatic valence",
            )
        else:
            need(valence in allowed[atom.element], "Invalid explicit fragment valence")
    mapping = p["mapping"]
    need(set(mapping) == {"heavy", "native_h", "caps"}, "Invalid fragment mapping")
    heavy = mapping["heavy"]
    native = mapping["native_h"]
    caps = mapping["caps"]
    need(
        all(type(i) is int and type(s) is int for i, s in heavy.items())
        and len(set(heavy.values())) == len(heavy),
        "Invalid local-heavy mapping",
    )
    need(
        set(native) == set(heavy) and set(caps) == {"head", "tail"},
        "Incomplete parent/cap mapping",
    )
    covered = list(heavy.values())
    for i, s in heavy.items():
        need(system.topology.sites[s].element != "H", "Heavy mapping contains H")
        need(
            system.topology.sites[s].metadata.get("source_repeat_atom_index") == i,
            "Heavy source provenance mismatch",
        )
        need(type(native[i]) is list, "Invalid native H group")
        for h in native[i]:
            need(
                type(h) is int
                and system.topology.sites[h].element == "H"
                and set(system.topology.neighbors(h)) == {s},
                "Invalid native H parent",
            )
            covered.append(h)
    for role, cap in caps.items():
        need(
            set(cap)
            == {
                "site_id",
                "parent_local_index",
                "dummy_source_index",
                "atom_map_number",
            },
            "Invalid cap fields",
        )
        h = cap["site_id"]
        i = cap["parent_local_index"]
        need(
            type(h) is int
            and type(i) is int
            and i in heavy
            and system.topology.sites[h].element == "H"
            and set(system.topology.neighbors(h)) == {heavy[i]},
            "Invalid cap parent",
        )
        need(
            type(cap["dummy_source_index"]) is int
            and type(cap["atom_map_number"]) is int,
            "Invalid attachment indices",
        )
        need(
            system.topology.sites[h].metadata.get("fragment_cap_role") == role
            and system.topology.sites[h].metadata.get("dummy_source_index")
            == cap["dummy_source_index"],
            "Cap source provenance mismatch",
        )
        covered.append(h)
    need(
        caps["head"]["atom_map_number"] == 1
        and caps["tail"]["atom_map_number"] == 2
        and caps["head"]["dummy_source_index"] != caps["tail"]["dummy_source_index"],
        "Explicit head/tail orientation required",
    )
    need(
        len(covered) == len(set(covered))
        and set(covered) == set(system.topology.sites),
        "Fragment coverage mismatch",
    )
    for b in system.topology.bonds.values():
        if "H" in (
            system.topology.sites[b.site1].element,
            system.topology.sites[b.site2].element,
        ):
            need(b.order == 1 and not b.aromatic, "Hydrogen attachment must be single")
    return {
        "graph_identity": identity(
            {
                "atoms": {
                    s: atom_identity(a)
                    for s, a in sorted(system.topology.sites.items())
                },
                "bonds": sorted(edges(system).items()),
            }
        ),
        "coordinate_identity": identity(p["system"]["coordinates"]),
        "cap_policy": CAP_POLICY,
        "formal_charge": 0,
        "atom_count": len(covered),
        **FLAGS,
    }


@dataclass(frozen=True)
class CappedFragment:
    json_text: str

    @boundary
    def validate_integrity(self):
        p = unpack(self.json_text)
        need(
            p["schema"] == "island_capped_fragment_v1"
            and pack(p["data"]) == pack(fragment_data(p)),
            "Fragment contradicts source",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


def charge_vector(q, ids):
    need(
        isinstance(q, dict) and all(type(i) is int for i in q) and set(q) == set(ids),
        "Exact integer charge-site coverage required",
    )
    need(all(number(v) for v in q.values()), "Finite nonboolean charges required")
    return {i: float(q[i]) for i in sorted(q)}


def charge_data(p):
    need(
        set(p) == {"schema", "fragment", "raw_charges", "backend", "evidence", "data"},
        "Invalid charge-result fields",
    )
    frag = CappedFragment(pack(p["fragment"])).payload
    q = charge_vector(p["raw_charges"], system_from(frag["system"]).topology.sites)
    backend = p["backend"]
    if backend["method"] == "provided":
        need(
            set(backend) == {"method", "source", "compatibility"}
            and type(backend["source"]) is str
            and backend["source"]
            and backend["compatibility"] == "unverified"
            and p["evidence"] == {},
            "Provided charges do not establish force-field compatibility",
        )
    elif backend["method"] == "am1bcc":
        from .amber import validate_evidence

        validate_evidence(frag, q, backend, p["evidence"])
    else:
        raise FragmentChargeError("Unsupported charge backend; RESP not implemented")
    return {
        "fragment_identity": identity(frag),
        "raw_total_e": fsum(q.values()),
        "formal_charge_e": 0,
        "raw_residual_e": fsum(q.values()),
        "raw_charge_sha256": identity(q),
        "unit": "elementary_charge",
        **FLAGS,
    }


@dataclass(frozen=True)
class FragmentChargeResult:
    json_text: str

    @boundary
    def validate_integrity(self):
        p = unpack(self.json_text)
        need(
            p["schema"] == "island_fragment_charges_v1"
            and pack(p["data"]) == pack(charge_data(p)),
            "Charge result contradicts evidence",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


def template_data(p):
    need(
        set(p)
        == {
            "schema",
            "calculation",
            "policy",
            "charge_tolerance_e",
            "hydrogen_policy",
            "data",
        },
        "Invalid template fields",
    )
    result = FragmentChargeResult(pack(p["calculation"])).payload
    need(p["hydrogen_policy"] == H_POLICY, "Unsupported hydrogen policy")
    tolerance = p["charge_tolerance_e"]
    need(number(tolerance) and tolerance > 0, "Invalid charge tolerance")
    q = result["raw_charges"]
    total = fsum(q.values())
    policy = p["policy"]
    need(
        policy in {"strict", "uniform_fragment_l2_v1"},
        "Explicit conservation policy required",
    )
    if policy == "strict":
        need(abs(total) <= tolerance, "Fragment outside strict total-charge tolerance")
    offset = -total / len(q) if policy == "uniform_fragment_l2_v1" else 0.0
    transformed = {i: q[i] + offset for i in sorted(q)}
    if policy != "strict":
        need(abs(fsum(transformed.values())) <= 1e-12, "Projection conservation failed")
    return {
        "calculation_identity": identity(result),
        "transformed_charges": transformed,
        "uniform_offset_e": offset,
        "transformed_total_e": fsum(transformed.values()),
        "transfer_policy": TRANSFER,
        "cache_key": identity({k: v for k, v in p.items() if k != "data"}),
        "unit": "elementary_charge",
        **FLAGS,
    }


@dataclass(frozen=True)
class FragmentChargeTemplate:
    json_text: str

    @boundary
    def validate_integrity(self):
        p = unpack(self.json_text)
        need(
            p["schema"] == "island_fragment_template_v1"
            and pack(p["data"]) == pack(template_data(p)),
            "Template contradicts calculation/policy",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)

    @property
    def cache_key(self):
        return self.payload["data"]["cache_key"]


@boundary
def create_fragment_template(
    calculation,
    *,
    conservation_policy,
    charge_tolerance=0.002,
    hydrogen_policy=H_POLICY,
):
    p = {
        "schema": "island_fragment_template_v1",
        "calculation": calculation.payload,
        "policy": conservation_policy,
        "charge_tolerance_e": charge_tolerance,
        "hydrogen_policy": hydrogen_policy,
        "data": {},
    }
    p["data"] = template_data(p)
    return FragmentChargeTemplate(pack(p))


@boundary
def validate_force_field_compatibility(template, force_field):
    """Explicit GAFF-family/data binding, not scientific suitability certification."""
    backend = template.payload["calculation"]["backend"]
    need(
        isinstance(force_field, dict) and set(force_field) == {"family", "data_sha256"},
        "Require specific family and force-field data checksum",
    )
    need(
        backend["method"] == "am1bcc",
        "Provided charges do not establish force-field compatibility",
    )
    need(
        force_field["family"] in {"gaff", "gaff2"}
        and force_field == backend["force_field"],
        "Incompatible force-field family/version; no OPLS/PCFF fallback",
    )
    return (
        "matched_charge_backend_family_and_data; scientific_suitability_not_established"
    )


def assignment_data(p):
    need(
        set(p)
        == {
            "schema",
            "template",
            "target_system",
            "force_field",
            "charge_tolerance_e",
            "data",
        },
        "Invalid assignment fields",
    )
    template = FragmentChargeTemplate(pack(p["template"]))
    t = template.payload
    compatibility = "unassessed; no force-field compatibility claim"
    if p["force_field"] is not None:
        compatibility = validate_force_field_compatibility(template, p["force_field"])
    system = system_from(p["target_system"])
    if system.metadata.get("force_field") is not None:
        need(
            system.metadata["force_field"] == p["force_field"],
            "Known target force-field identity must be checked explicitly",
        )
    top = system.topology
    poly = system.metadata["polymer"]
    f = t["calculation"]["fragment"]
    frag = system_from(f["system"])
    m = f["mapping"]
    dp = poly["degree_of_polymerization"]
    need(
        type(dp) is int and dp >= 1 and 0 < system.number_of_sites <= 1000,
        "Target DP/site ceiling unsupported",
    )
    need(len(top.connected_components()) == 1, "Connected chain required")
    need(
        poly["source_psmiles"] == f["definition"]
        and poly["architecture"] == "linear"
        and poly["polymer_type"] == "homopolymer",
        "Incompatible repeat definition/architecture",
    )
    need(
        poly["number_of_repeat_units"] == dp
        and poly["number_of_inter_repeat_unit_bonds"] == dp - 1,
        "Repeat counts mismatch",
    )
    sequence = poly["sequence"]
    need(len(sequence) == dp and len(set(sequence)) == 1, "Copolymer unsupported")
    label = sequence[0]
    need(
        poly["repeat_unit_definitions"] == {label: f["definition"]}
        and poly["composition_counts"] == {label: dp},
        "Ambiguous repeat definitions",
    )
    need(
        not any(k in poly for k in ("tacticity", "stereochemical_sequence")),
        "Assigned stereochemical provenance unsupported",
    )
    heavy = {}
    hydrogens = {}
    for site, atom in top.sites.items():
        need(type(site) is int and site == atom.id, "Invalid stable ID")
        atom_identity(atom)
        meta = atom.metadata
        r = meta["repeat_unit_index"]
        need(
            type(r) is int
            and 0 <= r < dp
            and meta["chain_id"] == poly["chain_id"]
            and meta["repeat_unit_type"] == label,
            "Invalid repeat provenance",
        )
        if atom.element == "H":
            ns = list(top.neighbors(site))
            need(len(ns) == 1, "Hydrogen has ambiguous parent")
            hydrogens.setdefault(ns[0], []).append(site)
        else:
            i = meta["source_repeat_atom_index"]
            need(
                type(i) is int and i in m["heavy"] and (r, i) not in heavy,
                "Missing/ambiguous local heavy identity",
            )
            need(
                atom_identity(atom)
                == atom_identity(frag.topology.sites[m["heavy"][i]]),
                "Repeat atomic identity mismatch",
            )
            heavy[r, i] = site
    need(
        set(heavy) == {(r, i) for r in range(dp) for i in m["heavy"]},
        "Heavy inventory mismatch",
    )
    need(
        poly["head_site_id"] == heavy[0, m["caps"]["head"]["parent_local_index"]]
        and poly["tail_site_id"]
        == heavy[dp - 1, m["caps"]["tail"]["parent_local_index"]],
        "Attachment orientation mismatch",
    )
    expected = {}
    inverse = {s: i for i, s in m["heavy"].items()}
    for r in range(dp):
        for (a, b), kind in edges(frag).items():
            if a in inverse and b in inverse:
                expected[
                    tuple(sorted((heavy[r, inverse[a]], heavy[r, inverse[b]])))
                ] = kind
    for r in range(dp - 1):
        expected[
            tuple(
                sorted(
                    (
                        heavy[r, m["caps"]["tail"]["parent_local_index"]],
                        heavy[r + 1, m["caps"]["head"]["parent_local_index"]],
                    )
                )
            )
        ] = (1.0, False)
    q = t["data"]["transformed_charges"]
    assigned = {}
    groups = []
    transfers = []
    for (r, i), site in sorted(heavy.items()):
        value = q[m["heavy"][i]]
        hvalues = [q[h] for h in m["native_h"][i]]
        retained = []
        removed = []
        for role, cap in m["caps"].items():
            if cap["parent_local_index"] != i:
                continue
            if (role == "head" and r == 0) or (role == "tail" and r == dp - 1):
                hvalues.append(q[cap["site_id"]])
                retained.append(role)
            else:
                value += q[cap["site_id"]]
                removed.append(role)
                transfers.append(
                    {
                        "repeat": r,
                        "role": role,
                        "parent_site_id": site,
                        "charge_e": q[cap["site_id"]],
                    }
                )
        hs = sorted(hydrogens.get(site, []))
        need(len(hs) == len(hvalues), "Native/terminal hydrogen count mismatch")
        assigned[site] = value
        for h in hs:
            need(
                top.sites[h].metadata["repeat_unit_index"] == r
                and atom_identity(top.sites[h])
                == atom_identity(frag.topology.sites[m["caps"]["head"]["site_id"]]),
                "Hydrogen identity/repeat mismatch",
            )
            expected[tuple(sorted((site, h)))] = (1.0, False)
            assigned[h] = fsum(hvalues) / len(hs)
        groups.append(
            {
                "repeat": r,
                "parent_site_id": site,
                "native_count": len(m["native_h"][i]),
                "retained_caps": retained,
                "removed_caps": removed,
                "target_hydrogen_ids": hs,
                "group_total_e": fsum(hvalues),
            }
        )
    need(
        edges(system) == expected and set(assigned) == set(top.sites),
        "Target graph/coverage differs from oriented fragment assembly",
    )
    assigned = charge_vector(assigned, top.sites)
    total = fsum(assigned.values())
    expected_total = dp * t["data"]["transformed_total_e"]
    need(abs(total - expected_total) <= 1e-12, "Cap transfer failed conservation")
    tol = p["charge_tolerance_e"]
    need(
        number(tol) and tol > 0 and abs(total) <= tol,
        "Assembled molecular charge exceeds declared tolerance; no repair",
    )
    return {
        "template_identity": template.identity,
        "template_cache_key": template.cache_key,
        "assignments": assigned,
        "total_e": total,
        "expected_fragment_sum_e": expected_total,
        "groups": groups,
        "removed_cap_transfers": transfers,
        "charge_method": t["calculation"]["backend"]["method"],
        "conservation_policy": t["policy"],
        "cap_policy": CAP_POLICY,
        "transfer_policy": TRANSFER,
        "hydrogen_policy": H_POLICY,
        "compatibility": compatibility,
        "force_field": p["force_field"],
        "unit": "elementary_charge",
        "target_graph_identity": identity(
            {"heavy_mapping": sorted(heavy.items()), "bonds": sorted(expected.items())}
        ),
        **FLAGS,
    }


@dataclass(frozen=True)
class FragmentChargeAssignment:
    json_text: str

    @boundary
    def validate_integrity(self):
        p = unpack(self.json_text)
        need(
            p["schema"] == "island_fragment_assignment_v1"
            and pack(p["data"]) == pack(assignment_data(p)),
            "Assignment contradicts source/template",
        )

    @property
    def payload(self):
        self.validate_integrity()
        return unpack(self.json_text)

    @property
    def identity(self):
        return identity(self.payload)


@boundary
def assign_fragment_charges(
    template, system, *, force_field=None, charge_tolerance=1e-12
):
    p = {
        "schema": "island_fragment_assignment_v1",
        "template": template.payload,
        "target_system": system_data(system),
        "force_field": force_field,
        "charge_tolerance_e": charge_tolerance,
        "data": {},
    }
    p["data"] = assignment_data(p)
    return FragmentChargeAssignment(pack(p))


@boundary
def save_fragment_record(record, path):
    need(
        isinstance(
            record,
            (
                CappedFragment,
                FragmentChargeResult,
                FragmentChargeTemplate,
                FragmentChargeAssignment,
            ),
        ),
        "Expected fragment record",
    )
    record.validate_integrity()
    storage.publish(Path(path), record.json_text.encode())


@boundary
def load_fragment_record(path, *, expected_cache_key=None):
    raw = Path(path).read_text()
    p = unpack(raw)
    kind = {
        "island_capped_fragment_v1": CappedFragment,
        "island_fragment_charges_v1": FragmentChargeResult,
        "island_fragment_template_v1": FragmentChargeTemplate,
        "island_fragment_assignment_v1": FragmentChargeAssignment,
    }[p["schema"]]
    result = kind(raw)
    result.validate_integrity()
    if expected_cache_key is not None:
        need(
            isinstance(result, FragmentChargeTemplate)
            and result.cache_key == expected_cache_key,
            "Cache identity/policy mismatch",
        )
    return result
