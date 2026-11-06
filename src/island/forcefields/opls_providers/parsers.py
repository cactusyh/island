"""Bounded data-only parsers. No includes, macros, wildcards or generated values."""

from island.forcefields.assignment_core.schema import identity, require
from island.forcefields.bundle_adapters.schema import keys, loads, text

from .schema import FORMS, number, row_key


def pointer(file, location):
    return {"path": file["path"], "sha256": file["sha256"], "pointer": location}


def canonical_json(raw, file, graph_identity, *, database=False):
    src = loads(raw)
    if database:
        keys(src, ("schema", "entries"), "database_snapshot_schema")
        require(
            src["schema"] == "island.oplsaa-db-snapshot.v1"
            and isinstance(src["entries"], list),
            "database_snapshot_schema",
        )
        hits = [
            (i, e)
            for i, e in enumerate(src["entries"])
            if isinstance(e, dict) and e.get("graph_identity") == graph_identity
        ]
        require(len(hits) == 1, "database_exact_graph_missing_or_ambiguous")
        index, src = hits[0]
        prefix = f"/entries/{index}"
    else:
        prefix = ""
    keys(
        src,
        (
            "schema",
            "family",
            "graph_identity",
            "forms",
            "atoms",
            "records",
            "pairs14",
            "exclusions",
        ),
        "source_exchange_schema",
    )
    require(
        src["schema"] == "island.opls-provider-records.v1"
        and src["family"] == "OPLS-AA",
        "OPLS_only_source_exchange",
    )
    require(
        src["graph_identity"] == graph_identity, "source_exact_graph_identity_mismatch"
    )
    require(
        isinstance(src["forms"], dict) and src["forms"], "structured_forms_required"
    )
    for name, form in src["forms"].items():
        require(
            name in FORMS and form == FORMS[name], "equation_unit_ordering_mismatch"
        )
    require(
        all(
            isinstance(src[k], list)
            for k in ("atoms", "records", "pairs14", "exclusions")
        ),
        "source_arrays_required",
    )
    atoms, records = [], []
    for i, atom in enumerate(src["atoms"]):
        keys(
            atom,
            ("site", "atom_type", "element", "mass", "location"),
            "source_atom_schema",
        )
        text(atom["location"], "exact_source_location_required")
        atoms.append(
            {k: v for k, v in atom.items() if k != "location"}
            | {"source": pointer(file, prefix + f"/atoms/{i}; " + atom["location"])}
        )
    for i, row in enumerate(src["records"]):
        keys(
            row,
            (
                "id",
                "family",
                "sites",
                "form",
                "values",
                "applicability",
                "reason",
                "location",
            ),
            "source_parameter_schema",
        )
        require(
            row["form"] is None or row["form"] in src["forms"],
            "source_form_declaration_missing",
        )
        text(row["location"], "exact_source_location_required")
        records.append(
            {k: v for k, v in row.items() if k != "location"}
            | {"source": pointer(file, prefix + f"/records/{i}; " + row["location"])}
        )
    return atoms, records, src["pairs14"], src["exclusions"]


def gromacs(raw, file, graph_identity, conventions):
    """PolyParGen/LigParGen self-contained .itp/.top subset, with source lines.

    Explicit numeric interaction rows only. Unsupported sections/functions block.
    An absent improper is *not* inferred to be zero or not applicable.
    """
    del graph_identity  # Identity is bound independently to the pinned graph file.
    section, sections = None, {}
    allowed = {
        "defaults",
        "atomtypes",
        "moleculetype",
        "atoms",
        "bonds",
        "angles",
        "dihedrals",
        "pairs",
        "exclusions",
    }
    for lineno, line in enumerate(raw.decode("utf-8").splitlines(), 1):
        stripped = line.split(";", 1)[0].strip()
        if not stripped:
            continue
        require(
            not stripped.startswith("#"),
            "preprocessor_include_or_macro_requires_explicit_flattened_source",
        )
        if stripped.startswith("["):
            require(stripped.endswith("]"), "malformed_section")
            section = stripped[1:-1].strip()
            require(section in allowed, "unsupported_gromacs_section:" + section)
            require(section not in sections, "duplicate_gromacs_section:" + section)
            sections[section] = []
        else:
            require(section is not None, "row_before_section")
            sections[section].append((lineno, stripped.split()))
    defaults = sections.get("defaults", [])
    require(len(defaults) == 1, "explicit_gromacs_defaults_required")
    _, d = defaults[0]
    require(
        len(d) == 5 and d[:3] == ["1", "3", "yes"],
        "unsupported_nonbonded_mixing_or_pair_generation",
    )
    require(
        number(d[3]) == number(conventions["scale14_lj"])
        and number(d[4]) == number(conventions["scale14_coulomb"]),
        "source_14_policy_mismatch",
    )
    molecules = sections.get("moleculetype", [])
    require(
        len(molecules) == 1 and len(molecules[0][1]) == 2 and molecules[0][1][1] == "3",
        "single_molecule_nrexcl3_required",
    )
    types, atoms, records = {}, [], []

    def add(family, sites, form, values, lineno, source_section=None):
        rid = row_key(family, sites)
        hit = next((r for r in records if r["id"] == rid), None)
        if hit is not None:
            require(
                family == "torsion" and form == hit["form"] == "torsion_periodic",
                "duplicate_or_conflicting_source_row:" + rid,
            )
            hit["values"].append(values)
            hit["source"]["pointer"] += f",line:{lineno}"
        else:
            records.append(
                {
                    "id": rid,
                    "family": family,
                    "sites": sites,
                    "form": form,
                    "values": [values],
                    "applicability": "required",
                    "reason": "",
                    "source": pointer(
                        file, f"line:{lineno};section:{source_section or section}"
                    ),
                }
            )

    section = "atomtypes"
    for lineno, fields in sections.get(section, []):
        # GROMACS permits optional bond_type / atomic_number columns. Last six
        # fields retain mass, charge, ptype, sigma, epsilon with atomic number
        # optional; require the unambiguous last-five layout and numeric mass.
        require(len(fields) in (6, 7, 8), "unsupported_atomtype_layout")
        name, mass, charge, ptype, sigma, epsilon = fields[0], *fields[-5:]
        require(ptype == "A" and name not in types, "duplicate_or_nonatom_type")
        number(mass)
        number(charge)
        number(sigma)
        number(epsilon)
        types[name] = (mass, sigma, epsilon, lineno)
    section = "atoms"
    for lineno, fields in sections.get(section, []):
        require(len(fields) == 8, "explicit_charge_mass_single_state_required")
        site, atom_type = int(fields[0]), fields[1]
        require(atom_type in types, "undefined_atom_type")
        mass, sigma, epsilon, type_line = types[atom_type]
        require(number(fields[7]) == number(mass), "source_mass_disagreement")
        # Element is supplied by an independently pinned exact graph; mass and
        # connectivity are checked later. Never infer chemical typing by mass.
        atoms.append(
            {
                "site": site,
                "atom_type": atom_type,
                "element": None,
                "mass": fields[7],
                "source": pointer(file, f"line:{lineno};section:atoms"),
            }
        )
        add("charge", [site], "charge", [fields[6]], lineno)
        add("nonbonded", [site], "lj12_6", [sigma, epsilon], type_line, "atomtypes")
    for section, family, arity, functions in (
        ("bonds", "bond", 2, {"1": ("bond_harmonic", 2)}),
        ("angles", "angle", 3, {"1": ("angle_harmonic", 2)}),
        (
            "dihedrals",
            "torsion",
            4,
            {
                "1": ("torsion_periodic", 3),
                "9": ("torsion_periodic", 3),
                "3": ("torsion_rb", 6),
                "5": ("torsion_opls", 4),
                "2": ("improper_harmonic", 2),
            },
        ),
    ):
        for lineno, fields in sections.get(section, []):
            require(
                len(fields) > arity and fields[arity] in functions,
                "unsupported_or_implicit_interaction",
            )
            form, count = functions[fields[arity]]
            require(
                len(fields) == arity + 1 + count, "missing_explicit_interaction_fields"
            )
            actual_family = "improper" if form.startswith("improper") else family
            add(
                actual_family,
                list(map(int, fields[:arity])),
                form,
                fields[arity + 1 :],
                lineno,
            )
    pairs = []
    for _, fields in sections.get("pairs", []):
        require(
            len(fields) == 3 and fields[2] == "1", "unsupported_explicit_pair_override"
        )
        pairs.append(list(map(int, fields[:2])))
    exclusions = []
    for _, fields in sections.get("exclusions", []):
        require(len(fields) >= 2, "invalid_explicit_exclusion")
        exclusions.extend([[int(fields[0]), int(other)] for other in fields[1:]])
    return atoms, records, pairs, exclusions


def claim_metadata(row):
    """A hash binds all raw values without serializing them into public reports."""
    return {
        "record_id": row["id"],
        "family": row["family"],
        "sites": row["sites"],
        "form": row["form"],
        "applicability": row["applicability"],
        "source": row["source"],
        "claim_sha256": identity(row),
    }
