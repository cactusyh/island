"""Additive source adjudication, independent of MM prediction and fitting.

A comparison contract describes a scalar row within its own source model. It
never supplies a missing PCFF typing or equilibrium-dependency contract. Unknown
metadata blocks numerical comparison; a digest is not source authority.
"""

from decimal import Decimal, localcontext

from island.pcff_completion import digest, require, signed

CATEGORIES = (
    "same_row_numeric_conflict",
    "unit_conflict",
    "equation_conflict",
    "phase_conflict",
    "prefactor_conflict",
    "resolution_conflict",
    "chemical_environment_conflict",
    "forcefield_family_conflict",
    "equilibrium_dependency_conflict",
    "source_authority_conflict",
    "insufficient_evidence",
)
AUTHORITY = {
    1: "exact native PCFF/FRC with complete convention",
    2: "exact author-primary artifact for the same contract",
    3: "exact published table for the same family, resolution and equation",
    4: "explicitly compatible non-native research contract",
    5: "analog or transferred value, unqualified for native use",
    6: "unexecuted QM-fit target",
}
# Fixed comparison precision only. Degree conversions carry this approximation;
# kJ/kcal and nm/angstrom factors are exact decimal/rational definitions.
PI = Decimal("3.1415926535897932384626433832795028841971693993751")


def normalized_cell(value, unit):
    """Known dimensional transformations only, including angular powers."""
    factor, target = Decimal(1), unit
    with localcontext() as ctx:
        ctx.prec = 48
        if unit == "degree":
            factor, target = PI / 180, "radian"
        elif unit == "nm":
            factor, target = Decimal(10), "angstrom"
        elif unit in ("kJ/mol", "kJ/(mol*radian^2)"):
            factor, target = Decimal(1) / Decimal("4.184"), unit.replace("kJ", "kcal")
        elif unit == "kcal/(mol*degree^2)":
            factor, target = (180 / PI) ** 2, "kcal/(mol*radian^2)"
        elif unit not in {
            "radian",
            "angstrom",
            "angstrom^-1",
            "dimensionless",
            "elementary_charge",
            "kcal/mol",
            "kcal*angstrom^6/mol",
            *{f"kcal/(mol*angstrom^{n})" for n in (2, 3, 4)},
            *{f"kcal/(mol*radian^{n})" for n in (2, 3, 4)},
            "kcal/(mol*angstrom)",
            "kcal/(mol*radian)",
            "kcal/(mol*angstrom*radian)",
        }:
            return {"value": None, "unit": unit, "factor": None, "resolved": False}
        return {
            "value": None if value is None else str(Decimal(value) * factor),
            "unit": target,
            "factor": str(factor),
            "resolved": value is not None,
        }


# Complete scalar comparison metadata. Charge/dispersion assumptions are also
# compared: equal bonded numbers do not establish interchangeability of models.
CONTRACT_KEYS = {
    "interaction_family": "equation_conflict",
    "ordered_sites": "chemical_environment_conflict",
    "chemical_environment": "chemical_environment_conflict",
    "resolution": "resolution_conflict",
    "forcefield_family": "forcefield_family_conflict",
    "energy_equation": "equation_conflict",
    "force_equation": "equation_conflict",
    "phase": "phase_conflict",
    "prefactor": "prefactor_conflict",
    "multiplicity": "phase_conflict",
    "orientation": "phase_conflict",
    "equilibrium_dependencies": "equilibrium_dependency_conflict",
    "charge_nonbonded_contract": "forcefield_family_conflict",
    "source_target": "source_authority_conflict",
}


def compare(left, right):
    """Classify structurally first; never subtract incompatible parameter values."""
    categories, missing = set(), []
    for key, category in CONTRACT_KEYS.items():
        a, b = left["contract"].get(key), right["contract"].get(key)
        if a is None or b is None:
            missing.append(key)
        elif a != b:
            categories.add(category)
    if left["fields"] != right["fields"]:
        categories.add("equation_conflict")
    a, b = left["comparison_values"], right["comparison_values"]
    if len(a) != len(b) or [x["unit"] for x in a] != [x["unit"] for x in b]:
        categories.add("unit_conflict")
    if any(not x["resolved"] for x in a + b):
        missing.append("value_or_unit")
    if missing:
        categories.add("insufficient_evidence")
    numeric = not categories
    delta, equal = None, None
    if numeric:
        with localcontext() as ctx:
            ctx.prec = 40
            delta = [
                str(Decimal(x["value"]) - Decimal(y["value"]))
                for x, y in zip(a, b, strict=True)
            ]
        # This allowance is solely for the declared pi approximation/roundoff,
        # never for rounding away disagreement in printed source decimals.
        equal = all(abs(Decimal(d)) <= Decimal("1e-36") for d in delta)
        if not equal:
            categories.add("same_row_numeric_conflict")
    return signed(
        {
            "schema": "island.source-comparison.v1",
            "alternatives": [left["record_id"], right["record_id"]],
            "comparison_identities": [left["identity"], right["identity"]],
            "categories": sorted(categories),
            "missing_metadata": sorted(set(missing)),
            "numerically_compared": numeric,
            "normalized_difference": delta,
            "equal_with_declared_arithmetic_allowance": equal,
            "authority_selected": None,
            "values_averaged": False,
            "native_PCFF_resolution": False,
        }
    )


def equation_metadata(row):
    required = (
        "energy_equation",
        "force_equation",
        "units",
        "phase",
        "prefactor",
        "multiplicity",
        "orientation",
        "equilibrium_dependencies",
    )
    metadata = row.get("equation_metadata", {})
    require(
        all(metadata.get(k) is not None and metadata[k] != "" for k in required),
        "complete_equation_metadata_required",
    )
    require(row.get("source_row_id") is not None, "exact_source_row_required")
    return True


def native_authority(row):
    """Hierarchy is necessary, never sufficient: require a complete exact contract."""
    require(
        row.get("selection_basis") == "source_evidence",
        "residual_or_MD_selection_forbidden",
    )
    require(type(row.get("authority_level")) is int, "authority_level_type")
    require(row["authority_level"] in AUTHORITY, "authority_level_unknown")
    require(row["authority_level"] <= 3, "analog_or_external_native_promotion")
    require(row.get("source_family") == "PCFF", "cross_family_native_overlay")
    require(row.get("resolution") == "all_atom", "resolution_native_overlay")
    require(
        row.get("exact_environment") is True
        and row.get("same_source_contract") is True,
        "analog_environment_not_exact",
    )
    require(
        row.get("complete_equation_conventions") is True,
        "incomplete_equation_conventions",
    )
    require(row.get("unresolved_conflicts") == [], "unresolved_source_conflict")
    require(
        row.get("zero_fill") is False and row.get("wildcard") is False,
        "zero_fill_wildcard_forbidden",
    )
    require(
        bool(row.get("source_sha256")) and bool(row.get("source_location")),
        "source_evidence_required",
    )
    equation_metadata(row)
    return True


def profile_gate(rows, targets, *, strict_conversion=False, static_closure=False):
    """Check exact identities and all gates; a count of 378 alone is insufficient."""
    require(
        len(targets) == 378 and len({t["id"] for t in targets}) == 378,
        "exact_378_target_inventory",
    )
    require(
        len(rows) == 378
        and {r["target_id"] for r in rows} == {t["id"] for t in targets},
        "partial_profile_forbidden",
    )
    lookup = {t["id"]: t for t in targets}
    contracts, families = set(), set()
    for row in rows:
        t = lookup[row["target_id"]]
        require(
            row.get("target_identity") == t["identity"]
            and row.get("stable_sites") == t["stable_sites"]
            and row.get("ordered_types") == t["ordered_types"],
            "exact_row_mapping_required",
        )
        require(
            row.get("exact") is True and row.get("complete") is True,
            "partial_convention_profile_forbidden",
        )
        require(
            row.get("zero_fill") is False and row.get("wildcard") is False,
            "zero_fill_wildcard_forbidden",
        )
        require(
            type(row.get("authority_level")) is int
            and row["authority_level"] in (1, 2, 3, 4),
            "analog_source_profile_forbidden",
        )
        require(
            row.get("selection_basis") == "source_evidence",
            "residual_or_MD_selection_forbidden",
        )
        require(
            row.get("resolution") == "all_atom"
            and row.get("unresolved_conflicts") == [],
            "profile_contract_unresolved",
        )
        require(
            bool(row.get("source_sha256"))
            and bool(row.get("source_location"))
            and bool(row.get("contract_id")),
            "profile_source_evidence_required",
        )
        require(
            row.get("exact_environment") is True
            and row.get("same_source_contract") is True
            and row.get("complete_equation_conventions") is True,
            "profile_exact_source_contract_required",
        )
        if row["authority_level"] <= 3:
            native_authority(row)
        equation_metadata(row)
        contracts.add(row["contract_id"])
        families.add(row.get("source_family"))
    require(
        len(contracts) == 1 and families == {"PCFF"}, "cross_family_profile_overlay"
    )
    require(
        strict_conversion is True and static_closure is True,
        "conversion_and_static_closure_required",
    )
    return True


def verify_signed(value):
    require(
        value.get("identity")
        == digest({k: v for k, v in value.items() if k != "identity"}),
        "adjudication_outer_identity",
    )


DEVIATION_SHA = "b1de0425346d6fbfe67f5eb66f2453750269da3cbbfbed6db705972627fb46dd"


def execution_deviation(root):
    """Reconstruct actual out-of-scope executions without rerunning them."""
    import hashlib
    import json
    import tarfile

    report_path = root / "scope-deviation/execution-deviation.json"
    require(
        hashlib.sha256(report_path.read_bytes()).hexdigest() == DEVIATION_SHA,
        "J72_deviation_report_mutation",
    )
    report = json.loads(report_path.read_text())
    archive = root / "scope-deviation" / report["archive"]
    require(
        hashlib.sha256(archive.read_bytes()).hexdigest() == report["archive_sha256"],
        "J72_deviation_archive_mutation",
    )
    raw = {}
    with tarfile.open(archive, "r:xz") as tar:
        for member in tar:
            require(
                member.isfile() and member.name in report["files"],
                "J72_quarantine_inventory",
            )
            data = tar.extractfile(member).read()
            require(
                hashlib.sha256(data).hexdigest()
                == report["files"][member.name]["sha256"]
                and len(data) == report["files"][member.name]["bytes"],
                "J72_quarantine_bytes",
            )
            raw[member.name] = data
    require(set(raw) == set(report["files"]), "J72_quarantine_missing_file")
    qm = [name for name in raw if name.endswith("/sqm.out")]
    require(
        len(qm) == 4 and all(b"Calculation Completed" in raw[n] for n in qm),
        "J72_actual_QM_count",
    )
    records = []
    for i in (1, 2, 3, 5, 6, 7):
        prefix = f"test_real_final_graph_source_c{i}/acceptance/"
        frames = {}
        for name, data in raw.items():
            if (
                name.startswith(prefix)
                and "/generation-" in name
                and name.endswith("/record.json")
            ):
                branch = name[len(prefix) :].split("/")[0]
                payload = json.loads(data)["payload"]
                frames.setdefault(branch, []).extend(
                    f["step"] for f in payload["frames"]
                )
        require(
            max(frames["full"]) == 4 and max(frames["split"]) == 2,
            "J72_actual_MD_frames",
        )
        resumed = max(frames["relocated"]) - 2
        require(resumed == (0 if i == 7 else 2), "J72_actual_resume_frames")
        records.append(
            {
                "case": i,
                "recorded_steps": 4 + 2 + resumed,
                "uncommitted_resume_steps_possible": 2 if i == 7 else 0,
            }
        )
    minimum = sum(r["recorded_steps"] for r in records)
    maximum = minimum + sum(r["uncommitted_resume_steps_possible"] for r in records)
    require(
        minimum == 46
        and maximum == 48
        and report["MD_steps_verified_minimum"] == minimum
        and report["MD_steps_possible_maximum"] == maximum,
        "J72_actual_MD_count",
    )
    return signed(
        {
            "schema": "island.j72-actual-execution-scope.v1",
            "zero_execution_scope_met": False,
            "authorized_new_QM": 0,
            "authorized_MD_steps": 0,
            "new_QM": len(qm),
            "new_HMDS_QM": 0,
            "QM_fits": 0,
            "MD_steps": None,
            "MD_steps_verified_minimum": minimum,
            "MD_steps_possible_maximum": maximum,
            "cases": records,
            "native_attempt_records": qm,
            "archive_sha256": report["archive_sha256"],
            "use_for_source_adjudication": False,
            "reserved_validation_activated": False,
        }
    )
