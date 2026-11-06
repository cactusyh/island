"""Local, non-production parameter fragments; pure serialization, never engines.

The embedded exchange payload preserves source literals and provenance. Readers
check the native fragment against the payload, not merely its comment checksum.
These are parameter/topology fragments, not coordinates or simulation scripts.
"""

import base64
import json
from decimal import Decimal
from pathlib import Path
from xml.etree import ElementTree as ET

from island.forcefields.assignment_core.schema import identity, require
from island.forcefields.bundle_adapters.schema import loads, sha

from .registry import OPLSProviderRegistry, fence
from .schema import FLAGS, FORMS, number, sealed, verify_seal

EXTENSIONS = {"gromacs": ".itp", "lammps": ".coeff", "openmm": ".xml"}


def _atoms(c):
    charge = {
        r["sites"][0]: r["values"][0][0]
        for r in c["records"]
        if r["family"] == "charge"
    }
    nb = {
        r["sites"][0]: r["values"][0]
        for r in c["records"]
        if r["family"] == "nonbonded"
    }
    return [(a, charge[a["site"]], *nb[a["site"]]) for a in c["atoms"]]


def _periodic(row):
    if row["form"] == "torsion_periodic":
        return row["values"]
    if row["form"] == "torsion_opls":
        return [
            ["0" if i % 2 else "180", str(number(v) / 2), str(i)]
            for i, v in enumerate(row["values"][0], 1)
        ]
    raise ValueError(
        "export_form_not_supported_no_torsion_reinterpretation:" + row["form"]
    )


def _gromacs(c):
    conv = c["conventions"]
    lines = [
        "[ defaults ]",
        f"1 3 yes {conv['scale14_lj']} {conv['scale14_coulomb']}",
        "[ atomtypes ]",
    ]
    for a, q, sigma, epsilon in _atoms(c):
        lines.append(f"S{a['site']} {a['mass']} {q} A {sigma} {epsilon}")
    lines.extend(["[ moleculetype ]", "ExactOligomer 3", "[ atoms ]"])
    for a, q, _, _ in _atoms(c):
        site = a["site"]
        lines.append(
            f"{site} S{site} 1 OLI {a['element']}{site} {site} {q} {a['mass']}"
        )
    for family, section in (
        ("bond", "bonds"),
        ("angle", "angles"),
        ("torsion", "dihedrals"),
        ("improper", "dihedrals"),
    ):
        if family != "improper":
            lines.append("[ " + section + " ]")
        for row in c["records"]:
            if row["family"] != family or row["applicability"] == "not_applicable":
                continue
            function = {
                "bond_harmonic": "1",
                "angle_harmonic": "1",
                "torsion_periodic": "9",
                "torsion_rb": "3",
                "torsion_opls": "5",
                "improper_harmonic": "2",
            }[row["form"]]
            for term in row["values"]:
                lines.append(" ".join([*map(str, row["sites"]), function, *term]))
    lines.append("[ pairs ]")
    lines.extend(f"{a} {b} 1" for a, b in c["pairs14"])
    lines.append("[ exclusions ]")
    lines.extend(f"{a} {b}" for a, b in c["exclusions"])
    return "\n".join(lines) + "\n"


def _decimal(value):
    return format(value, "f")


def _lammps(c):
    # Exactly declared algebraic unit/half-factor conversions. No fitting or
    # 9-6 conversion. Source numbers remain intact in the exchange envelope.
    lines = [
        "# Parameter fragment only; no coordinates, pair cutoff, integrator or run command.",
        "# units real; atom_style full; pair_style/cutoff compatibility requires a separate engine audit",
        "# source kJ/mol -> kcal/mol: divide by 4.184; nm -> Angstrom: multiply by 10",
        "pair_modify mix geometric",
        "special_bonds lj 0 0 "
        + c["conventions"]["scale14_lj"]
        + " coul 0 0 "
        + c["conventions"]["scale14_coulomb"],
        "bond_style harmonic",
        "angle_style harmonic",
        "dihedral_style fourier",
        "improper_style harmonic",
    ]
    scale = Decimal("4.184")
    for a, q, sigma, epsilon in _atoms(c):
        site = a["site"]
        lines.extend(
            [
                f"mass {site} {a['mass']}",
                f"set atom {site} charge {q}",
                f"pair_coeff {site} {site} {_decimal(number(epsilon) / scale)} {_decimal(number(sigma) * 10)}",
            ]
        )
    counters = {name: 0 for name in ("bond", "angle", "torsion", "improper")}
    for row in c["records"]:
        family = row["family"]
        if family not in counters or row["applicability"] == "not_applicable":
            continue
        counters[family] += 1
        index = counters[family]
        lines.append(
            "# "
            + family
            + "_type "
            + str(index)
            + " sites "
            + " ".join(map(str, row["sites"]))
        )
        if family == "torsion":
            terms = _periodic(row)
            values = " ".join(
                " ".join((_decimal(number(k) / scale), n, phase))
                for phase, k, n in terms
            )
            lines.append(f"dihedral_coeff {index} {len(terms)} {values}")
        else:
            x0, k = row["values"][0]
            factor = 100 if family == "bond" else 1
            kval = number(k) / (2 * scale * factor)
            equilibrium = number(x0) * 10 if family == "bond" else number(x0)
            lines.append(
                f"{family}_coeff {index} {_decimal(kval)} {_decimal(equilibrium)}"
            )
    return "\n".join(lines) + "\n"


def _openmm_exchange(c):
    # A data-only XML exchange, intentionally not a loadable OpenMM ForceField.
    # The consumer must implement geometric sigma and honor exact source order;
    # creating an OpenMM Context or registering a force field is out of scope.
    root = ET.Element(
        "IslandOPLSParameterExchange",
        engine="OpenMM",
        loadableForceField="false",
        productionValidated="false",
    )
    ET.SubElement(
        root,
        "NonbondedContract",
        form="12-6",
        mixing="geometric_sigma_epsilon",
        **{
            k: c["conventions"][k]
            for k in ("scale14_lj", "scale14_coulomb", "long_range", "cutoff", "tail")
        },
    )
    for name in sorted({r["form"] for r in c["records"] if r["form"]}):
        form = ET.SubElement(root, "FunctionalForm", id=name)
        form.text = json.dumps(FORMS[name], sort_keys=True)
    for a, q, sigma, epsilon in _atoms(c):
        ET.SubElement(
            root,
            "Atom",
            site=str(a["site"]),
            type=a["atom_type"],
            element=a["element"],
            mass=a["mass"],
            charge=q,
            sigma=sigma,
            epsilon=epsilon,
        )
    for row in c["records"]:
        item = ET.SubElement(
            root,
            "Interaction",
            id=row["id"],
            family=row["family"],
            sites=" ".join(map(str, row["sites"])),
            form=str(row["form"]),
            applicability=row["applicability"],
        )
        for term in row["values"]:
            ET.SubElement(item, "Values").text = " ".join(term)
    for name, pairs in (("Pair14", c["pairs14"]), ("Exclusion", c["exclusions"])):
        for a, b in pairs:
            ET.SubElement(root, name, site1=str(a), site2=str(b))
    return ET.tostring(root, encoding="unicode") + "\n"


def body(canonical, engine):
    require(engine in EXTENSIONS, "unsupported_export_engine")
    return {"gromacs": _gromacs, "lammps": _lammps, "openmm": _openmm_exchange}[engine](
        canonical
    )


def render(canonical, report, engine):
    require(
        report["decision"] in ("provisional_parameterization", "audited_exact_graph"),
        "incomplete_result_export_forbidden",
    )
    require(
        report["canonical_identity"] == identity(canonical),
        "export_graph_or_parameter_identity_mutation",
    )
    for k, v in FLAGS.items():
        require(report[k] == v and canonical[k] == v, "export_readiness_mutation")
    native = body(canonical, engine)
    receipt = sealed(
        {
            "schema": "island.opls-export-receipt.v1",
            "engine": engine,
            "result_identity": report["identity"],
            "canonical_identity": identity(canonical),
            "graph_identity": canonical["graph_identity"],
            "provider": canonical["provider"],
            "source_hashes": canonical["source_hashes"],
            "native_sha256": sha(native.encode()),
            "non_production": True,
            "scientifically_validated": False,
            "artifact_kind": "data_only_xml_exchange_not_OpenMM_ForceField"
            if engine == "openmm"
            else "parameter_topology_fragment_not_complete_simulation",
            **FLAGS,
        }
    )
    envelope = base64.b64encode(
        json.dumps(
            {"canonical": canonical, "report": report, "receipt": receipt},
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).decode()
    public_receipt = json.dumps(receipt, sort_keys=True, separators=(",", ":"))
    if engine == "openmm":
        return (
            "<!-- ISLAND_OPLS "
            + envelope
            + " -->\n<!-- PROVENANCE "
            + public_receipt
            + " -->\n"
            + native
        )
    prefix = ";" if engine == "gromacs" else "#"
    return (
        prefix
        + " ISLAND_OPLS "
        + envelope
        + "\n"
        + prefix
        + " PROVENANCE "
        + public_receipt
        + "\n"
        + prefix
        + " NON-PRODUCTION; NOT SCIENTIFICALLY VALIDATED\n"
        + native
    )


def read_export(
    path, *, expected_sha256, expected_source_hashes, expected_graph_identity
):
    """Offline round-trip. No engine imports; source pins remain caller supplied."""
    from .registry import pinned

    fence()
    raw = pinned(path, expected_sha256).decode()
    first = raw.splitlines()[0]
    require(" ISLAND_OPLS " in first, "export_provenance_receipt_missing")
    encoded = first.split(" ISLAND_OPLS ", 1)[1].removesuffix(" -->")
    try:
        value = loads(base64.b64decode(encoded, validate=True))
    except (ValueError, UnicodeError):
        raise ValueError("invalid_export_envelope") from None
    c, report, receipt = (value[k] for k in ("canonical", "report", "receipt"))
    verify_seal(report)
    verify_seal(receipt)
    require(
        c["source_hashes"] == sorted(expected_source_hashes),
        "export_source_pin_mismatch",
    )
    require(
        c["graph_identity"] == expected_graph_identity, "export_exact_graph_mismatch"
    )
    require(
        raw == render(c, report, receipt["engine"]), "export_body_or_metadata_mutation"
    )
    return value


def export_parameters(
    settings,
    *,
    engine,
    destination,
    provider,
    manifest,
    manifest_sha256,
    audit=None,
    audit_sha256=None,
):
    """Reopen originals before each export; write exclusively outside Git worktrees."""
    fence()
    require(
        settings["intended_export_engine"] == engine, "intended_export_engine_mismatch"
    )
    destination = Path(destination).absolute()
    resolved = destination.resolve()
    require(
        not any((parent / ".git").exists() for parent in (resolved, *resolved.parents)),
        "licensed_export_must_remain_outside_Git_checkout",
    )
    require(not destination.exists(), "export_destination_must_be_new")
    result = OPLSProviderRegistry().assign(
        settings,
        provider=provider,
        manifest=manifest,
        manifest_sha256=manifest_sha256,
        audit=audit,
        audit_sha256=audit_sha256,
    )
    require(result.canonical is not None, "source_bundle_required_no_export")
    raw = render(result.canonical.payload, result.report, engine).encode()
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("xb") as stream:
        stream.write(raw)
    return {
        "path": str(destination),
        "sha256": sha(raw),
        "source_hashes": result.report["source_hashes"],
        "graph_identity": result.report["graph_identity"],
        "decision": result.report["decision"],
        **FLAGS,
    }
