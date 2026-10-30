"""Data-bound external Class-II completion candidates; no profile or MD authority.

J54 source/type/charge reconstruction is mandatory. J55 numerical claims remain
independent. A functional form defines how a coefficient enters E, not its value.
"""

import json
from collections import Counter
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

from island.forcefields.pcff.class2 import FAMILIES, FIELDS
from island.forcefields.pcff.mixed_sip import DISTRIBUTION_PINS
from island.forcefields.selection.inventory import canonical, digest
from island.pcff_pdms import Class2Assessment, reconstruct

VERSION = "island.pcff-class2-completion.v1"
CLASSES = (
    "source_missing",
    "analytically_derivable",
    "requires_qm_fit",
    "requires_experimental_fit",
    "requires_author_confirmation",
    "not_identifiable",
)
SUPPORTED = (
    "bond-bond",
    "bond-angle",
    "angle-angle",
    "bond-bond_1_3",
    "end_bond-torsion_3",
    "middle_bond-torsion_3",
    "angle-torsion_3",
    "angle-angle-torsion_1",
    "torsion_3",
)
# Exact pinned distribution equation spans; never general PCFF defaults.
EQUATIONS = {
    "bond-bond": (3365, 3365),
    "bond-bond_1_3": (3618, 3618),
    "bond-angle": (3690, 3690),
    "angle-angle": (3932, 3932),
    "end_bond-torsion_3": (4210, 4210),
    "middle_bond-torsion_3": (4511, 4512),
    "angle-torsion_3": (4871, 4872),
    "angle-angle-torsion_1": (5185, 5185),
    "torsion_3": (2657, 2657),
}
BASE = "5fd06ec0f25532ed164ba534b04a52a7d7a0aedb"
BRANCH = "codex/phase-4j56-pcff-pdms-missing-term-fitting"


class CompletionError(ValueError):
    def __init__(self, reason, **details):
        self.reason, self.details = reason, details
        super().__init__(reason + ": " + canonical(details))


def require(ok, reason, **details):
    if not ok:
        raise CompletionError(reason, **details)


def pin(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signed(value):
    return {**value, "identity": digest(value)}


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def code_identity():
    return digest({p.name: pin(p) for p in sorted(Path(__file__).parent.glob("*.py"))})


def _verified(system, assessment):
    require(type(assessment) is Class2Assessment, "verified_J54_assessment_required")
    p = assessment.validate_integrity(system)
    require(p["family"] == "pcff", "wrong_profile_family")
    require(p["typing_status"] == "complete", "typing_evidence_required")
    require(p["charge_status"] == "complete", "charge_evidence_required")
    require(
        not p["normalization"] and not p["redistribution"], "charge_repair_forbidden"
    )
    return p


def equation(assessment, family):
    require(family in EQUATIONS, "unsupported_functional_form", family=family)
    assessment.assets.load()
    path = Path(assessment.assets.directory) / "pcff.frc"
    lines = path.read_text().splitlines()
    lo, hi = EQUATIONS[family]
    p = {
        "source_file": "pcff.frc",
        "source_sha256": DISTRIBUTION_PINS["pcff.frc"],
        "citation": "LAMMPS pcff.frc cff91 distribution, revision "
        "e891a3e10973c1a729e391a0aefaa02fd70f8c0f; internal source references retained by J54",
        "doi": None,
        "doi_status": "not assigned to this distribution",
        "source_lines": [lo, hi],
        "raw_text": "\n".join(lines[lo - 1 : hi]),
        "functional_form": family,
        "fields": FIELDS[family],
        "source_units": FAMILIES[family][2],
        "unit_symbols": {"E": "kcal/mol", "L": "angstrom", "A": "radian"},
        "unit_interpretation": "J54 pinned distribution converter; angular displacements are radians",
        "coefficient_value_equation": None,
        "reason": "energy equation contains independent coefficients; does not determine missing values",
    }
    return signed(p)


def classify(system, assessment):
    """Classify each exact missing site; retain directed source queries and chains."""
    p = _verified(system, assessment)
    equations = {
        f: equation(assessment, f)
        for f in sorted({r["family"] for r in p["missing_terms"]})
    }
    rows = []
    missing_ids = {r["id"] for r in p["missing_terms"]}
    for row in p["missing_terms"]:
        f = row["family"]
        torsion = "torsion" in f
        category = "requires_author_confirmation" if torsion else "requires_qm_fit"
        reason = (
            "linear_Phi_vs_cos_Phi_conflict_and_undefined_angular_power"
            if f == "angle-angle-torsion_1"
            else "source_torsion_geometry_not_resolved"
            if torsion
            else "independent_coupling_coefficient_requires_identifying_reference_data"
        )
        rows.append(
            {
                "id": row["id"],
                "family": f,
                "query": row["query"],
                "sites": row["sites"],
                "source_status": "source_missing",
                "classification": category,
                "reason": reason,
                "dependencies": row["dependencies"],
                "missing_dependencies": [
                    d for d in row["dependencies"] if d in missing_ids
                ],
                "equation_evidence": equations[f]["identity"],
                "analytical_derivation_authorized": False,
                "required_reference": "licensed, hash-bound distorted geometries with identifiable "
                "isolated/residual energies or forces; Hessian/vibrational data when available; "
                "source-bound typing/charges and domain-specific end-group evidence",
                "domain": "exact disiloxane control only; Si-H transfer forbidden"
                if p["domain"]["SiH_control"]
                else "exact linear TMS-capped methyl PDMS graph",
            }
        )
    return signed(
        {
            "schema": VERSION,
            "family": "pcff",
            "candidate_profile": p["candidate_profile"],
            "assessment_identity": p["identity"],
            "graph_identity": p["graph_identity"],
            "typing_identity": p["typing_identity"],
            "charge_identity": p["charge_identity"],
            "source_hashes": p["source_hashes"],
            "code_identity": code_identity(),
            "terms": sorted(rows, key=lambda r: r["id"]),
            "classification_counts": dict(
                sorted(Counter(r["classification"] for r in rows).items())
            ),
            "equations": equations,
            "missing_conventions": p["missing_conventions"],
            "parameter_completeness": "incomplete",
            "cross_engine_equivalence": "not_claimed",
            "scientific_validation": "not_established",
            "production_readiness": "not_established",
            "production_validated": False,
            "simulation_readiness": "not_established",
            "registered_profiles": [],
            "created_terms": [],
            "MD_authorized": False,
        }
    )


def derive(system, assessment, target_id):
    p = _verified(system, assessment)
    row = next((r for r in p["missing_terms"] if r["id"] == target_id), None)
    require(row is not None, "exact_missing_target_required")
    e = equation(assessment, row["family"])
    # Equilibrium lengths/angles cannot determine an independent mixed derivative.
    return signed(
        {
            "schema": VERSION,
            "operation": "derive",
            "status": "blocked",
            "target": row["id"],
            "query": row["query"],
            "graph_identity": p["graph_identity"],
            "typing_identity": p["typing_identity"],
            "charge_identity": p["charge_identity"],
            "reason": "no_source_equation_determines_missing_coefficient",
            "source_equation": e,
            "formula": None,
            "input_rows": row["dependencies"],
            "conversion": None,
            "implementation": VERSION,
            "code_identity": code_identity(),
            "uncertainty": {"status": "not_identifiable", "value": None},
            "independent_reconstruction": "no coefficient relation to reconstruct",
            "origin_on_success": "derived_external",
            "result": None,
            "parameter_completeness": "incomplete",
            "execution_authority": False,
        }
    )


def inventory(root):
    """Recursive payload inventory: no links, traversal, or mutable outside assets."""
    root = Path(root)
    require(root.is_dir() and not root.is_symlink(), "dataset_unavailable")
    entries = {}
    for p in sorted(root.rglob("*")):
        require(not p.is_symlink(), "dataset_symlink_forbidden", path=str(p))
        require(
            p.is_file() or p.is_dir(), "dataset_special_file_forbidden", path=str(p)
        )
        if p.is_file():
            entries[p.relative_to(root).as_posix()] = pin(p)
    require(bool(entries), "reference_data_unavailable")
    return entries


def read_pinned(root, name, hashes):
    name = Path(name)
    require(not name.is_absolute() and ".." not in name.parts, "unsafe_dataset_path")
    key = name.as_posix()
    require(key in hashes, "unlisted_dataset_file", file=key)
    p = Path(root) / name
    require(pin(p) == hashes[key], "dataset_hash_mutation", file=key)
    return json.loads(p.read_text())


def _finite(value):
    import numpy as np

    a = np.asarray(value, dtype=float)
    require(np.isfinite(a).all(), "nonfinite_reference_value")
    return a


def intake(system, assessment, request, data_root):
    """Validate a fitting contract; total QM data are never isolated by assertion.

    Contract intake accepts all reference channels. The bounded optimizer below
    supports independently documented isolated component energies for BB/BA/AA.
    Total energies, forces, Hessians and experiments remain visible pending a
    source-bound joint-residual/decomposition adapter. No arbitrary design matrix.
    """
    p = _verified(system, assessment)
    require(request.get("schema") == VERSION, "wrong_fitting_schema")
    require(request.get("family") == "pcff", "wrong_profile_family")
    require(request.get("source_hashes") == p["source_hashes"], "source_hash_mutation")
    for k in (
        "graph_identity",
        "typing_identity",
        "charge_identity",
        "atom_order",
        "charge_vector",
    ):
        require(request.get(k) == p[k], k + "_mutation")
    require(
        request.get("candidate_profile") == p["candidate_profile"], "profile_mutation"
    )
    require(request.get("domain") == p["domain"], "domain_or_terminal_group_mutation")
    target = next(
        (r for r in p["missing_terms"] if r["id"] == request.get("target_id")), None
    )
    require(target is not None, "exact_missing_target_required")
    require(request.get("query") == target["query"], "target_query_mutation")
    require(target["family"] in SUPPORTED, "unsupported_functional_form")
    eq = equation(assessment, target["family"])
    require(
        request.get("equation_identity") == eq["identity"], "functional_form_mutation"
    )
    hashes = inventory(data_root)
    require(request.get("files") == hashes, "dataset_inventory_or_hash_mutation")
    provenance = request.get("provenance", {})
    for k in ("citation", "license", "data_origin", "reference_version"):
        require(
            isinstance(provenance.get(k), str) and bool(provenance[k].strip()),
            "missing_" + k,
        )
    license_file = provenance.get("license_file")
    require(license_file in hashes, "license_evidence_not_pinned")
    require(
        (Path(data_root) / license_file).read_text().strip()
        == provenance["license"].strip(),
        "license_evidence_mutation",
    )
    require(
        provenance.get("permission_status") in ("licensed", "permission_granted"),
        "data_permission_unavailable",
    )
    require(
        provenance.get("kind") in ("qm", "experimental", "reference_engine"),
        "reference_kind_required",
    )
    if provenance["kind"] == "qm":
        require(
            all(provenance.get(k) for k in ("method", "basis", "software")),
            "qm_method_basis_required",
        )
    elif provenance["kind"] == "reference_engine":
        require(provenance.get("software"), "reference_engine_required")
    require(
        request.get("units")
        == {
            "coordinates": "angstrom",
            "energy": "kcal/mol",
            "forces": "kcal/mol/angstrom",
            "hessian": "kcal/mol/angstrom^2",
            "angles": "radian",
            "charge": "e",
        },
        "unresolved_or_mutated_units",
    )
    data = read_pinned(data_root, request.get("samples_file", ""), hashes)
    require(isinstance(data, list) and bool(data), "reference_samples_required")
    n = len(p["atom_order"])
    ids = []
    for s in data:
        require(
            isinstance(s.get("id"), str) and bool(s["id"]), "sample_identity_required"
        )
        ids.append(s["id"])
        require(
            _finite(s.get("coordinates", [])).shape == (n, 3),
            "geometry_or_atom_order_incomplete",
        )
        require(s.get("split") in ("train", "held_out"), "held_out_split_required")
        require(s.get("atom_order") == p["atom_order"], "sample_atom_order_mutation")
        require(
            any(
                s.get(k) is not None
                for k in ("energy", "forces", "hessian", "vibrations", "observables")
            ),
            "reference_observable_required",
        )
        for k, shape in (
            ("energy", ()),
            ("forces", (n, 3)),
            ("hessian", (3 * n, 3 * n)),
        ):
            if s.get(k) is not None:
                require(
                    _finite(s[k]).shape == shape, "reference_shape_mismatch", channel=k
                )
        if s.get("vibrations") is not None:
            v = s["vibrations"]
            require(
                v.get("units") == "cm^-1"
                and _finite(v.get("frequencies", [])).size > 0,
                "vibrational_units_or_values_required",
            )
        if s.get("observables") is not None:
            require(
                isinstance(s["observables"], list) and bool(s["observables"]),
                "experimental_observables_required",
            )
            for v in s["observables"]:
                require(
                    all(v.get(k) for k in ("name", "units", "method")),
                    "experimental_observable_metadata_required",
                )
                require(
                    _finite(v["uncertainty"]).shape == ()
                    and float(v["uncertainty"]) > 0,
                    "experimental_uncertainty_required",
                )
                _finite(v["value"])
    require(len(ids) == len(set(ids)), "duplicate_sample_identity")
    fit = request.get("fit", {})
    require(
        fit.get("objective") == "weighted_isolated_component_energy_least_squares",
        "objective_not_supported",
    )
    fields = FIELDS[target["family"]]
    require(fit.get("fields") == fields, "coefficient_order_mutation")
    for k in ("initial_values", "lower_bounds", "upper_bounds"):
        require(
            _finite(fit.get(k, [])).shape == (len(fields),),
            "fit_configuration_incomplete",
            field=k,
        )
    require(
        all(
            lo <= v <= hi
            for lo, v, hi in zip(
                fit["lower_bounds"],
                fit["initial_values"],
                fit["upper_bounds"],
                strict=True,
            )
        ),
        "invalid_bounds_or_initial_values",
    )
    require(
        fit.get("constraints")
        == ["source_equation", "fixed_verified_equilibria", "exact_graph_domain"],
        "constraints_required",
    )
    require(
        fit.get("optimizer") == "numpy.linalg.lstsq_SVD_v1", "optimizer_not_supported"
    )
    weights = fit.get("weights", {})
    require(
        set(weights) == set(ids)
        and all(_finite(v).shape == () and float(v) > 0 for v in weights.values()),
        "invalid_weights",
    )
    blockers = []
    if "torsion" in target["family"]:
        blockers.append(
            "angle_angle_torsion_form"
            if target["family"] == "angle-angle-torsion_1"
            else "torsion_geometry"
        )
    if not all(s.get("energy") is not None for s in data):
        blockers.append(
            "energy_channel_required_by_bounded_optimizer; force/Hessian/experimental_adapter_not_authorized"
        )
    if request.get("energy_scope") != "isolated_target_component":
        blockers.append(
            "total_reference_requires_verified_joint_residual_decomposition"
        )
    isolation_file = request.get("isolation_evidence_file")
    if isolation_file is None:
        blockers.append("independent_component_isolation_evidence_required")
    else:
        evidence = read_pinned(data_root, isolation_file, hashes)
        require(
            evidence.get("query") == target["query"]
            and evidence.get("target_id") == target["id"]
            and evidence.get("graph_identity") == p["graph_identity"],
            "isolation_target_mutation",
        )
        for k in (
            "method",
            "code_file",
            "raw_reference_file",
            "equation_identity",
            "reference_command",
        ):
            require(bool(evidence.get(k)), "isolation_evidence_incomplete", field=k)
        require(
            evidence["equation_identity"] == eq["identity"],
            "isolation_equation_mutation",
        )
        for k in ("code_file", "raw_reference_file"):
            require(evidence[k] in hashes, "isolation_asset_not_pinned", field=k)
        # Cross-check the reference engine's component values, not total E.
        raw = read_pinned(data_root, evidence["raw_reference_file"], hashes)
        require(
            raw.get("query") == target["query"]
            and raw.get("target_id") == target["id"]
            and raw.get("graph_identity") == p["graph_identity"],
            "raw_reference_identity_mutation",
        )
        require(raw.get("samples") == data, "raw_reference_samples_mutation")
    if len([s for s in data if s["split"] == "train"]) <= len(fields):
        blockers.append("insufficient_training_degrees_of_freedom")
    if not any(s["split"] == "held_out" for s in data):
        blockers.append("held_out_data_required")
    for d in target["dependencies"]:
        if any(r["id"] == d for r in p["missing_terms"]):
            blockers.append("missing_equilibrium_dependency:" + d)
    return signed(
        {
            "schema": VERSION,
            "operation": "fit_intake",
            "request": request,
            "dataset_identity": digest({"files": hashes, "provenance": provenance}),
            "target": target,
            "source_equation": eq,
            "assessment_identity": p["identity"],
            "status": "blocked" if blockers else "data_bound_candidate_not_executed",
            "blockers": sorted(set(blockers)),
            "code_identity": code_identity(),
            "origin_on_success": "fitted_external",
            "created_terms": [],
            "parameter_completeness": "incomplete",
            "execution_authority": False,
        }
    )


def _design(system, p, target, samples):
    """Independent K columns from actual Cartesian geometry; no supplied matrix."""
    from math import acos

    import numpy as np

    order = {i: j for j, i in enumerate(p["atom_order"])}
    inventory_rows = {r["id"]: r for r in p["interaction_inventory"]}

    def equilibrium(key):
        row = inventory_rows[key]["assignment"]
        require(
            row is not None and row["status"] == "assigned",
            "missing_equilibrium_dependency",
        )
        # J54 converted source units are kJ/mol, angstrom, radians.
        return float(row["values"][0])

    def rkey(i, j):
        return "quartic_bond:" + ",".join(map(str, sorted((i, j))))

    def akey(i, j, k):
        t = min((i, j, k), (k, j, i))
        return "quartic_angle:" + ",".join(map(str, t))

    def distance(x, i, j):
        return float(np.linalg.norm(x[order[i]] - x[order[j]]))

    def angle(x, i, j, k):
        a, b = x[order[i]] - x[order[j]], x[order[k]] - x[order[j]]
        den = np.linalg.norm(a) * np.linalg.norm(b)
        require(den > 0, "degenerate_geometry")
        return acos(float(np.clip(np.dot(a, b) / den, -1, 1)))

    f, ids = target["family"], target["sites"]
    design = []
    for s in samples:
        x = np.array(s["coordinates"], dtype=float)
        if f in ("bond-bond", "bond-angle"):
            i, j, k = ids
            a = distance(x, i, j) - equilibrium(rkey(i, j))
            b = distance(x, j, k) - equilibrium(rkey(j, k))
            if f == "bond-bond":
                v = [a * b]
            else:
                t = angle(x, i, j, k) - equilibrium(akey(i, j, k))
                v = [a * t, b * t]
        elif f == "bond-bond_1_3":
            i, j, k, l = ids
            v = [
                (distance(x, i, j) - equilibrium(rkey(i, j)))
                * (distance(x, k, l) - equilibrium(rkey(k, l)))
            ]
        elif f == "angle-angle":
            i, j, k, l = ids
            v = [
                (angle(x, i, j, k) - equilibrium(akey(i, j, k)))
                * (angle(x, k, j, l) - equilibrium(akey(k, j, l)))
            ]
        else:
            raise CompletionError(
                "source_torsion_convention_requires_author_confirmation"
            )
        design.append(v)
    return np.array(design)


def fit(system, assessment, request, data_root):
    """Fit a candidate only with identifiable, isolated, licensed reference data.

    The result never registers a profile or closes unrelated missing rows or
    conventions. Experimental/force/Hessian joint inversion is intake-only v1.
    """
    import numpy as np

    contract = intake(system, assessment, request, data_root)
    require(not contract["blockers"], "fitting_blocked", blockers=contract["blockers"])
    p = _verified(system, assessment)
    data = read_pinned(data_root, request["samples_file"], request["files"])
    x = _design(system, p, contract["target"], data)
    y = np.array([s["energy"] for s in data], dtype=float)
    w = np.array([request["fit"]["weights"][s["id"]] for s in data])
    train = np.array([s["split"] == "train" for s in data])
    a, b = x[train] * np.sqrt(w[train, None]), y[train] * np.sqrt(w[train])
    coeff, _, rank, singular = np.linalg.lstsq(a, b, rcond=None)
    require(rank == x.shape[1], "not_identifiable", rank=int(rank))
    require(np.isfinite(coeff).all(), "nonfinite_fit")
    cfg = request["fit"]
    require(
        np.all(coeff >= cfg["lower_bounds"]) and np.all(coeff <= cfg["upper_bounds"]),
        "unconstrained_optimum_outside_declared_bounds",
    )
    second = np.linalg.solve(a.T @ a, a.T @ b)
    require(
        np.allclose(coeff, second, rtol=1e-10, atol=1e-12),
        "independent_fit_reconstruction_failed",
    )
    residual = x @ coeff - y
    dof = int(train.sum() - len(coeff))
    objective = float(np.sum(w[train] * residual[train] ** 2))
    covariance = np.linalg.inv(a.T @ a) * objective / dof
    require(np.isfinite(covariance).all(), "uncertainty_not_identifiable")
    parameters = {k: repr(float(v)) for k, v in zip(cfg["fields"], coeff, strict=True)}
    return signed(
        {
            "schema": VERSION,
            "operation": "fit",
            "origin": "fitted_external",
            "status": "fitted_candidate_not_registered",
            "contract": contract,
            "dataset_identity": contract["dataset_identity"],
            "data_hashes": request["files"],
            "target_interaction": contract["target"],
            "functional_form": contract["source_equation"],
            "graph_identity": p["graph_identity"],
            "typing_identity": p["typing_identity"],
            "charge_identity": p["charge_identity"],
            "charge_origin": p["charge_origin"],
            "source_hashes": p["source_hashes"],
            "equilibrium_input_rows": {
                r["id"]: r["assignment"]
                for r in p["interaction_inventory"]
                if r["id"] in contract["target"]["dependencies"]
            },
            "parameters": parameters,
            "units": [
                u.replace("E", "kcal/mol")
                .replace("L", "angstrom")
                .replace("A", "radian")
                for u in contract["source_equation"]["source_units"]
            ],
            "conversion": "identity: kcal/mol, angstrom, radian; no charge or source-row modification",
            "objective": objective,
            "objective_function": cfg["objective"],
            "weights": cfg["weights"],
            "constraints": cfg["constraints"],
            "bounds": [cfg["lower_bounds"], cfg["upper_bounds"]],
            "initial_values": cfg["initial_values"],
            "optimizer": cfg["optimizer"],
            "convergence": {
                "status": "full_rank_global_linear_optimum",
                "rank": int(rank),
                "singular_values": singular.tolist(),
            },
            "uncertainty": {
                "method": "residual_based_coefficient_covariance",
                "covariance": covariance.tolist(),
                "degrees_of_freedom": dof,
                "limitation": "conditional statistical uncertainty; excludes source/QM/model error",
            },
            "held_out_validation": {
                "ids": [s["id"] for s in data if s["split"] == "held_out"],
                "residuals": residual[~train].tolist(),
                "weighted_objective": float(np.sum(w[~train] * residual[~train] ** 2)),
                "status": "recorded_not_scientifically_validated",
            },
            "independent_reconstruction": {
                "method": "normal_equations_separate_from_SVD",
                "values": second.tolist(),
            },
            "implementation": VERSION,
            "code_identity": code_identity(),
            "numpy_version": np.__version__,
            "final_parameter_checksum": digest(parameters),
            "evidence_checksum": digest(contract),
            "parameter_completeness": "incomplete",
            "execution_authority": False,
            "registered_profiles": [],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    )


@dataclass(frozen=True)
class CompletionRecord:
    json_text: str

    @property
    def payload(self):
        p = json.loads(self.json_text)
        require(
            p.get("identity")
            == digest({k: v for k, v in p.items() if k != "identity"}),
            "completion_checksum_mutation",
        )
        return p

    def validate_integrity(self, system, assessment, *, data_root=None):
        p = self.payload
        op = p["operation"]
        if op == "derive":
            expected = derive(system, assessment, p["target"])
        elif op == "fit_intake":
            expected = intake(system, assessment, p["request"], data_root)
        elif op == "fit":
            expected = fit(system, assessment, p["contract"]["request"], data_root)
        else:
            raise CompletionError("unsupported_completion_operation")
        require(p == expected, "completion_evidence_mutation")
        return p


def load_assessment(directory, expected_sha256, distribution):
    from island.pcff_pdms import PCFFPDMSAssets

    return reconstruct(directory, expected_sha256, assets=PCFFPDMSAssets(distribution))


def request_template(system, assessment, target_id):
    """An unfilled contract, with null initial values, is not reference data."""
    p = _verified(system, assessment)
    target = next((r for r in p["missing_terms"] if r["id"] == target_id), None)
    require(target is not None, "exact_missing_target_required")
    eq = equation(assessment, target["family"])
    count = len(FIELDS[target["family"]])
    return {
        "schema": VERSION,
        "family": "pcff",
        **{
            k: p[k]
            for k in (
                "candidate_profile",
                "source_hashes",
                "graph_identity",
                "typing_identity",
                "charge_identity",
                "atom_order",
                "charge_vector",
                "domain",
            )
        },
        "target_id": target_id,
        "query": target["query"],
        "equation_identity": eq["identity"],
        "files": {},
        "samples_file": None,
        "energy_scope": None,
        "isolation_evidence_file": None,
        "provenance": {
            "kind": None,
            "citation": None,
            "license": None,
            "license_file": None,
            "permission_status": "not_supplied",
            "data_origin": None,
            "reference_version": None,
            "method": None,
            "basis": None,
            "software": None,
        },
        "units": {
            "coordinates": "angstrom",
            "energy": "kcal/mol",
            "forces": "kcal/mol/angstrom",
            "hessian": "kcal/mol/angstrom^2",
            "angles": "radian",
            "charge": "e",
        },
        "fit": {
            "objective": "weighted_isolated_component_energy_least_squares",
            "fields": FIELDS[target["family"]],
            "initial_values": [None] * count,
            "lower_bounds": [None] * count,
            "upper_bounds": [None] * count,
            "constraints": [
                "source_equation",
                "fixed_verified_equilibria",
                "exact_graph_domain",
            ],
            "optimizer": "numpy.linalg.lstsq_SVD_v1",
            "weights": {},
        },
    }
