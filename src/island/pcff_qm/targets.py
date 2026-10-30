"""Exact-graph, single-coefficient QM targets with explicit proposed protocols."""

from island.pcff_completion import _verified, digest, require
from island.pcff_completion.targets import FIXTURE_GRAPHS
from island.pcff_completion.targets import targets as historical_targets

from .geometry import graph_record, target_scan
from .records import Record, attestations, code_identity, record

FIXTURES = tuple(FIXTURE_GRAPHS)
BINDINGS = (
    "graph_identity",
    "typing_identity",
    "charge_identity",
    "atom_order",
    "charge_vector",
    "source_hashes",
    "domain",
)
# These are transparent input-preparation choices, not author/source values and
# not assertions of QM adequacy. Changes create new immutable target identities.
PROTOCOL = {
    "authority": "J57 proposed reference-data preparation policy; not author/source authority",
    "method": "PBE0",
    "basis_set": "def2-SVP",
    "dispersion_treatment": "none",
    "scf_convergence": {
        "energy_hartree": 1e-10,
        "density": 1e-8,
        "max_iterations": 200,
    },
    "geometry_optimization": "none; frozen generated Cartesian geometry; gas-phase single point gradient",
    "boundary": "isolated nonperiodic molecule",
    "method_adequacy": "not_established; basis/method convergence study required before fit validation",
    "charge_state": 0,
    "spin_state": {"multiplicity": 1, "two_S": 0},
}


def targets(system, assessment, fixture):
    p = _verified(system, assessment)
    require(
        FIXTURE_GRAPHS.get(fixture) == p["graph_identity"], "fixture_graph_mismatch"
    )
    out = []
    for prior in historical_targets(system, assessment, fixture):
        t = prior.payload
        row = t["target_site"]
        scan = target_scan(row)
        value = {
            "schema": "island.pcff-qm-target.v2",
            "family": "pcff",
            "status": "qm_fit_required",
            "readiness": "qm_data_missing",
            "previous_target_identity": t["identity"],
            "query": t["query"],
            "interaction_type": t["interaction_type"],
            "atom_type_environment": t["atom_type_environment"],
            "molecular_fixture": fixture,
            "target_site": row,
            "applicable_sites": t["applicable_sites"],
            **{k: p[k] for k in BINDINGS},
            "graph": graph_record(system, p),
            "input_geometry_checksum": digest(
                [system.coordinates.get(i).tolist() for i in p["atom_order"]]
            ),
            "atom_labels": [
                system.topology.sites[i].element + str(i) for i in p["atom_order"]
            ],
            "charge_model": p["charge_origin"],
            "charge_state": 0,
            "spin_state": PROTOCOL["spin_state"],
            "functional_form": t["functional_form"],
            "convention_contract": [r.payload["identity"] for r in attestations()],
            "scan": scan,
            "required_energies": True,
            "required_forces": True,
            "hessian_frequency_requirement": {
                "policy": "full Cartesian Hessian required at central scan point for mixed-derivative identifiability",
                "sample_offsets": [0.0, 0.0],
                "frequencies_alone_sufficient": False,
            },
            "fitting_variables": t["functional_form"]["fields"],
            "parameter_bounds": {
                "policy": "explicit finite user bounds required before fit_ready",
                "numeric_bounds": "not_supplied; no coefficient bounds invented",
            },
            "physical_constraints": [
                "fixed source typing and charges",
                "source Class-II functional form unchanged",
                "positive stable full Hessian at an independently optimized minimum",
                "no silent zeros, charge redistribution or borrowed coefficients",
            ],
            "units": t["units"],
            "qm_protocol": PROTOCOL,
            "expected_output_schema": "island.pcff-qm-reference.v2",
            "held_out_validation": {
                "split_rule": "point index modulo five equals one is held_out; all other points train",
                "thresholds": "user must predeclare energy, force and mixed-Hessian tolerances before fitting",
                "checks": [
                    "all four signs of paired deformations",
                    "full rank residual Jacobian",
                    "no duplicate geometries",
                    "independent graph/domain reconstruction",
                    "QM basis/method convergence and model inadequacy assessment",
                ],
            },
            "uncertainty": {
                "status": "not_quantifiable_without_data",
                "requirements": [
                    "QM convergence",
                    "basis/method dependence",
                    "fit covariance",
                    "held_out residuals",
                    "Class-II model error",
                ],
                "expected_uncertainty": "numerical target to be declared with actual reference data",
            },
            "code_identity": code_identity(),
            "QM_executed": False,
            "parameter_completeness": "incomplete",
            "execution_authority": False,
        }
        out.append(record(value))
    return tuple(out)


def validate_target(system, assessment, target):
    require(type(target) is Record, "immutable_target_required")
    t = target.payload
    require(
        any(
            r.payload == t for r in targets(system, assessment, t["molecular_fixture"])
        ),
        "target_binding_mutation",
    )
    return t


def inventory(cases):
    """Select the smallest authorized fixture deterministically, never a name heuristic.

    Every candidate has already been matched to its pinned exact graph. Graph size
    minimizes computational cost; it confers no transferable domain authority.
    """
    grouped = {}
    for name in FIXTURES:
        s, a = cases[name]
        for t in targets(s, a, name):
            grouped.setdefault(t.payload["query"], []).append(t)
    require(len(grouped) == 19, "nineteen_target_inventory_required")
    selected = []
    for query, candidates in sorted(grouped.items()):
        chosen = min(
            candidates,
            key=lambda r: (
                len(r.payload["atom_order"]),
                r.payload["graph_identity"],
                r.payload["identity"],
            ),
        )
        selected.append(
            {
                "query": query,
                "target": chosen.payload,
                "fixture_bindings": [
                    {
                        "fixture": t.payload["molecular_fixture"],
                        "graph_identity": t.payload["graph_identity"],
                        "target_identity": t.payload["identity"],
                    }
                    for t in candidates
                ],
            }
        )
    return record(
        {
            "schema": "island.pcff-qm-target-inventory.v2",
            "targets": selected,
            "directed_queries": 19,
            "joint_fits_authorized": [],
            "coverage_note": "12 angle-angle and 7 bond-bond_1_3 directed queries, not 19 transferable independent coefficients",
            "inventory_checksum": digest(selected),
        }
    )
