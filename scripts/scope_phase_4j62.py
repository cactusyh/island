"""Map exact source-type environments; no Si-H to methyl transfer authority."""

import json
from collections import defaultdict
from pathlib import Path

from island.charge_references.records import unpack
from island.pcff_completion import pin, signed, write
from island.workflows.bundle import system_from

repo = Path.cwd()
root = repo / "docs/evidence/phase_4j62"
fixtures = {}
for name in ("hmds", "dp2"):
    directory = repo / "docs/evidence/phase_4j54/audit" / name
    report = json.loads((directory / "report.json").read_text())
    system = system_from(unpack((directory / "system.json").read_text()))
    neighbors = defaultdict(list)
    for b in system.topology.bonds.values():
        neighbors[b.site1].append(b.site2)
        neighbors[b.site2].append(b.site1)
    silicon = [
        {
            "atom_id": i,
            "source_type": report["typing"][str(i)],
            "neighbors": sorted(neighbors[i]),
            "neighbor_types": sorted(report["typing"][str(j)] for j in neighbors[i]),
        }
        for i, a in system.topology.sites.items()
        if a.element == "Si"
    ]
    missing = sorted(
        {
            r["family"] + ":" + ",".join(report["typing"][str(i)] for i in r["sites"])
            for r in report["missing_terms"]
        }
    )
    fixtures[name] = {
        "graph_identity": report["graph_identity"],
        "atom_count": len(system.topology.sites),
        "source_record_sha256": pin(directory / "report.json"),
        "silicon_environments": silicon,
        "directed_missing_queries": missing,
        "typing_status": report["typing_status"],
        "charge_status": report["charge_status"],
        "not_independent_coefficient_count": True,
    }
confirmation = (
    repo / "docs/evidence/phase_4j57/audit/author-attested-confirmations.json"
)
write(
    root / "chemical-scope-and-next-data.json",
    signed(
        {
            "schema": "island.j62-chemical-scope.v1",
            "SiH_control": "Exact nine-atom neutral-singlet disiloxane control only; hsi and sio/osi/hsi couplings cannot constrain c3/hc numerical coefficients. Qualitative pipeline and protocol lessons only.",
            "fixtures": fixtures,
            "methyl_numerical_values_transferred_from_SiH": False,
            "original_author_artifact_status": "not_yet_supplied",
            "user_attested_confirmations": {
                "count": 28,
                "record": str(confirmation.relative_to(repo)),
                "sha256": pin(confirmation),
                "evidence_kind": "user_attested_author_confirmation",
                "authority": "user_reported_author_confirmation",
                "supplies_missing_numbers": False,
            },
            "numerical_next_experiment_record": "HMDS-next-experiment-decision.json",
            "next_HMDS_experiment": [
                "Obtain a method-matched stationary candidate at a protocol justified by the measured basis/gradient sensitivity; the twelve-step SVP optimization is unfinished. The numerical next-experiment decision is authoritative after results.",
                "Verify the applicable methyl residual Hamiltonian and the 26 directed missing source queries, reducing only demonstrated reversal/type symmetries; account for all nuisance terms. No blind transfer of the Si-H twenty-column model.",
                "At that one fixed center, obtain two-step Si-C/C-H and methyl/Si-O torsional full Cartesian responses with matched SVP/TZVP or better reference controls; the present single-step/one-direction basis contrast is not enough for an identifiable fit.",
                "Select joint bond/angle/torsional deformations from the complete Jacobian and reserve distinct deformed conformers before fitting.",
            ],
            "HMDS_graph_bound_coordinate_pool": {
                "bonds": [[1, 2], [1, 10], [2, 5], [5, 6]],
                "angles": [[1, 2, 3], [1, 2, 5], [2, 5, 6], [2, 1, 10], [10, 1, 11]],
                "torsions": [[1, 2, 5, 6], [5, 2, 1, 10], [3, 2, 1, 10]],
                "joint_required": "C-Si-C/C-Si-O and Si-C-H/H-C-H angle pairs; first/third bond pairs and torsion-dependent bond/angle responses. Select independent rows only after binding the complete coupled Jacobian.",
                "status": "candidate source/graph mapping only; no additional J62 QM or future fit acceptance implied",
            },
            "subsequent_capped_DP2": "Separate 37-atom source inventory: silicon with two osi and two c3 neighbors introduces interior O-Si-O and Si-O-Si-O environments absent from HMDS end groups. Requires independent graph-bound source/type/charge and numerical targets.",
            "unsupported": [
                "cyclic siloxanes",
                "crosslinked networks",
                "Si-OH",
                "Si-O-C",
                "phosphorus",
                "automatic Si-H to Si-C transfer",
                "complete PCFF/COMPASS or .frc coverage",
            ],
            "ordinary_API": 'forcefield="pcff"',
            "new_registered_profiles": [],
            "production_validated": False,
            "simulation_readiness": "not_established",
        }
    ),
)
