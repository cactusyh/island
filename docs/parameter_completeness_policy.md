# Parameter completeness policy v1 (J55)

Parameter coverage, backend operation, cross-engine equivalence and scientific
validation are independent claims. Selecting a profile or assigning its parameters
never requires trajectory reproduction between engines. A numerical disagreement
cannot remove source-bound rows, change charges, or make an assigned row missing.

Use `island.forcefield_policy.prepare_polymer_forcefield(system, forcefield="opls")`
or `prepare_polymer(psmiles=..., forcefield="opls", ...)`. Profile and charge
selection remain internal to the canonical OPLS resolver. Administrators install
hash-bound assets once with `OPLSAssets` / `use_opls_assets` and, for PDMS source
conventions, `ConventionSources` / `use_convention_sources`.

Completeness requires exact graph/domain typing and charges, masses, all required
bond/angle/proper/improper/nonbonded/special-pair rows, all applicable Class-II
couplings and their dependency chains, evidence, source hashes, units, functional
forms, exclusions, separate Coulomb/LJ 1–4 scales, mixing, and declared cutoff and
long-range conventions. Missing requirements remain explicitly incomplete. A
source-declared absence can be inapplicable with evidence; it cannot be replaced
by an invented zero. No cross-family borrowing or charge repair is permitted.
Source-validated legacy converters keep their original scope and evidence.

The immutable J55 report owns separate `parameter_completeness`, `typing_status`,
`charge_status`, `lammps_export_status`, `lammps_runtime_status`,
`openmm_build_status`, `openmm_runtime_status`, `cross_engine_equivalence`,
`scientific_validation` and `production_readiness`. `parameter_complete` is the
Boolean alias for the first field. It includes the exact interaction inventory,
missing terms/conventions, original assignment, source identity, charge record,
selected external/native identity, backend choices, and historical source gaps.

Numerical statuses initialize to `not_run`; cross-engine equivalence initializes
to `not_claimed`. Constructing an evaluator records only build success. Evaluating
one static state records only that backend's runtime success. Export records only
export success. Each operation validates live source/graph/charge identity again.
These operations can succeed for an explicitly scoped diagnostic model while an
author-convention audit remains incomplete. Such a result cannot authorize MD.
Backend absence/failure cannot alter the parameter inventory.

Optional component observations can be attached through
`record_static_comparison`, with independently pinned raw evidence and exact
assignment/input identities. J49 energy/force/FD bounds remain unchanged. Failed
observations set equivalence to `failed`; passing bounded static observations
alone leave broad equivalence `not_claimed`. The importer labels these as external
observations, does not execute engines or certify their measurement provenance,
and never promotes scientific or production status. Existing J47–J53 diagnostic
drivers remain available with their original contracts. Their raw numerical
results remain historical evidence, not completeness gates.

`ClaimRequirements` and CLI flags independently require completeness, LAMMPS
runtime, OpenMM runtime, cross-engine equivalence or scientific validation. The
runtime flags perform a backend-specific static check, with explicit assets and
an independently pinned LAMMPS executable where needed. Equivalence/science flags
check those claims; they do not silently run MD or manufacture verification.
No current J55 adapter grants a broad equivalence or scientific-validation pass.

Saved bundles preserve observation snapshots under an independently pinned
manifest. Reconstruction rehashes relocated sources and regenerates parameters,
conventions and completeness. Observation snapshots do not restore an engine,
checkpoint, RNG state or MD authority. JSON does not grant execution rights.

GAFF/GAFF2 continue through the existing AmberTools adapter, retaining native
AM1-BCC, provided-charge and verified RESP/Psi4 provenance. No new charge model or
normalization is introduced. Their finite imported convention audit uses the
existing source exclusions, SCEE/SCNB, parameter units and converter definitions.
PCFF still requires the full exact Class-II inventory and explicit model/special
pair conventions; missing couplings reject completeness. Existing converter
interpretations remain labelled converter choices, with source discrepancies
preserved. No missing coupling is inserted by this policy.
