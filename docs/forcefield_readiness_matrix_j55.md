# Independent readiness matrix (J55)

| Claim | Sufficient evidence | Does not establish |
| --- | --- | --- |
| Parameter complete | Independently validated exact interaction and convention inventory, typing, charges, source integrity and deterministic serialization | Engine execution, agreement, equilibrium |
| LAMMPS export passed | Actual deterministic source-bound data/input generation | LAMMPS runtime or physical qualification |
| LAMMPS runtime passed | Pinned executable completes and returns valid static observables | OpenMM agreement, MD or pressure qualification |
| OpenMM build passed | Actual supported model construction | Static evaluation, MD or agreement |
| OpenMM runtime passed | Valid static evaluation with the constructed model | LAMMPS agreement or equilibrium |
| Cross-engine equivalence passed | Explicitly scoped, independently verified comparison contract | Universal force-field validity or production readiness |
| Scientific validation passed | Existing domain/state-point scientific gates independently satisfied | Universal polymer coverage |
| Production readiness | All applicable project-wide scientific/operational gates | Broader chemical-domain coverage |

A complete profile remains complete when optional comparisons fail. When no
comparison is requested its status is `not_claimed`, including after single-engine
runtime checks. Missing parameters or necessary physical conventions still produce
`incomplete` even if a diagnostic engine can build a system.

Current native ethane finite OPLS control: parameter complete; OpenMM build/static
passed; no native LAMMPS export adapter; cross-engine equivalence not claimed.
Existing analytical external XML tests also exercise a complete, explicitly
software-only domain with failed external comparison observations. They supply no
new physical literature parameter or registered PDMS profile.

Current HMDS/capped DP2/DP4/DP6: source-bound interaction rows, typing and charges
complete; **overall parameter/convention completeness incomplete** because source
`dielectric` and `periodic_image_policy` remain unresolved. LAMMPS diagnostic
export/static and OpenMM build/static are reported independently. The seven other
J52 gaps (Coulomb implementation constant, precision, neighbor list, pressure,
virial, table interpolation, erfc kernel) remain explicit backend/observable
qualifications. J55 neither assigns them source authority nor erases their old
MD blockers. Source-conflicting same-scope conventions still reject execution;
the archive DFT branch is not combined with the differently scoped SI table.

All current scientific/production states remain `not_established`,
`production_validated=False`, `simulation_readiness="not_established"`.
