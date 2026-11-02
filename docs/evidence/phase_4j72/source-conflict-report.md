# J72 source conflict adjudication

J71's 12 papers and three supporting documents (287 pages) are retained reading evidence, not 15 newly read documents in J72. All 467 extraction identities and 1,470 numerical cells are preserved. Original file hashes and archive lines are checked separately from offline reconstruction. No new parameter extraction or HMDS fitting/validation. The initial historical suite violated the zero-execution scope: four SQM AM1 minimizations and 46 recorded MD steps, with up to two uncommitted steps possible. These J31 CHO outputs are quarantined and unused in adjudication; execution-scope-report.json reconstructs the actual counts.

## Counts and scope

- J71_records_preserved: 467
- J71_numerical_cells_preserved: 1470
- comparison_pairs: 399
- numerically_compatible_pairs: 25
- equal_scalar_pairs: 23
- pair_conflict_category_counts_overlapping: {'same_row_numeric_conflict': 2, 'unit_conflict': 80, 'equation_conflict': 80, 'phase_conflict': 0, 'prefactor_conflict': 0, 'resolution_conflict': 79, 'chemical_environment_conflict': 371, 'forcefield_family_conflict': 233, 'equilibrium_dependency_conflict': 0, 'source_authority_conflict': 8, 'insufficient_evidence': 295}
- documented_unequal_source_assignment_conflicts: 2
- convention_obstruction_counts_overlapping: {'equilibrium_dependency_conflict': 1, 'insufficient_evidence': 5, 'equation_conflict': 3, 'source_authority_conflict': 2}
- unequal_source_conflicts_mathematically_resolved: 0
- real_duplicate_equivalences_replayed: 16
- synthetic_controls_count_as_real_resolutions: False
- exact_target_instances_resolved: 0
- native_PCFF_instances_resolved: 0
- external_complete_research_instances: 0
- non_native_reference_records: 362
- PCFF_publication_reference_records: 105
- analog_only_instances: 294
- instances_without_candidate_records: 24
- instances_unresolved: 378
- row_status_counts: {'convention_unresolved': 60, 'source_missing': 24, 'analog_unqualified': 294}
- instance_conflict_category_counts_overlapping: {'chemical_environment_conflict': 354, 'equation_conflict': 114, 'equilibrium_dependency_conflict': 354, 'forcefield_family_conflict': 354, 'insufficient_evidence': 378, 'phase_conflict': 60, 'source_authority_conflict': 354, 'unit_conflict': 60}
- rows_requiring_exact_primary_or_author_evidence: 378
- rows_definitely_requiring_QM_fitting: 0
- rows_explicitly_blocked_by_native_AAT_equation_phase_units: 60
- additional_order_and_dependency_blockers: all candidate mappings require complete native equilibrium/source typing; classifications overlap

Pair category counts overlap. They describe structural compatibility screening, not a count of incorrect publications. Numerical comparison occurs only after compatible scalar source contracts are established. The benchmark comparison scope is its nominal COMPASS labels and archive topology, never a native PCFF equivalence contract.

## Unresolved printed assignments

| Evidence | Alternatives | Result | Required evidence |
|---|---|---|---|
| benchmark_SI_archive_angle_c_si_o | cell-0237, cell-0345 | source_conflict_unresolved; neither selected | author correction of TablesS3 labels and topology/source-type mapping; full source/converter contract |
| benchmark_SI_archive_angle_c_si_c | cell-0246, cell-0336 | source_conflict_unresolved; neither selected | author correction of TablesS3 labels and topology/source-type mapping; full source/converter contract |

The C–Si–O SI row is 113.1855/36.2069/−20.3939/20.0172 while topology-bound archive type102 (O–Si–C) is 114.9060/23.0218/−31.3993/24.9814. C–Si–C has the opposite assignments (archive type29). The order reversal of a scalar angle is exact; it does not prove that swapping SI labels is an authoritative correction. No explicit revision is supplied.

## Mathematical equivalence and remaining conventions

Sixteen PolCA RB main/SI pairs agree exactly in their original kJ/mol cells and replay after E_kcal=E_kJ/4.184. These are duplicate confirmations, not repaired unequal conflicts or PCFF rows. Original alternatives and exact citation locations are retained in mathematical-resolutions.json. Degree/radian, coefficient angular-power and energy conversions have independent manufactured controls; those controls are not real-source corrections.

- SI COMPASS Si-C-H reference lengths are reversed relative to printed ordered sites; archive topology cannot by itself resolve converter orientation
- SI after improper_style class2 prints dihedral_coeff si_c_c_h and h_c_h_h; these are not an authoritative improper/AA mapping
- SI UA hybrid 9-6/12-6 description versus archive lj/cut/coul/cut 12-6
- PolCA Table3 OH-H 0.945 nm retained as printed; decimal/unit issue unresolved, no inferred correction
- native FRC literal Phi/AAT authority versus Sun1995 cosine AAT is an inherited unresolved convention, not silently reconciled

The UA printed lj/class2 command versus blanket LJ12-6 description is additionally blocked in comparison metadata, without changing J71 records. Unknown temperature, source target, charge/exclusion contract and equilibrium dependencies remain unknown. No equality or numeric disagreement is inferred across those gaps.

## Authority and missing instances

Native FRC > same-contract author-primary > same-contract published table > explicitly compatible external research > analog/transfer > unexecuted QM target. A publication's authority over its own printed cells does not establish exact c3/sio/osi/hc mapping. Native PCFF/CFF91, CFF93, COMPASS, Smith PDMS, HUANG, OPLS, GAFF, PolCA/UA, coarse-grained and QM corrections remain separate contracts. No cross-family overlay is authorized.

All 378 instances retain exact graph sites, original identities, source lookups, symmetry and equilibrium dependencies in row-adjudication.json. Element-only candidates are analogs. AAT additionally lacks native literal-Phi/cosine equation, phase and unit authority. Missing rows are not zeros. Source absence is not proof that a row requires QM fitting; exact J70 target templates remain unexecuted.

## Decision and capability

Decision: **source_conflict_unresolved**. Zero exact native or complete external instances resolved. All 378 remain unresolved. No execution profile, strict conversion, complete static closure or parameter assignment capability. Ordinary prepare_polymer(..., forcefield="pcff") remains unchanged and fails closed through its existing capability checks. All coefficients remain unaccepted. Chemical transfer, scientific validation and production readiness remain unestablished. Obtain exact primary rows and an authoritative benchmark correction before considering a separately scoped next phase.
