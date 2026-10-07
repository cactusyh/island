# Phase 4J19 — unified TypedGraph and assignment records

## Base and architecture

J18 was verified as merged before this phase: fetched `origin/main` is
`b28f262` (`Add immutable unified final graph transformations`), on top of the
J17 merge. The J19 branch starts at that fetched main.

J19 adds an additive `island.graph.unified` layer:

- `UnifiedForceFieldSource` (`island_unified_force_field_source_v1`) binds the
  family, exact source hash, source files/executable provenance and evidence.
- `UnifiedTypedGraph` (`island_unified_typed_graph_v1`) binds the exact
  `FinalChemicalGraph` identity, component/molecule identities, integer site
  IDs and labels, family/source/profile/method, provenance, evidence, and the
  automatic-perception state.
- `UnifiedGraphCharges` (`island_unified_graph_charges_v1`) binds final graph
  and typed-graph identities, family/source/method/origin, exact finite vector,
  units, component totals, total charge, provenance and evidence. It never
  normalizes or redistributes values.
- `UnifiedAssignmentDiagnostics`
  (`island_unified_assignment_diagnostics_v1`) separates typing, charge,
  native-charge availability, parameter-row, executable-potential, periodic
  backend, and scientific-verification states.
- `UnifiedParameterAssignment`
  (`island_unified_parameter_assignment_v1`) binds all three identities, family,
  source, native family payload, and the diagnostics record.

The `unify_pcff_*`, `unify_opls_*`, and `unify_amber_preparation` functions are
named adapters. They validate family-native records and preserve their source,
charge, parameter and executable conventions. They do not rewrite historical
records or fabricate automatic/native records. PCFF J17 records and identities
remain unchanged.

`ForceFieldRequest` accepts `final_graph`, `typed_graph`, and `graph_charges`.
The PCFF common path adapts unified records into the existing named PCFF
external-record APIs. OPLS-AA and GAFF/GAFF2 unified records reject before
backend execution when their current native backend has no externally assigned
input adapter; this is a precise no-fallback diagnostic. Existing automatic
OPLS and AmberTools APIs remain unchanged.

## Fixtures and validation

The compatibility fixtures are the J18 final graphs: PSMILES
`[*:1]CC[*:2]` at DP3 and DP10, the controlled nonperiodic crosslink graph,
and the periodic graph record. Unified records support both periodic and
nonperiodic neutral graphs. Current force-field evaluators remain finite
nonperiodic implementations; periodic preparation rejects explicitly.

The focused J19 suite checks exact IDs, graph/type/charge/source/family/method
identity, finite values and units, component totals, stale/tampered records,
wrong-family fallback, PCFF native/provided origins, and common PCFF preparation.
The ordinary suite, J17 PCFF regressions, OPLS regressions, AmberTools regressions,
Ruff and pip checks remain required gates.

No crosslink generation, packing, periodic force evaluation, MDAnalysis,
parameter fitting, wildcard fallback, zero coupling, source mixing, full PCFF
coverage, msi2lmp parity, or scientific production claim is added. Existing
historical bundles and workflows retain their records and identities.

`production_validated=False` and `simulation_readiness="not_established"`.

For the PCFF common path, callers select the existing expanded source profile
(`island_pcff_source_graph_v1` or a later authorized profile) and its explicit
resolution policy in `PCFFOptions`; this keeps the legacy acyclic automatic
profile guard unchanged. The unified records themselves remain the authoritative
final-graph/type/charge identity.

## Corrective integrity revision

The corrective J19 commit rejects the reproduced cross-family defect where a
PCFF unified source declared `source.family="GAFF"`. If a source family is
present it must exactly equal the requested family. Source hashes must be
lowercase hexadecimal SHA-256 strings of exactly 64 characters.

Typed graph validation now independently recomputes component and molecule
identities against the supplied `FinalChemicalGraph`. Charge validation checks
family/source equality, typed-graph identity, component-total keys and values,
finite nonboolean elementary-charge values, and total-charge equality. Parameter
assignment validation checks family/source equality across typed and charge
records, all graph/type/charge identities, and diagnostics schema/identity.

The PCFF common preparation path compares both unified source hashes to the
currently loaded pinned FRC before adapting into the legacy J17 record path.
Regression controls include valid-outer-identity mutations for component,
molecule, source, family, typed/charge, and assignment nesting. J17 and J18
serialized records remain unchanged.
