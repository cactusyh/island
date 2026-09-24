# Phase 4C1: charge-result and typing-result validation

Phase 4C1 fixes two validation defects in the graph-based Phase 4C workflow. It
does not change chemical topology, coordinates, numerical parameter forms, or
nonbonded policy semantics.

## Charge-result integrity

Reconstructing a valid `ChargeAssignmentResult` with `dataclasses.replace()` could
put `None` in an assignment, diagnostic, or coverage field, or put a string in
`tolerance`. Validation previously dereferenced these values or used them in
arithmetic first, leaking `AttributeError` or `TypeError`.

`validate_integrity()` now checks record types, stable site keys, diagnostic
reasons, coverage structure, component diagnostics, finite scalar fields, and
boolean status fields before evaluating charge sums. It raises
`InvalidChargeAssignmentResultError` with the affected field or site. The existing
`is_compatible_with()` method returns `False` for malformed result content.
`is_input_compatible_with()` remains a separate fingerprint check; matching input
fingerprints do not certify a reconstructed result's integrity.

Both `compose_parameterized_system()` and
`ParameterizedSystem.from_components()` validate the charge result before making
an owned snapshot. Failure leaves the chemical system and caller-owned results
unchanged. Valid incomplete results from diagnostic mode remain inspectable:
their integrity can be valid while `complete` is `False`.

## Shared precomputed typing validation

`ParameterAssignmentEngine` and `AtomTypeChargeEngine` now use the same
RDKit-independent `validate_complete_typing_result()` helper. It verifies graph,
ruleset, and combined signatures; exact stable-site coverage; assignment and
diagnostic record types and site IDs; selected/matched rule references; surviving
and eliminated partitions; and diagnostic consistency. Unknown rule IDs and a
diagnostic whose `site_id` disagrees with its mapping key raise
`InvalidTypingResultError` before either engine returns an assignment.

This is structural validation of a supplied typing result. It does not re-run
SMARTS matching or infer chemistry from coordinates. Library and charge-table
requirements remain checked by their respective engines.

## Compatibility and scope

No public signature schema or engine version changed: signed input and output
content has the same meaning. `validate_complete_typing_result()` is a reusable
public helper in `island.forcefields.typing`; existing builder and assignment
entry points retain their call signatures. Explicit zero partial charges,
per-component and whole-system formal-charge checks, and
`simulation_readiness="not_established"` remain unchanged.

ISLAND still has no production GAFF/GAFF2 library, QM charge calculation,
improper/Class II terms, energy evaluation, minimization, MD, or engine exporter.
