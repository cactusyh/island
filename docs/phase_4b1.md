# Phase 4B1: parameter-result integrity and copying fixes

Phase 4B1 is a focused correction to Phase 4B. It adds no charges, production
force fields, simulation engines, or new interaction forms.

## Reproduced issues and corrections

### Modified selections could retain apparent compatibility

Phase 4B's assignment signature contained graph, typing, library, engine, and
selected parameter IDs. It did not hash the complete `ParameterSelection` wrappers,
selected numerical records, diagnostics, coverage, or completeness state. As a
result, a caller could reconstruct a frozen result with a mapping key for one bond
and `selection.site_ids` for nonexistent sites while retaining apparently valid
input fingerprints.

Phase 4B1 separates two questions:

- `result.is_input_compatible_with(system, typing, library)` checks only whether
  current authoritative inputs match the stored input fingerprints.
- `result.validate_integrity(...)` validates the result's selected output content
  and raises `InvalidParameterAssignmentResultError` when malformed.
- `result.is_compatible_with(...)` now requires both input compatibility and valid
  result integrity. It returns `False`, rather than raising, for malformed result
  content.

Both `result.to_parameterized_system(system)` and
`ParameterizedSystem.from_assignment(system, result)` call integrity validation
before copying anything. They raise `InvalidParameterAssignmentResultError` for a
malformed result, leaving the input untouched.

Integrity validation checks:

- mapping keys, `selection.site_ids`, family, and arity;
- exact authoritative interaction inventory membership;
- selection types against stable-site type assignments;
- record pattern compatibility under complete reversal symmetry;
- wrapper parameter ID/source against the selected record;
- record family and library identity;
- selected record equality with the supplied library during compatibility checks;
- missing/ambiguous diagnostics and candidate IDs;
- coverage counts and completeness state against actual content;
- graph, typing, and library fingerprints;
- the versioned complete result-content signature.

The new `island_parameter_assignment_result_v2` signature hashes all selected
record values and units, wrappers, diagnostics, coverage, declared limitations,
and authoritative input fingerprints. It is deterministic and intended for stale
content detection and reuse checks, not authentication.

### Frozen mapping reconstruction failed

`ParameterAssignmentResult.__post_init__` previously called `deepcopy()` directly
on `MappingProxyType` metadata. `dataclasses.replace(result)` therefore failed
because mapping proxies are not pickleable.

Construction now first materializes ordinary dictionaries, deep-copies their
contents, and wraps assignment maps again. The supported copying contract is:

```python
from copy import deepcopy
from dataclasses import replace

same_content = replace(result)
owned_copy = deepcopy(result)
```

Both operations preserve protected assignment mappings and valid content
signatures. Metadata is deeply isolated from caller-owned inputs and from copied
results. Assignment maps remain read-only through mapping proxies.

### Explicit zero amplitudes were rejected

Lennard-Jones epsilon and periodic torsion amplitudes now accept finite values
greater than or equal to zero:

- `epsilon == 0` is a real LJ record with source and library provenance;
- `force_constant == 0` is a real periodic term inside a nonempty torsion record;
- negative, NaN, and infinite values remain invalid;
- values are never converted to absolute values.

The supported torsion convention remains
`k * [1 + cos(n*phi - phase)]` with non-negative `k`. An absent record remains a
missing parameter and is never replaced by a synthetic zero. Zero values participate
in coverage, library signatures, and result-content signatures; changing zero to a
nonzero value changes those signatures.

Lengths, harmonic force constants, equilibrium bond lengths, and equilibrium
angles retain their Phase 4B constraints. This patch does not generally relax
positive-value validation.

## Completeness terminology

- **Input compatibility** means graph, typing, and library fingerprints still
  match. It says nothing by itself about reconstructed result content.
- **Result integrity** means the output structure, selected records, diagnostics,
  coverage, and v2 signature are internally consistent.
- **Complete supported scope** means every required LJ, harmonic bond, harmonic
  angle, and proper periodic torsion has one selected record.
- **Production validation** remains absent.
- **Simulation readiness** remains `not_established`; charges remain `unassigned`.

These properties are intentionally distinct.

## Coordinate and parameter units

ISLAND coordinates currently use angstroms. Phase 4B parameter lengths use nm,
and angle/torsion records follow the units and formulas documented in
`docs/phase_4b.md`. Future energy evaluators and exporters must perform explicit,
functional-form-aware unit conversion. Phase 4B1 does not rescale coordinates or
silently convert stored parameter values.

## Still deferred

Phase 4B1 does not implement charge assignment, nonbonded mixing policies,
production parameter libraries, PCFF/Class II parsing, impropers, minimization,
MD, LAMMPS export, packing, or crosslinking.
