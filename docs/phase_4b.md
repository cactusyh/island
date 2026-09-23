# Phase 4B: typed numerical parameter assignment

Phase 4B assigns an explicitly limited set of numerical records after atom typing:

```text
authoritative chemical graph + compatible complete AtomTypingResult
    -> bond-graph interaction inventory
    -> exact atom-type parameter matching
    -> ParameterAssignmentResult
    -> owned ParameterizedSystem snapshot (only when complete)
```

Coordinates and coordinate provenance are not consulted. Chemical sites and their
metadata are not mutated. Atom typing and parameter assignment remain separate,
and charge assignment is explicitly absent.

## Public API

The RDKit-independent API is available from `island.forcefields` and
`island.forcefields.parameters`:

- `LennardJonesParameter`
- `HarmonicBondParameter`
- `HarmonicAngleParameter`
- `PeriodicTorsionTerm`
- `ProperTorsionParameter`
- `ParameterLibrary`
- `ParameterAssignmentEngine`
- `ParameterAssignmentResult`
- `ParameterSelection`
- `ParameterAssignmentDiagnostic`
- `FamilyCoverage`
- `derive_interaction_inventory()`
- `island_demo_parameters_v1()`

```python
from island.builders import build_linear_polymer
from island.forcefields import (
    ParameterAssignmentEngine,
    RDKitSmartsAtomTypingEngine,
    island_demo_parameters_v1,
    island_demo_v1_ruleset,
)

system = build_linear_polymer("[*]CC[*]", dp=5, generate_3d=False)
ruleset = island_demo_v1_ruleset()
typing = RDKitSmartsAtomTypingEngine().type_system(system, ruleset)
result = ParameterAssignmentEngine().assign(
    system, typing, ruleset, island_demo_parameters_v1()
)
snapshot = result.to_parameterized_system(system)
```

`ParameterizedSystem.from_assignment(system, result)` is the equivalent explicit
container entry point. Both paths copy the molecular system and assignment maps;
later changes to caller-owned coordinates, topology, or dictionaries do not alter
the validated snapshot. Existing direct `ParameterizedSystem(...)` construction
remains available.

## Functional forms and canonical units

Phase 4B accepts only the following forms and exact units. It rejects other unit
spellings or forms rather than guessing conversions.

### Lennard-Jones 12-6

```text
U(r) = 4 epsilon [(sigma/r)^12 - (sigma/r)^6]
epsilon: kJ/mol
sigma: nm
```

`sigma` is the zero-crossing distance. It is not the minimum-energy distance;
the latter is `2^(1/6) sigma`.

### Harmonic bond

```text
U(r) = 1/2 k (r - r0)^2
k: kJ/(mol*nm^2)
r0: nm
```

The factor `1/2` is part of the declared form.

### Harmonic angle

```text
U(theta) = 1/2 k (theta - theta0)^2
k: kJ/(mol*rad^2)
theta0 storage: degree
```

The stored equilibrium angle is in degrees. An evaluator must convert `theta` and
`theta0` to radians before applying the formula. The factor `1/2` is part of the
declared form.

### Proper periodic torsion

```text
U(phi) = sum_i k_i [1 + cos(n_i phi - delta_i)]
k_i: kJ/mol
n_i: positive integer
delta_i storage: degree, in [0, 360)
```

Angles are converted to radians for cosine evaluation. A multi-term torsion is one
parameter record containing multiple `PeriodicTorsionTerm` values. Its terms are
not competing records.

All numerical values must be finite and physically positive where required.
Missing terms are never filled with zero-valued records.

## Inventory and matching

Required bonds, angles, and proper torsions are derived afresh from authoritative
bonds. Cached `Topology.angles` and `Topology.dihedrals` may be empty or stale and
are ignored. Inventory keys use stable site IDs, remove duplicates, and canonicalize
complete reversal. A proper torsion is a simple length-three path containing four
distinct sites, which prevents false torsions in three-member rings.

Parameter matching uses exact atom-type tuples only. Bond `(A, B)` is equivalent
to `(B, A)`; angle `(A, B, C)` to `(C, B, A)`; and proper torsion
`(A, B, C, D)` to `(D, C, B, A)`. There are no wildcards, declaration-order
precedence, or specificity heuristics. Zero matches is missing; more than one
matching record is ambiguous even when its values happen to agree.

Impropers and class-II cross terms are unsupported. A topology containing required
impropers is rejected rather than silently dropping them.

## Input validation and diagnostics

The engine does not trust `AtomTypingResult.complete` by itself. It verifies graph,
ruleset, and combined signatures; exact site-key coverage; key/site identity;
untyped and ambiguous lists; selected and matched rules; and the partition between
surviving and eliminated rules. A stale or internally inconsistent result raises
`InvalidTypingResultError`.

Coverage and missing/ambiguous diagnostics are separate for site LJ, bonds, angles,
and proper torsions. Every unresolved diagnostic contains its family, stable site
IDs, exact type pattern, candidate parameter IDs, and failure reason. Strict mode
raises `IncompleteParameterAssignmentError` with the structured partial result at
`error.result`. `strict=False` returns that partial result without fabricating
parameters.

Results retain graph, typing, mutable typing-content, parameter-library, and final
assignment signatures. `result.is_compatible_with(system, typing, library)` ignores
coordinate-only changes but rejects relevant graph, typing-content, or library
changes. Parameter records, selections, and result mappings form an immutable
validated assignment snapshot.

## Synthetic demonstration library and boundaries

`island_demo_parameters_v1` is labeled:

> **SYNTHETIC SOFTWARE-TEST PARAMETERS — NOT FOR SCIENTIFIC SIMULATION.**

It covers the explicit-hydrogen polyethylene subset of `island_demo_v1`, including
terminal C-H environments. Values exist only to exercise the software. It is not
GAFF, GAFF2, OPLS, PCFF, CVFF, or a validated force field, and it deliberately
reports missing parameters for unsupported chemistry such as oxygen-containing
polymers.

`complete_supported_scope=True` means only that site LJ, harmonic bond, harmonic
angle, and proper periodic torsion records were uniquely selected for the declared
scope. The result separately states:

- `charges_status == "unassigned"`;
- `production_validated is False`;
- `simulation_readiness == "not_established"`.

It does not establish an MD-ready system. Phase 4B does not calculate charges,
import production libraries, assign impropers, mix nonbonded terms, minimize,
simulate, export LAMMPS input, pack systems, or crosslink polymers.

## Recommended next boundary

A later phase can add a separately versioned charge-assignment interface and
production-library adapters, followed by an explicit simulation-preparation
validator. Those steps should consume rather than mutate this result and should
not turn the chemical `MolecularSystem` into the parameter authority.
