# Phase 4C: charge assignment and explicit nonbonded policies

Phase 4C composes three independent, signed components without moving force-field
identity into the chemical graph:

```text
MolecularSystem chemical graph
    + complete ParameterAssignmentResult
    + complete ChargeAssignmentResult
    + explicit NonbondedPolicy
    -> owned ParameterizedSystem snapshot
```

Chemical formal charge remains `AtomSite.formal_charge`. Partial charges live only
in `ChargeAssignmentResult` and the composed parameterized snapshot. Coordinates
and conformation provenance are not charge or policy inputs in this phase.

## Charge APIs

Public interfaces include:

- `ChargeAssignmentEngine`
- `ProvidedChargeEngine`
- `AtomTypeChargeEngine`
- `AtomTypeChargeEntry` and `AtomTypeChargeTable`
- `ChargeAssignmentResult` and structured diagnostics
- `island_demo_charges_v1()`

All partial charges use `elementary_charge`. Positive, negative, and explicit zero
values are preserved. Booleans, NaN, infinity, missing provided sites, and unknown
provided sites are rejected. Absence is never interpreted as zero.

`ProvidedChargeEngine` takes an exact stable-site mapping and explicit source.
`AtomTypeChargeEngine` consumes a compatible complete `AtomTypingResult` and a
versioned exact type table. Zero matching entries are missing; multiple entries are
ambiguous. It never resolves conflicts by declaration order.

Both engines use stable summation and compare each connected component against the
sum of its authoritative formal charges, then compare the whole system. Thus a +1
and -1 error on disconnected ions cannot cancel globally. Charged molecules are
not forced neutral. Residual charge is never redistributed. The default absolute
tolerance is `1e-6` elementary charge. An explicit target charge must agree with
the graph's total formal charge within that tolerance.
The result records `total_charge_residual` and `total_within_tolerance`; the global
sum is checked even when every component separately meets the tolerance.

Charge results retain method/version, entry provenance, coverage, component and
total diagnostics, and deterministic input/content signatures. Their assignment
maps are read-only, `copy.deepcopy()` is supported, and metadata is caller-isolated.
Graph-only engines declare `requires_coordinates=False`; the abstract contract can
later support conformer-dependent methods without pretending they exist today.
For provided mappings, `ProvidedChargeEngine.input_signature(...)` computes the
current external-input fingerprint; pass it as `input_signature=` to result
compatibility checks when deciding whether a changed caller mapping can reuse an
existing result.

## LJ mixing policy

`NonbondedPolicy` is immutable, versioned, explicitly selected, and records units,
formulas, and source. It supports LJ 12-6 only:

### Lorentz–Berthelot

```text
sigma_ij   = (sigma_i + sigma_j) / 2
epsilon_ij = sqrt(epsilon_i * epsilon_j)
```

### Geometric

```text
sigma_ij   = sqrt(sigma_i * sigma_j)
epsilon_ij = sqrt(epsilon_i * epsilon_j)
```

Epsilon uses kJ/mol and sigma uses nm. Explicit zero epsilon remains zero after
mixing. Pair values are obtained on demand with `mix_lj()` or
`mixed_parameters_for_types()`; ISLAND does not allocate an atom-by-atom N×N
matrix. LJ 9-6, sixth-power mixing, and explicit pair overrides are rejected or
deferred. Policies are never inferred from library names.

LAMMPS documents related engine-level mixing and special-neighbor controls in
[`pair_modify`](https://docs.lammps.org/pair_modify.html) and
[`special_bonds`](https://docs.lammps.org/special_bonds.html). Phase 4C records an
engine-independent policy and does not export LAMMPS input.

## Bonded-neighbor scaling

LJ and Coulomb scale factors are independent for 1–2, 1–3, and 1–4 pairs. Each
factor must be finite and within `[0, 1]`. Classification uses the shortest path in
the authoritative bond graph:

- one bond: 1–2 settings;
- two bonds: 1–3 settings;
- three bonds: 1–4 settings;
- more distant or disconnected: full weight `(1, 1)`.

Shortest paths take precedence in rings, so an adjacent ring pair remains 1–2 even
when a longer alternative path exists. `scaling_for_pair()` queries one stable
unordered pair. `local_pair_scalings()` returns only local non-full exceptions; it
does not enumerate distant pairs. This is an explicit ISLAND policy convention,
not a claim that all force-field variants use the same exclusions.

## Validated composition

Use `compose_parameterized_system(...)` or
`ParameterizedSystem.from_components(...)`. Composition validates Phase 4B result
integrity, charge-result integrity, graph agreement, completeness, and policy type
before copying. The aggregate signature covers the unchanged Phase 4B assignment
signature, charge-result signature, policy signature, and graph signature.

The Phase 4B result is not rewritten to claim it assigned charges. The composed
snapshot reports separate statuses:

- supported parameter coverage complete;
- charge assignment complete;
- nonbonded policy available;
- production validation absent;
- simulation readiness `not_established`.

Synthetic completeness is not scientific validation and is not MD readiness.

## Units and deferred functionality

ISLAND coordinates currently use angstroms. Parameter lengths and LJ sigma use nm.
No conversion or coordinate rescaling occurs in Phase 4C; future evaluators and
exporters must explicitly convert units and functional forms.

Phase 4C does not implement RESP, AM1-BCC, Gasteiger, QM calculations, production
GAFF/GAFF2 or other libraries, PCFF/Class II terms, impropers, energy/force
evaluation, minimization, MD, LAMMPS export, packing, or crosslinking.
