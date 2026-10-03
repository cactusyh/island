# Phase 4G2 — source-resolved OPLS-AA single points

## Base and scope

Branch `codex/phase-4g2-oplsaa-parameters-energy` starts at main
`b61c96628d47e93417e809c175e2f77e5a7486e8`. Its tree exactly matches reviewed
4G1 `56a96f13b5036d7a21e1a8c9587160f303679f68` after squash merge. No old
branch was merged again. The separate Phase 4F5 research and historical
artifacts remain untouched.

This implements complete **source-defined** parameter coverage and finite,
nonperiodic, unconstrained NoCutoff single points. Coverage is not scientific
suitability. The source has no explicit improper section; this phase does not
invent impropers or establish that their absence is physically sufficient.
No minimization, dynamics, periodicity, solvent, switching, long-range
corrections, constraints, virtual sites or additional force families are enabled.
`production_validated=False`, `simulation_readiness="not_established"`.

The exact source and acquisition/license policy remain those of
[Phase 4G1](phase_4g1.md): Foyer revision
`dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd`, XML SHA256
`c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`.
No XML or parameter library is redistributed. Runtime never downloads one.
Historical typing/charge schemas, signatures and capability payloads are
unchanged; their old capability declarations describe the original API.

## API and ownership

```python
from island.forcefields.oplsaa import load_oplsaa_source, parameterize_oplsaa
from island.evaluation.oplsaa import OPLSSinglePointEvaluator

source = load_oplsaa_source('/installed/foyer/forcefields/xml/oplsaa.xml')
parameters = parameterize_oplsaa(system, source)
parameters.validate_integrity(system, source)
snapshot = parameters.to_parameterized_system(system, source)
evaluator = OPLSSinglePointEvaluator.from_parameterized_system(snapshot, source)
result = evaluator.evaluate()  # or complete stable-ID coordinates in angstrom
```

`OPLSParameterizationResult` is an owned JSON-envelope record with schema
`island_opls_parameters_v1`. `json_text` can be retained and reconstructed with
`OPLSParameterizationResult(text)`; validate against the supplied authoritative
system and exact source before use. `.payload` returns a separate owned value.
Validation recomputes every inventory, selected source record, conversion,
charge identity and policy. Rechecksummed missing or changed parameters are
rejected with `OPLSAssignmentError`. Checksums establish internal consistency,
not authenticity of an asserted upstream execution.

The explicit `ParameterizedSystem` snapshot carries native external charges,
site/interactions, geometric nonbonded policy and the complete parameter record.
It is mutable caller-owned data: evaluator binding validates actual contents,
including metadata, instead of trusting `aggregate_signature`. The evaluator
then owns a copy. Input systems, metadata, charges and prior evaluations are not
mutated. Binding rejects incomplete records before constructing an OpenMM model.
Coordinate errors use `EvaluationInputError`; backend errors use `EvaluationError`.

Typing/native-charge validation and new assignment require the pinned optional
Foyer stack. A new `resolution_pin.json` additionally pins OpenMM 8.6.1 and the
content of its matching methods. Offline parameter integrity checks require
core Python/NumPy and the source bytes, not Foyer, RDKit, ParmEd or an OpenMM
Context. Evaluation requires OpenMM; the builder example additionally requires
RDKit. Source or implementation changes require an explicitly reviewed new pin.

## Resolution and numerical conventions

The pinned XML has 307 bond, 964 angle, 1089 proper RB and 825 nonbonded records.
The implementation inspected Foyer `Forcefield.createSystem`, OpenMM 8.6.1
`HarmonicBondGenerator.createForce`, `HarmonicAngleGenerator.postprocessSystem`
and `RBTorsionGenerator.createForce`. Source matching supports types, classes,
empty wildcards and reversal. Angles choose the first match; RB selects the
first specific match, otherwise the first wildcard match. Upstream bond
candidates use a set: differing simultaneous matches are rejected as ambiguous;
identical duplicate rows retain all candidate IDs and a deterministic primary.
Every new assignment is checked against the pinned upstream force inventory and
numerical values. This prevents an unverified local matching choice from being
published as successful coverage.

Records carry zero-based section/record identifiers, source attributes,
interaction-order stable IDs, types/classes, raw units and converted values.
Bonds, angles and simple proper paths are regenerated from authoritative bonds;
stale graph-derived caches are ignored. Explicit zero source terms are retained
(the upstream model may omit a zero harmonic term). Missing terms are errors.

- Bonds: `0.5*k*(r-r0)^2`. Source nm → angstrom: `r0 *= 10`, `k /= 100`.
- Angles: `0.5*k*(theta-theta0)^2`. Stored equilibrium degrees are converted
  back to radians at model creation; the source radian spring coefficient stays
  unchanged.
- RB: `sum(c[n]*cos(theta-pi)**n, n=0..5)`, with theta=0 cis, theta=pi trans.
  All six coefficients, including the constant and negative/zero values, remain
  explicit. Reversal preserves the energy. No periodic-series approximation.
- LJ: `4*sqrt(epsilon_i*epsilon_j)*[(sqrt(sigma_i*sigma_j)/r)^12 -
  (sqrt(sigma_i*sigma_j)/r)^6]`.
- Graph shortest distances 1 and 2 are excluded; distance 3 is included once,
  with separate Coulomb and LJ factors 0.5. Rings and multiple proper paths do
  not multiply a pair. Longer-distance pairs are unscaled.

The evaluated model is constructed only from validated ISLAND records. A
standard `NonbondedForce` supplies Coulomb with LJ set to zero. A custom
nonbonded force supplies geometric-sigma/epsilon LJ; distance-3 LJ is a separate
custom bond term with the same mixing. Zero LJ never removes Coulomb.
Each call uses a private fresh Reference Context, installs the full coordinate
frame, requests energy/forces and components, and destroys resources. It never
steps the Integrator. Components are bond, angle, RB proper, Coulomb and LJ.
Outputs are kJ/mol and kJ/(mol*angstrom); OpenMM forces are multiplied by 0.1.

Primary implementation references:
[Foyer forcefield](https://github.com/mosdef-hub/foyer/blob/dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd/foyer/forcefield.py),
[OpenMM forcefield matching](https://github.com/openmm/openmm/blob/8.6.1/wrappers/python/openmm/app/forcefield.py),
[OpenMM RB force](https://docs.openmm.org/latest/api-python/generated/openmm.openmm.RBTorsionForce.html).
The installed source content hashes, rather than a floating documentation page,
are the implementation pin.

## Reproducible acceptance

```sh
python scripts/validate_oplsaa_energy.py --xml /installed/oplsaa.xml --output /new/report
# Inspect the immutable declaration, then execute precisely that plan:
python scripts/validate_oplsaa_energy.py --xml /installed/oplsaa.xml --output /new/report --execute-declared
python examples/oplsaa_singlepoint.py --xml /installed/oplsaa.xml
ISLAND_OPLSAA_XML=/installed/oplsaa.xml python -m pytest tests/test_oplsaa_parameters.py
```

The CLI distinguishes PS's expected native-charge rejection from an unexpected
failure in the other eight cases; any unexpected failure exits nonzero and keeps
diagnostics. Seeds, perturbations, settings and tolerances are published before
execution. The independent reference creates a separate ParmEd chemical graph,
calls pinned `Foyer.apply`, then `ParmEd.createSystem` with geometric combining.
It does not reuse ISLAND's resolution, conversion or pair-construction helpers.
Raw Foyer `createSystem()` alone has arithmetic-sigma ordinary LJ and is **not**
the energy reference. The acceptance code inspects actual custom LJ expression,
particle coefficients and scaled exception sigma values in the ParmEd model.

The first development acceptance retained failures caused by a reference
inspection bug (epsilon/sigma particle fields read in the wrong order). The
corrected run used a new directory, unchanged source, seeds and tolerances.
An RB fallback-precedence error found during source inspection was also fixed
and is covered by the first-wildcard regression. No tolerance was relaxed.

### Executed real-source results

The committed `phase_4g2_acceptance.json` records the declaration and actual
outcomes. Small molecules used coordinate seed 20261003; polymers used local
templates with template/assembly seed 2026. Seed 71823 generated perturbations
of 0.015, 0.03 and 0.05 angstrom (PE50: first only). These are nonsingular,
non-equilibrium frames, not optimized structures.

| Case | Sites | Bonds / angles / RB | max energy difference (kJ/mol) | max force difference (kJ/mol/A) |
|---|---:|---|---:|---:|
| ethane | 8 | 7 / 12 / 9 | 7.11e-15 | 5.68e-14 |
| butane | 14 | 13 / 24 / 27 | 1.42e-14 | 1.71e-13 |
| ethanol | 9 | 8 / 13 / 12 | 3.55e-15 | 5.68e-14 |
| dimethyl ether | 9 | 8 / 13 / 6 | 7.11e-15 | 5.68e-14 |
| benzene | 12 | 12 / 18 / 24 | 1.42e-13 | 3.41e-13 |
| PE DP3 | 20 | 19 / 36 / 45 | 2.84e-14 | 3.41e-13 |
| PEO DP3 | 23 | 22 / 39 / 42 | 2.84e-14 | 1.71e-13 |
| PE DP50 | 302 | 301 / 600 / 891 | 1.16e-10 | 4.66e-10 |

All component comparisons passed. Predetermined total/component and force
absolute tolerance was 1e-5, relative tolerance 2e-10 in the above units. Six
selected force components across butane, ethanol and benzene were checked with
central differences at 1e-4 and 5e-5 angstrom: maximum errors 4.54e-6 and
1.14e-6 kJ/(mol*A), below declared absolute 0.002 plus relative 2e-5 tolerance.
The approximately fourfold reduction supports the expected second-order error.
PE50's complete case, including construction and reference work, took 25.0 s;
this is a bounded smoke check, not a scalability claim.

PS DP3 remains rejected at native charge validation, with component charge
-0.23 e. No parameters or energy were claimed for it. Foyer warns that its
graph-theoretic impropers are unassigned: the pinned XML has no explicit improper
records and only RB propers. This remains a scientific coverage limitation.

Independent software checks include harmonic formulas, a six-coefficient RB
scan and force difference, unequal-sigma mixing, zero-LJ charged particles,
1–4 scaling, ring alternate paths, malformed/rechecksummed records, ownership,
stale caches and matching precedence. These analytical fixtures are explicitly
synthetic; their invented numerical parameters are not source-library evidence.

## Verification and remaining work

Verification counts and affected example outcomes are recorded with final
verification below. Ordinary tests do not require Foyer; the real-source
permutation test is opt-in. Amber archive/QM acceptance was not rerun: no Amber
behavior, parameterization or historical signature was changed.

Before OPLS minimization or dynamics: explicitly integrate bound-system
compatibility, fresh verification/session lifecycle and source-specific restart
identities into those workflows, and validate their numerical behavior. This
phase does not claim those capabilities. Agreement establishes implementation
fidelity to this pinned source, not universal polymer accuracy, thermal
readiness or missing scientific coverage.

### Final verification

- Ordinary suite: **1,035 passed, 10 skipped**. The existing nine opt-in
  archive/live checks plus the new explicitly enabled Foyer permutation test
  account for the skips.
- Pinned Foyer environment: **23 focused tests passed**, including actual
  remapped/reordered ethanol parameterization and energy/force agreement at
  absolute 1e-9 and relative 1e-12. Analytical fixture tests use real OpenMM;
  toy source pin substitution is labelled software evidence.
- Ruff and `python -m pip check`: passed (pip check in both environments).
- Executed examples: new `oplsaa_singlepoint.py`, historical
  `oplsaa_typing_charges.py`, `evaluate_singlepoint.py`, and retained-template
  `generic_fragment_charges.py` (PE cache
  `b679bba61e0263e17f1740f69eb4db29d28f1a14a1e2172ae6d661462e580e1c`).
  The latter required its documented template/cache arguments; its first
  argument-free invocation correctly exited with usage error.
- Actual numerical environment: Python 3.11; Foyer 1.2.0 at the pinned revision,
  OpenMM 8.6.1 Reference, ParmEd 4.3.1, NumPy 2.4.6, GMSO 0.17.0,
  NetworkX 3.6.1, Lark 1.3.1, lxml 6.1.3 and ele 0.2.0. The isolated
  environment received pytest for the opt-in test; project dependencies were
  not downgraded.

Retained developer runs are `island-validation/phase4g2-energy` (failed reference
inspection), `phase4g2-energy-verified` (eight numerical passes and expected PS
rejection), and `phase4g2-energy-final` (same declaration extended with exact
independent atom-type comparison, also passed). These paths are execution
locations, not public API requirements. All seeds and numerical thresholds
remained fixed. The final PE50 case took 90.1 s including redundant JSON decoding
inside its added type-check loop; that local bookkeeping has subsequently been
hoisted out of the loop. The original numerical acceptance took 25.0 s. Neither
measurement is an evaluator-only benchmark or a promised performance bound.
No omitted or failed scientific case was converted into a pass.
