# Phase 4J8 — authoritative operational source profile and aromatic vertical slice

## Base and purpose

Branch `codex/phase-4j8-authoritative-pcff-profile` starts at merged J7 main
`cb18a400e8e14e9c48e29c87aeab4f4652ce8702`, whose complete tree was verified equal
to reviewed `ef5939e877524e2c53069665f05b94d465dd9280`. No merge or main modification
is performed. Historical records retain their original schemas and identities.

J8 adds explicit operational authorization and completes one bounded aromatic
execution slice. The source remains exactly the established LAMMPS PCFF file:

- Revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
- `tools/msi2lmp/frc_files/pcff.frc`.
- SHA256 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
- Operational profile **`island_pcff_benzenoid_c6h6_v1`**.

The source's historical ingestion-profile name is retained unchanged. This new
operational profile is a separate ISLAND authorization identity. It does not
rename the source or claim universal commercial PCFF equivalence.

The selected domain is one connected, finite, neutral, explicit-H C6 aromatic
cycle with one ordinary H per carbon, standard masses, no isotopes, radicals or
assigned stereochemistry. Recognition uses authoritative graph structure,
bond orders/aromaticity, atomic inventories and chemistry, not molecule names,
coordinates, insertion order or a polymer whitelist.

This domain is chemically meaningful because it exercises active proper torsions,
zero-equilibrium Wilson out-of-plane terms, angle-angle multiplicity, ring paths,
nonbonded interactions and every required Class II coupling. Its complete source
coverage was already established in J1; J8 adds explicit source authorization and
the complete durable workflow slice. It does **not** claim newly discovered atom
types or repaired missing coefficients. Selecting more IFF inventory rows would
not supply the missing organic couplings or charge conventions.

## Operational contract and isolation

`PCFFOperationalProfile` is a JSON-owned, versioned
`island_pcff_operational_source_profile_v1` contract. It binds:

- Exact source hash, source profile and repository/revision/path provenance.
- Authorized labels `cp`, `hc`, and graph typing `island_pcff_source_graph_v1`.
- Existing source-native zero-base direct/ordinary bond-increment policy,
  1e-12 e charge tolerance, and all native connected-component checks.
- Explicit source families and implemented physical-model interpretation.
- Authorized typing, charges, model assembly and evaluation.
- Unsupported domains, automatic lower-order fallback forms, torsion-torsion,
  nonzero Wilson equilibrium and policy-derived parameter zeros.
- Conservative readiness flags and the retained reference evidence.

Only the named, exact registered contract is accepted. A user cannot authorize a
candidate by changing capability booleans, adding labels, replacing the hash,
editing provenance or recomputing an envelope checksum. Public profile validation
also invokes the existing full native model validator. Source rows, dependency
coverage, charges and numerical content are not replaced by `complete=True`.
Every active term must be source-derived, including explicit source zeros.
No policy-derived zero enters this selected slice.

All source files remain external. PCFF-IFF v1.5
`3ad5a1be7334c646ed6cb813b769d0e89a1aa03fe695013e940626b7df922693` and the LUNAR
v1.6 file `66a7798e8f7acd676c29c380ca050ad299a9b7c49c3fecf37821e67219c4e80b` remain
comparison-only. J7 fallback queries and `source_variant_only` results cannot
become an operational selection or a native model. The original native loader's
SHA authorization and all existing native-source safeguards remain intact.

## Public preparation and evaluation

```python
from island.forcefields import (
    ForceFieldRequest, PCFFOptions, prepare_forcefield, create_evaluator,
)
from island.forcefields.pcff import PCFFOperationalSelection

selected = PCFFOperationalSelection(
    "island_pcff_benzenoid_c6h6_v1",
    "e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c",
)
request = ForceFieldRequest("pcff", PCFFOptions(
    source_path="/explicit/local/pcff.frc",
    lj=(0, 0, 1), coulomb=(0, 0, 1),
    typing_profile="island_pcff_source_graph_v1",
    source_profile=selected,
))
prepared = prepare_forcefield(system, request)
evaluator = create_evaluator(system, prepared)
result = evaluator.evaluate_fresh()
```

An invalid selection/options combination fails before native typing. Chemistry
outside the selected graph domain fails before preparation/evaluation. The existing
native validators then check typing, component charges, complete parameters,
model interpretation and source consistency. The evaluator is the established
`PCFFSinglePointEvaluator`; fresh/session behavior, binding, physical constants,
units, numerical formulas and native fingerprints are unchanged.

Unprofiled `PCFFOptions` retains its previous behavior. Its facade remains
`island_prepared_forcefield_v1` with byte-identical description/identity rules.
Explicit profile selections create `island_prepared_forcefield_v2`, carrying the
complete authorization record. Facade identities change to bind that selection;
native model/evaluation identities remain the same for identical native inputs.
`adopt_forcefield(..., pcff_profile=profile)` offers the same guarded adoption path
for existing validated native results without rerunning preparation.

Owned profile access is `prepared.operational_profile`. Its data and all native
accessors are copies. Compatible coordinate replacement remains valid; changing
chemical graphs, masses, stable IDs or stereochemistry still requires a compatible
binding. AtomSite and caller metadata are not mutated to store types or charges.

## Bundles and workflows

A profiled facade uses the distinct
`island_prepared_forcefield_bundle_v2` bundle schema. It adds the canonical logical
artifact `pcff_profile` / `pcff-profile.json`. The original system, assignment and
model files retain their original native formats/signatures. Source resolution
remains explicit `PreparedForceFieldSources(pcff_frc=...)`.

Loading verifies the file manifest, profile's own envelope and exact registered
contract, then reconstructs/validates the native assignment/model and calls the
same guarded adoption path. The restored facade and native identities must match
the manifest. Removing the authorization, changing its source/policy or substituting
source records fails semantic validation even after checksums are recomputed.
Old v1 bundles retain their file sets and reconstruction behavior unchanged.

The existing prepared-bundle workflow consumes these bundles directly. Its
configuration, minimizer, velocity initialization, BAOAB integrator, checkpoint
schemas and semantic manifest validation are unchanged. Start/resume/read/status
all reconstruct the authorized preparation through the common loader. A completed
no-op resume still validates the complete committed state. Failed-attempt binding
and checkpoint authority inherited from corrected I3 remain in place.

No source-aware numerical engine, new optimizer, new thermostat or general plugin
framework is introduced. PCFF-IFF operational model assembly remains unauthorized.

## Declared acceptance and independent reference

[The declaration](evidence/phase_4j8_declaration.json) was published before
implementation. It fixes source/executable hashes, retained input coordinates,
negative controls, geometries, policies, numerical tolerances and all budgets.
The retained system is
`../island-validation/phase4j3-declared/matrix-final/benzene/system.json`, SHA256
`5ccbf3dbd7d5d6d4a0b99d1f1eb309672ab099ad5aa5c09f74e7b060720af8e1`.
Its seed-2026 coordinates are reused without rebuilding or re-embedding.

All required assignments are checked against the independently authored J5 raw
line/Decimal reader and source resolver. Its topology inventory is regenerated
from authoritative bonds independently of ISLAND's inventory. Expected labels are
authored from the domain/source descriptions (`cp` for C, `hc` for H), not derived
from production typing. Every selected row, oriented coefficient and normalized
value is compared. Native charges are independently summed from raw increments.

The executable reference uses the pinned msi2lmp converter's independently
constructed atom/bond/angle/proper/improper inventories and source assignments,
then LAMMPS Class II. It is a validation-only path, not runtime authority. No
`-ignore`, invented coefficient or assumed row is used. H5/J1's separately documented
angle-angle equilibrium-role correction is reconstructed from the reference's
own angle inventory; raw converter output and explicit overrides are retained.
The existing center-preserving AA/source limitations remain documented.

The five numerical configurations are original, asymmetric perturbation, and
trajectory steps 0, 2, 4 (step 0 is the accepted minimized geometry). Perturbation:
`0.015*sin(arange(3*N).reshape(N,3)+0.3)` angstrom. Reference parameters are never
exported from ISLAND's selected coefficient list. ISLAND constructs its own
OpenMM model from validated native records. LAMMPS uses finite all-pair cutoffs
covering the complete system, no periodicity, shifting or tail correction, and
the same independently selected special-pair weights.

Declared tolerances:

- Energy/components atol **1e-5 kJ/mol**.
- Forces atol **1e-5 kJ/(mol*angstrom)**.
- Reference rtol **2e-10**.
- Source coefficient conversion atol **1e-12**, rtol 0.
- Molecular/component charge tolerance **1e-12 e**.
- Split comparison atol **1e-10**, rtol **1e-12** in canonical units.
- Finite-difference displacements **1e-4, 1e-5, 1e-6 angstrom**.

Central finite differences of ISLAND total energy at the predeclared asymmetric geometry are
additional force checks. They do not replace independent source selection or
LAMMPS whole-system comparisons. No tolerances, parameters, seeds or geometry
were changed following a failure.

## Actual vertical-slice results

Preparation, exact source/profile, graph typing, native component charges,
complete source-resolved model, fresh/session execution and numerical comparisons
passed. The model contains 12 bonds, 18 angles, 24 proper torsions, 6 active Wilson
terms, 18 angle-angle couplings and every required adjacent/torsional cross term.
BB13 records are explicit source assignments; none is a missing-row zero.

Maximum measured LAMMPS discrepancy across total/grouped energies and all forces:

- Energy/component: **2.4913404672588513e-13 kJ/mol**.
- Force: **2.0649437372903274e-6 kJ/(mol*angstrom)**.

Finite-difference maximum force errors were **3.40698e-5**, **3.42988e-7** and
**4.63412e-8 kJ/(mol*angstrom)** at the three declared displacements, demonstrating
convergence before small-step rounding becomes relevant.

Minimization used the existing reusable session and generic minimizer:

| Metric | Result |
|---|---:|
| Declared iterations/evaluations | 5000 / 10000 |
| Consumed iterations/evaluations | 19 / 24 |
| Initial/final energy (kJ/mol) | 50.27027960957274 / 37.62672123337065 |
| Final maximum force | 0.09493220143247325 kJ/(mol*angstrom) |
| Final RMS force | 0.04417104527585573 kJ/(mol*angstrom) |
| Termination | force_converged |
| Independent fresh verification | passed |

BAOAB used 300 K, friction 5/ps, timestep 0.1 fs, seeds 78123/99181, four total
steps in two segments of two. Recording every step gives five deduplicated frames.
Each segment consumed its declared four evaluations; the uninterrupted comparison
consumed six. Start constructed five Contexts (minimization plus first segment),
child resume constructed three, and uninterrupted propagation constructed three.
Complete RNG states and trajectory origins match exactly. Coordinates, velocities,
energies, forces, components and times at all shared steps have measured maximum
discrepancy **0** on this host. This is not a cross-platform bitwise guarantee.

The comparison starts from the exact saved minimized coordinates and initialized
velocities. It does not repeat minimization or initialization. The bundle and
workflow directory were relocated before a genuinely separate Python-process
resume. Scientific imports were blocked during bundle reconstruction, then the
existing runtime dependencies were restored for propagation. Resume explicitly
forbids minimization and velocity initialization. No source/QM preparation occurs.

Copies of the real bundle were deliberately corrupted at the profile and model
source fields, with their inner/outer checksums recomputed. Both were rejected by
semantic validation. The original bundle remained byte-identical. Synthetic unit
fixtures are clearly labelled separately and are not real-source acceptance.

## Commands and retained artifacts

From repository root; every experiment requires a new output directory:

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
ROOT=../island-validation/phase4j8-declared
$PY scripts/validate_pcff_profile.py --source "$FRC" \
  --output "$ROOT/acceptance" --require-full-source
# Exit 1 is required while full-source coverage is incomplete.
$PY examples/pcff_operational_profile.py --source "$FRC"
$PY -m pytest -q tests/test_pcff_operational_profile.py \
  tests/test_forcefield_preparation.py tests/test_prepared_bundles.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

The original numerical/workflow report is immutable. Additional real-bundle
semantic-tampering and offline completed-state checks are linked in the final
receipt. The reusable CLI includes the same tamper gate before claiming vertical-
slice acceptance. Full-source success remains a separate, unmet requirement.

Exact commands, process IDs, tool build identity, artifact checksums, parameters,
source rows, forces, finite differences, counters and historical preservation are
in the J8 receipt and the external artifact root. Full source libraries, generated
scientific bundles and executables remain external. LAMMPS executable SHA256 is
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7` and reports
`30 Sep 2026`, `patch_30Sep2026-65-ge891a3e10`. The converter SHA256 is
`a6aa207a4e93f4ad7387a47f0b3230f4bbc6ef69645c98f75d9f10325075fe98`.

## Remaining source/chemical gates

All J6/J7 gaps remain explicit and unchanged:

- No source-backed `nh+` increments in the inspected libraries; pyridinium fails.
- Guanidinium's printed increments total 0.9999 e; no normalization convention
  is established or implemented.
- Sulfur/lactam required cross terms remain absent or depend on unresolved rows.
  IFF ordinary equivalences can resolve lactam base angles, not complete models.
- Conflicting automatic wildcard rows remain ambiguous under reversal; no numeric
  priority is introduced.
- Nonzero Wilson equilibrium signed-permutation semantics remain unsupported.
- Empty torsion-torsion source sections do not establish universal zero physics.
- PCFF-IFF candidate hashes retain no native assembly/evaluation authorization.

Global native graph-predicate coverage stays **86/133** with **47 unresolved labels**.
The new profile covers one explicit graph domain and two atom labels; it does not
complete the full source or authorize unrelated compounds. Periodicity, packing,
curing, crosslinked-network dynamics and new charge models are outside J8.

`production_validated=False`, `simulation_readiness="not_established"` throughout.
A converged local minimum and four-step trajectory establish a bounded software
vertical slice, not equilibration or chemical/property accuracy. No merge occurs.

Final verification: **1,622 passed, 10 skipped** in the ordinary suite
(1323.15 s), **64 passed** in focused source/profile/facade/bundle regressions
(117.12 s), and **2 supplemental controls passed** (0.88 s). Ruff and both
retained environment pip checks pass. The new completed workflow/bundle and six
historical bundles/two historical workflows passed offline inspection with
scientific imports blocked. **1,141 historical file hashes** and the source pin
remain unchanged.

The independently computed final LAMMPS force metrics are fmax
**0.0949322014295659** and atomic-norm RMS **0.04417104521520071**
kJ/(mol*angstrom), satisfying the same declared criterion and agreeing with the
fresh ISLAND result. Atomic RMS is `sqrt(mean_i(||F_i||^2))`.

Receipts: [vertical slice](evidence/phase_4j8.json),
[registered profile](evidence/phase_4j8_profile.json),
[independent source selection](evidence/phase_4j8_selection.json),
[tamper controls](evidence/phase_4j8_tamper.json),
[offline checks](evidence/phase_4j8_offline.json),
[preservation](evidence/phase_4j8_preservation.json), and
[verification](evidence/phase_4j8_verification.json). The standalone example
also ran using the actual pinned source. The reusable CLI's integrated tamper
helper was separately executed on the same retained real bundle after the
numerical/workflow run; no scientific calculation was repeated to add that gate.
