# Phase 4J12 — guarded native source resolution

## Outcome and boundary

J12 adds an executable, opt-in native resolver policy and validates it through a
PSMILES-built final graph, preparation, evaluation, minimization and durable
workflow. **The requested recovery of the six priority chemical cases is still
unmet.** Missing source rows and a nonconserving printed charge sum prevent those
models. This is not full-source completion or a new chemical typing tranche.

Branch: `codex/phase-4j12-pcff-authoritative-resolution`. Base:
`d00d44cd5ca6fa1be71de1480bbb3d99008a0e04`, fetched `origin/main`. It contains the
reviewed J11 correction `27782f5d7f3d79d13ea749fbc728c50e9011de88`; the only tree
additions relative to that reviewed commit are 20 notebook checkpoint files.
They are preserved. Main is not modified or merged.

`production_validated=False`; `simulation_readiness="not_established"`.

## Immutable source and predeclared experiment

The external FRC is unchanged:

- Repository: `https://github.com/lammps/lammps`.
- Revision: `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
- Path: `tools/msi2lmp/frc_files/pcff.frc`.
- SHA256: `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.

[The declaration](evidence/phase_4j12_declaration.json) preceded production edits.
It binds seven retained J10 graphs by checksum, source, policies, coordinate seed
2026, numerical checks and budgets. The executed positive case is phenylene DP2,
22 explicit sites, constructed from `[*:1]c1ccc([*:2])cc1`. Coordinates and graph
were reused without regeneration. This case was already source-complete under
historical rules; success verifies the new resolver, not expanded chemical scope.

Artifacts are new files under `../island-validation/phase4j12-declared`.
No retained experiment is overwritten. Full source libraries and compiled tools
remain external to Git. The original pre-change obligation audit exits 1 and
reproduces the source limitations before implementation.

## Audited converter behavior versus native authorization

Primary implementation references all use the pinned revision:

| Source routine | Audited behavior | J12 treatment |
|---|---|---|
| [ReadFrcFile.c:76–95](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/ReadFrcFile.c#L76) | Reads the configured base and cross-term tables | Namespace remains explicit; no invented automatic cross-term table |
| [SearchAndFill.c:168–242](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/SearchAndFill.c#L168) | Expands symmetric coefficient halves; substitutes a higher version for an identical ordered key | Preserve supported half expansion and highest decimal version per ordered key; reject conflicting equal-version candidates |
| [GetParameters.c:60–250](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L60) | Direct lookup before ordinary equivalences; separate base/coupling searches; directional bond-angle coefficients | Same tier order, family mapping and directional roles |
| [GetParameters.c:363–525](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L363) | Installs equilibrium dependencies, swaps torsional end blocks; initializes non-cp BB13 to zero | Preserve role/dependency semantics; new strict policy requires a source row instead of that initialization |
| [find_match/match_types:1055–1169](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L1055) | Exact pass before wildcard pass, forward/reverse per row, first matching file row | Same matching tiers, but all winning-tier physical coefficients must agree; incidental file order cannot authorize a conflict |
| [get_equivs/find_equiv_type:1241 onward](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/src/GetParameters.c#L1241) | Family-specific ordinary columns | Same nonbond/bond/angle/torsion/OOP columns; no universal equivalence mapping |

The pinned [README:140–149](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README#L140)
excludes automatic-equivalence supplementation and bond increments. msi2lmp
receives supplied types and charges. It cannot establish automatic typing,
missing increments, or an unprinted guanidinium correction. Automatic positional
lookup and native charges remain explicitly identified ISLAND source policies.

The compiled control extracts the actual pinned C matchers and ordinary
get-equivalence routines. GetParameters.c SHA256 is
`add196c115ed9066aa5eba6da5bd8036bd29a3619395385140276e108c769591`;
Forcefield.h is `4d6c8fa7772d235a3e0fa5ed7398da11b7ced28794170214d754f79d23c3c139`.
`source-audit-pins.json` records the other inspected file hashes and checkout HEAD.
Eight forward/reverse native/C checks and five ordinary-column checks passed.
Automatic-row C controls exercise leaf functions only, not full-converter support.
For conflicting wildcards the native result is ambiguous even when C chooses a
row. Reversing the synthetic table-order control can change C's selection; that
behavior remains comparison-only. No `-ignore` is used.

The known generic-reversal defect for angle-angle centers is not adopted. The
reference retains raw converter output and the independently reconstructed H5
angle-equilibrium-role override in `reference/aa-overrides.json`.

## Public native policy and identity

Select `island_pcff_msi_guarded_source_v3` explicitly:

```python
from island.builders import build_linear_polymer
from island.forcefields import (
    ForceFieldRequest, PCFFOptions, prepare_forcefield, create_evaluator,
)
from island.forcefields.pcff.fallbacks import MSI_POLICY

system = build_linear_polymer("[*:1]c1ccc([*:2])cc1", dp=2, random_seed=2026)
prepared = prepare_forcefield(system, ForceFieldRequest("pcff", PCFFOptions(
    "external/pcff.frc", (0, 0, 1), (0, 0, 1),
    typing_profile="island_pcff_source_graph_v4", resolution_policy=MSI_POLICY,
)))
evaluator = create_evaluator(system, prepared)
```

This uses the existing native pipeline: final graph typing → oriented native
increments → Class II assignment → model validation → facade. No CAR/MDF,
converter, LAMMPS or new QM is required for runtime. CAR/MDF are written only by
the independent acceptance reference. Optional imports remain lazy.

Resolution is deterministic:

1. In the ordinary family, try supplied labels: exact, then wildcard.
2. If missing, use that family's ordinary equivalence column: exact, then wildcard.
3. Only missing supported base terms enter their declared automatic lower-order
   namespace and positional end/apex/center columns. The existing automatic
   nonbond column can query the same-source 9–6 family. No automatic cross terms.
4. Keep the highest version of each identical ordered pattern. At the selected
   exact/wildcard tier, compare all legal oriented coefficients. Differing values
   remain ambiguous, irrespective of specificity, numeric wildcard suffix or
   file order. An ambiguous higher tier never falls through to another tier.
5. Exact/direct success does not consult irrelevant lower-priority equivalence
   conflicts. AA center/shared arm and Wilson center are fixed; proper reversal
   swaps directional coefficient blocks and the existing equilibrium roles.
6. Every active coupling, including BB13, needs a selected source row and all
   equilibrium dependencies. An explicit source zero is retained. No derived
   zero is authorized by this new policy.

Each new assignment retains source lines/IDs, supplied and transformed labels,
legal permutations, equivalence rows, candidate versions, rejection/tier reasons,
selected values and stable sites. Semantic integrity recomputes this trace.
The model's versioned compatibility profile binds the policy and strict BB13
interpretation; source/model/native/facade identities remain distinct.

Historical exact/ordinary and positional v1/v2 policies are unchanged, including
their documented BB13 compatibility zeros. The new policy is additive within
explicitly versioned native records; it does not reinterpret signed old data.
No evaluator equation, charge coefficient, source pin or numerical tolerance was
changed. Shared offline PCFF identity derivation is reused. Bundles and workflows
reconstruct these same native records and run their validators. No new registry
chemical profile, schema migration or alternate numerical engine is introduced.

## Priority cases and independent source checks

The independent Decimal reader selects its own raw rows, permutations and
conversions. The seven-case check compares source selections and partial charge
vectors even on failure; it never feeds incomplete charges into parameterization.
[The gap receipt](evidence/phase_4j12_gaps.json) identifies every request, site,
label sequence, candidate row, search path and equilibrium dependency.

| Retained final graph | Typing | Native charges | New strict model | Remaining blocker |
|---|---|---|---|---|
| Phenylene DP2 | complete | complete | complete | none within this graph |
| PEO DP3 | complete | complete | incomplete | 42 missing BB13 rows |
| Amine chain DP2 | complete | incomplete | not publishable | na–hn2 increments; structural inspection also finds 5 ambiguous base angles and 100 missing coupling requests |
| Thioether DP2 | complete | complete | incomplete | 27 missing BB13 rows |
| Chlorinated chain DP2 | complete | complete | incomplete | 93 missing BB13/BB/BA/AA/EBT/MBT/AT/AAT requests |
| Pyridinium | complete | incomplete | not publishable | 3 bonds lack nh+ increments |
| Guanidinium | complete | incomplete | not publishable | printed 0.9999 e versus +1; 12 BB13 rows also absent under strict source requirements |

Pyridinium's nh+–cp and nh+–hn searches remain absent through direct, ordinary and
automatic columns. No increment contains nh+. Amine na–hn2 searches remain absent;
hn is not substituted. Guanidinium uses line 528 (`c+ nr`, 0.2653, 0.0680) and line
808 (`h* nr`, 0.4068, −0.4068). Independent Decimal summation is exactly 0.9999 e;
−0.0001 e residual fails the unchanged 1e-12 gate. Missing-charge partial totals
are not full molecular charge predictions.

Neither the converter nor the pinned FRC provides the missing scientific data.
Required next evidence is specific: charged-pyridine and amino-H increment rows,
an authoritative guanidinium convention or higher precision data, missing
couplings and reversal-invariant authority for conflicting automatic angles.

| Measure | Before | J12 | Denominator |
|---|---:|---:|---:|
| Declared graph-predicate labels | 86 | 86 | 133 |
| Unresolved labels | 47 | 47 | 133 |
| Typed declared graphs | 7 | 7 | 7 |
| Complete native charge graphs | 4 | 4 | 7 |
| Complete strict source-row models | 1 | 1 | 7 |
| New resolver whole-system numerical/workflow verification | 0 | 1 | 7 |

Historical operational model counts are not strict-source counts: historical PEO
and thioether models can use the explicitly documented non-cp zero convention.
Those records remain valid under their original profiles. They are not evidence
of new source rows. The fresh full row ledger inventories 5,319 versioned records,
133 atom types, 134 ordinary equivalences, 108 automatic equivalences and 564
increments; queryability is not executable coverage. `--require-full-source`
continues to exit **1**.

## Context-sensitive domains

[The context contract](evidence/phase_4j12_contexts.json) lists exact unresolved
source labels/rows and needed information. It is diagnostic, not an authorization
through caller-supplied JSON. Metals require phase, oxidation/charge, coordination
and an applicable source charge model. Ionic labels require explicit component
charge and compensation. Zeolite/surface labels require framework substitution,
bridging/protonation, boundary and coordination roles. Uppercase Br/Cl remain
halogen-family obligations with their ionic source descriptions preserved.

Isotopes require explicit isotope metadata, atomic number, mass and a justified
parent environment; there is no inferred mass substitution. Existing D2O/dw
behavior remains covered by historical tests. Nonzero Wilson row
`wilson_out_of_plane:cff91:3231` (`az ob hb sz`, chi0=3.8934°) is still unsupported
until signed/peripheral semantics are established. Empty torsion–torsion data is
not converted into universal source zeros; existing named scope limitations stay
explicit. No source profile for metals/surfaces/coordination is authorized here.
Native typing and public preparation reject the declared Na+, Cu2+ and Al controls
before parameter assignment or evaluation. Unrelated metadata cannot authorize
an organic fallback.

## Executed real numerical and durability acceptance

The positive graph has 513 independently checked source assignments: 23 bonds,
36 angles, 52 propers, 12 Wilson terms, 36 AA couplings, all required cross terms,
and 22 nonbonded records. All selected coefficients have source rows, including
explicit zeros. The raw reader rebuilds inventories from the graph and resolves
rows independently. LAMMPS receives independently specified cp/hc types and raw
increment charges, not ISLAND-selected coefficients.

LAMMPS reports `30 Sep 2026`, `patch_30Sep2026-65-ge891a3e10`; executable SHA256
`db56822e75ec1f61af453e7727d6de04d351bb91d0465d8e0116011cafa19bd7`.
msi2lmp SHA256 is `a6aa207a4e93f4ad7387a47f0b3230f4bbc6ef69645c98f75d9f10325075fe98`.
Commands, raw converter data, LAMMPS inputs/logs/forces and hashes are retained.

Initial, asymmetrically perturbed, minimized and steps 2/4 comparisons passed
`atol=1e-5` kJ/mol and kJ/(mol·angstrom), `rtol=2e-10`. Maximum grouped-energy error:
**2.8421709430404007e-13 kJ/mol**; maximum force-component error:
**3.885958221871988e-11 kJ/(mol·angstrom)**. No cutoff/shift/tail/periodic changes.
Finite differences at 1e-4, 1e-5 and 1e-6 Å gave maximum force errors
3.84298e-5, 3.84472e-7 and 1.37385e-7 respectively.

Minimization used 78 iterations and 83 evaluations, independently fresh-verified.
Energy: **536.144955544251 → 506.03072523648075 kJ/mol**. Maximum atomic force:
**0.07660584226766147**, RMS force **0.04374344769618155 kJ/(mol·angstrom)**.
Budgets remained 5,000 iterations/10,000 evaluations, force criterion 0.1.

Four BAOAB steps used 0.1 fs, 300 K, friction 5/ps, velocity seed 78123 and
thermostat seed 99181. Two 2-step segments were compared with one 4-step run from
exactly the saved initialization. Bundles and workflow were relocated; a new
Python process reconstructed native records and resumed without preparation,
minimization or velocity initialization. All five common frames had measured
maximum differences **0.0** in coordinates, velocities, energies, forces,
components and time; full PCG64 and origin matched exactly. Declared split
tolerances were 1e-10 absolute/1e-12 relative, not a bitwise portability guarantee.
Uninterrupted evaluations: 6; split cumulative evaluations: 8. Context counts:
start/minimum/first segment 5, uninterrupted 3, child continuation 3.

A separate offline child then inspected completed status, frames and no-op resume
with OpenMM/RDKit/SciPy/Foyer/ParmEd imports blocked. Native identity
`ee8f2666b9b6ef7a43183a3d72b12dac45ddb5e32e3af41ffb5f41725aea4e4c` and facade identity
`16d87f13c6cf4e219289e77b200cddf8cc5703024174b28b90351f8ca3661572` were preserved.
Rechecksummed source/policy tamper controls were rejected semantically.

This checks force-model fidelity and same-integrator continuation consistency.
It is not independent integrator/thermostat validation or chemical accuracy.

## Reproduction

```sh
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
ROOT=../island-validation/phase4j12-declared
$PY scripts/inspect_pcff_msi_resolution.py --source "$FRC" \
  --output /new/native-inspection --require-full-source
# Exit 1; complete diagnostic receipt, incomplete full-source acceptance.
$PY scripts/audit_pcff_c_matcher.py --source "$FRC" \
  --msi2lmp-source ../island-validation/phase4h4-upstream/lammps/tools/msi2lmp/src \
  --output /new/c-native-equivalence --compare-native
$PY scripts/validate_pcff_profile.py --source "$FRC" \
  --declaration docs/evidence/phase_4j12_declaration.json --output /new/vertical-slice
$PY scripts/audit_pcff_source_rows.py --source "$FRC" \
  --output /new/source-rows.json --require-full-source
# Exit 1; no unsupported label has been promoted.
$PY examples/pcff_source_resolution.py --source "$FRC"
$PY -m pytest -q tests/test_pcff*.py tests/test_prepared_workflow.py \
  tests/test_prepared_bundles.py tests/test_psmiles.py
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

Outputs are exclusive. Runtime reconstruction requires the explicitly resolved
FRC and core NumPy; evaluation needs OpenMM, optimization SciPy, initial PSMILES
construction RDKit. The independent oracle additionally uses external compiled
LAMMPS/msi2lmp. No downloads or external converters run during native preparation,
loading or offline inspection.

## Verification and preservation

Focused historical PCFF/prepared/PSMILES tests: **536 passed**, 218.70 s. The final
resolver controls separately passed **28 tests**, including the added direct-tier
precedence control. Synthetic fixtures test versions, exact/wildcard conflicts,
reversal, positional equivalence, repeated labels, explicit zeros, missing BB13,
semantic trace tampering and refusal to adopt incomplete models. Initial synthetic
harness namespace/end-marker mistakes are retained in logs; they are not chemical
or numerical experiments.

Ruff, both environment pip checks and the end-to-end PSMILES example pass.
J8/J9/J10 retained bundles/workflows passed read-only child reconstruction and
completed inspection with unchanged native/facade/profile identities.
**2,692 snapshotted paths (2,354 unique historical files)** are byte-identical.
No historical receipt, source file, workflow, checkpoint or notebook checkpoint
was rewritten. The ordinary suite passed **1,679 tests with 10 skips**, 329.56 s. Final
command/log hashes are recorded in `phase_4j12_verification.json`.

The implemented advancement is a guarded, traceable native resolution policy
with executable end-to-end evidence. No new atom label, charge coefficient or
coupling is invented. The six requested chemical recoveries, nonzero Wilson
semantics and full-source completion remain explicit unmet gates.
