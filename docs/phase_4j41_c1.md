# Phase 4J41-C1 — mixed-silicon coverage and phosphorus state evidence

C1 preserves J41 finite neutral silane and keeps **every additional retained
Si/P fixture blocked**. The pinned source has no mixed-Si bond–bond-1–3 or
end-bond–torsion rows that can complete these fixtures under the declared
strict policy. Phosphorus chemical environments are distinguished, but the
distribution does not establish an unambiguous executable atom-type mapping.
No parameters, charges, zeros or successful dynamics are fabricated.
The tranche acceptance gate remains **1**.
`production_validated=False`; `simulation_readiness="not_established"`.

The base is merged J41, `885fa039d4755b56d2c0aac66633b2e4e8f6d5b0`.
The frozen tested implementation is
`077c48b9ef98d9b571efeccc561e2f5958dbefbf` on
`codex/phase-4j41-c1-mixed-silicon-phosphorus`. Only additive C1 modules,
scripts, tests and evidence are introduced. Historical typing, source resolution,
profile registries and numerical behavior are unchanged. No changes to `main`
or historical documentation/evidence are made.

## Independent source inventory

The ASCII scanner retains source section, namespace, line number, version,
original text and values. For every Si–C/H/O/N/S/F/Cl/Br/I/Si/P environment,
`mixed-source-inventory.json` lists atom definitions, ordinary and automatic
equivalences, bond increments, ordinary and automatic bonds/angles/torsions,
Wilson out-of-plane rows, all Class-II couplings, and nonbonded parameters.
The empty torsion–torsion family remains explicit. Companion nonbonded and
equivalence records are retained; membership in a candidate list never proves
chemical applicability or source completeness. Ordinary and automatic aliases
are followed for one source-declared step, with no transitive chemical fallback.
Wildcard candidates are separate, not evidence for missing exact Si/P rows.

Selected inventory counts below include source versions and aliases, not
deduplicated executable coefficients. Automatic quadratic rows are shown
separately from ordinary quartic rows. A matching lower-order row cannot supply
a missing cross term.

| Source environment | Increments | Quartic bonds | Automatic quadratic bonds | Quartic angles | Bond–bond-1–3 | End-bond–torsion |
|---|---:|---:|---:|---:|---:|---:|
| Si–C | 13 | 3 | 10 | 11 | 0 | 0 |
| Si–H | 4 | 2 | 1 | 12 | 0 | 0 |
| Si–O | 8 | 5 | 3 | 20 | 0 | 0 |
| Si–N | 9 | 0 | 10 | 0 | 0 | 0 |
| Si–S | 3 | 0 | 5 | 0 | 0 | 0 |
| Si–F | 1 | 0 | 1 | 0 | 0 | 0 |
| Si–Cl | 1 | 0 | 2 | 0 | 0 | 0 |
| Si–Br | 1 | 0 | 1 | 0 | 0 | 0 |
| Si–I | 1 | 0 | 1 | 0 | 0 | 0 |
| Si–P | 1 | 0 | 1 | 0 | 0 | 0 |

The audit compares **4,359 numerical catalog records** with raw ASCII values
and executes **347 queries** through a newly compiled pinned msi2lmp C leaf
matcher. Of these, 337 assign existing exact source rows and 10 deliberately
query absent direct mixed-Si keys. The graph audit separately follows ordinary
equivalences before deciding whether a required coupling is missing.
Independent raw matching checks selected
versions, orientations and normalized coefficients. Existing compiled
ordinary-equivalence and conflicting-wildcard controls are also rerun. These
are selection tests, not claims of complete molecular parameterization.

Distribution SHA-256:

| File | SHA-256 |
|---|---|
| `pcff.frc` | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| `pcff.rlb` | `82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36` |
| `pcff_templates.dat` | `dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093` |

The companion RLB contains only `VERSION` / `elib`, not an executable chemical
typing rule set. `msi2lmp` consumes CAR/MDF atom labels and charges; it does not
perform chemical perception, resolve phosphorus oxidation states or assign
native charges. Its successful assembly cannot authorize a guessed input label.
The source code, compiler command/output, executable bytes/hashes and actual
converter outputs are retained in the portable evidence.

## Exact graph assessments and blocked stages

All systems have explicit hydrogens and stable integer atom IDs, seed 2026 and
recorded full coordinates. Construction applies the same declared
`0.015*sin(arange(3*N)+0.3)` Å diagnostic perturbation as J41. The assessments
bind final graph, atom order, coordinates, source distribution, chemical state,
profile, resolution policy, charge origin, native component totals and all
applicable interaction requests. Graph connectivity determines component
membership; source-query status never substitutes for a complete parameter set.

| Fixture / SMILES | Atoms | Result / exact blocking stage |
|---|---:|---|
| Silane `[SiH4]` | 5 | J41 model unchanged; fresh finite reference, minimization and continuation replay passes |
| Methylsilane `C[SiH3]` | 8 | 9 bond–bond-1–3 and 9 end-bond–torsion requests missing; converter exits 17 |
| Disiloxane `[SiH3]O[SiH3]` | 9 | 54 requests missing; converter exits 17 |
| Hexamethyldisiloxane `C[Si](C)(C)O[Si](C)(C)C` | 27 | J41's 378 missing requests reproduced; converter exits 17 |
| Sila-oligomer `C[SiH2]C[SiH2]C` | 17 | 98 missing requests, 1 ambiguous oriented angle request, 6 additional wildcard-only torsions rejected; converter exits 13 |
| Si–N `N[SiH3]`, Si–S `S[SiH3]` | 7 / 6 | Complete chemical type mapping not independently qualified; no guessed converter input |
| Si–F/Cl/Br/I `X[SiH3]` | 5 each | Same typing boundary; increments/automatic rows do not authorize missing Class-II couplings |
| Phosphate ester `COP(=O)(OC)OC` | 17 | P(V) oxo template/source connections conflict |
| Phosphonate `CP(=O)(OC)OC` | 16 | P(V) oxo template/source connections conflict |
| Phosphine `P` | 4 | Trivalent P has no degree-three P source atom type |
| Phosphine oxide `CP(C)(C)=O` | 14 | Distinct oxo P environment; same source connections conflict |
| PF5 `FP(F)(F)(F)F` | 6 | Degree-five all-single-bond P(V) is not an established phosphazene mapping |
| Mixed C/H/O/Si/P `CO[Si](C)(C)OP(=O)(OC)OC` | 27 | P typing blocked; no model, charge vector or dynamics published |

The 54 disiloxane requests are 24 angle–angle plus 6 each of
bond–bond-1–3, end-bond–torsion, middle-bond–torsion, angle–torsion and
angle–angle–torsion. The 378 requests remain 78 angle–angle plus 60 each of
those five torsional coupling families. Exact site tuples, source type tuples,
equivalence paths, missing dependencies and candidate rows are stored for each
request. Converter assembly is executed without `-ignore`; failure output is
retained even when guarded native selection rejects earlier. No incomplete
fixture proceeds to energy evaluation, minimization or MD.

The converter's first failed requests are `hc c3 si hsi` (methylsilane),
`hsi sio osi sio` (disiloxane), and `hc c3 sio c3` (hexamethyldisiloxane),
all in end-bond–torsion assembly. The sila-oligomer fails angle assembly for
`si c2 si`. Native independent fallback inspection exposes conflicting automatic
quadratic-angle rows 1891 and 1892: `* c_ *` has `(109.5°, 60.0)`, whereas
`*1 c_ si_` has `(112.3°, 34.6)` in source real units. Choosing file order or
wildcard specificity would conceal this disagreement; C1 rejects it.

Independent Decimal reconstruction gives component residual `0.0000 e` for
silane, methylsilane, disiloxane, hexamethyldisiloxane and the sila-oligomer.
Native increment charges agree within the unchanged `1e-12 e` bound. These
charge-complete records do not repair missing or ambiguous parameter terms.
The remaining fixtures receive no guessed charge vector.

For neutral oxo P, FRC line 163 describes general `p` with four connections;
line 164 describes phosphazene `p=` with five. The template block starting at
line 1562 requests
`p=` for `P(=O/S/N)(-*)(-*)(-*)`, which has four neighbors and bond-order sum
five. FRC `c=` and `n=` establish that connections elsewhere means neighbor
count, not bond-order sum. Reinterpreting the P field globally is unjustified.
The available `o p` and `o p=` bond parameters differ and are not interchangeable.
C1 reports phosphate, phosphonate, other oxo-P, trivalent phosphine-like and
unsupported connectivity separately. Formal charge, bond orders, neighbors
and conventional formal oxidation numbers are persisted; these descriptions
are not evidence that the source supports an oxidation state. No P type is
selected, and no neighboring oxygen type or charge is guessed.

## Opt-in API and strict publication boundary

The additive API lives in `island.forcefields.pcff.mixed_sip`:
`load_pcff_sip_distribution`, `assess_pcff_mixed_sip`,
`load_pcff_mixed_sip_assessment`, and `prepare_pcff_mixed_sip`.
The distribution pins all three original files. Versioned profiles are:

- `island_pcff_mixed_silane_c1_v1`
- `island_pcff_siloxane_c1_v1`
- `island_pcff_organosilicon_polymer_c1_v1`
- `island_pcff_phosphorus_audit_c1_v1`
- `island_pcff_mixed_sip_audit_c1_v1`

These are explicit assessment profiles, not new default execution permissions.
Only the already qualified finite silane can be prepared through C1, by
delegating to the unchanged J41 preparation contract. Missing/ambiguous
couplings, wildcard-only Si/P terms, source-byte changes, atom order, graph,
coordinates, charge vectors, component totals and policy mutations invalidate
the record. Even hypothetical complete source coverage requires independent
fixture qualification before broader execution. This tranche establishes no
new mixed-Si or P executable profile.

Charge origin remains `PCFF/source_native_bond_increments`; no normalization,
redistribution, provided charges, RESP/Psi4 or AM1-BCC enters this tranche.
An explicit numeric zero in an authoritative row differs from an absent row;
only the former can be used. Native ambiguous candidates remain blocked even
when a converter's incidental file order chooses a row.
The converter's historical non-`cp` bond–bond-1–3 zero convention is not an
authoritative numerical source row and does not satisfy C1's explicit
no-missing-coupling rule.

## Finite control, reproducibility and verification

The finite-silane model identity remains
`a231146e09dbe71bd56baa4e3006cd97410aee8e297b6f8cd69f8dd15021486b`.
The complete J41 study is freshly executed into C1's own `baseline/` directory;
no old receipt is copied as new execution evidence. It repeats independent
coefficient/component/force checks, central differences at `1e-3`, `1e-4` and
`1e-5 Å`, minimization, NVE and BAOAB 200-step finite software controls,
relocation, and fresh-process reconstruction. Reference bounds remain
`atol=1e-5`, `rtol=2e-10`; minimization requires `fmax <= 0.1 kJ/(mol Å)`.
Uninterrupted and split final coordinates, velocities and full RNG match
exactly. Checkpoint bytes match for identical segmented lineage in-process
versus relocated fresh-process continuation. Whole-run checkpoint lineage and
evaluation counters intentionally differ; whole/split checkpoint-file equality
is not claimed.

The converter is rebuilt from the pinned local source (3.9.11 / 6 Sep 2024,
SHA-256 `59f3718d4e537d584011315b76d0bf80ca105b43468f87b384cfa4ce19d76123`).
LAMMPS SHA-256 is
`6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191`.
The actual numerical path remains finite `lj/class2/coul/cut`, Class-II bonded
styles, sixth-power LJ mixing, and `special_bonds lj/coul 0 0 1`.
OpenMM uses Reference double precision; exact versions and hardware are recorded
in the environment artifact. A single silane has no active intramolecular
nonbonded pairs under these exclusions, so intermolecular Si interactions
are not numerically qualified.

Reproduce in the declared supported environment:

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41_c1.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --converter /tmp/j41-c1-validation/converter/msi2lmp \
  --matcher /tmp/j41-c1-validation/matcher/leaf.so \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --output /tmp/j41-c1-validation/final-study
```

The expected study exit is **1**, because source completeness fails for the
additional fixtures. The evidence archives retain compiler commands, native
outputs, raw sources, exact graph/mapping/charge records, independent audits,
environment/package inventories, verification commands and historical hashes.
After extracting `execution.tar.xz`, use `--reconstruct /path/to/study` for
fresh-process validation. The publication validator is:

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j41_c1_receipt.py docs/evidence/phase_4j41_c1
```

Final verification: **37 focused tests**, **1,182 J22–J41 regressions**, and
**3,187 full-suite passes with 10 existing optional skips**. Ruff, compileall,
all five environment pip checks and diff checks pass. Commands, JUnit reports
and exit statuses are retained in the C1 receipt and validation archive.
Eight additional file-level mutations are rejected, including raw source,
source inventory, pressure-state promotion, policy, atom order, component total,
P type selection and compiled coefficient values. The original J41 receipt
validator still passes, and the new C1 receipt/manifest validator passes.
Fresh-process reconstruction from another directory
leaves all 184 retained study files byte-identical.
All 2,503 tracked historical documentation/evidence byte hashes, including
J41's evidence archives, are checked against the base snapshot.

No periodic scientific acceptance changes: `ewald_unvalidated`,
`pressure_unavailable` and long-NPT blockers remain. No universal PCFF,
COMPASS, msi2lmp, mixed-Si/P, charged/ionic/metallic coverage, experimental
agreement, equilibration or production readiness is claimed.
