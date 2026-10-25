# Phase 4J41-C3 — authorized PCFF/COMPASS source variants

No additional mixed-silicon or phosphorus domain is qualified. The two
available authorized local sources still lack required Class-II couplings or
an independently justified chemical mapping. The acceptance gate remains **1**.
The historical finite neutral PCFF silane remains valid with model identity
`a231146e09dbe71bd56baa4e3006cd97410aee8e297b6f8cd69f8dd15021486b`.

`production_validated=False`; `simulation_readiness="not_established"`.

The base is merged C2 on `origin/main`,
`369759b7a886aa9623652afb3e281771efbb474d`.
The working branch is `codex/phase-4j41-c3-authorized-source-variants`.
The receipt identifies the frozen implementation tested before evidence
publication. Historical receipts, identities and archives are preserved.

## Authorized inventory and source isolation

The audit reads only the user-authorized local LAMMPS `tools/msi2lmp/frc_files`
distribution. It selects these complete source variants:

| Source | Printed version | SHA-256 |
|---|---|---|
| `pcff.frc` | 4.0, 8 October 2013; all earlier declarations retained | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| `compass_published.frc` | 1.1, 30 June 2009; references 1–10 | `28156c9c95346a9efde17d0ae0aff8053410e1de86b05f8d284ef781c0bbadfe` |

The distribution README states that its shipped force-field files are openly
available; it also distinguishes proprietary force fields that are not
shipped. This records local provenance and the authorized audit scope, without
asserting access to any proprietary COMPASS distribution. No source is
retrieved from the network. CFF91, CVFF, ClayFF and OPLS-AA files are inventoried
as directory entries and excluded as different source families.

PCFF retains its own `pcff.rlb` and `pcff_templates.dat`, with the unchanged J41
hashes. COMPASS retains its embedded templates and has no declared companion
files. PCFF's `#define cff91`, connectivity column, automatic equivalences and
terminating `#end` are preserved. The published COMPASS file instead declares
`#define compass`, omits atom connection counts and automatic equivalences,
and ends after reference 10 without `#end`. These differences are recorded;
the original bytes are not adapted into PCFF or silently supplemented.

The additive `sip_variants` API requires an explicit source hash, family,
version/provenance and authorization basis. Each immutable inventory binds the
complete FRC and its own declared companions. It records every numerical row,
atom declaration, ordinary and automatic equivalence, increment, template,
Si/P alias and wildcard candidate, including empty coupling families. Its
identity confers no chemical typing or model-execution authority. Reloading an
inventory reconstructs it from the original bytes; changing a serialized row,
source label, template or authority flag is rejected. Explicitly supplied
additional sources can be inventoried through `load_sip_source_variant` with
their own hash and provenance; they gain no execution authority from ingestion.

## Executed fixture matrix

Each source reconstructs the exact retained final graphs and coordinates for
silane, methylsilane, disiloxane, hexamethyldisiloxane, phosphate ester and the
mixed Si/P probe. A neutral cyclic H-substituted phosphazene,
`[PH2]1=N[PH2]=N[PH2]=N1`, adds a probe of COMPASS's explicitly described P
scope. The retained PCFF assessments are re-created and checked against C2.
All 17 historical C2 decisions are separately reconstructed from the unchanged
C2 raw archive, including candidates outside the new comparison matrix.

| Fixture | Pinned PCFF | Published COMPASS |
|---|---|---|
| Silane | Historical source-complete finite control; converter exit 0 | Source-complete finite diagnostic; converter exit 0; no new execution profile |
| Methylsilane | 18 missing interaction contexts; converter exit 17 | 27 missing contexts; converter exit 17 |
| Disiloxane | 54 missing contexts; converter exit 17 | 48 missing contexts; converter exit 17 |
| Hexamethyldisiloxane | 378 missing contexts; converter exit 17 | 252 missing contexts; converter exit 17 |
| Phosphate ester | Original P typing ambiguity; conversion blocked | No state-specific phosphate mapping; conversion blocked |
| Phosphazene probe | Unqualified P mapping; conversion blocked | Source-supported trial labels; 129 missing contexts; converter exit 15 |
| Mixed Si/P | Original P typing ambiguity; conversion blocked | No complete Si/P mapping or interaction set; conversion blocked |

The raw reconstruction files preserve every missing family, ordered site tuple,
supplied type tuple, search table hash, equivalence path and candidate source
row. The counts above describe required unresolved contexts, not zero-valued
terms. COMPASS's entire `bond-bond_1_3` table is empty. Its methylsilane also
lacks `end_bond-torsion_3` and `angle-torsion_3` contexts. No ordinary bond or
angle supplies a missing cross term. File-order wildcard selection is observed
in the compiled matcher but never grants chemical or completeness authority.

The COMPASS siloxane trial uses the printed `si4`/`si4c`/`o2z` definitions and
ordinary equivalence `si4c -> si4` for bonded families. These labels are scoped
diagnostic candidates, with per-atom source/template evidence. They are not a
new general chemical type-assignment profile. PCFF's existing `si`/`sio` and
hydrogen distinctions remain unchanged.

PCFF's general `p` and phosphazene `p=` declarations retain the original
connection/template conflict. COMPASS's `p4=` template has four neighbors,
while its template comment and reference 4 explicitly describe phosphazenes.
The broad template does not independently establish phosphate, phosphonate,
phosphine-like, charged P or other oxidation-state coverage. Only the neutral
phosphazene probe receives trial labels; missing couplings still block it.
No P parameters are transferred between these states, and no guessed labels
or charges are supplied for the blocked phosphate and mixed Si/P graphs.

## Independent converter and numerical checks

The original local `msi2lmp` C sources are freshly compiled. Source bytes,
compiler command/version, executable bytes/hash, matcher build and execution
commands are retained. The converter consumes CAR/MDF labels and charges;
it does not independently perceive chemical atom types. Full conversion is
attempted only for graphs with source-supported diagnostic labels and complete
neutral increment vectors, without `-ignore`. Partial failure logs remain
failure evidence and do not establish coefficient or force completeness.

The independent ASCII oracle checks **5,319 PCFF rows** and **820 COMPASS rows**,
including **4,923** and **713** numerical rows respectively. It reconstructs the
same-source Si-partner inventories for C/H/O/N/S/halogens/Si/P. The compiled
C leaf matcher executes **374** distinct source searches: **144 assigned** and
**230 missing**, across separately bound source tables. Searches retain
ordinary/automatic paths and exact original values; nothing is mixed between
files. The actual converter outcomes and first missing keys are preserved.

Both finite silane diagnostics compare source-derived masses, atom/type
ordering, charge vectors, interaction counts, nonbonded coefficients and all
applicable Class-II coefficient blocks against compiled converter output.
An independent raw-row energy oracle reconstructs bond, angle, bond–bond,
bond–angle and angle–angle contributions and compares them with LAMMPS.
Central energy differences use `h = 1e-3, 1e-4, 1e-5 Å`; the final maximum force
error is approximately `1.63e-7 kJ/(mol Å)` for each source. The unchanged bounds
are `atol=1e-5` and `rtol=2e-10`. At this graph, both files contain the same
applicable Si–H numerical coefficients; this does not resolve any previously
missing mixed-environment coupling.

The control is isolated and nonperiodic, using Class-II styles,
`lj/class2/coul/cut`, sixth-power mixing and `special_bonds lj/coul 0 0 1`.
There are no active intramolecular nonbonded pairs or proper dihedrals.
The graph-inapplicable degree-four Wilson term is distinguished from a
required missing row. No missing coupling is assigned zero. The new COMPASS
silane diagnostic does not publish a new transferable execution profile.
No new mixed Si/P domain is promoted, and no blocked-domain dynamics runs.

## Reproduction and verification

```bash
PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41_c3.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --converter /tmp/j41-c3-validation/converter/msi2lmp \
  --matcher /tmp/j41-c3-validation/matcher/leaf.so \
  --lammps ../island-validation/phase4j22/lmp-kspace \
  --output /tmp/j41-c3-validation/study
```

A completed study exits **1** for retained scientific/source blockers. A
Python exception is not a completed result. Raw converter inputs/outputs,
source inventories, numerical arrays, search results and decisions must also
pass independent reconstruction. Extract `execution.tar.xz` into a new location
and run `scripts/phase_4j41_c3.py --reconstruct <location>/study`. This checks
both variants, the original C2 study and retained finite-silane checkpoint/RNG
continuation in a fresh process.

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j41_c3_receipt.py docs/evidence/phase_4j41_c3
```

The receipt/manifest validator checks raw archive membership and byte hashes,
reconstructs all decisions and numerical controls, checks JUnit results, and
verifies historical and implementation file hashes. Exact verification
commands, disjoint supported test inventories and results are retained in
`validation.tar.xz`. Notebook checkpoint copies are excluded using pytest's
normal hidden-directory discovery rules.

Periodic `ewald_unvalidated`, `pressure_unavailable` and blocked long-NPT
states remain unchanged. No universal PCFF/COMPASS coverage, ionic/metallic
support, AM1-BCC support, experimental agreement, scientific equilibrium or
production readiness is established.
