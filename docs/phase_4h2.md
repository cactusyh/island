# Phase 4H2 — bounded automatic PCFF typing

## Base and scope

Branch `codex/phase-4h2-pcff-automatic-typing` starts at reviewed 4H1
`e125512bc4a9fd64fa8a3926f4fbb6b0ca0a8017`. At preparation, `origin/main` was
`98902ebf1f235326632bd7359206efddaa35fdff`; **4H1 remains an unmerged dependency**.
No merge, main modification, source edit, historical-record rewrite, or changes
to unrelated notebook checkpoint directories were made.

This adds the **ISLAND `island_pcff_acyclic_cho_v1` profile**, not universal PCFF
or a claim of equivalence to commercial PCFF typing. It types an authoritative
graph and bridges complete decisions to 4H1's unchanged native bond-increment
resolver. Bonded/Class II assignment and a PCFF `ParameterizedSystem` remain
unimplemented. `production_validated=False`, `simulation_readiness="not_established"`.

Source is unchanged:

- LAMMPS revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
- `tools/msi2lmp/frc_files/pcff.frc`.
- SHA256 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
- [Pinned FRC](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/frc_files/pcff.frc).

The source remains an explicitly supplied local file. Nothing is downloaded at
runtime. See [4H1](phase_4h1.md) for source/license handling and charge semantics.

## Chemical domain and rule table

Support is determined from the final graph, including ends and inter-repeat
bonds: neutral individual C/H/O atoms, saturated single bonds, acyclic connected
components, explicit H, finite nonperiodic atomistic representation. Branching
is supported within that domain; topology need not be a named polymer.

The component exclusion pass precedes **all** matching. Aromaticity,
unsaturation, rings, water, peroxides, carbons with two or more oxygen neighbors
(including acetals/hemiacetals), other elements, charge, radical/isotopic states,
nonstandard masses and incomplete explicit H/valence are diagnosed. An offending
component receives no assignments; an unaffected component can have partial
coverage, but the combined result cannot proceed to charge assignment.

Standard mass evidence is conservative: C 12.011/12.01115, H 1.008/1.00797,
O 15.999/15.9994 dalton, absolute tolerance `1e-6` dalton. This also rejects
isotope masses when an older construction adapter did not preserve isotope
metadata. Other user mass conventions need a separately reviewed profile;
there is no mass replacement. Explicit isotope/radical metadata also excludes
a site. Degree must be C=4, O=2, H=1; hydrogen parents must be C or O.

The following predicates run **only after these guards**:

| Rule ID | Environment | Type | Explicit precedence | Pinned source evidence |
|---|---|---|---|---|
| `C.sp3` | saturated C in the accepted domain | `c` | overridden by H-count specifics | line 67: generic sp3 C |
| `C.H1` | C with 1 H, 3 heavy neighbors | `c1` | overrides `C.sp3` | line 70 |
| `C.H2` | C with 2 H, 2 heavy neighbors | `c2` | overrides `C.sp3` | line 71 |
| `C.H3` | C with 3 H, 1 heavy neighbor | `c3` | overrides `C.sp3` | line 72 |
| `H.carbon` | H attached to supported C | `hc` | unique | line 102 |
| `H.alcohol` | H on ordinary supported C–O–H | `ho` | unique, declared convention below | lines 107, 249, 810 |
| `O.alcohol` | O with one H and one C neighbor | `oh` | unique | lines 156, 298 |
| `O.ether` | O with two C neighbors | `oc` | unique | lines 154, 296 |

`c` is deliberate for quaternary C (and methane), not an unknown-chemistry
fallback. All matching rule IDs, suppressed rule IDs, the selected rule/version,
source atom-type record and neighborhood evidence are retained. Multiple
surviving rules produce an ambiguous diagnostic rather than order-based selection.
The full signed profile lives in `forcefields/pcff/typing_profile.json`, including
predicates, overrides, source lines, mass policy and external audit file hashes.

## Hydroxyl audit and intentional external disagreement

Audited external implementation:
[LUNAR PCFF.py at 67dabeda9e6bd3cc8f968aa0c6a88730886138ef](https://github.com/CMMRLab/LUNAR/blob/67dabeda9e6bd3cc8f968aa0c6a88730886138ef/src/atom_typing/typing/PCFF.py),
plus `typing_functions.py` at that same revision. Both hashes are pinned in the
profile. Files remain external in `../island-validation/phase4h2-upstream`.
The GPL-3.0 external implementation was inspected and optionally executed;
its code is not copied or bundled in ISLAND.

The discrepancy is real:

- FRC line 107 describes `ho` as oxygen-bound hydrogen; line 108 calls `ho2`
  hydroxyl hydrogen. `oh` is oxygen bonded to H; `oc` is ether/acetal sp3 O.
- Ordinary equivalence maps `ho -> h*`, `oh/oc -> o`, but `ho2 -> ho2`.
  Auto-equivalence's **bond-increment** fields make the same mappings. They
  do not turn `ho2` into `ho` or `h*`.
- The FRC has `h*–o` increments (+0.4241, -0.4241 e), line 810. Its `ho2`
  increment rows are `ho2–o_2` and `ho2–oz`, lines 817–818, associated with the
  separately typed ester/carbonate oxygen contexts. There is no `ho2–oh` or
  `ho2–o` increment in this pin.
- LUNAR's earlier rule, lines 513–516, selects `ho2` for O-bound H with a
  second-neighbor C (or several other elements). It precedes the general
  `ho` rule. Ordinary ethanol therefore becomes `oh/ho2` there. The external
  oxygen rules still select `oh` for ordinary alcohol and `oc` for ether.

**ISLAND convention:** use `oh/ho` for ordinary saturated alcohol, based on the
literal `ho`/`oh` source descriptions and their explicit bond equivalences.
Acids, esters, carbonates, multiple-O carbons and their alternate types are
outside this profile. This is a named chemical convention, not a claim that
`ho2` is erroneous in every PCFF implementation. Charge neutrality is not the
justification. The lack of `ho2–o` coverage is recorded, not repaired.

The actual external comparison confirms the declared difference at ethanol's
hydroxyl H and the terminal hydroxyl H in PEO DP1 and DP3. No rename, changed
source row, enabled auto-equivalence fallback, or retrospective change to 4H1's
manual ethanol record is made. A future `ho2` convention would need distinct
chemical justification and parameter coverage; it cannot silently reuse this
profile's charge result.

## API, compatibility and records

```python
from island.forcefields.pcff import (
    load_pcff_source, type_pcff_atoms, assign_automatic_pcff_charges,
    save_pcff_automatic_record, load_pcff_automatic_record,
)

source = load_pcff_source("/external/pcff.frc")
typing = type_pcff_atoms(system, source, profile="island_pcff_acyclic_cho_v1")
print(typing.payload["coverage"], typing.payload["diagnostics"])
if typing.complete:
    result = assign_automatic_pcff_charges(system, typing)
    print(result.complete)  # typing completeness does not guarantee charge coverage
    if result.complete:
        charges = result.charges  # owned external stable-ID mapping, unit e
    save_pcff_automatic_record(result, "automatic-charges.json")
    restored = load_pcff_automatic_record("automatic-charges.json", source, system=system)
```

Schemas are new: `island_pcff_automatic_typing_v1` and
`island_pcff_automatic_charges_v1`. Historical manual schemas/signatures are
unchanged. Automatic typing binds the exact source, entire profile/rule content,
canonical chemical graph, per-site decisions, coverage and diagnostics. Integrity
recomputes decisions and the native charge result, rejecting contradictory
assignments/evidence even with a new envelope checksum. Data access returns owned
copies. These are internal consistency checks, not computational authentication.

The chemical graph includes stable IDs, element/atomic number, formal charge,
mass, bond endpoints/order/aromaticity and normalized atom metadata:
`aromatic`, `isotope`, `radical_electrons`, `chiral_tag`, `cip_label`. Assigned
stereo metadata is preserved and bound; stereo does not change these scalar
charge rules. No stereo assignment or correctness certification is performed.
Coordinates, atom names, residue names, repeat indices, builder seeds, coordinate
history, unrelated metadata and stale derived angle/torsion caches are not inputs.
Absent metadata uses the profile's documented neutral/unspecified defaults.
The actual bonded H inventory, not a cached H-count annotation, controls matching.

The bridge builds an owned **canonical chemical snapshot** and invokes the
unchanged explicit-type/native-charge resolver on it. Its provenance explicitly
names automatic typing and embeds that record's identity. A new wrapper retains
the full automatic evidence and existing native record; it does not pretend the
labels were manually validated. Wrapper compatibility checks use the original
chemical graph definition. The embedded historical-format record is bound to the
canonical snapshot, not to unrelated caller metadata.

`typing.assignments` rejects incomplete results. Diagnostic payloads retain
partial assignments and every site's status. The charge bridge rejects incomplete
typing; complete typing with missing increments returns an incomplete charge
result whose `.charges` accessor rejects use. Component and molecular checks,
`1e-12 e` tolerance, zero base charges, oriented contributions and blocked
4H1 auto-equivalence remain unchanged. There is no normalization or fallback.

Malformed structures, invalid numeric/field types, source/profile mismatch and
contradictory records raise `PCFFError`; well-formed unsupported chemistry returns
a validated diagnostic record. Saving uses strict checksummed data-only JSON,
exclusive atomic publication and validation before writing. Loading requires the
verified source and optionally the external system. No executable objects or
live backends are serialized.

Typing, charging and offline validation require core ISLAND/NumPy only. Tests
block RDKit, OpenMM, ParmEd, Foyer and SciPy during reconstruction. Construction
and the separate reference oracle can require RDKit. The SMILES CLI uses the
existing seeded `from_smiles()` embedding API during **construction**; the typer
never embeds. The PSMILES example uses the existing builder's 2D construction.
An already authoritative graph may have no coordinates at all.

## Commands and actual acceptance

```sh
# Unchanged source inspector:
python -m island.forcefields.pcff /external/pcff.frc

# Construct, type and charge, with separate structured coverage fields:
python -m island.forcefields.pcff.typing_cli --source /external/pcff.frc --smiles CCO
python -m island.forcefields.pcff.typing_cli --source /external/pcff.frc --psmiles '[*:1]CCO[*:2]' --dp 3 --output /new/pcff-result
python examples/pcff_automatic_charges.py --source /external/pcff.frc --dp 3

# Exclusive experiment; declaration written before any case executes:
python scripts/validate_pcff_automatic.py --source /external/pcff.frc --lunar-dir /external/pinned-lunar-files --output /new/acceptance
```

The typing CLI exits nonzero for incomplete typing or charge coverage and refuses
an existing output directory. The acceptance CLI distinguishes required supported
cases from expected unsupported diagnostics, retains failures and exits nonzero
for unmet gates. External comparison is explicitly requested with `--lunar-dir`;
without it, that comparison is reported unavailable rather than implied executed.
No runtime source download occurs.

Executed output: `../island-validation/phase4h2-acceptance`. Declaration, every
typing/charge record, all external assignments/information/tallies (including the
external assumed-type table), unsupported diagnostics and report are retained.
Hashes and summarized outcomes: [evidence/phase_4h2.json](evidence/phase_4h2.json).

The predeclared oracle combines manually authored per-case type counts with
separately authored RDKit SMARTS for the expected environments. Expected labels
are not obtained from the production typer. Native charges use independent
Decimal summation of audited FRC lines 493, 505, 519 and 810 and fixed ordinary
bond-equivalence rows, without production selection/summation helpers. Tolerance
was fixed at `1e-12 e` for both charge comparisons and totals before execution.
SMILES construction used seed 2026; polymer construction used existing 2D mode.

| Required case | Sites | Main environments | Typing + charge result | External difference |
|---|---:|---|---|---|
| Ethane | 8 | methyl | passed | none |
| Propane | 11 | methyl/methylene | passed | none |
| Isobutane | 14 | tertiary `c1` | passed | none |
| Neopentane | 17 | quaternary `c` | passed | none |
| Ethanol | 9 | alcohol | passed | hydroxyl `ho` / `ho2` |
| Dimethyl ether | 9 | ether | passed | none |
| PE DP1 | 8 | actual H-terminated ends | passed | none |
| PE DP3 | 20 | end/interior C | passed | none |
| PEO DP1 | 9 | actual terminal OH | passed | terminal H `ho` / `ho2` |
| PEO DP3 | 23 | internal ether, terminal OH | passed | terminal H `ho` / `ho2` |
| PE DP50 | 302 | 2 methyl, 98 methylene, 202 H | passed | none |

All five declared unsupported graphs (cyclohexane, water, peroxide, acetal and
unsaturated aldehyde) returned incomplete diagnostics as expected. The maximum
per-site charge discrepancy was **1.1102230246251565e-16 e**; maximum molecular
residual was **4.85722573273506e-17 e**. PE DP50 typing, including its public
integrity validation, took **0.20664 s** in this run. This is a smoke timing,
not a repeated benchmark; it excludes construction, charging and reference work.

External execution called the exact pinned `PCFF.nta` through an independently
constructed BFS neighborhood adapter for these acyclic graphs, not the full LUNAR
bond/ring-perception pipeline. All 11 cases reported zero assumed and zero failed
assignments. The three expected `ho/ho2` differences were retained; all other
labels agreed. Upstream agreement is profile comparison, not a typing truth set
or scientific accuracy validation. No external charge repair or charge resolver
was used to hide the hydroxyl discrepancy.

Environment: Python 3.11.16, NumPy 2.4.6, RDKit 2026.03.6. No new QM,
parameterization, energy evaluation or dynamics ran. All five files from the
4H1 acceptance manifest retained their hashes; its three charge records reloaded
with their exact original identities.

## Verification and remaining work

Focused regressions cover specificity, quaternary/tertiary C, actual chain ends,
ID remapping and orientation, coordinate/repeat-metadata independence, ownership,
isotope/radical/unsupported chemistry, malformed inputs, partial charge coverage,
profile/source mismatch, rechecksummed tampering, strict persistence, optional
import isolation and CLI gates. Synthetic source fixtures are explicitly labeled;
they are separate from the real-source execution above.

The supported coverage gates are met under the declared ISLAND convention.
Unmet scientific/extension gates include an independently curated chemical typing
truth set, resolution of the external hydroxyl convention for broader chemistry,
ring/unsaturation/other-element rules, auto-equivalence-only charges, full Class II
parameter coverage and an independently validated PCFF energy model. This phase
does not enable simulation or claim compatibility with all commercial PCFF releases.

### Final verification

- Full ordinary suite: **1,163 passed, 10 skipped**, 96.08 s. Existing opt-in
  backend/archive skips remain; no automatic-PCFF regression was skipped.
- Focused automatic-PCFF suite: **28 passed**, 8.08 s, including synthetic
  missing-parameter/conflicting-rule injections and real builder graph inputs.
- Ruff and pip check: passed.
- Automatic PEO DP3 example and ethanol typing CLI: executed successfully.
- Real pinned-source acceptance: 11 supported passes, five expected rejections;
  exact pinned external typing executed on all 11 supported cases.
- Historical 4H1 artifact hashes and identities: unchanged.

Executed commands used `../island-validation/phase4e2-env/bin/python`:

```sh
../island-validation/phase4e2-env/bin/python -m pytest -q
ruff check .
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4e2-env/bin/python examples/pcff_automatic_charges.py --source ../island-validation/phase4h1-sources/lammps/pcff.frc --dp 3
../island-validation/phase4e2-env/bin/python -m island.forcefields.pcff.typing_cli --source ../island-validation/phase4h1-sources/lammps/pcff.frc --smiles CCO
../island-validation/phase4e2-env/bin/python scripts/validate_pcff_automatic.py --source ../island-validation/phase4h1-sources/lammps/pcff.frc --lunar-dir ../island-validation/phase4h2-upstream --output ../island-validation/phase4h2-acceptance
```

Use a new exclusive output directory for repeat execution. The experiment's
source, rules, expected labels, seed and tolerances were not adjusted after its
results to obtain agreement.
