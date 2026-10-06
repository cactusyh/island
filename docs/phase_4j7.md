# Phase 4J7 — explicit PCFF source variants

## Repository and scope

Branch `codex/phase-4j7-pcff-source-variants` starts from merged J6
`origin/main` (`13c096d`; full SHA in the declaration). Its complete tree was
verified identical to reviewed `5ccf0cd00a2f3c858ee3323cd21055c8966e23e3`.
No prerequisite was merged and main was not modified.

J7 adds owned, versioned source comparisons, explicit hash/profile selection,
source-bound fallback **queries**, physical coefficient reversal checks, and
cross-source diagnostics for the six J6 targets. Native assignment remains bound
to the audited source. A comparison-only candidate is not an operational profile.
No equations, parameters, charges, historical schemas or readiness flags change.

The declaration preceded implementation:
[phase_4j7_declaration.json](evidence/phase_4j7_declaration.json).
It fixes source bytes/provenance, retained systems, type labels, the J2 resolver,
controls, child-process budgets and acceptance gates. The separate
[conversion declaration](evidence/phase_4j7_conversion_declaration.json) preceded
numerical coefficient comparison. It permits 1e-12 absolute conversion roundoff
between `degrees*pi/180` and `degrees*(pi/180)`; source selection is exact and
physical reversal compares role-transformed coefficients exactly.

## Sources and acquisition

Four declared files have three distinct byte hashes. Sources are acquired only
by explicit validation setup, then supplied locally. Runtime APIs perform no
network access, source search, installation, converter call or force evaluation.

| File | Immutable repository revision | SHA256 |
|---|---|---|
| LAMMPS `tools/msi2lmp/frc_files/pcff.frc` | `e891a3e10973c1a729e391a0aefaa02fd70f8c0f` | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| LUNAR `frc_files/pcff.frc` | `67dabeda9e6bd3cc8f968aa0c6a88730886138ef` | same as LAMMPS |
| Heinz `pcff_interface_v1_5.frc` | `584179265906d93aa40f7d5b275143985871b654` | `3ad5a1be7334c646ed6cb813b769d0e89a1aa03fe695013e940626b7df922693` |
| LUNAR `frc_files/pcff_interface_v1_6mBN.frc` | `67dabeda9e6bd3cc8f968aa0c6a88730886138ef` | `66a7798e8f7acd676c29c380ca050ad299a9b7c49c3fecf37821e67219c4e80b` |

The [pinned Heinz repository README](https://github.com/hendrikheinz/INTERFACE-force-field-and-surface-models/blob/584179265906d93aa40f7d5b275143985871b654/README.md)
identifies PCFF-INTERFACE and its extensions, citing Heinz et al., Langmuir 2013,
29, 1754. The [author's distribution page](https://bionanostructures.com/interface-md/)
is additional acquisition context. PCFF-IFF is a related, explicitly separate
source family; this does not establish equivalence to commercial PCFF releases.
The v1.6 name is the supplied filename, not an invented semantic version. The
record preserves actual `#version` declarations, including the older organic
parameter history and IFF additions.

The pinned LUNAR source documents GPL-3 distribution. The author repository has
no separate LICENSE file at the pinned revision. Redistribution of full libraries
is not assumed: full files and generated full-catalog comparison records remain
external. Git contains synthetic software fixtures, manifests, hashes and compact
receipts. The pinned LUNAR path `pcff_iff_v1_5_CNT_poly_solv.frc` was unavailable
(HTTP 404); the repository's actual file listing identified the separately named
v1.6 file. It did not replace the original PCFF source.

## Inventory and comparison

`load_pcff_source_variant(path, expected_sha256=..., provenance=...)` owns immutable
bytes and encoded provenance. Required fields are repository, revision, source
path, family, license handling and citation. Every supplied file is represented,
including byte-identical mirrors with distinct acquisition provenance.

`compare_pcff_sources(variants)` produces
`island_pcff_source_comparison_v1`. Its inventories preserve:

- Exact labels, atom descriptions, ordinary and automatic family/position maps.
- Oriented increment endpoint values, reference numbers and decimal versions.
- All supported numerical Class II and automatic record forms, original values,
  unit conversions, source row IDs, source lines, duplicate candidates and namespaces.
- Empty sections and raw unsupported sections/reference notes.
- Canonical content differences between every pair, independent of source line
  movement and numeric formatting; full-source hashes still distinguish byte changes.

Known catalog construction is factored into one data interpreter used by both
paths. The original `inspect_pcff_full_source()` still requires the original
native source authorization. Candidate inventory does not weaken this guard.
Unknown source features remain raw or explicitly unsupported query semantics.

| Source | Atom labels/records | Ordinary equivalences | Automatic equivalences | Bond increments |
|---|---:|---:|---:|---:|
| LAMMPS/LUNAR PCFF | 133 | 134 | 108 | 564 |
| PCFF-IFF v1.5 | 214 | 215 | 189 | 566 |
| LUNAR PCFF-IFF v1.6 | 240 | 241 | 193 | 566 |

These are source inventory counts, not chemical recognition or executable-model
counts. The original 86/133 graph-predicate count and 47 unresolved labels remain
unchanged. Added IFF labels are not silently promoted to ISLAND automatic typing.
All family counts and pairwise row differences are in the separate J7 receipt.

## Strict source selection and fallback

```python
from island.forcefields.pcff import (
    PCFFVariantSelection, PCFFVariantPolicy,
    resolve_pcff_variant_record, select_native_pcff_source,
)

selected = PCFFVariantSelection(expected_sha256, declared_profile)
policy = PCFFVariantPolicy(primary=selected)  # no other source searched
query = resolve_pcff_variant_record(
    variants, policy=policy, family="bond_increments", namespace="cff91_auto",
    types=["cp", "nh+"], sites=[11, 29],
)
print(query.payload["classification"])
```

An opt-in fallback is a tuple of separately specified `PCFFVariantSelection`
objects. Every declared hash/profile must be available and verified, including
unused fallback entries. Fallback is attempted only for an absent row; an
ambiguous primary source is not silently bypassed. Successful alternate queries
retain `source_variant_only`, every attempted source, original ordered types/sites,
equivalence evidence, candidate IDs, legal permutations and selected coefficients.
All candidate data stay bound to one selected source. No mixed-source native model
is assembled; query records carry `model_assembly_authorized=False`.

`select_native_pcff_source(variants, selection)` permits the existing native
preparation path only for the audited bytes/profile. It returns the original
`PCFFSource` identity unchanged, including for a verified equal-byte mirror.
Other hashes fail with the explicit missing native-semantics authorization.
They are inspectable candidates, not successful native preparations. Existing
`prepare_forcefield`, bundles, evaluator/session and workflow contracts retain
that native validation boundary.

Profiles are `lammps_pcff_explicit_cho_v1` for the exact established hash and
`island_pcff_frc_candidate_v1` for comparison-only variants. The latter names an
ISLAND inspection contract, not scientific compatibility or force-field version.
Source fallback cannot graft one successful row into an existing signed model.

## Reversal and source-state contracts

The new query contract reuses J2's exact/ordinary/automatic positional resolution
and highest-version/conflicting-candidate policy. It does not assign numeric
priority to wildcard suffixes. Forward and legally reversed queries are both
resolved; directional coefficient blocks are transformed back to original roles
before exact coefficient comparison. Angle centers, improper centers and shared
angle-angle arms retain their required positions.

Conflicting equally specific wildcard choices or inconsistent reversed values
are `wildcard_ambiguous` and cannot pass. Synthetic controls use asymmetric
bond-angle blocks, opposing wildcard endpoint patterns, identical-coefficient
positive controls and illegal center movement. Nonzero Wilson equilibrium
retains `model_incomplete`: coefficient presence cannot establish its unresolved
signed-angle/permutation interpretation. Missing rows do not pass a physical
coefficient-invariance gate merely because both searches found nothing.

`island_pcff_variant_resolution_v1` and
`island_pcff_source_variant_audit_v1` identify the requested states:

| State | Meaning |
|---|---|
| `source_row_present` | A supported, source-bound lookup has consistent coefficients; no model completeness claim |
| `source_row_absent` | No row after declared searches; structurally inapplicable requests are separately flagged |
| `source_variant_only` | Original source lacks the row but an explicitly identified other source contains it |
| `wildcard_ambiguous` | Conflicting selection or failed physical reversal check |
| `equilibrium_dependency_missing` | A coupling has a row but lacks a required source-bound equilibrium value |
| `policy_derived_zero` | Existing named original-profile rule, with dependencies; never a source row |
| `native_charge_incomplete` | Missing/conflicting increments, incompatible atom inventory or failed component formal-charge check |
| `model_incomplete` | Native model or interpretation gate unmet |

The original non-`cp` BB13 zero is applied only under its existing pinned profile.
It is not inherited by IFF candidates. More apparent gaps in an IFF summary can
therefore reflect an unadopted compatibility policy, not different raw rows.
No website recommendation to use converter `-ignore` establishes applicability
for these organic targets. That option is not used as acceptance or runtime authority.

## Cross-source audits and independent verification

`audit_pcff_source_variants(original_assessment, variants)` consumes an integrity-
checked J6 assessment. It checks every original ordered interaction, recomputes
equilibrium dependencies separately in each source, preserves the original
chemical graph/native identities and compares bond increments without alias
substitution or normalization.

Atom-label element/connectivity evidence is retained. Cross-variant use of original
PCFF labels does not establish equivalent chemical typing rules for an alternate
library. Charge results are explicitly comparison calculations; candidates cannot
be represented as validated native charge/preparation results.

The acceptance script independently scans raw source bytes with a separately
authored line/Decimal reader. It checks every source row and section, then uses
the independently authored J5 raw resolver to compare row selections, candidate
sets, coefficient orientation and normalized values. Its Decimal increment
summation uses authoritative bonds and connected components. It does not use
production selection, catalog, charge or inventory helpers as its oracle.

No new energy equations or complete target model are introduced. The unchanged
J6 isolated LAMMPS/finite-difference evidence is retained; it does not solve
wildcard priority. Whole-system force, minimization and workflow acceptance for
these incomplete models remains unmet. No QM or additional trajectory is run.

## Persistence, relocation and integrity

Save/load APIs are explicit:

- `save_pcff_source_comparison` / `load_pcff_source_comparison`.
- `save_pcff_variant_resolution` / `load_pcff_variant_resolution`.
- `save_pcff_source_variant_audit` / `load_pcff_source_variant_audit`.

All use strict data-only envelopes and existing exclusive atomic publication.
Loaders require explicit verified external variant objects, plus the original
assessment for an audit. They recompute inventories, differences, source decisions,
reversal outcomes, charges, dependencies and gates. Rechecksummed contradictions
are rejected with `PCFFError`; checks establish consistency, not authentication
of caller-provided scientific provenance. Candidate source files are never
selected from paths embedded in historical records.

The acceptance relocates saved comparison/audit/system/assessment records and
loads them in a genuinely separate process. Source libraries are supplied by an
explicit current declaration, remain external and are hash-checked. Scientific
imports (OpenMM, RDKit, SciPy, ParmEd, Foyer) are blocked during reconstruction.
No typing execution, parameterization, Context, force evaluation or converter is
needed. Original assessment bytes and systems are copied without re-signing.

## Reproduction

From repository root, use new output directories on every repeat:

```sh
PY=../island-validation/phase4e2-env/bin/python
ROOT=../island-validation/phase4j7-declared
DECL=docs/evidence/phase_4j7_declaration.json
$PY -m island.forcefields.pcff.variants_cli --declaration "$DECL" \
  --output "$ROOT/inspection" --require-full-source
# Expected exit 1: an inventory does not establish full-source completion.
$PY examples/pcff_source_variants.py --sources "$DECL" \
  --sha256 e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c \
  --profile lammps_pcff_explicit_cho_v1 --types cp nh+ --output "$ROOT/example"
$PY scripts/validate_pcff_variants.py --declaration "$DECL" \
  --tolerances docs/evidence/phase_4j7_conversion_declaration.json \
  --assessments ../island-validation/phase4j6-declared/assessment/relocated \
  --output "$ROOT/acceptance" --require-full-source
# Expected exit 1 even when the separately reported comparison gate succeeds.
$PY -m pytest -q
ruff check .
$PY -m pip check
../island-validation/phase4g1-env/bin/python -m pip check
```

The generic CLI's exit 0 means validated source inspection/query publication.
The bounded acceptance separately requires all four declared files, all six
required diagnoses and successful child reconstruction. Neither gate can make
the full-source or operational-model gate true. Errors and incomplete cases
remain in new receipts; prior J1–J6 records and directories are preserved.

Exact executed results, source differences, blocker rows, reconstruction receipts,
commands, artifact hashes and preserved historical identities appear in the
separate J7 evidence files. All readiness flags remain
`production_validated=False`, `simulation_readiness="not_established"`.

## Actual findings and correction receipt

Both IFF variants provide a different **ordinary equivalence** route for the
lactam angle. Three-membered `c3h–c_1–n3n` and four-membered `c4h–c_1–n4n` map to
`c–c_1–n`. PCFF-IFF v1.5 selects quartic-angle row **2583**; the LUNAR v1.6 file
selects row **2640**. Source values are **116.9257 degrees, 39.4193, −10.9945,
−8.7733** (kcal-based angular polynomial coefficients). The c3h/c4h and n3n/n4n
ordinary equivalence rows explicitly supply these mappings. This is a different
source policy as printed in that file, not permission to change the old library.
The raw c–c_1–n parameter is also available in the old catalog; its old ring-type
equivalences do not reach it. The improvement is an alternative source's lookup
path, not an invented new interaction coefficient.

This clears the two base-angle ambiguities and their 8/12 dependent equilibrium
requests under the alternate queries. Other required couplings remain absent.
Direct automatic quadratic-angle requests still contain the conflicting wildcard
patterns. The separate predeclared wildcard controls query all four environments
in every distinct source and execute the unchanged pinned LUNAR matcher. An
ordinary route resolving a base angle does not establish wildcard priority.

| Target | Original pinned structural blockers | IFF v1.5 / v1.6 structural blockers | Charge result in all three sources |
|---|---:|---:|---|
| Thioformaldehyde | 9 | 9 / 9 | complete zero total |
| Thioacetone | 33 | 45 / 45 | complete zero total |
| Pyridinium | 0 | 0 / 0 | three missing nh+ increments; total remains partial |
| Guanidinium | 0 | 12 / 12 | exactly 0.9999 e, residual −0.0001 e |
| Alpha lactam | 14 | 18 / 18 | complete zero total |
| Beta lactam | 20 | 32 / 32 | complete zero total |

IFF counts include source-missing BB13 rows for which its operational compatibility
policy has not been registered. The original 62 conditional non-cp policy zeros
remain available only under the established original profile. All six models
remain incomplete in every inspected source. These findings are bounded to these
exact files and searches, not all commercial variants or unpublished conventions.

The initial reconstruction hit the declared **600-second** timeout. Its original
records, manifest, inputs and failure are retained at `acceptance/`. The driver
had fully validated each immutable record through its native loader, then invoked
its public `identity` property and replayed the same expensive validation.
The correction computes content identity directly from the immutable result
**just returned by the validated loader**. Public validators, identity formulas,
records and scientific settings are unchanged. Software regressions verify loader
invocation and failure propagation; no check is replaced by a shallow checksum.

A separate [reconstruction declaration](evidence/phase_4j7_reconstruction_declaration.json)
fixes byte-identical copied records and the same 600-second budget. Corrected
separate-process reconstruction passed in **269.525 s**, with scientific imports
blocked and every original/copy hash and identity preserved. The independently
repeated raw-source comparison also passes. The corrected comparison gate is true;
the original timed-out experiment remains failed. Operational and full-source
gates remain false. This is inspection/persistence acceptance, not successful
preparation or force evaluation under PCFF-IFF.

Reconstruction can be repeated without any new scientific preparation by copying
the checked `acceptance/relocated` records to a new output root, then running:

```sh
$PY scripts/validate_pcff_variants.py --child /new/root/relocated \
  --declaration "$DECL"
```

The caller must supply the current verified source paths in the declaration.
Relative source paths are resolved against the process working directory;
absolute current paths can be supplied for a different working directory. They
are not recovered from historical provenance. The bounded correction invocation,
copy hashes, command and timing are retained in `reconstruction-corrected/`.

## Executed verification and limitations

The ordinary suite passed **1,613 tests**, with **10 skips** (1144.53 s). The
focused source/PCFF/workflow suite passed **78 tests** (526.17 s). Two additional
inspection-driver regressions passed after the consumer correction (0.45 s);
public scientific/native implementation and validators were unchanged by that
correction. Ruff and pip checks in both retained environments passed. The public
source-variant example and inspection CLI ran; `--require-full-source` exited 1.

Independent checks covered four source files and **2,991 ordered requests** across
six cases and three distinct variants, including the 20 structurally inapplicable
requests per variant (2,931 applicable requests). Eighteen case/variant charge
calculations were checked independently in Decimal arithmetic. Twelve direct
automatic-angle ambiguity checks, six positive/negative controls, and twelve
case/source invocations of the pinned external matcher passed their declared
control expectations. None established wildcard priority.

Six historical bundles and two completed workflows passed offline status/frame/
no-op-resume checks with scientific imports blocked, matching J6 identities and
hashes exactly. **1,033 historical files** and all declared source hashes remain
unchanged. Original J1–J6 failures, libraries, bundles and checkpoint directories
were not rewritten. J7's original timeout is preserved independently of the
successful corrected reconstruction.

Receipts:

- [Corrected comparison gate](evidence/phase_4j7.json).
- [Source inventories](evidence/phase_4j7_sources.json) and
  [difference manifests](evidence/phase_4j7_differences.json).
- [Interaction gaps and alternate-source resolutions](evidence/phase_4j7_gaps.json).
- [Wildcard controls](evidence/phase_4j7_wildcards.json).
- [Verification](evidence/phase_4j7_verification.json),
  [historical inspections](evidence/phase_4j7_historical_offline.json), and
  [preservation hashes](evidence/phase_4j7_preservation.json).

Remaining requirements are concrete: source-backed nh+ increments; a justified
charged-guanidinium conservation convention; missing sulfur/lactam Class II
couplings; a scientifically reviewed whole-source IFF operational profile,
including type meaning, applicability, special-pair policies and numerical
verification; and authority for any wildcard precedence rule. A new source's
ordinary equivalence is not grafted into the old signed model. The known source
profile's historical `_cho_v1` name is preserved as an identity contract; it does
not redefine the project's full chemical-coverage target as CHO.

No newly complete native model, whole-system force calculation, minimization or
workflow is claimed. The comparison gate does not establish chemical accuracy,
production readiness or full-source completion. No merge occurred.
