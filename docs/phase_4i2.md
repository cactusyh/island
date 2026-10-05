# Phase 4I2 — durable prepared force-field bundles

## Repository and scope

Branch `codex/phase-4i2-prepared-forcefield-bundles` starts at
`af34439ae8b1caedc8f53325d056a3a846acd1ec`, updated main. Its content was compared
with reviewed I1 `fb62b931c6622aecb0e715659c45d66fb375ae8c` and matched after squash
merge. No prerequisite branch was merged locally.

This is additive persistence for I1's facade and existing native records. No
force-field mathematics, source pin, charge policy, evaluator, historical record
schema or signature is changed. It does not generalize the Amber-only workflow,
provide a dynamics restart, or add preparation during reconstruction.

## Public contract

Imports are from `island.forcefields`:

```python
from island.forcefields import (
    AmberBundleArtifacts, PreparedForceFieldSources,
    save_prepared_forcefield, load_prepared_forcefield, create_evaluator,
)

# PCFF example: the caller explicitly resolves its installed, pinned FRC.
sources = PreparedForceFieldSources(pcff_frc="/installed/pcff.frc")
save_prepared_forcefield(original_system, prepared, "/existing/parent/new-bundle",
                         sources=sources)
# The bundle directory can now be moved to another location/process.
loaded = load_prepared_forcefield("/relocated/new-bundle", sources=sources)
system, prepared = loaded.system, loaded.prepared
result = create_evaluator(system, prepared).evaluate_fresh()
```

`save_prepared_forcefield()` returns the destination `Path`. The destination must
not exist, and its parent must exist. There is deliberately no overwrite option.
It requires the original system retained by the facade, including coordinates and
metadata. Saving a coordinate-replaced system under the old preparation is
rejected; the original Amber evidence is never edited to make it fit.

`load_prepared_forcefield()` returns `LoadedPreparedForceField`. Its `.system`,
`.prepared`, and `.manifest` properties return independent owned copies of the
original `MolecularSystem`, validated `PreparedForceField`, and manifest envelope.
No Context or evaluator is retained or serialized. Subsequent compatible coordinate
replacement uses the unchanged `create_evaluator()` binding checks. Changes to
graph, stable IDs, masses, elements, assigned stereochemistry or relevant native
policy remain subject to the existing validators.

All malformed input, missing/incompatible source/artifact, native reconstruction
or publication failures cross the `PreparedBundleError` boundary (a
`ForceFieldError`). Original exceptions are chained, preserving native diagnostics.
No error triggers fallback, retyping, parameter regeneration or source search.

### Required inputs and native reconstruction

| Family | Additional save input | Save/load source resolution | Reconstruction |
|---|---|---|---|
| GAFF / GAFF2 | `artifacts=AmberBundleArtifacts(original_prmtop)` explicitly authorizes copying that file | None; OPLS/PCFF source options are rejected | ParmEd reparse of verified original prmtop; restore original import provenance/signatures; native preparation validation |
| OPLS-AA | None | `PreparedForceFieldSources(opls_xml=local_xml)` | Existing `OPLSParameterizationResult` from its unchanged native JSON, validated against the exact pinned XML and system |
| PCFF | None | `PreparedForceFieldSources(pcff_frc=local_frc)` | Existing H3 assignment and H4 model records, validated with their external pinned source and authoritative system |

For Amber, a minimal save is:

```python
save_prepared_forcefield(
    original_system, prepared, "/existing/parent/amber-bundle",
    artifacts=AmberBundleArtifacts("/retained/job/result.prmtop"),
)
loaded = load_prepared_forcefield("/relocated/amber-bundle")
```

The original prmtop must match the signed imported-result source hash. It supplies
the complete supported numerical parameter inventory and original index mapping;
it is never regenerated. Preparation/import signatures, native charges, lineage,
charge method, synthetic/user source labels and historical tool/artifact provenance
remain unchanged. Other historical logs, MOL2 and restart files are not necessary
for this reconstruction and are not copied. Their original hashes and paths remain
recorded, without claiming their files were rechecked by each bundle load.

Loading Amber snapshots the checked prmtop bytes into a private temporary file for
ParmEd. Historical absolute paths are not used for operational reads. It reuses
`preparation_from()` and the native validators, then `adopt_forcefield()`. No new
successful signature is substituted for a mismatching historical signature.

OPLS source/typing/native-charge evidence, RB terms and geometric mixing stay in
native records. PCFF retains the full source-derived H3 assignment, native charge
and typing evidence, H4 policy-derived zeros, explicit LJ/Coulomb policy and Class
II/9–6 representation. Neither is made into an Amber import or generic harmonic
parameter container. The native record constructors and integrity validators are
the same ones used by the existing native loaders; fixed classes are selected by
a static family allowlist, not by arbitrary module/class names in the files.

Source libraries remain external in v1. There is no library-copy or download mode.
Explicit resolver paths can change after relocation as long as their pinned bytes
and native identities match. A source path is never persisted as an operational
requirement. GAFF source libraries/executables are not needed to reparse resolved
prmtop data. The original native source/version identities remain in provenance.

## Format and validation

Directory manifest schema: **`island_prepared_forcefield_bundle_v1`**.
`manifest.json` is an envelope `{payload, sha256}`; its SHA256 covers canonical
strict JSON of the payload. The payload contains exactly:

- `schema` and allowlisted `family`;
- the original I1 `prepared` description and `prepared_identity`;
- `files`, mapping logical names to fixed relative `path` and complete-byte `sha256`.

File inventories are fixed by family:

| Family | Logical artifacts → filenames |
|---|---|
| GAFF / GAFF2 | system → `system.json`; preparation → `preparation.json`; prmtop → `result.prmtop` |
| OPLS-AA | system → `system.json`; parameters → `parameters.json` |
| PCFF | system → `system.json`; assignment → `assignment.json`; model → `model.json` |

`system.json` and Amber `preparation.json` use the established explicit container
tags (`map`, `list`, `tuple`) to preserve integer keys and tuple/list distinctions.
This is data-only JSON, not Python object deserialization. The system encoder
preserves sites, authoritative bonds, explicit impropers, original coordinates,
metadata and coordinate provenance. Stale derived interaction caches are not
promoted to authoritative inventories. Units continue to be the existing system
angstrom/native-record conventions.

OPLS/PCFF native signed JSON text is stored unchanged. Amber reconstruction data
uses the existing `{record, record_signature, import_source, import_provenance}`
payload plus original prmtop bytes. All native validators run before exposing a
loaded facade. The reconstructed facade description and identity must exactly
match the original manifest. Checksums alone do not establish those relationships.

Strict JSON rejects duplicate keys and nonfinite values. Typed container decoding,
canonical lossless system round trips and native validators reject malformed
mappings, boolean-coordinate coercion, duplicate IDs, truncated records and
contradictory data. Fixed artifact paths and symlink rejection prevent traversal
or substitution through bundle artifact references. Extra unreferenced files are
not read or treated as scientific inputs.

The existing nonperiodic Amber bounding-box metadata is preserved with finite
positive lengths and three explicitly false periodic flags. This is not a periodic
cell implementation. Existing OPLS/PCFF no-box restrictions remain in force.

Validation establishes internal consistency, not computational authenticity.
Checksums cannot prove that arbitrary caller-supplied evidence came from the
claimed calculation. Nor does persistence turn a scientifically unsuitable model
into a suitable one.

## Transactional publication and recovery

1. Validate the original system/facade, explicit source resolver and required
   artifact before publication. Read the original prmtop once and verify its hash.
2. Write/fsync a private sibling `.name.candidate-*` directory. Native files and
   system data precede the manifest.
3. Reconstruct the entire candidate through the same semantic validation used by
   loading. No scientific execution occurs during this check.
4. Exclusively reserve the new destination with `mkdir`. An existing destination
   or another writer wins unchanged.
5. Atomically rename the complete candidate over **only this writer's empty
   reservation**, then fsync the parent directory.

This uses the existing same-host POSIX filesystem publication model. It never
replaces a pre-existing bundle, including an empty user directory. Exceptions
before the rename remove the writer's temporary candidate/reservation where
possible; cleanup failures do not mask the original reconstruction/rename error.

A process interruption before publication can leave an unpublished hidden candidate
or empty reservation. Neither is an authoritative bundle at the requested target;
loading an empty target fails. After confirming no writer is active, these leftovers
may be inspected and removed explicitly. Do not automatically select a candidate
based on its filename. An interruption after rename leaves the complete validated
bundle; failure of the final directory fsync means durability was not confirmed
and is reported as an error. It does not expose a partially populated bundle or
replace an earlier valid target.

Loading snapshots the manifest and checked file bytes. Immutable bundle files are
never rewritten during loading. The verified snapshot is reconstructed without
re-reading scientific JSON through uncontrolled historical paths.

## Dependencies

Importing the save/load APIs requires no optional scientific stack. They do not
import OpenMM, RDKit, ParmEd, Foyer or SciPy on import. Dependencies during actual
reconstruction are:

- **Amber:** core NumPy plus ParmEd. A real retained GAFF bundle was also loaded
  in a fresh Python process with OpenMM, RDKit, Foyer and SciPy imports blocked.
  ParmEd may import installed optional packages transitively when not blocked;
  this is distinct from requiring them for the demonstrated reconstruction.
- **OPLS / PCFF:** core/NumPy and the explicit pinned local XML/FRC. No Foyer typing,
  OpenMM, RDKit, AmberTools or SciPy execution is needed. Existing offline source
  and assignment consistency arithmetic still runs.
- **Evaluation afterward:** the existing OpenMM dependency and native Reference
  evaluator contract. A loaded preparation alone is not an evaluation.

The compatible separate OPLS environment remains available, but the actual I2
OPLS reconstruction/evaluation succeeded in the main environment without importing
Foyer. Both environments passed `pip check`. No dependency was installed or
project environment downgraded for this phase. Native validators may reject data
if an incompatible parser/software change alters its original identity; no
cross-platform bitwise guarantee is made.

## Actual bounded acceptance

`validate_prepared_bundles.py` consumes only manifest-listed, checksum-verified I1
artifacts. It does not construct molecules, parameterize, retype or run QM.
The I1 declaration stored ordinary JSON coordinate mappings; the acceptance
reader restores only their schema-defined canonical integer site-ID keys and
preserves the recorded metadata without guessing other types. I2 bundles themselves
use the lossless tagged format above.

For each retained PE DP3 preparation it verifies the I1 outcome and exact facade
identity, publishes a declaration, evaluates the pre-save state, saves, relocates,
and launches a genuinely separate process from the new working directory. During
child loading, explicit guards reject subprocess execution, preparation/typing
entry points and OpenMM Context construction. After loading, it compares a fresh
evaluation, a session evaluation and a session fresh evaluation with the pre-save
state. Named components and all stable-ID force components are checked, alongside
native parameter/model/coordinate/evaluation identities and facade identity.

Predeclared tolerances are unchanged from I1: **atol=1e-10, rtol=1e-12**, in kJ/mol
and kJ/(mol·angstrom). These comparisons establish persistence/adapter fidelity,
not independent force-field accuracy. Same-host equality is not a portability
promise. GAFF/GAFF2 retain their explicitly synthetic alternating ±0.01 e labels;
OPLS and PCFF retain their source-native charge policies.

All four retained preparations passed actual separate-process reconstruction and
all three evaluation modes (fresh, session, session-fresh):

| Selector | Status | Max energy error (kJ/mol) | Max named component error (kJ/mol) | Max force error (kJ/(mol·Å)) | Seconds |
|---|---|---:|---:|---:|---:|
| gaff | passed | 0 | 0 | 0 | 1.21 |
| gaff2 | passed | 0 | 0 | 0 | 1.41 |
| oplsaa | passed | 0 | 0 | 0 | 6.09 |
| pcff | passed | 0 | 0 | 0 | 196.46 |

Facade/native identities and complete saved system data matched exactly. All I1
manifest-listed input hashes were rechecked unchanged after execution. Outputs are
retained in `../island-validation/phase4i2/{gaff,gaff2,oplsaa,pcff}`. The
[machine-readable evidence](evidence/phase_4i2.json) records source pins, input/output
hashes, bundle manifests, commands, parent/child PIDs, outcomes and environment.

All cases used the existing main Python 3.11.16 environment (NumPy 2.4.6,
OpenMM 8.6.1, ParmEd 4.3.1). No Foyer import was observed in any child. OPLS/PCFF
children also did not import RDKit. Installed ParmEd transitively imported RDKit
in the unrestricted Amber children; a separate real-record offline test blocked
that import and still passed.

Supplementary fresh-process offline checks passed for all four bundles with
OpenMM, RDKit, Foyer and SciPy imports blocked; OPLS/PCFF additionally blocked
ParmEd. These checks loaded/validated records without evaluation. The actual
numerical acceptance imported OpenMM only for the guarded load check and subsequent
evaluation. No typing, parameterization, QM, minimization or dynamics was rerun.

Elapsed times include saving, full candidate reconstruction, relocation, process
startup, native integrity checks and numerical comparisons. They exclude the
initial retained-record audit. PCFF's deeper validation dominates its measured
cost; no speedup or new cache is claimed.

### Reproduction

```bash
PY=../island-validation/phase4e2-env/bin/python
$PY scripts/validate_prepared_bundles.py --family gaff \
  --retained ../island-validation/phase4i1-final/gaff --output /new/gaff-check
$PY scripts/validate_prepared_bundles.py --family gaff2 \
  --retained ../island-validation/phase4i1-final/gaff2 --output /new/gaff2-check
$PY scripts/validate_prepared_bundles.py --family oplsaa \
  --retained ../island-validation/phase4i1-installed/oplsaa \
  --xml ../island-validation/foyer-4g1-upstream/foyer/forcefields/xml/oplsaa.xml \
  --output /new/opls-check
$PY scripts/validate_prepared_bundles.py --family pcff \
  --retained ../island-validation/phase4i1/pcff \
  --frc ../island-validation/phase4h1-sources/lammps/pcff.frc \
  --output /new/pcff-check

# Public API example: load an existing bundle, save a copy, relocate, spawn child.
$PY examples/prepared_forcefield_bundle.py \
  --bundle /completed/gaff-check/relocated --output /new/example
# Add --xml or --frc for the corresponding source-dependent bundle.
```

Each acceptance directory is exclusive. Missing retained inputs or failed child
checks produce nonzero exit status and a durable failure outcome, never replacement
scientific inputs. Per-case declarations, commands, PIDs, baseline records, bundle
manifests, results, hashes and timings remain external. Full libraries and generated
prmtop data are not committed to Git.

## Verification and remaining limits

Executed checks:

- Complete ordinary suite: **1,338 passed, 10 skipped**, 317.40 s.
- Focused bundle regressions: **32 passed**, 25.51 s.
- Ruff, whitespace checks, and `pip check` in both the main and existing separate
  OPLS environments passed.
- The new example actually saved, relocated and evaluated the GAFF bundle in a
  child process; its synthetic provided-charge provenance was preserved.
- All four retained real-backend reconstruction/evaluation gates passed, with
  zero measured errors on this host. No retained artifact or required environment
  was unavailable. The ordinary suite's existing opt-in/dependency skips were not
  substituted for this separately executed acceptance.

Software regressions separately exercise native round trips, coordinate replacement,
mass/stereo/site/charge/model tampering with recomputed envelopes, malformed JSON,
missing files, checksum failures, source pins/relocation, typed metadata, boolean
coordinates, source/artifact mismatches, overwrite protection, injected rename and
cleanup failures, and concurrent destination ownership. These injected failures
are software evidence, not additional live scientific calculations.

No new automatic chemistry, numerical engine, charge policy, minimizer, integrator,
workflow migration, periodic interaction, packing, crosslinking or MLIP support was
introduced. PCFF native integrity reconstruction remains relatively expensive;
public validation is not bypassed or cached. Missing libraries and Amber prmtop
artifacts cannot be regenerated by this API.

These bundles preserve preparation inputs and native models. They are not a
trajectory checkpoint and contain neither velocities nor thermostat RNG state.
The next force-field-independent single-chain workflow can build on this layer,
but is not implemented here.

`production_validated=False`, `simulation_readiness="not_established"`.
