# Phase 4D2.1: AmberTools correction and live validation

This is a correction to the [Phase 4D2 backend](phase_4d2.md), not a new
force-field implementation. The authoritative ISLAND graph and stable site IDs
remain unchanged. The backend still delegates numerical parsing and nonbonded
checks to `import_amber_prmtop()`, and successful preparation still reports
`production_validated=False` and `simulation_readiness="not_established"`.

## Reproduced defects and contracts

1. **Atom names.** The old `A000`/`A001` names were rewritten by real
   antechamber's `adjustatomname()` (for example, `A000` on carbon became
   `C000`), breaking the exact-name lineage check. Generated names now begin
   with the element symbol and use a global base-36 index within the four-byte
   atom-name limit: `C000`, `O001`, `Cl02`, `Br03`. This supports the backend's
   100-site maximum, including one- and two-letter elements. The input, typed
   MOL2 and final prmtop must retain a unique bijection; reordering is allowed,
   renaming/truncation/collision is not. There is no atom-order or coordinate
   proximity fallback. The rule follows the reviewed AmberClassic
   [`adjustatomname()` implementation](https://github.com/Amber-MD/AmberClassic/blob/e483b09f60872db28c2349dbef0a290a4c64e02a/src/antechamber/common.c)
   and its [correction table](https://github.com/Amber-MD/AmberClassic/blob/e483b09f60872db28c2349dbef0a290a4c64e02a/dat/antechamber/CORR_NAME_TYPE.DAT).
2. **Preparation integrity.** The old wrapper accepted a rehashed outer
   record claiming GAFF/AM1-BCC while its signed imported result claimed
   GAFF2/provided charges. Schema `island_ambertools_preparation_v2` and engine
   version `2` require a real `ImportedAmberResult`, exact agreement between
   the outer record and its signed inner provenance, and structural/semantic
   checks for charge method, force-field family, coordinate digest, atom maps,
   source/artifact checksums, toolchain, and all three successful stages. The
   outer-only `imported_result_signature` is excluded from the inner copy to
   avoid circular hashing. A matching hash is necessary but cannot make a
   contradictory record valid. `validate_integrity()` and
   `to_parameterized_system()` raise `InvalidAmberImportResultError` before
   snapshot construction. Valid `dataclasses.replace`, `deepcopy`, and owned
   snapshots remain supported. This version change intentionally invalidates
   old v1 preparation wrappers; it does not alter legacy imported-result or
   Phase 4B/4C signatures.
3. **Restart coordinates.** Previously an achiral molecule bypassed the CIP
   path entirely, so a readable restart with a later NaN was accepted. Final
   lineage validation now requires exactly `N × 3` numeric, finite coordinates
   for *every* molecule, before any optional stereochemistry check. Finite
   source charges are checked as well. A bad restart produces a `tleap`-stage
   error with retained job artifacts; it never changes the input system.

The preparation record includes the input coordinates in angstroms. Parameter
lengths remain in nm after import; this correction does not rescale source
coordinates or establish an energy evaluator. The default charge-total
tolerance remains `1e-4 e`; the real AM1-BCC phenol runs explicitly requested
`0.002 e` because output serialization left approximately `-0.001 e` total.
No residual was redistributed. Provided per-site values remain checked at
`1e-5 e` across typed MOL2 and prmtop.

## Live reference execution

Five cases were run without mocking `_run_stage`, lineage validation, or
`import_amber_prmtop()`. The read-only existing Conda installation was
AmberTools **24.8**, build `cuda_None_nompi_py310h834fefc_101` from
`conda-forge/linux-64`; its package metadata SHA-256 was
`cf35bc9959f0485ed11c9c2f5521bf71cea8a0199379a6af63ad9ac02d5fc376`.
Individual executable version banners were unavailable from their help output;
the preparation records instead include executable checksums and the package
metadata. The loaded GAFF data header said **1.81** (SHA-256
`21ce262ab4af254c2f2b4bf95f4413c2a6593c67faa31ffc7108b7dd07a63d09`);
GAFF2 said **2.2.20** (SHA-256
`14ad62c8e532c47e2e400e2ca6ad8052b33bda4c4b9bc6f49b6a528e527512be`).
The matching `leaprc.gaff`/`leaprc.gaff2` checksums were respectively
`901104629d8c7566a1179cea5f4499d69ac7d5e914d5ea30f80021393c8bcde7`
and `4f8f713bf7cc38d72ee0c4f2fa4178cef2a87f2a3b37abef1f77ed39e0ff94ca`.

| Input and method | Periodic impropers | Source 1–4 pairs | Zero-LJ sites | Outcome |
| --- | ---: | ---: | ---: | --- |
| Phenol, GAFF, AM1-BCC | 6 | 23 | 1 | Passed |
| Phenol, GAFF2, AM1-BCC | 6 | 23 | 0 | Passed |
| Phenol, GAFF2, nonuniform provided charges | 6 | 23 | 0 | Passed; no SQM file |
| Capped PE DP=3, GAFF2, provided charges | 0 | 45 | 0 | Passed; local-template input |
| Assigned `F[C@H](Cl)Br`, GAFF2, provided charges | 0 | 0 | 0 | Passed; `Cl`/`Br` names and R center preserved |

The reference generator independently compared source Amber and converted
ISLAND bond, angle, proper/improper torsion, LJ, and nonexcluded Coulomb
energies at explicit non-equilibrium geometries. The ordered improper uses a
dihedral computed from four 3D points. All applicable comparisons passed; a
family absent from a source is recorded as absent, not fabricated. It also
records mapping, exclusions, 1–4 counts, commands, checksums, SQM markers,
charge settings, and result signatures. This is conversion agreement, **not**
scientific suitability validation. The older pinned upstream phenol fixture
still has unknown exact force-field and charge-method provenance.

Reproduce with the project Python, ParmEd/RDKit installed, and AmberTools
executables on `PATH` with `AMBERHOME` set:

```bash
python scripts/generate_ambertools_references.py --output /tmp/island-ambertools-references
ISLAND_RUN_AMBERTOOLS_INTEGRATION=1 python -m pytest -q tests/test_ambertools_integration.py
```

The generated `references.json` and unique per-job directories contain the
exact source files and logs. The original reviewed run at
`/tmp/island-4d21-final-RpoirW` is historical, temporary evidence. The completed
2026-09-29 rerun is retained on project storage at:

```
/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-validation/phase-4d2-1-20260929-final
```

Its [committed provenance manifest](references/phase_4d2_1/manifest.json) has
SHA-256 `eb4281fdd9275fc4be0a8b22a72e65860e51ed003b9e7eb6ac56e79f2b81edb2`.
It includes all five cases, exact commands, source inputs and maps, executable
and data checksums, Python/package versions, generator and production-source
checksums, and checksums for all 118 retained files. Each job contains the
original ISLAND `input_system.json` as well as input MOL2, charge inputs where
applicable, output topology/restart and logs. `artifact_subdirectory` locates a
job relative to the archive root after relocation; absolute paths in signed
records describe the original run and must not be rewritten.

Redistribution rights for files generated from the mixed-license AmberTools
data installation have not been separately reviewed, so the raw real outputs
are retained locally rather than published. The committed manifest contains
provenance and numerical comparison results, not parameter tables. Offline
regressions generate real parsed synthetic prmtop/restart fixtures and retain
the separately provenanced upstream phenol fixture. Live tests remain opt-in.

To keep a new live test's outputs on persistent storage, choose a fresh output
path (the generator refuses to overwrite an existing manifest):

```bash
ISLAND_RUN_AMBERTOOLS_INTEGRATION=1 \
ISLAND_AMBERTOOLS_REFERENCE_OUTPUT=/path/to/new-reference-directory \
python -m pytest -q -rs tests/test_ambertools_integration.py
```

## Follow-up integrity acceptance (2026-09-29)

The audit of `b35586e` reproduced ten additional failures using actual validated
synthetic imports with both preparation copies re-signed: nontext force-field,
charge-method, CIP and command fields leaked `TypeError`; contradictory duplicate
command flags passed; missing/inconsistent serialization tolerances passed;
and a provided-charge outcome could claim an SQM completion marker.

Both public result methods now reject these through
`InvalidAmberImportResultError`. All three stage commands use the same validation
of text option/value pairs, unique flags, successful return status and required
settings. Four further regressions bind input/output filenames, the authoritative
formal charge and singlet multiplicity to the actual preparation command. The per-site serialization tolerance is the engine's fixed `1e-5 e`
and the provided-charge outcome must agree. The AM1-BCC completion marker must
match the same expression used for actual SQM output; fabricated text is rejected.
Canonical JSON digests compare the shared payload, excluding only the outer
`imported_result_signature`. Additional tests cover missing required fields,
invalid completion markers, valid reconstruction/deepcopy and owned snapshots.

Schema v2 and engine version 2 remain: this closes validation holes without
changing the emitted preparation format. Legitimate v2 records remain valid;
previously accepted malformed records are intentionally rejected. These are
content hashes, not cryptographic attestations of execution; live retained
artifacts provide the execution evidence.

Verification in the existing `miniforge3/envs/island` Python environment:

| Check | Result |
| --- | --- |
| Complete ordinary suite | 416 passed, 1 skipped (explicit opt-in live test) |
| Real AmberTools suite | 1 passed, exercising all five cases above; no mocks |
| ParmEd deliberately unavailable in an isolated subprocess | 328 passed, 7 skipped: three ParmEd-dependent modules, three backend cases, and the opt-in live test |
| Ruff | All checks passed |
| `python -m pip check` | No broken requirements |
| All 12 example scripts | Passed, including actual short-PE preparation |

AmberTools was read from the existing `miniforge3/envs/md_rg` installation;
no environment was installed into or changed. All five live cases also create
validated snapshots and check input non-mutation and readiness flags. No real
integration acceptance case remains unexecuted. Raw-output public redistribution
and scientific suitability remain outside this software acceptance.

Remote verification found `origin/main` at `fc2d93e`, Phase 4D1 at `2d346b3`,
and unmerged Phase 4D2 at `b35586e` before this follow-up. Corrections therefore
stay on `codex/phase-4d2-ambertools-backend`. Merge Phase 4D1 (including its
optional-ParmEd fix) first, then Phase 4D2 including both correction commits.

Remaining scope limits are the 100-site, single connected closed-shell input;
whole-chain AM1-BCC scalability; unvalidated scientific charge/parameter
quality; unsupported Amber energy terms and pair overrides; and no MD,
minimization, export, packing, or production-readiness claim. Phase 4D2.1
does not add repeat-unit charge transfer.
