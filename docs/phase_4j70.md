# Phase 4J70 — authoritative PCFF missing-row intake

The single decision is **`parameter_incomplete`**. The declared local sources do not supply authoritative numerical values for any of the 378 missing applicable Class-II interaction instances in the exact neutral-singlet, explicit-H, 27-atom HMDS graph. J70 executes the source intake, exact interaction mapping, two strict converter checks and offline static reconstruction. It creates no source-complete bundle, execution profile or automatic assignment capability. No values are fabricated or inferred from a favorable QM residual.

The branch is `codex/phase-4j70-authoritative-pcff-missing-row-intake`. Fetched `origin/main` contains required J69 ancestor `9f379f85f58807d9ae6502da1db712d8fc5680c9`; the actual base is fetched main `efacd4b798a1f345ea05f928905be88ad49acf53`. No main merge or push is performed. Unrelated working-tree and untracked files are preserved.

`docs/evidence/phase_4j70/plan.json` was frozen before source parsing, SHA-256 `162cc2d03c8619c9b753152d198aba342dc6c9a15e8162d32ea52219d6a4cd00`. It binds the source inventory, J69 receipt, base, exact-only intake policy, 26 static fixtures, one converter worker, a 4 GiB allocation, 60-second converter timeout and two-attempt cap. The allocation is declared; the subprocess timeout is enforced. J70 performs **zero new QM attempts, zero MD steps, zero production fits and zero research solver calls**. All three reserved validation configurations remain untouched. Si–H and capped-DP2 calculations remain inactive.

`source-inventory.json`, SHA-256 `b1a83f81b75506c5c3669ffcbcdf4d4590fbdee6f01490d291a60930f396dbd9`, freezes 89 original local files before the parser sees them. Discovery includes the mounted `island_ffref` literature/archive directory, declared phase4h1/h2/h3/h4 and phase4j7 source directories, the pinned msi2lmp distribution and converter sources, and the declared author-MD directory's README. The inventory records original absolute paths, hashes, sizes, source categories, version/provenance and licensing status. `source-review.json` adds archive/member hashes, exact prior J44/J54 citation, bibliography, license and row-location reviews with their original identities. Reusing those reviews is explicitly recorded; J70 does not claim a new PDF transcription. No network download or unrelated force-field substitution occurs. Large raw literature archives and production trajectories are not committed.

The sources remain separated:

| Source category | Discovered assets | Intake outcome |
| --- | ---: | --- |
| Native PCFF distribution | 3 | Same pinned FRC bytes; all 378 applicable rows absent |
| Mounted literature/external material | 62 | Exact existing reviews and archive identities retained; no authoritative missing rows supplied |
| External comparison assets | 11 | Other force-field files and IFF/LUNAR material remain separate |
| Converter C/header assets | 12 | Pinned msi2lmp provenance retained |
| Research-derived author-MD README | 1 | Benchmark/simulation material is not an original author parameter artifact |

Five PCFF-labelled FRC files are independently parsed and queried without pooling their records. The three native copies share SHA-256 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`; the IFF/interface candidates have separate hashes and provenance. Every one returns 378 missing queries. The engineering `LAMMPS-pcff-main.zip` contains a byte-identical native FRC and adds no rows. The benchmark archive uses COMPASS/OPLS/other models; its values are not native PCFF authority. No native, author-primary, literature-external or research-derived coefficient is adopted. Original author artifacts remain `not_yet_supplied`; the 28 confirmations remain user-attested.

The immutable versioned `island.pcff-source-row.v1` records bind graph identity, stable sites and ordered types, permitted symmetry mapping, namespace and source lookup, original raw values/location/hash, normalized values and conversion factors, coefficient units, scalar energy and force equation, prefactor, multiplicity, orientation, phase and equilibrium dependencies. Missing coefficients are JSON `null`. They are never explicit zeros. Supplied rows must be reparsed from their original bytes; a matching outer identity alone is insufficient. Unit conversion is stored separately from original source cells. Frozen serialized records return copied payloads. Full graph/dependency validation occurs against the unchanged retained model. Cross-family borrowing, wildcard or lower-family fallback, duplicate/partial bundles and silent native/external overlay are rejected. Research values cannot grant native source authority or execution permission.

`completeness-matrix.json` enumerates every missing interaction instance, rather than reporting only symmetry classes:

| Family | Instances | Current row status |
| --- | ---: | --- |
| Angle-angle | 78 | `qm_fit_required`, inherited unexecuted target requirements |
| Bond-bond 1–3 | 60 | `missing` |
| End-bond-torsion | 60 | `missing` |
| Middle-bond-torsion | 60 | `missing` |
| Angle-torsion | 60 | `missing` |
| Angle-angle-torsion | 60 | `missing` |
| Total | 378 | Incomplete |

Enumeration is checked against the entire unchanged graph interaction inventory, with unique stable instance IDs and exact family counts. Each row preserves its source-equilibrium dependency IDs, values and units, lookup evidence and required missing authority. The historical 56 symmetry-reduced coupling unknowns are a separate research representation; with the two stretching corrections they form the existing 58-variable research Hamiltonian. Neither a missing native row nor a symmetry-reduced fitted coefficient is silently relabelled as the other.

`missing-row-and-convention-report.json` records all 378 source-unavailable instances. Separately, 240 torsion-dependent instances retain unresolved native Phi conventions; angle-angle-torsion additionally has unresolved literal linear-Phi versus cosine angular units/equation authority. The research cosine AAT form does not resolve the native FRC convention. Empty source sections do not establish numerical zeros. Exact coefficients without complete conventions would remain ineligible. `qm-fit-targets.json` freezes 78 exact site/type/dependency target records and explicit input templates without coefficient guesses. Their producer, QM protocol and fitting contract remain unset, so they are labelled non-executable. J70 authorizes neither their QM execution nor fitting; a later fitted row would remain research-derived without independent native authority.

`parameter-bundle-schema.json` records the immutable row contract and completeness gates. `parameter-bundle.json` references every instance identity and remains `parameter_incomplete`, with zero supplied coefficients, no profile and no execution authority. Completeness also requires closure of the full Class-II inventory, resolved conventions, deterministic complete strict conversion and independent full static closure. Parsing a coefficient or passing a synthetic table check cannot satisfy these gates. Future source-supported conventions require an explicit versioned contract rather than relaxing the present unresolved declarations.

Two actual, non-QM converter attempts used the unchanged hash-verified J69 msi2lmp executable, v3.9.11 / 6 Sep 2024, built from pinned LAMMPS revision `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`. Each used a separate fresh relocated directory, exact retained HMDS CAR/MDF inputs, stable types, original charges and native FRC, with `-class II -print 3 -nocenter` and **without `-ignore`**. Both returned code 17:

```
Unable to find endbonddihedral data for hc c3 sio c3
```

The executable, FRC, inputs, commands, runtimes, return codes, stdout and stderr hashes are retained under `converter/`. Stdout is byte-identical across both attempts: deterministic failure is established. Deterministic complete coefficient emission and complete conversion are not established. Both permitted attempts are consumed; no hidden retry or permissive conversion follows. The software test for a complete converter-table transcript is a manufactured control, not an actual successful HMDS conversion.

`static-closure.json` records fresh offline replay of J69's available-source energy/force audit over all 26 retained configurations, with original producer identity. Maximum available-family energy difference is approximately 1.78e-15 kcal/mol and force-component difference approximately 4.75e-14 kcal/mol/angstrom. This closes the existing declared available subset. No new complete conservative potential can be reconstructed for the missing coefficients; full-bundle static closure is explicitly blocked. The original native energies/gradients and PBE0/def2-TZVP dense-grid protocol are retained unchanged. No static replay creates a new QM execution identity.

`canonical-capability.json` distinguishes the existing bounded family-only inspection route, `island.pcff_pdms.prepare_polymer_forcefield`, from the ordinary public facade, `island.forcefield_policy.prepare_polymer_forcefield`. The former reports the exact 378 missing stable interaction IDs for the HMDS fixture. The latter retains its existing fail-closed `typing_ambiguous` diagnostics. Tests also invoke the existing family-only `prepare_polymer(..., forcefield="pcff", ...)` route on the exact 27-atom HMDS graph and check that requiring completeness fails. No API behavior changes or bounded domain registration are activated. No manual research profile or charge-method name is introduced. A complete bounded capability does not exist, and the ordinary selector cannot automatically assign incomplete rows or fall back across force-field families.

J68's bounded and bound-free incompatibility, the protected original22 result, the four J67 observations and every historical failure remain unchanged. Supplying particular coefficients for already free research directions cannot remove the checked bound-free incompatibility of the unchanged span and baseline. Source completeness and reproduction of the declared QM target remain independent questions. J70 does not attribute that incompatibility to an incorrect source coefficient, and does not claim that filling the 378 instances will fix it.

The next dependency is exact, hash-bound authoritative numerical source rows with complete native Phi/AAT, orientation, equilibrium and pair-convention authority. The matrix identifies every blocked site and its evidence requirement. When such authority becomes available, it can be ingested into its proper source namespace, reparsed, converted strictly and independently checked statically before any bounded registration decision. No generic additional sampling, QM campaign, fitting allocation or automatic bound revision follows from this intake.

The implementation inventory is explicit: `source_intake.py`, the J70 intake/replay/receipt scripts, J70 tests and this document. Historical manifests are not broadened with future-module globs. Tests cover source SHA/row mutations, missing and duplicate detection, ordered/reversed sites, units and normalization, equations/prefactors, phases/multiplicity, equilibrium dependencies, source-domain separation, zero-fill/wildcard/partial rejection, actual strict converter failure, manufactured complete-table gating, available-subset closure, canonical fail-closed diagnostics, re-signed scientific declarations and receipt/manifest schema mutations. Fresh-process relocation copies both code and evidence and reconstructs the retained audit while blocking Psi4 and solver imports. The final verifier independently reconstructs every static summary metric and complete capability declaration, including the ordinary public facade and every missing interaction. It also checks parsed-source membership against the frozen inventory and cross-checks both converter sidecars, exact flags, input/source hashes and raw logs. Fifteen additional re-signed nested-report mutation controls reject false metrics, completed exports, profile registration, provenance overlays and inconsistent converter requests. These controls were actually executed in J70; they are not claimed as reproductions of a historical attack. Synthetic source cells and converter transcripts remain explicitly separate from actual HMDS evidence.

Preliminary failures and corrections are retained separately from terminal successful checks. Initial manufactured FRC inputs lacked mandatory parser sections; they were corrected as synthetic fixtures. Initial focused checks preceded completed capability reports. Report preparation exposed different existing public/bounded status schemas; it was corrected to read `source_assignment.status` and `parameter_completeness`. Intake equivalence lookup was cached without changing scientific records. These are software/intake corrections, not QM attempts or successful source discoveries. The first focused run had 20 failures and 39 passes; later overlapping focused runs passed. They add no unique terminal coverage.

Terminal verification results, exact collection, scope skips and original failures are retained in `validation/`. The supported full suite uses the existing J69 zero-QM/zero-MD policy and verified historical native fixtures. Tests requiring forbidden native calculations or dynamics are explicit scope skips, including manufactured integrator steps; static checks and manufactured conic controls remain enabled. The original full run completed 5,552 cases: 5,139 passed and 413 skipped. After that run, the frozen verifier strengthening passed all 85 final J70 cases. Raw full XML, logs and implementation snapshots are retained; its overlapping 70 J70 cases are replaced once by the final 85 in terminal accounting. Unique terminal coverage is 5,567 cases: 5,154 passed and 413 skipped, including 1,012 J54–J69 regression cases and 85 J70 cases. Focused/prototype repetitions add no unique coverage. All 8,606 base files remain byte-identical. J54–J69 receipt validators, Ruff, compileall, execution-environment pip checks, diff checks and historical byte preservation are recorded separately. The final source hashes are checked again against the original 89-file inventory. The verifier revision was frozen against 21 scientific artifact hashes; none changes. Exact terminal counts and all skip reasons are supplied in the receipt and `validation/terminal-test-counts.json`.

Offline replay needs the repository and retained evidence, not the original external source paths or new native execution:

```bash
PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/replay_phase_4j70.py

PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/replay_phase_4j70.py --deep-parent

PYTHONPATH=.:src OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4g1-env/bin/python scripts/validate_phase_4j70_receipt.py \
  docs/evidence/phase_4j70
```

Source completeness is blocked; native/external provenance remains separate; available-subset source-to-evaluator correctness passes; complete converter emission and parameter assignment capability are blocked. QM-fit targets are unexecuted requirements. Chemical transfer, scientific validation and production readiness are not established. No native PCFF equivalence, complete PDMS coverage, capped-DP2 scope, experimental or condensed-phase accuracy is claimed.

```
accepted_coefficients = 0
automatic_assignment = False
production_validated = False
simulation_readiness = "not_established"
```
