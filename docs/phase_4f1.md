# Phase 4F1 — bounded larger-chain provided-charge preparation

## Scope, branch and size policy

The prerequisite remains unmerged. Branch
`codex/phase-4f1-long-chain-provided-charges` starts from verified remote feature
tip `1d216234df071768e2673500892e348f4a5c8f62`, depending on Phase 4E8/4E8.1.
Fetched main was `a563d7f27451c17659a68bac8686bdae95f2c475`. Merge the prerequisite
before this branch; no merge, main modification or force push was performed.

`AmberToolsOptions.max_atoms` and `WorkflowConfig.max_atoms` expose one execution
budget. `ambertools.policy.validate_size_policy()` is shared by options, input
serialization, engine preflight, workflow inspection/prepared startup, signed
preparation validation and durable bundle reconstruction.

| Charge mode | Default | Allowed configured maximum |
| --- | ---: | ---: |
| Explicit provided charges | 100 | 1,000 |
| Whole-molecule AM1-BCC | 100 | 100 |

Budgets must be positive Python integers; bool, floating-point, null, negative,
zero and out-of-range values fail. The **actual authoritative explicit-site
inventory**, including hydrogen, must fit the configured budget. Diagnostics give
actual count, configured count and charge mode. A larger ceiling is an opt-in
resource policy, not a promise that every chemistry parameterizes successfully.

Oversized systems are rejected before external tools. Workflow preflight checks
installation files without invoking executable version probes; probing occurs
only in the engine after actual size/chemistry checks. The connected, finite,
nonperiodic, closed-shell explicit-hydrogen model, supported elements, valence,
bond semantics and assigned stereochemistry restrictions remain unchanged.

```python
from island.forcefields import AmberToolsOptions, AmberToolsParameterizationEngine

options = AmberToolsOptions(
    "gaff2", "provided", provided_charges=charges,
    max_atoms=1000,
    charge_source="User charge file / citation and method description",
    retain_success_artifacts=True,
)
prepared = AmberToolsParameterizationEngine().parameterize(system, options)
```

Provided mode retains antechamber `-c rc -cf charges.txt`, then parmchk2, tleap
and the existing validated Amber importer. It never requests SQM. There is no
zero-charge fallback, repeat transfer, renormalization or charge-method switch.
Exact stable-ID coverage, component/formal totals, configured charge tolerance
and existing per-site serialization tolerance **1e-5 e** are enforced.
Four-character element-prefixed base-36 names are unchanged: the two-character
Cl/Br suffix capacity (1,296 indices) covers all 1,000 permitted sites. Exact
name/element/connectivity lineage checks remain mandatory.

## Signed preparation and configuration contracts

New preparations use **`island_ambertools_preparation_v3`, engine version `3`**.
The v3 additions are inside both identical shared preparation copies and their
existing signed import/outer-record relationships:

- `size_policy`: schema `island_ambertools_size_policy_v1`, configured `max_atoms`,
  `actual_atoms`, and `charge_method`. Integrity checks bind these to the real
  inventory and mode; contradictory reconstructed records fail even after a
  caller recomputes their content hashes.
- `provided_charge_input` is null for AM1-BCC. Provided mode records the complete
  original string-ID charge mapping, its canonical `digest()` SHA-256, a source
  label and `method_provenance="unverified_user_supplied"`. The default source is
  `user-supplied; charge method unknown`. A caller's descriptive source/citation
  is retained without claiming that ISLAND scientifically verified the method.
- The exact ten-decimal `charges.txt` bytes have their own SHA-256, also required
  in `artifact_sha256`. Integrity validation checks serialization, coverage,
  totals and agreement with imported charges. Mathematical charge content and
  serialized bytes are distinct checksums.

`InvalidAmberImportResultError` remains the malformed-preparation boundary.
Existing **v2 / engine 2** records retain their original 100-site contract and
signatures. They are read as v2, never upgraded or re-signed; v3 policy fields
cannot be smuggled into a v2 record. Historical test fixtures explicitly retain
the v2 schema. The archive loader accepts both versions and restores only the
schema-defined integer mapping keys as before.

New `WorkflowConfig` serialization uses **`island_single_chain_config_v2`** with
`max_atoms` and optional `charge_source`. Historical config v1 remains readable
and serializes with its original fields/default 100; size/source fields require
v2. Prepared startup and bundle reconstruction require the configured limit to
agree with the signed preparation's selected limit (historical v2 means 100).
They also check an explicitly configured source label against saved provenance.
Manifest, bundle, checkpoint, dynamics and analysis schemas remain unchanged.

### Model versus audit identity

The exhaustive imported `content_signature()` and public parameter fingerprint
continue to include **all signed provenance**, including the execution budget.
No historical signatures are rewritten. New `model_content_signature()` gives
the same result for every historical import; for v3 preparation only, it omits
`size_policy` from the model-identity input. OpenMM model identity uses that
signature. Thus changing only a v3 execution-size budget does not change the
physical model fingerprint, force construction or results, while its exhaustive
parameter/audit binding remains distinct. Full import validation still runs;
this is not permission to replace parameters during continuation. No integrator,
checkpoint numerical rule, force-field term or physical setting changed.

## Explicit charge-file example

First inspect exactly the chain that will be built, using fixed default template
and assembly seeds 2026:

```sh
python examples/parameterize_provided_chain.py inspect \
  --psmiles '[*]CC[*]' --dp 50 --max-atoms 1000 --output pe50-sites.json
```

This outputs the actual count and every site's integer ID, element and formal
charge. It does not generate a charge model. Supply a strict JSON file:

```json
{"source":"Describe the actual source, method and citation (or null if unknown)",
 "charges":{"1":-0.10,"2":0.10}}
```

That two-entry illustration is **not** sufficient for PE50: the file must contain
every inspected stable ID exactly once with finite charges and valid component
and total formal charge. Keys must be canonical integer strings. Duplicate JSON
keys and nonfinite constants are rejected. Source may be JSON null to retain the
explicit unknown-method label. The command neither fills missing values nor
copies charges to other repeats.

```sh
python examples/parameterize_provided_chain.py prepare \
  --psmiles '[*]CC[*]' --dp 50 --max-atoms 1000 --force-field gaff2 \
  --charges pe50-charges.json --amberhome "$AMBERHOME" --output new-pe50-preparation
```

The example selects that installation's bin directory in its own process PATH.
It retains input-system data, original signed preparation/import provenance and
all successful parameterization artifacts in the new user output directory.
Use `WorkflowConfig(..., charge_method="provided", provided_charges=charges,
max_atoms=1000, charge_source=source)` for the existing workflow. Keep
AM1-BCC at <=100; future transferable-charge work is not implemented here.

## Predeclared actual execution

`scripts/validate_long_chain_provided.py --output NEW_DIRECTORY --amberhome PATH`
writes its declared matrix/settings before execution and records all failures.
The actual existing AmberTools **24.8**, conda-forge build
`cuda_None_nompi_py310h834fefc_101`, was used from the `md_rg` environment without
altering that environment. Loaded GAFF 1.81 and GAFF2 2.2.20 data headers, data and
executable hashes, package metadata, commands, mappings and artifact checksums
are retained in each preparation record. OpenMM **8.6.1 Reference** evaluates the
resolved records; the independent path uses OpenMM's original-prmtop reader,
NoCutoff, no constraints, switching, dispersion correction or COM removal.

No traceable long-chain per-site charge references were available. All cases
therefore used explicitly labelled **synthetic alternating +0.01/-0.01 e** on
ascending stable IDs (nonzero, neutral, deterministic). These test command,
lineage, serialization, import and numerical consistency only. They do not
validate a scientific charge model. All actual per-site import charge errors
were **0 e**, and no case produced SQM output.

The fixed matrix was PE DP50/GAFF, PE DP100/GAFF2, and oxygen-containing
`[*]CCO[*]` DP20/GAFF2. Both template and assembly seeds were 2026. The configured
limit was 1,000, charge tolerance 1e-4 e, and timeout 600 seconds **per external
stage**. Numerical acceptance was declared before execution:
`atol=2e-5`, `rtol=2e-6`, applied to total/component energies in kJ/mol and forces
in kJ/(mol*angstrom). None was relaxed after observing results.

| Actual case | Sites | Total-energy absolute error | Maximum force-component absolute error |
| --- | ---: | ---: | ---: |
| PE DP50, GAFF | 302 | 6.519e-9 | 1.485e-8 |
| PE DP100, GAFF2 | 602 | 1.132e-6 | 1.490e-8 |
| PEO DP20, GAFF2 | 142 | 4.435e-8 | 1.339e-8 |

All passed total energy, bonds, angles, combined proper/improper torsions,
nonbonded energy and per-site forces. Bond/angle/torsion discrepancies were at
most 6.83e-13 kJ/mol. Total-energy relative errors were 9.10e-15, 1.29e-13 and
1.49e-13 respectively. Reported maximum componentwise force relative errors
(using a 1e-12 denominator floor) were 2.80e-7, 1.86e-5 and 9.34e-9; near-zero
components are governed by the predeclared **combined absolute/relative** test.
These statistics and every component discrepancy are retained in the report.

Initial energies were approximately 7.16e5, 8.77e6 and 2.98e5 kJ/mol, with very
large residual forces. Local-template coordinates are finite and nonsingular
but can contain severe nonbonded strain. Numerical agreement does not establish
that these are suitable MD starting states.

### Workflow smoke: acceptance remains unmet

The declared PEO20 workflow used fmax=0.1 kJ/(mol*angstrom), 500 iterations,
1,000 evaluations, default bounded line search, and otherwise unchanged stopping
criteria. Planned dynamics were 300 K, friction 5/ps, 0.1 fs, four steps in two
segments, four evaluator calls and three retained frames per segment; velocity
seed 78123 and thermostat seed 99181.

The existing workflow reached minimization successfully through its expanded
configuration and durable preparation reconstruction. It returned
`maximum_iterations`: **500 iterations, 512 evaluations**, independently verified
final evaluation. Energy decreased from **297652.8080 to 438.0768 kJ/mol**; fmax
from **1910787.2506 to 2.14017**, and RMS atomic force from **222024.2027 to
0.915389 kJ/(mol*angstrom)**. This fails the unchanged fmax=0.1 requirement.

The workflow correctly published `stage_failed`, accepted step 0, and the valid
unconverged minimization diagnostic. It did not initialize velocities or start
propagation. There are **no accepted dynamics frames to analyze**, so the
above-100-site end-to-end dynamics/analysis acceptance is explicitly **unmet**.
No budget or force tolerance was increased to relabel the run successful.
Parameterization-only success is not complete long-chain MD support.

## Measurements and the next bottleneck

Wall times in seconds from the final declared matrix:

| Measurement | PE50 GAFF | PE100 GAFF2 | PEO20 GAFF2 |
| --- | ---: | ---: | ---: |
| Local-template construction | 0.404 | 4.134 | 0.151 |
| antechamber | 0.388 | 0.434 | 0.475 |
| parmchk2 | 9.195 | 129.014 | 0.789 |
| tleap | 0.072 | 0.340 | 0.079 |
| Complete preparation | 11.824 | 134.897 | 2.140 |
| Separate reimport including validation | 0.762 | 2.502 | 0.205 |
| Separate preparation-integrity recheck | 0.265 | 0.911 | 0.078 |
| OpenMM binding | 0.519 | 1.307 | 0.181 |
| Fresh single point, including Context | 0.134 | 0.373 | 0.047 |

Stage times are subprocess wall times; preparation total also includes discovery,
serialization, chemistry/lineage, import and repeated integrity checks. Reimport
and integrity entries are separately timed post-preparation operations, not an
additive decomposition of total preparation. Binding includes validation/model
construction; evaluation includes input/output checks, all components and Context
cleanup. The shared host was not an isolated benchmarking machine. Cases differ
in both chemistry and force-field family, so these observations do not establish
a pure asymptotic size exponent.

Artifacts occupied **424857 / 811603 / 206586 bytes**. Linux process high-water RSS
was **118528 / 138272 / 151244 KiB** cumulatively across this sequential script;
child-process high-water readings were **80492 / 118528 / 138272 KiB**. These are
coarse process/fork high-water values, not isolated per-tool allocations or live
memory increments. The workflow attempt took **24.33 s**.

The main measured preparation bottleneck is **parmchk2 on the 602-site GAFF2
case** (129 of 135 seconds). No replacement, cache, or validation bypass was
introduced. Workflow progress is instead limited by the minimization budget and
large starting strain. Those are the next concrete engineering questions; no
long-chain equilibration or automatic charge transfer was added.

## Verification and retained evidence

Final execution files are outside the repository in
`island-validation/phase4f1-final/`: per-case input/charge files, signed preparation
records and tool directories, declared settings, comparisons and workflow failure
artifacts. Earlier attempts remain separately retained: the first exposed a PATH
setup issue; the next completed parameter generation but hit an acceptance-harness
MappingProxy/dict argument mismatch during separate reimport. Only the harness
was corrected, with unchanged charges, seeds, budgets and acceptance tolerances.
The final matrix above actually ran the tools and comparisons successfully.

Ordinary tests include default/opt-in/AM1-BCC limits, invalid values, rejection
before external discovery, noncontiguous stable-ID charge/MOL2 round trips,
Cl/Br names throughout 0..999 and actual halogen MOL2 round trips, contradictory
reconstructed policies/charges, input ownership/failure retention, workflow
propagation and v1/v2/v3 historical contracts. Injected timeouts are software
failure tests, separate from real execution. Existing process-group timeout,
stereo, improper/exclusion/1-4, numerical, checkpoint, analysis and optional-import
regressions remain in the full suite.

Delivery verification:

- Complete ordinary suite: **847 passed, 9 skipped** (70.23 s), including 19
  new size/policy regressions. Skips are the existing opt-in live/archive tests.
- Separately enabled historical acceptance: **4 passed** (18.89 s): five original
  AmberTools reference cases at three configurations each, archived PE workflow
  continuation, archived checkpoint continuation, and both retained PE analyses.
  This verifies historical v2 preparation, config v1, checkpoints and reports
  under their original signatures and tolerances.
- Actual larger-chain script: all three parameterization and independent
  single-point cases passed; the declared larger-chain workflow/analysis gate
  remains unmet as detailed above. This execution is separate from mocked
  failure tests and skipped pytest opt-ins.
- Ruff passed; `python -m pip check` found no broken requirements.
- New example: inspection reported 302 sites; preparation accepted the explicit
  saved nonzero charge file and completed a real PE50/GAFF run. The existing
  short-polymer provided-charge example also ran actual AmberTools successfully.
- Evaluation-session and cross-process checkpoint examples passed. Checkpoint
  split discrepancies remained zero for NVE and Langevin, with matching RNG
  state and 21 split versus 17 uninterrupted calls.

Environment: Python 3.11.16, NumPy 2.4.6, RDKit 2026.03.6, ParmEd 4.3.1,
OpenMM 8.6.1 Reference, SciPy 1.17.1 and AmberTools 24.8. Tool help output did
not expose individual executable version banners; package metadata and executable
hashes are retained instead. GAFF data SHA-256 is
`21ce262ab4af254c2f2b4bf95f4413c2a6593c67faa31ffc7108b7dd07a63d09`;
GAFF2 data SHA-256 is
`14ad62c8e532c47e2e400e2ca6ad8052b33bda4c4b9bc6f49b6a528e527512be`.
No optional scientific package became a mandatory core dependency.

`production_validated=False` and `simulation_readiness="not_established"` remain.
No OPLS/PCFF, transferable-charge templates, periodic cells, packing, crosslinking,
MLIP integration or new dynamics backend is included. Actual parameterization
acceptance reaches 602 sites; 1,000 is a bounded allowed ceiling with name-capacity
coverage, not a claim of live acceptance for every 1,000-site molecule.

### Acceptance CLI correction (Phase 4F2 prerequisite)

The CLI now exits nonzero unless every declared parameterization/numerical
case passes. `--require-workflow` additionally requires completed four-step
propagation and validated five-frame analysis. `report.json` records both
aggregate gates and the requested gate. A later workflow failure preserves
successful parameterization evidence. Missing dependencies cannot pass a gate.
Injected-outcome CLI tests are software tests, not AmberTools acceptance.
