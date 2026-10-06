# Phase 4J10 — source-row accounting, final-graph parity and bounded validation

## Repository and acceptance status

Branch: `codex/phase-4j10-pcff-msi2lmp-parity`.
Base: merged J9 squash `0f8f7e1` on origin/main, with tree identical to reviewed
`c20d5ac2063cfe3786fdefec217d7b6a0f5dbe7e`. Main is not changed or merged.

**The requested exhaustive operational acceptance is unmet.** This phase does
not make the full-source CLI return zero. The pinned file lacks required rows
for some declared domains, has unresolved aliases/permutation semantics, and
prints a guanidinium total of 0.9999 e. Manufacturing coefficients, treating
missing couplings as zero, or normalizing that charge would contradict the
requested strict contract. Those cases remain explicit failures.

Implemented work advances executable fused-polymer coverage, traces every source
record, and removes the repeated semantic/JSON workload exposed by J9. It does
not relabel bounded coverage as all-source completion.

## Exact source and primary implementation audit

Operational source:

- LAMMPS `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`.
- `tools/msi2lmp/frc_files/pcff.frc`.
- SHA256 `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.

Source libraries/executables remain external. Neither runtime assignment nor
bundle/workflow reconstruction needs CAR, MDF, msi2lmp or LAMMPS. Runtime source
selection still requires the exact native hash. No PCFF-IFF promotion or source
mixing is introduced.

`audit_pcff_source_rows.py` uses an independently authored flat line/Decimal
reader and checks row IDs against the semantic/numerical native catalogs. It
accounts for **5,319** versioned records, including 21 force-field definition
records, 133 atom types, 134 ordinary equivalences, 108 automatic equivalences,
564 bond increments, and all numerical families/namespaces. Reference/comment
text is retained in the source inventory rather than counted as coefficients.
Empty sections, including torsion-torsion, remain explicit.

`phase_4j10_source_rows.json` binds every row to source hash, section, namespace,
line, raw-line hash, version/reference, parser status, rule availability,
charge/matching route, numerical implementation, bounded observations and
remaining reason. Row presence is never global independent verification.
The native type ledger still has 86 implemented labels and 47 unresolved labels.

Pinned primary implementation evidence:

- `tools/msi2lmp/README:18–28`: CAR/MDF supply coordinates, topology and atom
  types. The converter does not establish chemical atom types from PSMILES.
- `README:140–149`: no automatic-equivalence supplementation or bond increments.
- `src/GetParameters.c:1055–1084`, `find_match`: exact records are scanned first,
  then wildcard records; the first matching file-order candidate is returned.
- `GetParameters.c:1087–1170`, `match_types`: each candidate is tried forward,
  then reversed unless reversal was disabled; any pattern beginning with `*`
  matches. A numeric suffix is not an authoritative chemical priority.
- `GetParameters.c:1241–1277`, `get_equivs`: distinct columns for nonbonded,
  bonds, angles, proper and out-of-plane families.
- `GetParameters.c:1171+`, `get_r0`, and callers at 212–213/374–376: equilibrium
  lengths are resolved separately for Class II dependencies.
- `GetParameters.c:735–787` and `find_angleangle_data`: the existing generic
  reversal can move the angle-angle center; it is not adopted as an ISLAND rule.
- `ReadMdfFile.c:388–396`: optional connection bond-order decoding.
- `WriteDataFile.c:391`: converter output rounds coordinates to nine decimals
  and charges to six. Independent numerical inputs install declared coordinates
  explicitly; raw converter files and the H5 angle-equilibrium overrides remain.

These source locations are at the pinned revision, not a moving develop URL.
msi2lmp's file-order behavior is reproducible converter behavior, not a resolution
of physically contradictory wildcard candidates. Native conflict/reversal
checks remain fail-closed. No `-ignore` is used.

## PSMILES-first native path and source-complete extension

The J9 final-graph workflow remains:

`PSMILES -> final MolecularSystem -> native perception/charges -> Class II
assignment/model -> evaluator -> minimum/bundle/workflow`.

Stable IDs, repeat/chain/end roles, explicit ends, branch/crosslink metadata,
chemical state and coordinate provenance remain owned and persistent. Interaction
inventories and equilibrium dependencies come from final authoritative bonds.
Topology edits invalidate prior binding; the J9 local update API retypes changed
neighborhoods and validates/rebuilds the full final graph. Its tests and historical
prepared-workflow attempt-authority tests remain in the ordinary suite.

New operational profile: `island_pcff_fused_benzenoid_graph_v1`. It authorizes
neutral explicit standard-mass C/H benzenoids where every aromatic carbon is in
a smallest six-membered aromatic cycle, including fused rings and linked repeats,
plus the already justified methyl-end environment. Five-ring, heteroatom,
charged/radical/isotopic and other saturated environments remain rejected.
Native cp/hc/c3 rules and source equations are unchanged; no atom-label count is
inflated by this graph-domain authorization.

The declared naphthylene PSMILES is
`[*:1]c1ccc2cc([*:2])ccc2c1`, DP2, seed2026, 34 explicit sites. It completes
native component charges and every active source-backed coupling, including
proper, Wilson-zero-equilibrium and angle-angle terms. No policy zero enters
the strict preparation. The unchanged J8/J9 profiles retain their narrower
behavior and identities.

## Practical integrity validation

J9 profiling exposed 134 million function calls for a profiled cold load plus
one repeated check: 58.76 seconds cold and 16.24 seconds for one hot check.
1,036 FRC parses and repeated full native recomputation/JSON work dominated.

The process-local validation cache contains only successful validation of
immutable JSON/source bytes, typed frozen chemical graph data and original
content identities. It never contains a live system, Context or external handle.
Every entry is keyed by record type, complete content, nested assignment,
source bytes/hash and trusted rule/constants/function context. Failed validations
are never cached. Rechecksummed changes miss the cache and rerun the original
semantic recomputation. Every caller system is checked anew for topology,
stable IDs, elements, masses, isotope/radical/stereo state; coordinate changes
remain compatible. Rule changes invalidate reuse. A compact owned authorization
summary avoids repeatedly decoding the complete source catalog.

Limits: 16 entries, 64 MiB, LRU eviction. Oversized records remain validated but
uncached. Source parse reuse is separately bounded to eight immutable sources;
public inventories are copied. No disk cache is trusted during reconstruction.
Historical signatures, numerical/evaluator fingerprints and file schemas are
unchanged. Public cache inspection/clear APIs are data-only and dependencies
remain lazy.

The declared unprofiled budget is cold load <=60 seconds and 20 repeated checks
<=5 seconds. Measured cold load was 2.98 seconds; 20 checks took 0.714 seconds.
Profiling overhead measurements are reported separately rather than presented
as a directly comparable unprofiled speedup.

## Broad diagnosis and retained failures

The predeclared matrix includes phenylene, fused naphthylene, PEO, amine,
thioether, chloro-substituted chains, DP20 benzenoid, pyridinium and guanidinium.
J9's end-group, branch, crosslink and mixed-sequence receipts remain unchanged.

Actual matrix outcomes:

- Phenylene and fused naphthylene: complete strict source-backed preparations.
- PEO and thioether: typing/native charges complete; required BB13 rows missing.
- Chloro-substituted chain: typing/charges complete; BB/BA and other couplings
  remain missing after the declared searches.
- Amine: na–hn2 increments still incomplete; no hn substitution.
- Pyridinium: nh+ charge coverage still incomplete; no nh substitution.
- Guanidinium: printed source total 0.9999 e; formal +1 check fails at 1e-12.
- Original DP20 ETKDG construction failed; that result is preserved. A separately
  declared local-template construction with fixed template/assembly seeds2026
  yielded 202 sites and a strict source-backed model. This is graph/inventory
  scalability evidence, not long-chain minimization/dynamics or charge accuracy.

`phase_4j10_parity.json` retains exact interaction requests, site IDs, candidates,
equivalence paths and charge/component evidence. A structurally complete diagnostic
model with documented historical policy zeros is not a strict preparation.
Diagnostic mode returns no executable facade. All requested major-domain minima
and workflows cannot be accepted while their required native inputs are incomplete.

## Independent fused-domain acceptance

A separate declaration fixes input checksum, raw-source reader, pinned reference
executables, asymmetric perturbation, FD sizes and budgets before execution.
The reference chooses raw FRC records independently and constructs its own
converter topology; ISLAND's compiled coefficient list is not the oracle.
It retains the H5 independently reconstructed angle-angle equilibrium correction.

Actual DP2 native model:
`40dd86aed3c5a482ed43a37988ca060592c173ad51d3a8593c426308d5897ff0`.
Profile identity:
`078d850a2beaff135f850e4a4b27ff7291f582fe86fe94ef2b7616ed0ff2bd06`.
Prepared identity:
`64f8e4eb9dd037c27e6f8fefd6d99ed55de790af13efe46e9d6f44ed657d6f03`.

Independent initial/asymmetric/final energy and force comparisons, FD convergence,
source selections and native charges passed. Minimization converged in 148
iterations/158 evaluations: 1479.5061521726466 ->1377.4828656514414 kJ/mol,
maximum force 0.09847393671189145, RMS 0.052162012712675274 kJ/(mol*angstrom),
independent fresh verification passed. Four-step split/uninterrupted continuation
relocated bundle/run directories and resumed in a separate process. All measured
common-frame differences were zero; complete PCG64 state matched exactly.
Completed offline status/frame reading/no-op resume and semantic tamper gates pass.
This establishes model and continuation fidelity within this domain, not universal
PCFF science or independent thermostat/integrator validation.

Unchanged tolerances: energy/force atol1e-5 in canonical units, rtol2e-10;
source coefficients/native charges1e-12; FD displacements1e-4/1e-5/1e-6 angstrom;
split atol1e-10/rtol1e-12. Minimum budgets5000 iterations/10000 evaluations,
force criterion0.1. Four-step BAOAB, two per segment, 0.1fs,300K,5/ps,
velocity seed78123, thermostat seed99181; evaluations4 per segment/6 whole.

## Reproduction

Primary environment: `../island-validation/phase4e2-env/bin/python`.
FRC: `../island-validation/phase4h1-sources/lammps/pcff.frc`.
External artifact root: `../island-validation/phase4j10-declared`.

```sh
python scripts/audit_pcff_source_rows.py --source /verified/pcff.frc --output /new/rows.json --require-full-source
python scripts/inspect_pcff_parity_matrix.py --source /verified/pcff.frc --output /new/parity --require-full-source
python scripts/validate_pcff_profile.py --source /verified/pcff.frc --declaration docs/evidence/phase_4j10_fused_declaration.json --output /new/fused-vertical
python -m pytest -q
ruff check .
python -m pip check
```

The first two commands intentionally return nonzero while required all-source
coverage remains incomplete. The successful fused vertical-slice command is
separate from that gate. Full-source completion, major-domain minima/workflows,
remaining source aliases, metals/ions/zeolites, missing increments/couplings,
nonzero Wilson interpretation and empty torsion-torsion applicability remain
explicitly unmet. No missing parameter is invented or normalized.

`production_validated=False`; `simulation_readiness="not_established"`.

## Executed verification and preserved evidence

Final ordinary suite: **1,652 passed, 10 skipped**, 290.17 seconds. Focused
PCFF/cache/prepared-workflow checks:60 passed; combined release checks28 passed;
cache-specific final checks5 passed. Ruff and both primary/isolated pip checks
pass. The existing PSMILES single-point example ran against the exact real source.

Matched cProfile measurements after the final authorization-summary change:
5.30 seconds cold and0.0865 seconds for one hot check, compared with the original
58.76/16.24 profiled measurements. The separately declared unprofiled budget
passed at2.98 seconds cold and0.714 seconds for20 checks. Cache limits and live
binding/tamper checks were tested; no scientific tolerance was relaxed.

Fused numerical maximum grouped energy error was1.13687e-12 kJ/mol and maximum
force-component error3.79217e-11 kJ/(mol*angstrom). FD errors for the declared
step sizes were3.62349e-5,4.16946e-7 and5.76351e-7 in force units. All used the
unchanged declared tolerances. Source-row ledger observations remain bounded to
this case, not global row/type/domain completion.

Independent raw-source charge calculations separately confirm the pyridinium,
guanidinium and amine failures; the exact missing searches, endpoint rows and
Decimal component totals are retained. No source-charge correction was applied.

All2,112 files in the pre-execution snapshot remain present and byte-identical.
J8/J9 relocated bundles and completed workflows revalidated read-only with
scientific imports blocked, preserving native/profile/facade identities and five
frames. Historical hashes, schemas, failed receipts and checkpoints are not
rewritten. Generated systems, models, source libraries and reference outputs
remain external. Exact receipts/commands/checksums are in
`phase_4j10.json`, `phase_4j10_verification.json` and the external artifact root.

Required full-source and all-major-domain operational gates remain **false**.
This is useful reviewable progress, not an exhaustive PCFF completion claim.
