# Phase 4J79 — PDMS literature reconciliation

The decision is **literature_reconciliation_partial**. The supplied Fritz–Hofmann
claims have a deterministic evidence record and a complete J78 crosswalk. The
primary paper was **not available or read**. No reported parameter is promoted
to a source-verified value, a resolved extension term, or a native PCFF row.

## Base, scope and preservation

Branch: `codex/phase-4j79-pdms-literature-reconciliation`.
Exact base: `3c9f9fb5066ac2132d8a9aa46f228a3e53a77ddd`.
Fetched `origin/main`: `594593ee656aace09133e0eba99648225c303ca3`.
The explicit requested base is used; main is untouched. The receipt identifies
the tested implementation commit separately from the final evidence commit.

The frozen inventory covers **10,322 historical tracked files**, including all
J1–J78 evidence, native PCFF implementation/registry, historical receipts and
the C1 quarantine. No historical globs, source hashes, parameter values,
readiness flags or public API implementation change. Unrelated untracked files
and local source assets remain untouched. J72's four unintended SQM attempts
and 46 verified / 48 possible MD steps remain historical violations.

This phase performs metadata access, evidence parsing, unit/charge arithmetic,
schema validation and offline matching only. QM, physical fitting, MD,
converters, scientific evaluators, parameter generation and reserved validation
are not executed. Verification uses the inherited J72-C1 seccomp/ptrace fence
and J73 evidence-only policy before evidence access. Synthetic executables and
imports test denial without invoking a real scientific backend. Real scientific
controls are explicit skips, never passed scientific checks.

## Source access and authority

Lydia Fritz and Dieter Hofmann, “Molecular dynamics simulations of the transport
of water-ethanol mixtures through polydimethylsiloxane membranes,” *Polymer* 38
(1997), 1035–1045. DOI: `10.1016/S0032-3861(96)00600-3`;
PII: `S0032-3861(96)00600-3`.

The unchanged J75 inventory contains 27 files. Their retained citations do not
include this paper. In particular, `main.pdf` is Patil and Mark (2000), a
different work. No other local paper is substituted. Crossref verifies the
bibliographic identity. The exact ScienceDirect PDF endpoint returned HTTP 403.
The publisher XML endpoint returned HTTP 200 **with core metadata only**, no
full text or parameter tables; its plain-text endpoint returned HTTP 400.
Response hashes, exact URLs and outcomes are retained in `source-access.json`.
No original paper, FRC file, archive or other copyrighted source asset is added.
The Crossref TDM license link does not establish general redistribution rights.

Accordingly, the source PDF hash, printed pages/tables/equation locations,
table-cell text and source coefficient units are **null**. Bibliographic page
range is not a parameter-table location. User-reported numbers and equations
are explicitly labeled `user_reported_transcription_not_primary_extraction`.
Their record hashes identify the claims; they are not fabricated paper hashes.
The evidence record cannot be upgraded by supplying an arbitrary hash or
re-signing scientific metadata. A later authorized intake must supply and
independently verify the original artifact and its exact locations.

## Numerical reconciliation

All **20 reported rows / 56 numerical cells** are retained as exact decimal
strings: four charge rows, three quartic bond rows, five quartic angle rows,
four torsion rows and four 9–6 van der Waals rows. None is source-extracted.
The unusual reported C–Si–O `H4=-924.9814` is preserved and flagged for original
table verification; it is not repaired using the C–Si–C value or another force
field. Printed zero claims remain reported claims, never missing-term fills.

Original reported kJ/mol and angstrom quantities remain separate from derived
fields. Angstrom-to-nm uses exactly `0.1`; kJ-to-kcal uses exactly `1/4.184`,
with deterministic Decimal arithmetic. Coefficient energy conversions are
**conditional on verifying kJ/mol in the original table**. Coordinate powers
are retained. These derived numbers are not assignable coefficients.

Equilibrium angles are reported in degrees. Whether the polynomial coefficients
multiply degree or radian deviations remains unresolved. The documented
conditional transformation is `K_rad,n = K_deg,n * (180/pi)^n`; it is not
applied without source evidence. J78's proposed no-half quartic equations do
not establish the paper's equation or prefactor.

Negative torsion amplitudes are retained, including `-1.3513`, `-0.1300` and
`-0.0580`. The paper's complete three-term expression, sign, prefactor,
dihedral zero, phases, multiplicities and reversal convention are unavailable.
No phase is inferred from an amplitude sign, and no Fourier translation is
performed. J78's proposed `Vn*(1-cos(n*phi-phase_n))` is only a comparison target.

The reported charges give exactly zero for `Si O C2 H6`. Applying those same
element charges to each frozen capped HMDS/DP2/DP4/DP6 graph instead yields
**+0.240 e**. Repeat neutrality does not prove cap neutrality. No charge
renormalization or end-group transfer is performed.

## Nonbonded audit

The reported function is preserved exactly:

`Uij = epsilon_ij * [2*(sigma_ij/r_ij)^9 - 3*(sigma_ij/r_ij)^6]`.

Under this equation, the named `sigma` denotes the minimum distance. It must
not be passed as the zero-crossing sigma of a 12–6 LJ model. Source verification
of the reported equation remains necessary.

| ISLAND layer | Existing contract | Reconciliation |
| --- | --- | --- |
| `forcefields.parameters.models.LennardJonesParameter` | `4*epsilon*((sigma/r)^12-(sigma/r)^6)`, kJ/mol and nm | Incompatible; no conversion/coercion into this model |
| Isolated `forcefields.pcff.model` / `terms` | 9–6 minimum-distance equation; sixth-power mixing | Symbolic powers exist, but paper mixing and policies are unverified; no evaluator call |
| J78 extension inventory | Proposed 9–6 target with explicit pair weights | Symbolic form compatible with the reported expression; source conventions and qualification blocked |

The isolated PCFF code's sixth-power rule is
`rij=((ri^6+rj^6)/2)^(1/6)` and
`epsij=2*sqrt(epsi*epsj)*ri^3*rj^3/(ri^6+rj^6)`.
This is recorded as **ISLAND's existing rule**, not attributed to the paper.
The paper's exact mixing rule, dielectric, cutoff, tail treatment, electrostatic
method, exclusions and 1–4 scaling remain null with explicit evidence requests.
No convention is inferred from COMPASS, CFF91, OPLS, GAFF or another PDMS model.

## Complete J78 crosswalk and remaining gaps

The crosswalk covers all **129 unique terms and 5,135 site instances**, retaining
the original parameter identities, source order, stable site IDs, graph
identities and equilibrium-dependency references. Matching an element sequence
in either direction only identifies a candidate; it proves neither source atom
typing nor exact chemical environment, orientation equivalence or cap transfer.

| Crosswalk status | Unique terms |
| --- | ---: |
| `literature_candidate_exact_form` | 0 |
| `literature_candidate_convention_blocked` | 20 |
| `end_group_unresolved` | 21 |
| `cross_term_unresolved` | 87 |
| `not_covered_by_source` | 1 |

Status precedence is cross term, end group, convention block, then uncovered.
Orthogonal flags preserve overlap: **69 terms involve end groups**, including
48 of the 87 cross terms and the 21 non-cross end-group terms. All 20 reported
source rows find an element-sequence candidate somewhere in the graphs;
**zero exact chemical-environment matches are source-verified**. The uncovered
non-cap/non-cross term is H–C–H. There is no automatic source-type aliasing.

The seven required class groups are each explicitly unresolved: bond–bond,
bond–angle, angle–angle, torsion–bond, torsion–angle, improper/out-of-plane and
torsion–torsion. They are absent from the supplied claims; without the full
paper, their absence from the paper itself cannot be asserted. No missing
coefficient is assigned zero. J78's structural improper nonapplicability and
unresolved torsion–torsion policy remain unchanged.

All **129 terms, 365 scalar parameters and 11 policies remain unresolved**.
The 423 J78 QM requests remain plans with zero completions. Native PCFF's
378 outstanding requests and zero exact resolutions remain unchanged.
The candidate retains its proposed `profile_origin="qm_fitted_extension"`;
literature claims are not QM-fitted observations or coefficients.

## Artifacts, validation and replay

Evidence is under `docs/evidence/phase_4j79/`:

- `literature-evidence.json`: exact reported strings and explicit missing source fields.
- `j78-literature-crosswalk.json`: all unique terms and instance references.
- `convention-compatibility.json`: nonbonded, Fourier, angular and policy audit.
- `unresolved-terms.json`: complete remaining term IDs and evidence requests.
- `source-access.json`, `local-source-search.json`: bibliographic access outcomes.
- `j78-input.json`, `model-inputs.json`, `frozen-inputs.json`: hash-pinned replay inputs.
- `historical-files.json`: explicit historical inventory.
- `reconciliation-report.json`, `receipt.json`, `manifest.json`: deterministic results and seal.
- `validation/`: terminal JUnit/collection, check logs, every fenced attempt and tested Git blobs.

The explicit implementation inventory contains this document, the J79 builder,
strict receipt validator and focused tests only. No library/runtime module is
changed. Replay needs no network, source PDF, PDF parser or scientific backend.
The validator pins frozen inputs and regenerates semantics, so re-signing an
invented source, a page/table, a value/unit or a readiness change fails closed.

Use a fresh output directory for each fenced invocation:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j79-replay -- -m scripts.phase_4j79_literature replay
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j79-validate -- -m scripts.validate_phase_4j79_receipt
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python -m scripts.run_phase_4j72_c1 --output /tmp/j79-tests -- -m pytest -p scripts.phase_4j73_scope tests/test_phase_4j73.py tests/test_phase_4j74.py tests/test_phase_4j75.py tests/test_phase_4j76.py tests/test_phase_4j78.py tests/test_phase_4j79.py
```

Focused tests cover exact transcription, provenance and semantic mutation,
repeat/cap charge arithmetic, 9–6/12–6 separation, units, unresolved Fourier
metadata, end groups, cross terms, complete identity coverage, runtime
non-registration, child-process denial and relocated fresh-process replay.
J73–J78 evidence-only regression tests run under the same fence. Historical
J78/J76/C1 receipts are independently validated; Ruff, compileall, supported
environment pip check and diff checks are recorded. The receipt provides exact
unique terminal counts and skips rather than adding overlapping development
runs. Earlier software probes (missing venv Ruff module, initial lint findings
and Decimal reciprocal-rounding failure) are retained separately with their
outcomes. The rounding issue was corrected by direct Decimal division before
terminal verification; none of these attempts launched science.

Final flags remain `runtime_registered=false`,
`automatic_assignment_enabled=false`, `accepted_coefficients=0`,
`fitting_attempted=false`, `production_validated=false`,
`simulation_readiness="not_established"`. All J79 prohibited execution counters
are zero. Source reconciliation, source authority, native completeness,
chemical transfer, scientific validation and production readiness remain
separate conclusions.
