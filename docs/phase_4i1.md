# Phase 4I1 — explicit unified force-field preparation

## Base and scope

`codex/phase-4i1-unified-forcefield-api` starts from
`d0fc10ad45b47de955777064817d060aa70ffa79` (updated main). Its tree was verified
identical to reviewed H7 `2f5d3eae0d0075d3f36e92b3ffedaa76ee974547` after
squash merge. No prerequisite was merged locally.

This additive facade unifies invocation, ownership and validation. Existing
GAFF/GAFF2, OPLS and PCFF numerical code and historical contracts are unchanged.
It does not extend the Amber-only high-level workflow bundle, migrate saved runs,
or alter dynamics, checkpoints, minimization or analysis schemas.

## Public API

All names are exported from `island.forcefields`:

```python
from island.forcefields import (
    ForceFieldRequest, PCFFOptions, prepare_forcefield, create_evaluator,
)

request = ForceFieldRequest(
    "pcff",
    PCFFOptions("/installed/pcff.frc", lj=(0, 0, 1), coulomb=(0, 0, 1)),
)
prepared = prepare_forcefield(system, request)
evaluator = create_evaluator(system, prepared)
with evaluator.open_session() as session:
    result = session.evaluate()
    independent = session.evaluate_fresh()
```

The system must already be constructed and chemically valid. The example
`examples/forcefield_singlepoint.py` demonstrates PSMILES construction with local
templates and explicit seeds followed by this path. No charge method is selected
implicitly.

| Selector | Required options | Charges and policy |
|---|---|---|
| `gaff` | `AmberToolsOptions("gaff", explicit_method, ...)` | Explicit `provided` mapping or explicitly requested `am1bcc`; existing Amber policy |
| `gaff2` | `AmberToolsOptions("gaff2", explicit_method, ...)` | Same explicit choice; selected GAFF2 data |
| `oplsaa` | `OPLSOptions(local_xml_path)` | Pinned source-native library charges; geometric LJ, native RB, source 1–4 scaling |
| `pcff` | `PCFFOptions(local_frc_path, lj=..., coulomb=...)` | Native bond increments; Class II/9–6 model and mandatory caller pair policy |

`ForceFieldRequest` is a small discriminated record, not a keyword forwarding
interface. Cross-family option objects, unknown selectors, GAFF/GAFF2 mismatch,
and invalid reconstructed options fail before scientific execution. Options own
provided mappings and policy sequences. Backend-specific native validators retain
size, timeout, source, charge, chemistry and coverage checks. There is no fallback.

Amber keeps its default 100-site ceiling, explicit provided-charge opt-in up to
1,000, and whole-molecule AM1-BCC limit of 100. Source files are explicitly local;
no facade operation downloads or updates a source. Every evaluator selected here
uses its existing Reference-platform finite nonperiodic contract.

### Source/profile distinctions

- Amber resolves the selected installed data and records its existing tool/data
  provenance. Facade profiles are `island_ambertools_gaff_adapter_v1` and
  `island_ambertools_gaff2_adapter_v1`; these are adapter versions, not invented
  semantic releases of GAFF data. Native signed records identify the actual data.
- OPLS uses `pinned_foyer_oplsaa_parameters_v1`, Foyer revision
  `dd2f6eaa0ec271432ccd0c5c729f17a3ea3364bd`, XML SHA256
  `c78ccb763cda33e3456a3f8b4ed8a8f0361b10c92c33aa02f2a9abf618c163e4`.
  It remains the selected Foyer OPLS-AA implementation. Native-charge failures,
  including the historical PS DP3 case, are not normalized or repaired.
- PCFF retains `island_pcff_acyclic_cho_v1` typing and
  `island_lammps_pcff_acyclic_cho_v1` operational model, LAMMPS revision
  `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`, FRC SHA256
  `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
  Its ordinary-alcohol oh/ho convention, policy-derived zeros, coefficient roles,
  native increments and restricted neutral saturated acyclic explicit-H C/H/O
  domain remain unchanged. Explicit special-pair weights are caller choices,
  not universal commercial PCFF defaults.

OPLS library charges, PCFF increments and GAFF provided/AM1-BCC charges are
separate native contracts. Provided charges alone do not establish compatibility
with another force field.

## Adoption and integrity

```python
from island.forcefields import adopt_forcefield

# Amber: original preparation system and complete validated preparation.
prepared = adopt_forcefield(original_system, "gaff2", amber_preparation)
# OPLS: existing native resolved parameters plus pinned owned source.
prepared = adopt_forcefield(system, "oplsaa", opls_parameters, source=opls_source)
# PCFF: existing H4 specification carrying its H3 assignment and source.
prepared = adopt_forcefield(system, "pcff", pcff_specification)
```

Adoption invokes native integrity checks, never an external command, typing job or
new parameterization. It does not accept a bare type/charge result as a complete
potential. PCFF model-definition incompleteness fails with native diagnostics.
Raw H3 source-row coverage and complete H4 operational-model definition remain
distinct; H4's documented policy-derived zeros do not rewrite H3 missing rows.

`PreparedForceField` owns a copy of the original system, native result and source.
`native_result`, `source` and `metadata` return owned data. Its public
`validate_integrity(system=None)` revalidates native content and recomputes all
facade metadata. The canonical `island_prepared_forcefield_v1` description binds
adapter, selector/profile, source, charge method, nonbonded policy, native identity,
chemical binding, status and limitations. Reconstructed contradictory metadata or
native data is rejected. `identity` hashes this description and must never be
compared as though it were a native parameter/model fingerprint.

The status `prepared` means the native preparation/model contract passed. It is
not a scientific validation flag. Changing caller objects or copies returned by
accessors cannot change the accepted preparation or an already bound evaluator.

For Amber, the saved preparation evidence is validated against its original
coordinates. `create_evaluator(new_coordinate_system, prepared)` validates the
unchanged imported parameter record against the supplied system, preserving its
existing coordinate-replacement contract. This does not re-sign old preparation
coordinates. OPLS/PCFF coordinate-independent records continue to permit coordinate
replacement without reparameterization. Graph, element, atomic number, mass,
stable IDs and assigned stereo changes are rejected. Native checks additionally
retain their backend-specific compatibility restrictions.

There is no new facade restart/persistence format. Persist native records through
their established contracts, reconstruct compatible external inputs, then adopt.
The small facade description is inspection metadata, not sufficient standalone
restart data. No `ImportedAmberResult` is fabricated for OPLS/PCFF and no Class II
term is coerced into a harmonic/LJ12–6 `ParameterizedSystem`. The existing abstract
`ForceFieldBackend.parameterize()` contract returns that historical container;
this additive API deliberately does not change that contract or pretend every
native model has that representation.

## Errors and dependencies

- `ForceFieldRequestError`: invalid selector or backend options, before execution.
- Native backend errors propagate from actual preparation, preserving missing
  dependencies/files, unsupported chemistry, site IDs and coverage diagnostics.
- `PreparedForceFieldError`: invalid/incompatible adoption or reconstructed facade,
  with the original native exception chained when applicable.
- Native evaluator errors retain their existing types and causes.

Importing the facade or inspecting requests imports no RDKit, Foyer, ParmEd,
OpenMM or SciPy and executes no tools. Actual preparation requires only the chosen
backend's dependencies. PSMILES construction separately requires its existing
chemistry stack. Native reconstruction retains its existing validation dependencies;
no dependency requirement is bypassed by adoption. Evaluation requires OpenMM.
Foyer preparation uses its separate compatible environment; it is not installed
by downgrading the main environment.

Native integrity checks can be expensive, particularly for PCFF. Public facade
access/dispatch deliberately retains those checks; it does not add a validation
cache or skip catalog validation. Numerical evaluations use the existing compiled
owned evaluator/session and do not revalidate the facade on each force call.

## Bounded acceptance and reproduction

`scripts/validate_forcefield_facade.py` exclusively creates an external output
directory, publishes its declaration, prepares PE DP3 at fixed local-template
seeds 2026/2026, adopts the native result, and compares both facade evaluators with
a direct native evaluator. It compares energy, components, every force, evaluation,
parameter and model identities, plus fresh/session behavior and input nonmutation.
Tolerances were declared as atol=1e-10 and rtol=1e-12; no scientific accuracy
threshold is inferred. Failures are recorded and exit nonzero; no automatic retry
or altered scientific settings occurs.

```bash
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
$PY scripts/validate_forcefield_facade.py --family pcff --source "$FRC" \
  --output /new/pcff-acceptance
$PY examples/forcefield_singlepoint.py --family pcff --source "$FRC" \
  --lj 0 0 1 --coulomb 0 0 1

# GAFF/GAFF2 require all three existing executables on PATH, not just AMBERHOME.
# Use absolute installed paths; retain parameterization artifacts outside the repo.
PATH=/installed/ambertools/bin:$PATH $PY scripts/validate_forcefield_facade.py \
  --family gaff2 --amberhome /installed/ambertools --output /new/gaff2-acceptance
# Repeat with --family gaff for the GAFF selector.

../island-validation/phase4g1-env/bin/python scripts/validate_forcefield_facade.py \
  --family oplsaa \
  --source ../island-validation/foyer-4g1-upstream/foyer/forcefields/xml/oplsaa.xml \
  --output /new/opls-acceptance
```

Actual execution (all four fresh preparation and adoption paths passed):

| Selector | Status | Energy (kJ/mol) | Max facade/native discrepancy | Seconds |
|---|---|---:|---:|---:|
| gaff | passed | 66.7403330303 | 0.0 | 1.64 |
| gaff2 | passed | 57.1122362678 | 0.0 | 1.62 |
| oplsaa | passed | 65.9376756825 | 0.0 | 11.10 |
| pcff | passed | 6.8473697140 | 0.0 | 127.23 |

All fresh/session energy, component and force comparisons were exactly equal on
this host, including native parameter/model/evaluation identities. This measures
adapter fidelity, not independent physical-model accuracy. GAFF/GAFF2 used explicit
alternating ±0.01 e synthetic charges solely as a software fixture; no SQM command
was executed. OPLS and PCFF used their native source charges.

Successful outputs are retained under `../island-validation/phase4i1-final/{gaff,gaff2}`,
`phase4i1-installed/oplsaa`, and `phase4i1/pcff`. Source identities, declarations,
original artifacts and result hashes are recorded in
[the evidence manifest](evidence/phase_4i1.json).

The main environment used Python 3.11.16, NumPy 2.4.6, OpenMM 8.6.1, ParmEd 4.3.1
and RDKit 2026.03.6. The separate OPLS environment additionally supplied the pinned
Foyer 1.2.0 implementation. AmberTools was the existing installation at
`../miniforge3/envs/md_rg`; signed preparation files retain exact executable hashes,
package/version responses and GAFF/GAFF2 data identities. Both environments passed
`pip check`.

Timings include preparation, repeated native integrity checks, adoption, direct
and facade evaluations, and fresh/session comparisons. PCFF's deeper catalog
validation is included; these timings do not establish a numerical-engine speedup.

Early setup failures are retained at `../island-validation/phase4i1/gaff`, `gaff2`
(missing executables on PATH) and `oplsaa` (incorrect local source path). They did
not execute scientific preparation. The first installed Amber executions in
`phase4i1-installed/gaff` and `gaff2` generated valid native preparations, then the
acceptance script failed while trying to encode the complete dataclass including
tuple-keyed interaction maps with the workflow metadata encoder. The script was
corrected to use the established preparation reconstruction payload (record,
original signatures and import provenance); no scientific engine changed. Those
failed outputs remain, and new bounded executions used new directories and the
same charges/settings. This was an explicitly retained tooling correction, not
charge-method fallback or normalization.

AM1-BCC selection/routing was verified using labelled software injection only;
no new QM was launched. Synthetic unit fixtures also cover typed option rejection,
source mismatch, native-charge/unsupported-chemistry failure propagation,
GAFF/GAFF2 mismatch, reconstruction tampering, ownership, changed chemical binding,
coordinate replacement, and optional import isolation. Real interoperability is
reported separately from those software tests.

Final verification:

- Complete ordinary suite: **1,306 passed, 10 skipped**, 288.18 s. The existing
  opt-in/dependency skips remain separate from the real four-backend execution.
- Focused facade tests: **25 passed**, 16.97 s.
- Ruff and whitespace checks passed; `pip check` passed in both the main and
  separate OPLS environments.
- The new PCFF PSMILES-to-single-point example executed successfully.
- All four real fresh-preparation/adoption/session gates passed; AM1-BCC routing
  remains an explicitly mocked test, not new QM acceptance.
- Acceptance artifact checksums were rechecked before recording this report.

A preliminary suite exposed a test-fixture error that assumed site ID 0; the
fixture uses noncontiguous IDs. The test was corrected to select an actual site,
and the full suite above was rerun. No native engine correction was required.

No requested bounded backend interoperability gate remains unavailable in this
environment. Scientific validity, expanded chemistry, AM1-BCC research and unified
high-level workflow bundles remain outside this phase.

## Limits

This milestone is orchestration fidelity for four existing native implementations.
It adds no source, chemistry, force-field mathematics or evaluator. It does not
claim all polymers are covered, extend the Amber workflow bundle, or establish
scientific accuracy, equilibration or production readiness.

`production_validated=False`, `simulation_readiness="not_established"`.
