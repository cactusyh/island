# Phase 4E1: nonperiodic single-point energy and forces

The branch `codex/phase-4e1-singlepoint-evaluation` starts at current main
`051458036fddac93c6a0aeed666855e8975b981f`: Phase 4D1 was merged in PR #7,
then Phase 4D2 including reviewed `bc929fd` in PR #8. The reviewed feature
commit and this main have identical trees. Main and checkpoint directories
were not edited.

## Contract and ownership

Install `pip install '.[evaluation,amber]'` to run the example and import an
Amber topology. OpenMM is optional and imported only on backend construction.
Importing `island.evaluation` also works without OpenMM, RDKit, or ParmEd.
The parameter binding and evaluation paths do not import ParmEd or read prmtop
files. ParmEd is needed only when initially importing the external topology.

```python
from island.evaluation import OpenMMSinglePointEvaluator

potential = OpenMMSinglePointEvaluator(system, imported_result, platform="Reference")
first = potential.evaluate()  # original system coordinates, in angstroms
second = potential.evaluate(new_positions_by_site_id, coordinate_unit="angstrom")
# Or check the graph as well as coordinates of a new MolecularSystem:
third = potential.evaluate_system(new_system)
```

`PotentialEvaluator` is a small structural protocol; `EvaluationResult` is an
immutable backend-independent result. It contains potential energy, components,
forces keyed by stable ID, explicit units, coordinate/parameter/model/evaluation
fingerprints, backend version, actual platform and calculation settings. It does
not establish a second parameter-assignment framework or implement MLIP.

The evaluator validates the actual `ImportedAmberResult` signature and graph,
then reconstructs parameter and torsion records to check their numerical values,
units and functional forms even after malicious re-signing. It validates the
nonbonded policy and source pair inventory, and owns copies of the authoritative
system and parameters. OpenMM is built directly from those records. Each frame
uses a private Context reconstructed from the owned model; no mutable Context
is exposed or reused across callers. A Context requires an integrator object,
but its `step()` method is never invoked.

`from_preparation(original_system, preparation_result)` first validates the
preparation against its original coordinates. Afterwards new evaluation frames
are independent of that historical coordinate digest. AmberTools is never rerun.
`from_parameterized_system(snapshot, imported_result)` requires the signed import
and compares all actual assignments and policy to it. An aggregate signature
alone is insufficient. Snapshot metadata is not trusted as parameter evidence.

`evaluate_system()` checks graph identity, stable IDs, masses and recorded
stereochemical labels. Coordinate-only mappings intentionally refer to the bound
graph. Mutating caller parameters cannot modify an existing model; changed
parameters require constructing a new evaluator. Fingerprints bind canonical
sorted IDs and explicit coordinate units; moving coordinates changes the
coordinate and evaluation hashes while preserving the model hash. The parameter
fingerprint includes signed source provenance, not just numerical coefficients.
These are content hashes, not attestations of scientific correctness.

## Physical model and units

The implementation uses native harmonic bonds/angles, separate periodic proper
and ordered improper forces, and a combined LJ/Coulomb `NonbondedForce`.
It uses nonperiodic `NoCutoff`, Lorentz–Berthelot mixing, fixed charges, and
exactly one exception per unique source exclusion pair. Source 1–4 pairs use
separate LJ and Coulomb scaling; 1–2/1–3 pairs are excluded. Exceptions are not
created per torsion term, so multi-term torsions and rings do not multiply them.
Improper order is preserved even when its central atom is second. Zero LJ does
not remove electrostatics.

No constraints, implicit solvent, switching, dispersion correction, PME,
restraints, HMR, minimization or integration steps are enabled. Periodic systems,
virtual sites, pair overrides, unsupported record families and mixing policies
are rejected. The snapshot entry point also rejects extra interaction families.

| Quantity | Public/storage convention | OpenMM conversion |
| --- | --- | --- |
| Coordinates | angstrom | multiply by 0.1 to nm |
| Energy | kJ/mol | unchanged |
| Force | kJ/(mol*angstrom), **F = -dE/dR** | multiply kJ/(mol*nm) by 0.1 |
| Bond equilibrium length | stored nm | unchanged |
| Bond force constant | stored kJ/(mol*nm²) | unchanged; energy is `0.5*k*(r-r0)^2` |
| Angle equilibrium and torsion phase | stored degrees | multiply by pi/180 |
| Angle force constant | stored kJ/(mol*rad²) | unchanged; energy is `0.5*k*(theta-theta0)^2` |
| Periodic torsion | kJ/mol | sum `k*(1+cos(n*phi-phase))` over every term once |

The Amber importer has already applied the Amber-to-ISLAND factors. Applying
those factors again in the evaluator would be an error. OpenMM's conventions
are documented in its [standard force definitions](https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html).
Bare arrays and non-angstrom unit labels are rejected. Mappings must have exactly
the bound stable IDs and finite numeric N×3 coordinates. Coincident sites within
`1e-10` angstrom (including excluded pairs), and assigned angular/torsion planes
with normalized cross-product magnitude at most `1e-12`, are conservatively
rejected. Extreme arithmetic and nonfinite or malformed backend outputs raise
`EvaluationError`; invalid inputs use `EvaluationInputError`. Missing OpenMM or
platform availability uses `EvaluationUnavailableError`. No successful result
contains NaN/Inf. Collinear configurations are outside this first evaluator even
where a special analytic limiting force could be defined.

## Independent acceptance

`tests/test_singlepoint.py` constructs path A from resolved ISLAND records and
path B from **OpenMM's original-prmtop reader**, with the same nonperiodic and
unconstrained settings. Path B is never built from path A's parameters. It covers:

- Three perturbed frames per synthetic chain/improper case and external phenol.
- Nonzero charges, nonzero phases, multi-term torsions, separate 1–4 factors,
  ring exclusions, a charged zero-LJ site, and source-index/insertion reordering.
- Bond, angle, combined proper+improper, nonbonded and total energy, and each
  atom's force. OpenMM's reader combines torsions; separate analytical improper
  checks preserve a central atom at both the second and third positions.
- Independent source formulas for harmonic terms, torsions, isolated LJ/Coulomb
  pairs, analytic pair forces, and LJ isolation using an uncharged source.
- Central finite differences at 1e-4 and 5e-5 angstrom, rotation/translation
  covariance, zero net force, frame/result isolation and immutable outputs.
- Re-signed malformed parameters, mutated snapshots, unsupported policies,
  missing/malformed coordinates, singular geometries and invalid backend output.

Whole-system tolerances: energy/components `rtol=2e-10, atol=2e-7 kJ/mol`;
forces `rtol=2e-9, atol=2e-6 kJ/(mol*angstrom)`. Analytical LJ comparisons allow
1e-9 kJ/mol for Amber coefficient serialization. Finite differences use
`rtol=2e-6, atol=2e-5 kJ/(mol*angstrom)` and must agree at both step sizes.
The older phenol fixture retains its unknown exact force-field/charge provenance.

The separately retained five AmberTools 24.8 cases from Phase 4D2.1 were also
compared, at three seeded perturbed frames each (15 whole-system comparisons):
phenol GAFF/AM1-BCC, phenol GAFF2/AM1-BCC, phenol GAFF2/provided nonuniform charges,
local-template capped PE DP=3 and assigned `F[C@H](Cl)Br`. Reproduce with:

```bash
python scripts/validate_singlepoint_references.py \
  --preparation-manifest /path/to/phase-4d2-1-20260929-final/references.json \
  --output /path/to/new-acceptance.json
```

The script verifies all retained artifact checksums and their agreement with
the signed preparation. Before hashing a JSON-loaded v2 preparation, it restores
integer keys only in `typed_mol2_index_to_site_id`, `prmtop_index_to_site_id`, and
`expected_cip_by_site`; coordinate and atom-name keys remain strings. It then
reimports the original topology with the recorded source identity, force field,
charge method and tolerance, attaches the unchanged shared preparation payload,
and verifies both historical preparation and imported-result signatures. It does
not compute replacement signatures or rewrite historical paths. Incompatible
parser/converter versions fail signature verification; use the original recorded
environment instead of re-signing. A relocated archive is supported through
relative `artifact_subdirectory` entries.

The evaluator is bound through `from_preparation()` on the original system
before evaluating new perturbed frames. The acceptance output records the two
verified historical signatures, actual charge evidence, exact evaluation
coordinates/fingerprints/settings and comparisons of both paths. Offline
synthetic archives test round trips, relocation, changed records, rehashed
contradictory provenance, malformed keys and changed artifact contents. It needs the retained source directory, not a new AmberTools run. Raw
AmberTools data remain locally retained under the prior redistribution policy.
The ordinary test suite uses network-free generated synthetic fixtures and the
pinned external phenol; no live AmberTools installation is needed. The three
validation paths are explicitly separate:

- Offline OpenMM numerics: `tests/test_singlepoint.py`.
- Archived real references: `tests/test_singlepoint_archived_references.py`,
  marker `archived_amber_reference`; opt in using
  `ISLAND_AMBERTOOLS_REFERENCE_MANIFEST=/path/to/references.json`. An unset or
  unavailable archive path produces an explicit skip reason.
- Live AmberTools regeneration: `tests/test_ambertools_integration.py`, marker
  `ambertools_integration`; opt in with `ISLAND_RUN_AMBERTOOLS_INTEGRATION=1`.
  Regenerate into a new directory; preserve historical manifests.

For example, archived acceptance alone is reproducible with:

```bash
ISLAND_AMBERTOOLS_REFERENCE_MANIFEST=/path/to/references.json \
python -m pytest -q -rs tests/test_singlepoint_archived_references.py
```

Both phenol AM1-BCC cases retain their explicitly requested **0.002 e** tolerance.
Their final total charge/residual is **-0.001000000000000112 e**; the verified
typed MOL2 already totals **-0.0010000000000000148 e**. This shows the residual
predates final prmtop import; it does not establish its precise cause or a
particular originating serialization stage. No charges were changed or
neutralized, and default-tolerance acceptance is not claimed. The other three
provided-charge cases have zero residual at their recorded `1e-4 e` tolerance.
These measurements are retained under `charge_report` in every acceptance row.

The validation environment is separate project storage:
`../island-validation/phase4e1-env`, created using the existing ISLAND Python
with `venv --system-site-packages`, then installing OpenMM **8.6.1** into that
venv only. The tested interpreter is Python **3.11.16**, with NumPy **2.4.6**
and ParmEd **4.3.1**. Numerical acceptance uses **Reference / double**. No GPU performance
claim is made. See the retained [acceptance manifest](references/phase_4e1/acceptance.json).

Run `python -m pytest -q -rs`, `python -m ruff check .`, `python -m pip check`,
and `python examples/evaluate_singlepoint.py`. When OpenMM or ParmEd is absent,
only the corresponding numerical test module skips; the lazy-import contract
still runs with OpenMM, RDKit and ParmEd blocked.

This validates software energy/force conversion on the listed cases. It does
not establish production suitability, force-field quality, MD readiness or
long-chain performance. Results retain `production_validated=False` and
`simulation_readiness="not_established"`.


## Verification on 2026-09-29

| Check | Result |
| --- | --- |
| Full suite with OpenMM 8.6.1 | 454 passed, 2 skipped (opt-in AmberTools execution and archived-reference suite) |
| Full suite without OpenMM, existing ISLAND environment | 425 passed, 3 skipped (OpenMM numerical module and both opt-in suites) |
| Explicit archived-reference suite | 1 passed; 15 independent comparisons, both historical signatures verified for every case |
| Ruff | All checks passed |
| pip check in evaluation environment | No broken requirements |
| All 13 examples, including new single-point and existing actual short-polymer preparation | Passed |

The five-case maximum absolute total-energy error was 2.867e-08 kJ/mol;
the maximum absolute force-component error was 7.795e-09
kJ/(mol*angstrom). These figures measure agreement between independently built
systems with the same OpenMM numerical engine; analytical and finite-difference
tests supply separate checks of formulas, force sign and unit conversion.
The full numerical acceptance ran without mocking either calculation path.
A separate negative test injects invalid backend output solely to check rejection.

The archived-reference report additionally records 30 selected finite-difference
checks (three Cartesian components in one frame of each of five cases, at two
step sizes). Maximum absolute force errors were 8.931e-06 at 1e-4 angstrom
and 2.233e-06 at 5e-5 angstrom, in kJ/(mol*angstrom). Both force agreement
and step-size agreement satisfy `rtol=2e-6, atol=2e-5`.

The completion correction also removes the stale README claim that all
GAFF/GAFF2 parameterization is unimplemented. No scientific acceptance or
production-readiness flag changed. There are no remaining unmet Phase 4E1
software acceptance conditions in the tested environment; the documented
physical-model and redistribution limitations remain.
