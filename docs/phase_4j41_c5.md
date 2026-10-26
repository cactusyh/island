# Phase 4J41-C5 — external PDMS Class-II bundle boundary

No authorized external PDMS parameter bundle was supplied or identified in the
inspected local assets. C5 adds an immutable import/schema validation path and
executes the unavailable-asset and incomplete-record controls. It publishes no
PDMS execution profile, coefficients, charges, or backend adapter. The reserved
schema profile name is `island_pcff_external_pdms_class2_v1`; its presence in an
input template does not register an executable profile.

`production_validated=False`; `simulation_readiness="not_established"`.

The branch is `codex/phase-4j41-c5-external-pdms-parameter-bundle`. After
`git fetch origin`, the verified base was
`cc3794db61f84fa14bedfeca5388945526d43aab`, merging C4 in PR #67. C4's
`734195bb5cb9661be87da9ac929d76e2e8103079` is an ancestor. Main is unchanged.
The tested implementation is
`21d51e5f4b2fa7e2f049bdcbd99381e3551bc84d`.

## Public input and immutable record

The public input remains PSMILES or an explicit final graph. No user CAR/MDF is
required. Internally generated CAR/MDF is used solely for the independent
converter oracle. `island.forcefields.external_pdms` supplies:

- `requirements(system)`: exact graph, stable sorted atom IDs, component and
  molecule identities, all applicable Class-II interaction keys, and explicit
  equilibrium dependencies.
- `bundle_template(system)`: an intentionally incomplete input template with
  no assigned types, charge vector, or numerical coefficients.
- `ExternalParameterBundle`: a frozen, JSON-owned checksummed record. Accessing
  its payload returns a separate decoded object.
- `load_external_pdms_bundle(path, system, expected_sha256=..., assets=...)`:
  loads a user-supplied path and compares its identity with an out-of-band pin.
  Missing paths return `asset_unavailable`; invalid or incomplete records return
  `parameter_incomplete`, with all detected stages and missing interactions.

The historical generic final-graph validator does not support silicon. C5 uses
a distinct `island_external_pdms_final_graph_v1` contract, without changing any
historical graph identity. It validates neutral, saturated explicit-H C/H/O/Si
graphs, Si atomic number and four-coordinate valence, C/O/H valence, single
bonds, and the presence of a Si–O bond. Charged, radical, isotope, aromatic,
phosphorus and unrelated domains are rejected. The graph includes masses,
chemical metadata, component and molecule membership, and polymer/end-group
provenance. Coordinates remain in the separate system record.

The bundle schema is `island_external_pdms_class2_bundle_v1`. It contains:

| Group | Required records |
|---|---|
| Identity | Exact graph and graph hash, atom order, schema/profile versions |
| Atoms | Explicit type per atom, element, mass, source-row evidence |
| Sources | One `external:` source family, primary row-table identity, all file hashes, provenance and evidence references |
| Executables | Hash-bound files plus version, build, revision and role |
| Charges | Exact decimal vector in atom order, explicit origin, evidence for each atom; complete oriented increments when requested |
| Nonbonded | Per-atom 9–6 records, sixth-power mixing and explicit 1–4 LJ/Coulomb scaling `0 0 1`, evidenced policy row |
| Bonded | Quartic bonds/angles, proper torsions, applicable Wilson impropers |
| Couplings | Bond–bond, bond–angle, bond–bond 1–3, angle–angle, angle–torsion, end/middle-bond–torsion and angle–angle–torsion |
| Dependencies | Every required equilibrium bond/angle reference, including exact ordered sites |
| Boundary | Fixed units/conversion policy, no charge transformations, no execution authority or scientific-readiness promotion |

The external supplier's primary JSON document has schema
`island_external_class2_source_rows_v1`, `source_family`, `provenance`, and a
`rows` mapping. Each row has `family`, exact ordered `types`, decimal-string
`values`, `units`, `provenance`, and `citations`. Each citation has a relative
`file`, one-based `line`, exact `raw` text, and zero-based numeric `columns`.
All cited files must be hash-checked assets. Reconstruction checks the literal
source values, including zeros. A primary row document that merely contains a
zero without independently reconstructible raw evidence is rejected. No
wildcard resolution, interpolation, native-source fallback, or zero filling is
performed. The original source unit vectors are those in `FAMILIES`; they use
kcal/mol, angstroms, degree equilibrium angles and radian angular displacements.
No numerical conversion occurs during schema validation. A future qualified
adapter must implement the recorded conversion policy explicitly.

Charge vectors are never normalized or redistributed. Exact decimal charge
totals must be zero for this neutral domain, including each connected component.
The accepted data-origin labels are `external_explicit_vector` and
`external_source_bond_increments`. Neither asserts native PCFF, published
COMPASS, OPLS-AA, GAFF/GAFF2, or RESP/Psi4 provenance. Missing increments block
increment-derived charges. PCFF and COMPASS source families cannot be combined
inside a profile. A changed bundle cannot authorize itself by updating its own
checksum: the caller supplies the expected bundle hash separately.

Checksums prove consistency, not authorship or physical validity. A hypothetical
complete import would receive only `parameter_bundle_validated`, with
`execution_authorized=False`, `execution_profile=None`, and a
`scientifically_unqualified` backend. This tranche has no complete authorized
external bundle with which to qualify an execution adapter. Neither the import
record nor an embedded receipt can grant that authority.

## Executed graph and source matrix

All five chemically validated final graphs were constructed and retained:

| Fixture | Atoms | External result | Original pinned-PCFF result |
|---|---:|---|---|
| Short Si–O–Si, `[SiH3]O[SiH3]` | 9 | `asset_unavailable` | `source_incomplete`; 54 unresolved coupling requests |
| Hexamethyldisiloxane | 27 | `asset_unavailable` | `source_incomplete`; 378 unresolved coupling requests |
| Explicit-H capped PDMS repeat, DP1 | 16 | `asset_unavailable` | `typing_blocked`; outside native siloxane profile |
| PSMILES PDMS chain, DP2 | 26 | `asset_unavailable` | `typing_blocked`; methoxy/hydroxyl terminal environments unresolved |
| PDMS chain with explicit end groups, DP4 | 46 | `asset_unavailable` | `typing_blocked`; methoxy/hydroxyl terminal environments unresolved |

The polymer fixtures use `[*:1]O[Si](C)(C)[*:2]` with head `[*]C` and tail
`[*]O`. They are finite, hydrogen-completed methoxy/hydroxyl-ended graphs, not
uncapped periodic repeats. The archive retains actual valence, end groups,
atom IDs, masses, coordinates and every interaction/equilibrium requirement.

The original and diagnostic compiled converters both ran on the two typed,
charge-complete native controls, without `-ignore`. Both exit 17:

- Short Si–O–Si: missing end-bond–torsion key `hsi sio osi sio`.
- Hexamethyldisiloxane: missing end-bond–torsion key `hc c3 sio c3`.

The independent raw-source reconstruction checks all 74 observed converter
lookup calls. Original and diagnostic stdout/stderr agree after removing trace
lines. The three typing-blocked graphs are not passed to the converter, and
are not misreported as unavailable executable assets. No original failure or
converter initialization becomes a coefficient. There are no new PDMS energy,
force, minimization, dynamics, checkpoint or RNG results because no complete
external model exists. Retained finite-silane evidence remains unchanged.

The native source hashes are:

| Asset | SHA-256 |
|---|---|
| `pcff.frc` | `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c` |
| `pcff.rlb` | `82087ba7e9117d3689ee93870b35cd21ddbc2db6b16142cdcf7696b5af7dcc36` |
| `pcff_templates.dat` | `dbde1d51aae2ce3a38d6bcb3e70c2c4ca2116a69904ea8c1b8ede58f06a60093` |
| Original compiled msi2lmp | `59f3718d4e537d584011315b76d0bf80ca105b43468f87b384cfa4ce19d76123` |
| Retained LAMMPS executable | `6659313ef24e7a8477008c98d19c9410cd0a04e47cff0a9d91216f1040bf9191` |

The archive retains the converter's source hashes, clean local source diff,
revision `36f7c81ced8dd456430b13647e48f0c7d7bd8faa`, compiler/version,
`-O2 -lm` build, executable hashes, runtime versions and hardware. The existing
C4 compiled builds are reused, and original converters are reexecuted after
relocation. No external PDMS source or charge hash is asserted when the asset
is absent. No force-field row was downloaded or borrowed from another family.

## Validation and reproduction

The focused C5 suite passed **51** tests. The supported full suite passed
**3,345**, with **10 existing optional skips** and no failures. Its exact
J22–J41-C4 subset contains **1,326 passing tests**. Ruff, compileall, all five
applicable pip checks, and the diff check passed. J41 and C1–C4 receipt/manifest
validators passed. Independent comparison against the fetched Git base
preserved **2,523 historical paths**, including the archive symlink and its
target bytes. Fresh-process relocation preserved **83 study files** exactly.

The receipt and archives contain exact commands, JUnit results, source
reconstructions, independently executed converter failures, relocated bundle
validation, and historical byte hashes. Schema controls reconstruct literal
pinned-PCFF rows only to test the row contract; these are never published as
external PDMS parameters. Explicit source zeros and inserted zeros have
separate controls.

Negative controls cover missing torsions/couplings, bond–bond and bond–angle,
increment requests, graph/order/atom identity, source hashes and families,
elements/valence, molecule identity, fabricated zero rows, changed charges or
charge transformations, units, mixing, 1–4 scaling, authority, equilibrium
dependencies, evidence and receipt/manifest/archive tampering. A related
Si–C–O graph separately confirms that its coupling requirements cannot be
dropped. Incomplete templates are validated after relocation in a fresh process;
all missing interactions remain blocked. The study reconstruction independently
reconstructs raw native rows and reexecutes original converter failures.

```bash
git fetch origin
git rev-parse origin/main
git merge-base --is-ancestor 734195bb5cb9661be87da9ac929d76e2e8103079 origin/main

# Extract the unchanged C4 oracle builds into a fresh directory.
mkdir -p /tmp/island-c5-oracle
tar -xJf docs/evidence/phase_4j41_c4/execution.tar.xz -C /tmp/island-c5-oracle

PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python scripts/phase_4j41_c5.py \
  --distribution ../lammps/lammps/tools/msi2lmp/frc_files \
  --oracle-tools /tmp/island-c5-oracle/study/tools \
  --output /tmp/island-c5-study
```

The completed study exits **1**, reporting `asset_unavailable`, native PCFF
acceptance **1**, and external software acceptance **1**. An audit error exits
2 and is not accepted as a completed blocked study. The unchanged C4
full-source command also remains nonzero. Supplying a candidate requires both
`--bundle /path/to/bundle.json` and its independent `--bundle-sha256`; a single
bundle must not be reused to authorize a different graph in the matrix.

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/phase_4j41_c5.py --reconstruct /path/to/relocated/study

PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  scripts/validate_phase_4j41_c5_receipt.py docs/evidence/phase_4j41_c5

PYTHONPATH=.:src OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1 \
  ../island-validation/phase4e2-env/bin/python -m pytest -q tests/test_phase_4j41_c5.py
```

Supported full-suite environment variables, disjoint test-shard commands,
J22–J41-C4 regression derivation, all applicable pip checks, historical receipt
validators, and fresh-process reconstruction commands are retained in
`validation.tar.xz`. Regression results are an exact checked subset of the
supported full-suite execution, not a second execution claimed from metadata.
No J41/C1–C4 evidence or historical identity is rewritten.

Polymer Ewald/PME remains `ewald_unvalidated`, backend pressure remains
`pressure_unavailable`, and cutoff recovery and long NPT remain blocked.
There is no periodic scientific closure, equilibrium, experimental agreement,
or production-readiness claim.
