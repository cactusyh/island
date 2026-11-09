# J92 — audited LAMMPS OPLS bundle intake

J92 adds an offline, strict two-file LigParGen adapter for one explicitly requested
finite PE five-mer. It builds on J90 graph/obligation identities and canonical
containers, J91 external-path and pinned-byte guards, and the inherited C1
execution fence. It changes no historical files or runtime registries.

The requested external directory is
`/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-opls-bundles/pe5`.
It was absent during this implementation. The delivered real-source decision is
`source_bundle_required`: zero source files/rows inspected, zero graphs qualified,
zero accepted coefficients. Provider version, graph identity, interactions and
conventions are **not assessed**. No alternative directory is searched. Synthetic
controls exercise software contracts only; their metadata explicitly disclaims
scientific validation.

## External contract

Keep all six files outside every Git worktree. No network lookup, discovery,
engine command, energy/force evaluation, minimization, conversion, fitting or
simulation is performed. The parser never executes LAMMPS settings. Recognized
`run 0` and energy-print commands are recorded as inert source text.

`manifest.json` has exactly `schema="island.j92.lammps-bundle.v1"`, `bundle_id`,
`family="OPLS-AA"`, and `files`. The latter lists exactly `Lig.lammps.lmp`,
`Lig.lammps.in`, `request.json`, `provider_metadata.json`, and
`independent_audit.json`, each with `path`, SHA-256 `sha256`, and positive byte
`size`. Supply the manifest SHA-256 separately using `--manifest-sha256` or
`ISLAND_OPLS_MANIFEST_SHA256`. This external pin binds all six files without an
impossible manifest self-hash. Missing, modified, extra or symlinked files block
intake. All five source/control byte sequences are checked before parsing any
of them. The inspected-source counter excludes the manifest control itself;
`source_hashes` includes all six. Replay reopens and rehashes every original.

`request.json` has exactly:

- `schema="island.j92.request.v1"`, `target="finite_PE_five_mer"`, `family="OPLS-AA"`;
- `settings`: inherited finite-chain fields `psmiles`, `repeat_count` (exactly 5),
  two explicit `end_groups`, `protonation_state="as_encoded"`, integer
  `charge_state`, `stereochemistry`, `tacticity`, `stereo_seed`, and
  `intended_export_engine`;
- `backends`: a nonempty unique list drawn from `lammps`, `openmm`;
- positive `atom_count`, exact `graph_identity`, and a complete bijective
  `atom_mapping` list of `[source_atom_id, graph_site_id]` pairs.

The repeat graph must be ethylene; caps are never supplied by the adapter.
Construction must yield one connected, explicit-H, radical-free C/H single-bond
molecule, with the declared charge and stereochemistry. The source bond graph,
atom inventory, mapped elemental masses and charge sum must match exactly.
Coordinates are not used to guess bonds, caps or stereochemistry.

`provider_metadata.json` has exactly `schema="island.j92.provider.v1"`,
`provider="LigParGen"`, `version`, `family="OPLS-AA"`, `source_authority`,
`license={"identifier":...,"authorization":"authorized_local_use"}`,
`provenance`, `atom_typing_origin`, `charge_origin`, boolean `synthetic`, and
`improper_policy`. The version must match the data header. The improper policy
requires `center_position` (zero-based position in source four-site rows) and
`omissions`, each with four source `sites` and an evidential `reason`. No absent
improper is assumed unnecessary: each graph obligation needs an explicit record
or an explicit not-applicable declaration bound by the independent audit.

`independent_audit.json` uses `schema="island.j92.audit.v1"`, `status`,
`independent`, `reviewer`, `authority`, `source_binding`, `graph_identity`,
`settings_identity`, `expansion_hash`, sorted `claim_hashes`, and
`improper_policy_identity`. A pending audit has `status="pending"`,
`independent=false`, and all other fields null. An approved audit requires
`status="approved"`, `independent=true`, reviewer/authority provenance, and
all identities matching freshly reconstructed claims. `source_binding` hashes
the sorted `(path,sha256,size)` manifest rows excluding the audit itself;
the separately pinned manifest binds the audit. This avoids circular hashes.
An audit claim is a source-bound reviewer assertion, not an engine test or
independent physical validation performed by island.

## Conventions and completeness

The narrow settings grammar requires explicit `units real`, `atom_style full`,
finite `boundary m m m`, `pair_style lj/cut/coul/cut 40 40`,
`pair_modify mix geometric tail no`, harmonic bond/angle styles, OPLS dihedrals,
CVFF impropers, and the exact data filename. Declarations needed by `read_data`
must precede it. Unknown/dynamic commands, overrides and unsupported styles
block intake. The complete observed command text, line and raw file hash are
reported. On any mismatch the validated convention object is withheld; expected
values are never reported as observed values.

The adapter records kcal/mol, Angstrom and elementary-charge units without any
conversion. It retains harmonic normalization, OPLS four-coefficient alternating
cosine signs and factor 1/2, CVFF order/sign/periodicity, and 12–6 LJ sigma as the
zero crossing and epsilon as positive well depth. Mixing is geometric. Both
cutoffs are explicitly 40 Angstrom, and Coulomb is truncated (no k-space or
long-range correction). No 9–6 coercion or torsion reinterpretation occurs.

`special_bonds lj/coul 0 0 0.5` or separate explicit LJ and Coulomb triples must
provide both 1–4 scales of 0.5 and 1–2/1–3 exclusions. A connected source graph
is traversed by shortest bond distance. Each unordered distance-3 pair is emitted
once as pair14; distance-1/2 pairs become explicit exclusions. The expansion
records the exact command/location, requested graph identity, indexed bond-graph
hash, algorithm version `island.j92.shortest-bond-distance-1-2-3.v1`, and expansion
hash. Its complete pair inventory is checked independently against inherited
graph obligations. No coefficients are generated by this policy expansion.

Data parsing checks all counts, type IDs/references, masses, charges, pair/bond/
angle/dihedral/improper coefficient widths, finite numerical tokens and row
locations. Every atom, bond, angle, torsion, improper, pair14, topology and
exclusion obligation appears in per-graph reports. Failures retain blockers and
unverified required terms; they never zero-fill missing parameters. Coefficients
remain ephemeral and unconverted in a native-real canonical container. Persisted
reports contain only hashes, source pointers, graph/policy and audit metadata.

Complete source/graph/convention evidence with a pending independent audit yields
`provisional_parameterization`. Only a matching approved independent audit and
complete obligations yield `audited_exact_graph`. Invalid or incomplete supplied
bundles yield `parameter_bundle_incomplete`; absent directory/manifest pin yields
`source_bundle_required`. Synthetic controls never increment scientific graph
counts. Accepted coefficients remain zero in this implementation even if a
future real audit passes.

## Replay and validation

Run from the repository with the inherited environment. Each fence output must
be fresh. The default replay reconstructs the delivered missing-source result:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j92-replay-new -- \
  -m scripts.phase_4j92_intake replay
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j92-strict-new -- \
  -m scripts.validate_phase_4j92_receipt
```

For an explicitly supplied future source, set `ISLAND_OPLS_SOURCE_DIR` and
`ISLAND_OPLS_MANIFEST_SHA256`, then use `ingest --directory /tmp/j92-intake-new`
under the same fence. Retain the returned snapshot SHA-256. Replay with
`replay --directory /tmp/j92-intake-new --snapshot-sha256 PIN`, the same manifest
pin and explicitly supplied originals. `--source-dir` may name a relocated
external copy; content and reconstructed reports must match. Metadata-only
replay cannot substitute for missing external originals.

`docs/evidence/phase_4j92/` contains the reports, sealed receipt, historical
SHA-256 inventory, tested implementation identity and compressed verification
logs. The verification inventory runs the inherited full portable J73–J91 and
graph/typing/assignment regression suite plus J92, under the execution fence.
Scientific tests outside that scope are explicitly skipped, not claimed passed.
It also checks inherited J81–J91 receipts/replays, deterministic and relocated
J92 replay, Ruff, compileall and installed dependencies. Focused controls test
external relocation, mutation, conventions, graph identity and per-instance
completeness; relocated missing-source replay alone proves no real coverage.

Runtime registration, automatic assignment and production validation remain
false. Simulation readiness and chemical-family transfer remain
`not_established`. Qualification never transfers to another cap, repeat count,
stereoisomer, network or infinite chain. Only OpenMM and LAMMPS are in J92 scope;
J92 emits no engine artifact and invokes no engine.
