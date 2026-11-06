# J91 — first real OPLS provider intake and exact-polymer qualification

Base: `7b1a332c83f16a40b28c8f482524eb60404d7ad2`.
Branch: `codex/phase-4j91-real-opls-provider-intake`.

The delivered decision is **source_bundle_required**. No external directory or
control pins were supplied. No provider files were inspected, no graph was
requested or qualified, and no scientific OPLS coverage is established. Existing
untracked literature assets were not searched. The software controls are synthetic
and cannot establish scientific qualification. The optional real integration test
skips explicitly when `ISLAND_OPLS_SOURCE_DIR` is absent.

## Intake and exact scope

`intake91.py` extends the unchanged J90 provider parsers, graph builder, Class-I
inventory, row validator, audit validator, and metadata-preserving serializers.
J81–J90 source/evidence files remain unchanged. This adds no runtime registry hook,
force-field family ranking, provider discovery, downloading, or scientific engine.
The C1 fence and J73 import restrictions must be installed before intake/replay.

Supply one external directory using `ISLAND_OPLS_SOURCE_DIR` or `--source-dir`.
It must be outside every Git worktree, with no symlink paths, and contain:

- `request.json` with schema `island.j91.request.v1`;
- `manifest.json` with the existing `island.opls-provider-bundle.v1` schema;
- `audit.json` with schema `island.j91.audit.v1`;
- exactly the source files in the manifest, including every opaque attachment.

Pin the three control files separately with `ISLAND_OPLS_REQUEST_SHA256`,
`ISLAND_OPLS_MANIFEST_SHA256`, and `ISLAND_OPLS_AUDIT_SHA256` (or the corresponding
`--request-sha256`, `--manifest-sha256`, `--audit-sha256` options). Missing inputs,
authorization, or complete metadata produce `source_bundle_required` before any
parameter/graph/opaque source bytes are opened. Reading pinned control metadata is
separate from source inspection. Unlisted, missing, changed, symlinked, or unsafe
source paths block. All source bytes are hashed before any are parsed. Every
intake, export and replay repeats the hash checks; no source cache is trusted.

The request fields are exactly:

- `schema`, `provider` (`polypargen` or `ligpargen`), `source_authority`;
- `settings`, `graph_identity`, `self_contained_finite_molecule=true`;
- `expected_conventions`, `expected_forms`, `nonbonded_definitions`.

`settings` retains J90's complete explicit finite construction contract: `psmiles`,
positive integer `repeat_count`, two cap SMILES in `end_groups`,
`protonation_state=as_encoded`, integer `charge_state`, `stereochemistry`,
`tacticity`, `stereo_seed`, and `intended_export_engine`. No caps, protonation or
stereochemistry are guessed. Unsupported multiple components and radicals block.
J91 restricts intake to C/H/O finite chains; prefer PE or PEO/PEG. It does not
qualify PDMS, repeat-unit fragments, infinite chains, other caps, other repeat
counts, stereoisomers, crosslinked networks or chemical-family transfer.

The manifest contains one exact ordered graph JSON and at most one self-contained
parameter file: canonical OPLS exchange JSON or a bounded GROMACS ITP/TOP from
PolyParGen/LigParGen. Includes, macros and implicit parameter lookup block. A
separate explicitly sourced improper-omission file is supported by J90; an absent
improper is never inferred to be zero. Provider identity/version, source version,
family, file hashes/sizes/roles, license authorization, provenance and conventions
are mandatory. LigParGen must describe the entire finite molecule. Separate
providers cannot be merged. OPLSAA-DB and other provider families are not selected
by the J91 entry point.

## Convention and audit contract

`expected_conventions` must contain the complete J90 `CONVENTIONS` schema and
match the manifest literally, including separate LJ/Coulomb 1–4 scaling, charge
origin, typing origin, exclusions, long-range/cutoff/tail policy, phases and
normalization. `expected_forms` must map all six families (`bond`, `angle`,
`torsion`, `improper`, `charge`, `nonbonded`) to complete structured definitions
from J90 `FORMS`. Malformed schemas produce structured blockers. Source forms
must match these declarations, including units and atom ordering.

`nonbonded_definitions` must equal the public `DEFINITIONS` contract: sigma is the
zero of the 12–6 potential; epsilon is the positive well depth; the minimum is at
`2^(1/6)*sigma`. Supported intake is kJ/mol, nm, charge e, geometric sigma/epsilon,
explicit 12–6 LJ, graph-distance-1/2 exclusions, ordered center-second impropers,
and explicit torsion forms. Arithmetic mixing and other unsupported conventions
are reported as incompatible; they are never relabeled as geometric. There is no
intake unit conversion, 9–6 coercion, RB/periodic reinterpretation or analog transfer.

The independent audit has exactly these fields:

- `schema=island.j91.audit.v1`, `request_sha256`;
- `source_audit`: J90 independent audit with reviewer, authority basis, independent
  flag, manifest/graph/convention identities, record and atom claim hashes;
- `pair_claim_hashes`, `exclusion_claim_hashes`: hashes of `{sites, source}` records
  with exact JSON pointers or physical file lines generated by `pair_claims`;
- `export_conversions`: empty, or the exact `CONVERSION` declaration for LAMMPS.

The audit binds the request and all original files through the manifest. Hashes
prove byte identity; the user supplies the reviewer's authority and authorization.
An audit is not generated automatically for a real source. Source locations and
claim hashes are retained, while numeric coefficients stay only in memory or in
external exports. No raw tables, coefficient arrays or opaque attachments enter
GitHub. `source_rows_inspected` counts parsed atom/parameter/pair/exclusion rows;
control files are separately identified by their pins.

Every required interaction instance is reported. Source exclusions must cover
every required pair explicitly; J90's graph-derived exclusion completion cannot
qualify a J91 result. Incomplete interactions or conventions block the graph.
Only a complete real source plus exact independent audit can set `qualified_graph`.
Synthetic fixtures always leave it null and scientific graph coverage zero, even
when the software contract returns `audited_exact_graph`.

Result vocabulary is preserved: missing required inputs → `source_bundle_required`;
a source-only bundle → `provenance_only`; failed source/term/convention gates →
`parameter_bundle_incomplete`; complete independently audited contract →
`audited_exact_graph`. `provisional_parameterization` remains a supported export
class, but J91 does not downgrade a missing audit into provisional qualification.

## Export boundaries

Exports require a fresh successful intake and an external destination outside both
Git and the source bundle. GROMACS is an ITP fragment and LAMMPS is a coefficient
fragment, neither a complete runnable system. OpenMM is data-only XML with
`loadableForceField=false`; no ForceField loader or Context is invoked. Embedded
receipts retain hashes, graph/provider identities and non-production flags.

LAMMPS's kJ/mol → kcal/mol, nm → Angstrom and harmonic half-factor transformations
require the exact independently audited `CONVERSION` record, including formulae,
source/target units, version and conversion hash. Without it, export blocks.
Only periodic torsions can use this LAMMPS conversion; RB and OPLS-Fourier
reinterpretation is blocked. GROMACS and data-only XML preserve source literals.
The exported source envelope and native fragment are both checked by J90's
`read_export`, without executing any engine.

## Replay and validation

Use the validation Python environment and unique fence directories on each run:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j91-replay-unique -- \
  -m scripts.phase_4j91_intake replay

PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j91-validate-unique -- \
  -m scripts.validate_phase_4j91_receipt
```

For supplied external originals, set the four environment variables above and run
`-m scripts.phase_4j91_intake ingest --directory /outside/checkout/new-report`.
Record the returned snapshot SHA-256. Replay with the same pins and
`replay --directory /outside/checkout/new-report --snapshot-sha256 SHA`.
Relocation changes only `--source-dir`; relative filenames and bytes stay exact.
Metadata-only replay fails explicitly if the original bundle or pins are missing.
Exports use `export --engine gromacs --destination /outside/checkout/oligomer.itp`
(or `lammps` / `openmm` with the request's intended engine).

The delivered sealed receipt is `docs/evidence/phase_4j91/receipt.json`.
`result.json`, `per-graph-report.json` and `missing-terms-and-blockers.json` hold
reconstructable metadata. Strict validation checks all artifact hashes, replay
semantics, exact tested-commit blobs, terminal JUnit/collection inventory, archived
fence events and every base-tracked historical file. Validation logs include
failed software attempts and explicit scope skips. Relocated external replay is
a software control; relocated code/evidence replay checks the delivered result.

All scientific execution counters, accepted coefficients, runtime registration,
automatic assignment and production validation remain zero/false.
`simulation_readiness` and `chemical_family_transfer` remain `not_established`.
No scientific or production validation is inferred from these software tests.
