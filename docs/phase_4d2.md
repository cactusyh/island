# Phase 4D2: auditable AmberTools GAFF/GAFF2 preparation


Current scope update: [Phase 4F1](phase_4f1.md) preserves the default 100-site
limit and adds an explicit up-to-1,000-site provided-charge budget. AM1-BCC
remains limited to 100; the original phase scope below is historical.

`AmberToolsParameterizationEngine` is a small *external preparation* adapter. It
does not implement GAFF SMARTS typing, parse prmtop itself, or rewrite ISLAND's
chemical graph. It runs `antechamber`, `parmchk2`, and `tleap`, verifies atom
lineage and intermediate chemistry, then delegates the resolved topology to
`import_amber_prmtop()` (Phase 4D1.1). Install optional Python dependencies
with `pip install '.[chemistry,amber]'`; install AmberTools separately and set
`AMBERHOME` and `PATH`. Importing the core package does not load RDKit or ParmEd.

```python
from island.builders import build_linear_polymer
from island.forcefields import (
    AmberToolsOptions, AmberToolsParameterizationEngine,
)

system = build_linear_polymer(
    "[*]CC[*]", dp=3, coordinate_method="local_templates",
    template_seed=2026, assembly_seed=2026,
)
charges = {site_id: 0.0 for site_id in system.topology.sites}  # test input only
result = AmberToolsParameterizationEngine().parameterize(
    system,
    AmberToolsOptions("gaff2", "provided", charges,
                      retain_success_artifacts=True),
)
snapshot = result.to_parameterized_system(system)
```

The result is an `AmberToolsPreparationResult` containing a signed preparation
record and the existing signed `ImportedAmberResult`; it is **not** a competing
numerical parameter assignment framework. The imported result also signs the
preparation record in its provenance. `to_parameterized_system()` validates
both results and creates an owned snapshot. Input topology, site IDs, repeat
provenance, stereochemistry, coordinates, and coordinate-source metadata are
unchanged on success and failure. `production_validated=False` and
`simulation_readiness="not_established"` remain.

## Supported scope and preflight

The first backend accepts one connected finite closed-shell molecule of at most
100 explicit sites, including a fully capped short polymer. Supported elements
are H, C, N, O, F, Cl, Br, and S. It validates authoritative atomic numbers,
RDKit valence/sanitization, radical count, even-electron singlet assumption,
integer formal charges, exact explicit-hydrogen sites, finite coordinates,
plausible bond lengths, and nondegenerate tetrahedral geometry. Known 2D
coordinate provenance is rejected. It does not embed or optimize a full polymer;
the short-PE example uses existing local 3D templates. Genuine 3D input is a
precondition, including for AM1-BCC. The charge method may internally change
its working geometry; source and antechamber coordinates are recorded, while
ISLAND coordinates are never replaced. Whole-chain AM1-BCC cost and
convergence at useful polymer sizes are **not established**.

Two charge modes are explicit and mutually exclusive:

- `provided`: exact stable-site mapping validated against total formal charge;
  serialized to `charges.txt`, passed with `antechamber -c rc -cf`, and checked
  at typed MOL2 and final prmtop within `1e-5 e` per site. It must not run SQM.
- `am1bcc`: `antechamber -c bcc` for the supported closed-shell case. A nonempty
`sqm.out` with a final `Calculation Completed` marker and no failure marker is
  required. The actual SQM output checksum and warnings are recorded. Final
  charges must pass the importer's formal-charge checks. No fallback or
  redistribution is performed.

The graph/charge total check uses the caller-selected absolute
`charge_tolerance` (default `1e-4 e`). Per-site provided-charge preservation
uses the separate `1e-5 e` serialization tolerance. An absent charge is never
zero-filled; initial MOL2 formal charges in AM1-BCC mode are *inputs* for that
calculation, not a claimed final charge model.

## Atom lineage and external stages

Sorted stable IDs receive unique four-character element-prefixed generated names
(`C000`, `O001`, `Cl02`, `Br03`, ...)
and a sidecar map. Typed MOL2 names must be unique and exactly preserved,
including after atom-index reordering. The backend checks mapped elements and
all bonds; RDKit aromaticity perception accepts equivalent aromatic/Kekule
forms but not changed nonaromatic bond orders. Explicit assigned R/S centers
must agree with coordinates before external execution, after antechamber, and
in the LEaP restart. Final prmtop atom names map bijectively to stable IDs;
the Phase 4D1 importer checks elements, connectivity, numerical parameters,
charges, exclusions, and source 1–4 interactions. Neither nearest-coordinate
matching nor unverified symmetric-graph isomorphism is used. If names are
renamed, truncated, duplicated, or lost, the job fails rather than guessing.

The explicit family is sent consistently to `antechamber -at gaff|gaff2`,
`parmchk2 -s 1|2` (GAFF|GAFF2), and the matching `leaprc.gaff|gaff2` in `tleap`.
Antechamber also receives explicit singlet multiplicity, full atom/bond typing,
duplicate-name checking, verbose status, and intermediate-file retention.
The selected leaprc must explicitly name the corresponding data file;
its header and both file checksums are recorded. All commands are argument-list
processes with isolated unique job directories, captured logs, timeouts, and
process-group termination on timeout. Empty or missing MOL2, frcmod, prmtop,
restart, and LEaP log artifacts fail. `ATTN, need revision`, missing-parameter
messages, or fatal tool diagnostics fail; legitimate numerical zero terms do
not. Parmchk analogy/estimate lines are preserved as provenance, not mislabeled
as unresolved. Failures raise `AmberToolsStageError` with stage, log details,
and a retained job directory. Successful artifacts are removed by default;
set `retain_success_artifacts=True` to keep them. A supplied `work_root`
controls the parent directory; filenames remain isolated per job.

All three executables must resolve under the selected `AMBERHOME/bin`; mixing
tools and data from different installations is rejected. Executable checksums
are retained alongside best-effort version probes.

The version-2 signed preparation record includes input chemistry/coordinate digest,
all three atom-index maps, commands, settings, executable version probes (or
explicit unavailability), source data/leaprc headers and checksums, MOL2,
frcmod, prmtop and restart checksums, charge outcome, warnings, and imported
result signature. The outer preparation record must exactly match the copy
inside the signed imported-result provenance, except for the outer
`imported_result_signature` field (excluded from the inner copy to avoid a
circular hash). Both records undergo structural and semantic validation before
snapshot creation, including consistent fixed per-site charge tolerances,
method-specific SQM outcome fields and unique command options. Malformed
reconstructed fields raise `InvalidAmberImportResultError` from either public
result method. Valid v2 records remain compatible. This is not a proof of scientific parameter quality. In
particular, successful execution does not establish that every analogy is
appropriate, that a charge model is transferable, or that an MD simulation
is ready.

## Real references and verification boundary

Run `python scripts/generate_ambertools_references.py --output PATH` with a
working AmberTools installation. The script refuses to overwrite an existing
manifest and runs separately named phenol GAFF/GAFF2 AM1-BCC cases, a phenol
GAFF2 provided-charge case, and a capped PE DP=3 local-template GAFF2 case, and an assigned Cl/Br stereocenter
case.
It stores the exact source prmtops, logs, lineage files, tool/data versions,
commands, checksums, and result signatures in unique job directories plus
`references.json`. The manifest also contains independent source-versus-
converted energies for representative non-equilibrium bonds, angles, proper
and ordered improper torsions, LJ and nonexcluded Coulomb pairs where present.
Zero-LJ site counts are reported from each actual selected library, not
assumed in advance. Phenol cases require periodic impropers. These outputs
must be independently reviewed before being checked in as reference fixtures.
The previously pinned ParmEd phenol file retains its **unknown** exact GAFF
version and charge-method provenance; this script never relabels it.

Normal `pytest` uses network-free unit/mock tests. The separate
`tests/test_ambertools_integration.py` invokes *real* tools only when
`ISLAND_RUN_AMBERTOOLS_INTEGRATION=1` and the executables are available. When
AmberTools is absent, it skips explicitly; mocked tests are not counted as
real integration acceptance. Phase 4D2.1 subsequently exercised all five
real-tool cases with AmberTools 24.8; see [the correction and reference report](phase_4d2_1.md)
for exact provenance, checksums, durable artifact location, the committed
five-case manifest, follow-up verification results and remaining scientific limits.

This phase does not add RESP, repeat-unit charge transfer, OPLS, PCFF/Class II,
MLIP, periodic packing, crosslinking, general minimization, MD, or export.

Published Amber antechamber option semantics: [Amber antechamber manual](https://ambermd.org/antechamber/ac.html).
The numeric `parmchk2 -s 1|2` mapping is also documented in
[AmberClassic's parmchk2 source](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/parmchk2.c);
its [GAFF](https://github.com/Amber-MD/AmberClassic/blob/main/dat/leap/cmd/leaprc.gaff)
and [GAFF2](https://github.com/Amber-MD/AmberClassic/blob/main/dat/leap/cmd/leaprc.gaff2)
leaprc files explicitly load `gaff.dat` and `gaff2.dat`, respectively.
Amber topology semantics: [Amber file formats](https://ambermd.org/FileFormats.php)
and [ParmEd Amber API](https://parmed.github.io/ParmEd/html/amber.html).
