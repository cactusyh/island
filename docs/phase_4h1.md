# Phase 4H1 — PCFF source records and explicit-type charges

## Base and capability boundary

Branch: `codex/phase-4h1-pcff-source-charges`.
Base: `98902ebf1f235326632bd7359206efddaa35fdff` (`origin/main`). Its tree
`4d7d3b44b8241ba1a4c1dfba444d84162a1ccc92` exactly matches reviewed 4G4
`c0698050b1fe3c611b738ce04193379a65f9f0e2`; no prerequisite merge was needed.
Historical GAFF/OPLS records and notebook checkpoint directories are unchanged.

Implemented: local source inspection, explicit type binding, native bond-increment
charges for a deliberately bounded domain, and validated data-only persistence.
Not implemented: automatic PCFF typing, complete parameter assignment, a
`ParameterizedSystem`, energy evaluation, export to LAMMPS, or dynamics.
`production_validated=False`, `simulation_readiness="not_established"`.

## Exact source and audit

The selected source is **LAMMPS-distributed pcff.frc**, not a claim of equivalence
to every commercial PCFF release:

- Repository: <https://github.com/lammps/lammps>
- Revision: `e891a3e10973c1a729e391a0aefaa02fd70f8c0f`
- Path: `tools/msi2lmp/frc_files/pcff.frc`
- Complete-byte SHA256: `e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`
- Header: BIOSYM forcefield dialect 1; latest declared file version 4.0,
  8 October 2013. Earlier version declarations remain preserved.
- Parser: `island_frc_v1`; assignment profile: `lammps_pcff_explicit_cho_v1`.

[The pinned source](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/frc_files/pcff.frc)
contains Biosym, Shenghua Shi, Huai Sun, Joerg-R. Hill, Behnam Vessal and LAMMPS
reference blocks (13 in total). These and all original comments are retained in
inspection. The repository's inspected root license is GPL-2.0; the FRC itself
has historical author notices, not a separate unambiguous data-redistribution
grant. The full file is therefore **external**, not bundled with ISLAND. Obtain
it explicitly, inspect the notices, and provide it locally. Import, loading and
assignment never access the network.

The [pinned msi2lmp README](https://github.com/lammps/lammps/blob/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/README)
explicitly says bond increments and auto-equivalence are unsupported. It cannot
serve as a charge oracle. The independent implementation audit used
[CMMRLab/LUNAR revision 67dabeda9e6bd3cc8f968aa0c6a88730886138ef](https://github.com/CMMRLab/LUNAR/tree/67dabeda9e6bd3cc8f968aa0c6a88730886138ef):
`src/all2lmp/ff_functions.py` (charge summation, lines 217–306) and
`src/all2lmp/read_frc.py` (version selection, lines 575–633). File hashes are in
`src/island/forcefields/pcff/pin.json` and [the evidence manifest](evidence/phase_4h1.json).
LUNAR's inspected license is GPL-3.0. ISLAND does not import, bundle or copy its
implementation; its own implementation uses explicit source records and
independently written bond traversal. The primary methodological publication is
[LUNAR, J. Chem. Inf. Model. (2024)](https://doi.org/10.1021/acs.jcim.4c00730).

### Established semantics and deliberate exclusions

- Each atom starts at zero; each authoritative bond contributes the two recorded
  endpoint increments. Reverse lookup swaps endpoint contributions, not signs.
- Lookup: direct supplied type pair, then ordinary **bond-family equivalence**.
  This is the supported lookup in the audited implementation. Nonbond, angle,
  torsion and out-of-plane equivalences remain separate.
- Auto-equivalence has a distinct `bond_increment` column, also preserved. Its
  apparent fallback branch in this pinned LUNAR routine follows an exhaustive
  `elif not types_flag` branch and cannot execute. We do **not** claim runtime
  validation of that fallback. Auto-equivalence-only coverage is blocked and
  reported with the available source evidence. No wildcard inference is added.
- `highest_version_conflict_reject_v1` selects the highest **decimal** record
  version per ordered key and namespace. Identical highest-version candidates
  retain all IDs; conflicting candidates fail rather than taking the last line.
  This deliberately tightens the upstream last-line handling of equal versions.
  Conflicting reversed rows also fail. Identical endpoint labels with unequal
  increments are ambiguous and fail. No endpoint-conservation assumption is used.
- Missing increments are incomplete diagnostics, never zero parameters. Explicit
  zero rows remain successful source selections. Every component and the whole
  molecule must satisfy `abs(sum(q) - Q_formal) <= 1e-12 e`, using `math.fsum`.
  This fixed profile tolerance tests accumulation of the published coefficients;
  it does not claim unrounded physical charges. Nothing is normalized.

## Parser and supported chemical domain

The pin contains 133 atom-type, 134 equivalence, 108 auto-equivalence and 564
bond-increment records. The parser preserves section namespaces, exact labels
(including `h*`, punctuation and case), version/reference tokens, line numbers,
raw rows, all duplicate candidates, and every unsupported section's original
lines. Semantic masses are in dalton and increments in elementary-charge units.
Unsupported energy sections keep their original numeric text, formula comments
and directives; their units are **not** guessed or converted into ISLAND values.
Malformed/truncated/nonfinite supported rows and a missing `#end` are errors.

`load_pcff_source(path)` verifies the exact pinned hash. An explicit
`expected_sha256=` allows structural inspection of another local file; such a
source has no audited assignment profile and cannot assign charges. A parsed
source is not necessarily an assignable source. No custom source was supplied
for this phase.

Charge assignment is bounded to nonperiodic atomistic saturated C/H/O graphs,
neutral individual atoms, explicit H, ordinary single bonds, positive finite
masses, and no isotope/radical/aromatic state. Carbon degree must be 4, oxygen 2,
and hydrogen 1 with a C/O parent. Accepted source types are `c`, `c1`, `c2`, `c3`,
`h`, `hc`, `h*`, `ho`, `o`, `oc`, `oh`. Explicit type element/connection counts,
carbon H multiplicities, hydrogen parents and `oh`/`oc` H environments are checked.
Other parsed types are explicitly outside this initial assignment domain.
Disconnected molecules are allowed only when **each** component passes.

This is a consistency check on **supplied manual typing**, not proof that these
labels are chemically optimal. There is no polymer-name dispatch, graph type
inference, inherited GAFF/OPLS label interpretation, or charge fallback.
Coordinates are neither needed nor read; no embedding occurs. Atom, graph and
system metadata (including assigned stereo) are copied into the bound identity;
coordinates, names and stale angle/torsion caches do not define that identity.

## APIs, ownership and persistence

```python
from island.forcefields.pcff import (
    load_pcff_source, assign_pcff_types, assign_pcff_charges,
    save_pcff_record, load_pcff_record,
)

source = load_pcff_source("/external/pcff.frc")
typing = assign_pcff_types(system, source, stable_id_to_type,
                           provenance="Manual assignment; describe its justification")
result = assign_pcff_charges(system, typing)
result.validate_integrity(system)
if result.complete:
    charges = result.charges  # owned stable-ID -> e mapping, external to AtomSite
    save_pcff_record(result, "charges.json")
else:
    print(result.payload["diagnostics"])
restored = load_pcff_record("charges.json", source, system=system)
```

Schemas: `island_pcff_explicit_typing_v1` and
`island_pcff_bond_increment_charges_v1`, using the existing strict checksummed
JSON envelope and tagged integer-key mappings. Records contain graph/source/type
identities, explicit typing provenance, selected source rows and candidates,
equivalence paths, oriented contributions, component totals and completeness.
A diagnostic record stores **partial_charges** for explanation; `.charges`
rejects it. A complete record uses the same field with `complete=True`.

Frozen records own a JSON string and immutable source bytes. Every payload/charge
access validates and returns a newly decoded owned value. Integrity recomputes
typing and all contributions from the supplied checksum-verified local source;
rechecksummed contradictory values do not pass. Loading requires that source,
and can additionally check the authoritative system. Source bytes, live backend
objects and arbitrary Python objects are not serialized. Signatures establish
internal consistency, **not computational authenticity or correct manual typing**.

Saving validates before publication and refuses overwrites. It uses the existing
atomic exclusive persistence helper. Loading rejects duplicate JSON keys,
nonfinite constants, unsupported schemas and malformed/reconstructed evidence
with `PCFFError`. Parsing, assignment, loading and comparison require only the
core ISLAND/NumPy environment; no RDKit, ParmEd, Foyer, SciPy, OpenMM or AmberTools.
A charge result cannot construct a complete PCFF parameterized snapshot.

## Reproducible commands and actual acceptance

Explicit source acquisition (not done by the library):

```sh
curl -L https://raw.githubusercontent.com/lammps/lammps/e891a3e10973c1a729e391a0aefaa02fd70f8c0f/tools/msi2lmp/frc_files/pcff.frc -o /external/pcff.frc
sha256sum /external/pcff.frc
python -m island.forcefields.pcff /external/pcff.frc
python -m island.forcefields.pcff /external/pcff.frc --full
python examples/pcff_charges.py /external/pcff.frc
python scripts/validate_pcff_charges.py --source /external/pcff.frc --output /new/acceptance
```

The CLI publishes a declaration before assignment and separate immutable records
and a machine-readable report. All three cases must pass for exit 0. Missing
source/dependencies/coverage or numerical gates produce nonzero exit with durable
case failures. Existing output directories are refused.

Executed in `../island-validation/phase4e2-env/bin/python`; no charge-generation
backend was used. Source audit files are in `../island-validation/phase4h1-sources`.
Final results are in `../island-validation/phase4h1-acceptance-final`; an earlier
successful execution remains in `phase4h1-acceptance`. The example's deliberately
incomplete type mapping was also rejected as expected. Source hashes remained
unchanged. The declaration fixed `1e-12 e` charge/comparison tolerance in advance.

The independent reference uses fixed audited FRC lines **493, 505, 519, 810**
(`c-c`, `c-h`, `c-o`, `h*-o`), fixed ordinary bond-equivalence rows, and Decimal
summation. It does not call the production parser, selection or charge helpers.
Types are illustrative manual assignments justified from source comments:
`c3`/`c2` for alkane C, `hc` for C-bound H, `oh`/`ho` for alcohol, `oc` for ether.
No trustworthy upstream typed fixture was established; that gate is unavailable.

| Actual case | Sites / bonds | Maximum per-site error vs independent rows (e) | Molecular sum (e) |
|---|---:|---:|---:|
| Ethane | 8 / 7 | 0 | -1.38778e-17 |
| Ethanol | 9 / 8 | 1.11022e-16 | 4.85723e-17 |
| Dimethyl ether | 9 / 8 | 1.04083e-17 | 0 |

All three passed exact selected-row coverage, per-component and total charge,
reversed bonds/insertion order, noncontiguous IDs, nonmutation and persistence.
Ethanol's C/C/O/O-H charges were -0.159, 0.027, -0.5571 and 0.4241 e; its C-bound
hydrogens were 0.053 e. These are selected-source charges, not synthetic values.
Upstream typing correctness and scientific accuracy are not inferred.

## Class II inventory and next requirements

The source includes quadratic/quartic bonds and angles, one-/three-term torsions,
Wilson out-of-plane, 9–6 nonbonded parameters and sixth-power combination
rules. Cross-term sections are bond–bond, bond–bond 1–3, bond–angle, angle–angle,
end-bond–torsion, middle-bond–torsion, angle–torsion, angle–angle–torsion and
an empty torsion–torsion section. All remain `raw_only`.
The quartic bond/angle and Wilson energy conventions and torsion phases are
preserved as source comments. No conversion into harmonic/LJ 12–6, OPLS mixing,
or OPLS 1–4 policy is attempted. Electrostatic pair exclusions/scaling and complete
Class II evaluation still require their own source audit and reference validation.

A next automatic-typing phase needs independently vetted typed graph fixtures and
explicit environmental priority rules: H count and carbon substitution, C/O
functional groups and ring contexts, charged/aromatic/heteroatom restrictions,
end-group environments, and ambiguity/override behavior. Source parameter
**equivalence is not a typing rule**. Extending beyond the initial labels needs
chemical evidence, not just a new whitelist. Auto-equivalence charge fallback
needs an independently verified implementation/reference before enabling it.
No further GAFF/OPLS study is prerequisite to this work.

## Verification

Synthetic regressions cover independent analytic charges, orientations, IDs,
explicit zero vs missing, version precedence/conflicts, unsupported auto lookup,
component failures hidden by molecular cancellation, invalid types/parents/H
coverage, persistence and rechecksummed tampering, ownership, CLI failure gates,
and assignment/loading with optional scientific imports blocked. Real-source
execution above is separate from these tests. Full-suite counts and final checks
are recorded below after execution.

Final executed checks:

- Full ordinary suite: **1,135 passed, 10 skipped**, 88.09 s. The skips are the
  existing opt-in backend/archive checks; no PCFF regression was skipped.
- Focused synthetic PCFF suite: **27 passed**. These include injected source
  conflicts/malformed records, not actual quantum or force-field calculations.
- `ruff check .`: passed; `python -m pip check`: no broken requirements.
- Inspection CLI and `examples/pcff_charges.py`: passed with the real pinned file.
- Real-source acceptance: **3/3 passed**, retained separate from synthetic tests.
- Environment: Python 3.11.16, NumPy 2.4.6; no optional chemical backend used by
  PCFF parsing, assignment or reconstruction.

Exact validation commands in this checkout:

```sh
../island-validation/phase4e2-env/bin/python -m pytest -q
ruff check .
../island-validation/phase4e2-env/bin/python -m pip check
../island-validation/phase4e2-env/bin/python -m island.forcefields.pcff ../island-validation/phase4h1-sources/lammps/pcff.frc
../island-validation/phase4e2-env/bin/python examples/pcff_charges.py ../island-validation/phase4h1-sources/lammps/pcff.frc
../island-validation/phase4e2-env/bin/python scripts/validate_pcff_charges.py --source ../island-validation/phase4h1-sources/lammps/pcff.frc --output ../island-validation/phase4h1-acceptance-final
```

Use a new output directory when repeating acceptance. All original successful
experiments remain intact. The unmet upstream-typed-fixture and auto-equivalence
coverage gates are not converted into scientific acceptance by passing tests.
