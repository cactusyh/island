# Phase 4H3 — PCFF Class II source records and assignment

## Base and capability boundary

Branch `codex/phase-4h3-pcff-class2-parameters` starts at updated main
`a3a480b0a5ba587ade88a36d16c3ad9120178d18`. Its tree is exactly the reviewed
4H2 tree at `fc030d4ebd9bac4c03e8c5d68b009179eb683116`:
`6c76e781faf3c6c7f4158960f6ca05d4205ddef8`. Both prerequisites are present.
No prerequisite merge, main modification, or historical-record rewrite is needed.

This phase resolves **source-native numerical parameter records**, with explicit
coverage, for the unchanged `island_pcff_acyclic_cho_v1` automatic profile.
It does not construct a PCFF `ParameterizedSystem`, energy evaluator, exporter,
or dynamics model. Complete typing and charges do not imply parameter coverage.
Numerical row coverage does not resolve every physical convention.

The PCFF source remains the external local file from LAMMPS revision
`e891a3e10973c1a729e391a0aefaa02fd70f8c0f`, path
`tools/msi2lmp/frc_files/pcff.frc`, SHA256
`e9c9da613d663c9613fba7823839c71cd07b8679363f629d8ad83c1b1349fd4c`.
No source data or upstream implementation is bundled. See [4H1](phase_4h1.md)
for licensing/source notices and [4H2](phase_4h2.md) for chemical restrictions.
Ordinary alcohols retain **oh/ho**. Auto-equivalence remains disabled.

`production_validated=False`, `simulation_readiness="not_established"`.

## APIs and persistence

```python
from island.forcefields.pcff import (
    load_pcff_source, type_pcff_atoms, assign_automatic_pcff_charges,
    assign_pcff_parameters, save_pcff_parameters, load_pcff_parameters,
)
source = load_pcff_source("/local/pinned/pcff.frc")
typing = type_pcff_atoms(system, source)
charges = assign_automatic_pcff_charges(system, typing)
parameters = assign_pcff_parameters(system, typing, charges)
report = parameters.payload
print(report["coverage"], report["physical_model_complete"])
save_pcff_parameters(parameters, "/new/parameters.json")
loaded = load_pcff_parameters("/new/parameters.json", source, system=system)
loaded.validate_integrity(system)
```

`inspect_pcff_class2(source)` adds semantic interpretation without changing
`PCFFSource.inventory`, its identity, or old signed records.
`PCFFClass2Result` uses schema `island_pcff_class2_assignment_v1`. It embeds the
validated automatic charge/typing evidence, graph identity, source identity,
source catalog, conventions, authoritative inventories, ordered assignments,
equilibrium dependencies and family coverage. Catalog entries contain original
text/line/namespace/version/reference, field names, original and expanded values,
native dimensions, conversion factors and normalized values.

The existing strict JSON envelope preserves integer stable-ID keys and rejects
nonfinite data/duplicate keys. Loading requires the compatible local source;
it does not need original coordinates, external tools or original artifact paths.
Validation reconstructs typing/charges, inventories, selection, dependencies,
conversions and completeness. Recomputing an envelope checksum does not allow
contradictory contents. These are consistency checks, not authentication.
Saving is exclusive; an existing file is preserved. Returned payloads are owned
copies. Coordinates, stale angle/dihedral/improper caches, repeat names and
coordinate-generation history do not enter selection. Assigned chemical stereo,
masses and the authoritative graph retain H2 compatibility checks.

Core assignment/loading require the project's core NumPy dependency, but not
RDKit, ParmEd, Foyer, OpenMM, SciPy or AmberTools. SMILES/PSMILES construction
has its existing RDKit dependency. No typing backend, QM or embedding is invoked
by the parameter API or saved-record CLI.

## Functional forms and units

Normalized numerical units are **kJ/mol, angstrom, radians**. `E`, `L`, `A` in
native dimension labels mean kcal/mol, angstrom and radians respectively.
Every energy-bearing coefficient is multiplied by exactly **4.184**; lengths
remain unchanged. Only equilibrium angles/phases multiply by pi/180.
Angular-displacement coefficients are **not** multiplied by degree/radian factors.
No one-half prefactor is introduced. Negative coefficients and explicit zeros
are preserved.

| Source section | Fields in source order | Meaning / dependencies |
|---|---|---|
| quartic_bond | r0, K2, K3, K4 | K2 dr² + K3 dr³ + K4 dr⁴; energy/Å², /Å³, /Å⁴ |
| quartic_angle | theta0, K2, K3, K4 | K2 dtheta² + K3 dtheta³ + K4 dtheta⁴; displacement radians |
| torsion_3 | V1, phase1, V2, phase2, V3, phase3 | Native three-term values retained; formula discrepancy below |
| wilson_out_of_plane | Kchi, chi0 | Three-connected center; source Kchi(chi-chi0)² |
| nonbond(9-6) | r_min, epsilon | epsilon[2(r_min/r)^9−3(r_min/r)^6] |
| bond-bond | Kbb | Adjacent bond-length displacements; both r0 dependencies |
| bond-angle | Kleft, Kright | Left/right dr times dtheta; two r0 and theta0 |
| bond-bond_1_3 | Kbb13 | Outer bonds of proper; two r0 dependencies |
| end_bond-torsion_3 | three left, three right | Outer dr times cos(n phi); two r0 dependencies |
| middle_bond-torsion_3 | three coefficients | Central dr times cos(n phi); central r0 |
| angle-torsion_3 | three left, three right | Adjacent dtheta times cos(n phi); two theta0 |
| angle-angle-torsion_1 | Kaat | Native scalar retained; angular power/form unresolved below |
| angle-angle | Kaa | Product of two angle displacements; two theta0 dependencies |

For bond–angle coefficients the dimension is energy/(Å rad); bond–bond and
bond–bond 1–3 are energy/Å²; bond–torsion is energy/Å; angle–torsion is energy/rad;
angle–angle is energy/rad². All equilibrium dependencies include the referenced
assignment and selected source row identities. A missing/ambiguous equilibrium
parameter prevents the dependent assignment from being complete.

The 9–6 distance is the **minimum-energy distance**, not 12–6 zero-crossing sigma.
The source declares sixth-power mixing:

```
r_ij = ((r_i**6 + r_j**6)/2)**(1/6)
epsilon_ij = 2*sqrt(epsilon_i*epsilon_j)*r_i**3*r_j**3/(r_i**6+r_j**6)
```

No pair energy is evaluated in this phase.

### Audited references and unresolved physical conventions

Primary documentation: [bond class2](https://docs.lammps.org/bond_class2.html),
[angle class2](https://docs.lammps.org/angle_class2.html),
[dihedral class2](https://docs.lammps.org/dihedral_class2.html),
[improper class2](https://docs.lammps.org/improper_class2.html), and
[pair class2](https://docs.lammps.org/pair_class2.html).
Actual implementation inspection is pinned to the **same LAMMPS revision as the
FRC**, not a moving documentation release. Downloaded source/documentation
hashes are in `docs/evidence/phase_4h3.json`.

- FRC `torsion_3` says **1+cos(n phi−phase)**; pinned LAMMPS documents
  **1−cos(n phi−phase)**. `GetParameters.c` copies the six values unchanged.
  This phase retains source values/phases and flags the discrepancy; it does not
  assert a resolved cis/trans energy convention or change coefficients/phases.
- FRC `angle-angle-torsion_1` writes a linear Phi displacement; LAMMPS uses
  cos(Phi). No Phi equilibrium is supplied in the FRC row. The record therefore
  labels its angular power **unresolved**, retains the scalar with the energy
  conversion only, and marks its functional form unresolved against the reference.
- Wilson source comments specify one out-of-plane displacement; LAMMPS averages
  three out-of-plane angles before squaring. This domain has no applicable
  three-connected center, so no Wilson energy convention is exercised.
- The FRC does not establish LJ and Coulomb special-pair weights. Generic LAMMPS
  `special_bonds` defaults are not evidence of this source's intended policy.
  Neither OPLS nor GAFF exclusions/scaling are inherited.
- The source declares an empty `torsion-torsion_1` section. Its empty inventory
  is retained as raw source information; its absence/applicability is unresolved,
  not treated as a zero interaction.
- The pinned converter restricts bond–bond 1–3 lookup to `cp`-containing propers
  and initializes other coefficients to zero. We do **not** adopt that as a
  justified source-native zero for saturated CHO. Missing rows stay missing.

These gates keep `physical_model_complete=False` for every result. There is no
energy-ready escape hatch, even if `parameter_coverage_complete` is true for a
smaller graph. A future evaluator needs explicit resolution of these conventions,
analytical term checks, and independent energy/force comparisons first.

## Selection, orientation and multiplicity

The versioned policy is
`exact_then_ordinary_family_equivalence_highest_version_v1`:

1. Try source labels directly, in the permitted equivalent orientations.
2. Only if no exact candidates exist, use the **specific family's** ordinary
   equivalences: nonbond, bond, angle, torsion or out_of_plane. Cross terms use
   angle/torsion/out_of_plane as in the audited converter.
3. Select the highest Decimal version. Preserve every candidate and identical
   highest-version selection. Conflicting oriented numerical values are
   `ambiguous`; file order cannot settle the conflict.
4. No wildcard, generic auto-equivalence, quadratic/lower-order or approximate
   fallback. Missing rows remain `missing` with exact sites and labels.

`SearchAndFill.c` documents symmetric shorthand. Only bond–angle (one value)
and end-bond–torsion/angle–torsion (three values) are expanded into two identical
blocks. Arbitrary truncated rows are rejected rather than zero-filled.

Bonds/angles/propers derive from authoritative bonds. Reversal exchanges the
left/right coefficient blocks and associated physical dependencies. Symmetric
repeated labels with unequal directional coefficients are ambiguous.

For angle–angle, `(I,J,K,L)` means the pair `(I,J,K)` and `(K,J,L)`: **J is
central and K is the shared arm**. Only I/L exchange preserves that row's
meaning. Every unordered neighbor triple produces three distinct angle pairs.
A tetrahedral C has four neighbor triples and twelve couplings. This follows the
pinned `MakeLists.c`/`find_angleangle_data` role logic; it does not sort an improper
as four interchangeable atoms. Tetrahedral Wilson candidates are explicitly
`not_applicable`; their angle–angle terms remain required.

`assigned` means numerical rows and dependencies were resolved, including source
zeros. `not_applicable` carries a structural reason. `missing` and `ambiguous`
prevent numerical completeness. Functional-form uncertainties are separately
recorded in source entries/conventions and prevent physical-model completeness.

## Acceptance and reproduction

Commands, with the existing environment and external source:

```bash
PY=../island-validation/phase4e2-env/bin/python
FRC=../island-validation/phase4h1-sources/lammps/pcff.frc
$PY scripts/validate_pcff_class2.py --source "$FRC" \
  --output ../island-validation/phase4h3-acceptance-final
$PY examples/pcff_class2_parameters.py --source "$FRC"
$PY -m island.forcefields.pcff.class2_cli --source "$FRC" \
  --charges /path/to/saved-automatic-charges.json --output /new/parameters.json
```

Outputs are exclusive: rerunning requires a new directory/file. The acceptance
CLI publishes its declaration **before** the cases and persists each outcome.
Default exit 0 means inventory/independent-row/integrity gates passed; it does
not claim complete coverage. `--require-complete` additionally requires numerical
parameter completeness and returns nonzero for the present matrix. The saved
record CLI has the same explicit coverage gate.

The independent oracle uses a separate raw-line scanner, explicit family
column/role specification and Decimal kcal conversion. It does not import the
production parser/selection/conversion helpers. It checks all selected records,
versions/orientations and normalized values, as well as independently declared
inventory counts. Tolerances were declared before execution: coefficient
absolute `1e-10`, relative `1e-13`, in each recorded normalized unit. These are
**parameter comparisons**, not energy/force comparisons or scientific accuracy.
No external parameter-assignment binary was executed; inspected upstream source
is additional semantic evidence, not a claimed numerical oracle run.

The first audit is preserved at `../island-validation/phase4h3-acceptance`.
It correctly exited 1 because the hand-declared dimethyl-ether proper count was
12 instead of 6. The independent graph count is two C–O bonds, each contributing
(4−1)(2−1)=3. The final declaration corrects that analytical expectation. Source,
chemical choices, seeds and numerical tolerances are unchanged. Final records
also make field names and the unresolved angle–angle–torsion unit explicit.
No failed historical H1/H2 records were altered.

Final measured results and verification are recorded below and in the evidence
manifest. Unmet source coverage and scientific limitations are retained even
when the software gate passes.

### Final measured matrix

| Case | Sites | Proper torsions / missing BB1–3 | AA couplings | Checked assignments | Assignment + validation (s) | Total audit (s) |
|---|---:|---:|---:|---:|---:|---:|
| butane | 14 | 27 | 48 | 282 | 4.791 | 14.593 |
| neopentane | 17 | 36 | 60 | 363 | 4.525 | 16.074 |
| ethanol | 9 | 12 | 24 | 140 | 4.684 | 13.451 |
| dimethyl_ether | 9 | 6 | 24 | 110 | 3.985 | 12.874 |
| PE_DP3 | 20 | 45 | 72 | 444 | 5.812 | 17.529 |
| PEO_DP3 | 23 | 42 | 72 | 444 | 5.103 | 16.153 |
| PE_DP50 | 302 | 891 | 1200 | 8058 | 35.920 | 133.973 |

All seven final inventory/source-comparison gates passed. Every listed family
other than bond–bond 1–3 resolved all required numerical rows. Wilson candidates
were structurally inapplicable; angle–angle couplings remained present. Every
case has `parameter_coverage_complete=False` and `physical_model_complete=False`.

Maximum absolute normalized coefficient discrepancy across the matrix: **9.094947017729282e-13**
(in each coefficient’s declared unit), below the predeclared tolerances. This
is Decimal-versus-binary conversion roundoff, not a measured energy error.

The 302-site record is 29,578,572 bytes. Assignment time includes reconstruction
and integrity checking; total audit also includes construction, typing, charge
validation, repeated result validation, oracle checks and serialization. These
are single measured executions, not median benchmarks. Repeated validation and
the intentionally detailed catalog/evidence dominate this smoke test; no checks
were disabled and no broad cache was introduced.

### Verification

- Full ordinary suite: **1,180 passed, 10 expected optional skips**.
- **17 synthetic Class II tests**, including asymmetric reversal, shared-arm
  permutations, repeated labels/conflicts, signed/zero terms, missing dependencies,
  arbitrary ID remapping, stale caches, ownership, source/charge mismatch,
  rechecksummed tampering, optional dependency isolation and exclusive CLI gates.
- Ruff and `python -m pip check`: passed.
- New Class II, existing automatic-charge and manual-charge examples: executed,
  exit 0. Saved-automatic-charge CLI: executed, diagnostic persisted, exit 1
  under `--require-complete` as expected for missing coverage.
- All **5 H1 and 42 H2 retained files** still match their committed hashes;
  pinned source bytes remain unchanged.
- Real pinned-source assignment audit: executed, seven software gates passed.
  No external assignment binary, energy/force evaluation, QM, optimization,
  dynamics or scientific acceptance is claimed.

The next bounded phase must settle BB1–3 applicability, torsion/angle–angle–torsion
conventions, torsion–torsion absence and separate LJ/Coulomb special-pair policy
before declaring a complete physical model. Then validate each native energy
form analytically and against a pinned independent implementation. No PCFF
simulation-ready snapshot is exposed here.
