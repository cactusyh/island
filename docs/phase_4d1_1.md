# Phase 4D1.1: Amber import semantics and external-reference validation

This correction keeps the imported-result path separate from SMARTS typing and
leaves Phase 4A–4C result schemas unchanged. The resolved Amber adapter version
is now `2`, and the imported content schema is `island_amber_resolved_v2`:
zero-LJ and improper-center meaning changed, and snapshots now expose canonical
interaction-family keys. Old imported results must be regenerated from their
source prmtop; legacy native assignment signatures are unaffected.

## Three corrections

1. A prmtop may contain `LENNARD_JONES_CCOEF`, which adds an r^-4 12-6-4
   contribution. The old A/B-only check accepted a saved-and-reloaded
   two-component topology while losing that energy. The adapter now validates
   CCOEF length and every value, rejects any nonzero or nonfinite coefficient
   with a flag-specific `UnsupportedAmberFeatureError`, and records a correctly
   sized all-zero CCOEF as inert. ParmEd may reject a malformed section before
   the adapter sees it; that parser error is wrapped as `AmberImportError` with
   the source flag in its message. Nearby energetic extensions are audited:
   alternate LJ arrays, active hydrogen-bond/solvation coefficients,
   polarizability, CMAP, and other detected energy families are rejected.
   Radius/screen fields are metadata, not an activated implicit-solvent energy.
2. The external phenol source contains an `ho` atom with epsilon=0 and
   `Rmin/2=0`. The old adapter rejected it. A no-LJ site is now represented
   canonically as `(epsilon=0, sigma=0)` in kJ/mol and nm. Here sigma=0 is a
   **sentinel convention, not a recovered physical radius**. The adapter audits
   every source A/B coefficient involving its numerical LJ type, including
   unlike pairs; any nonzero coefficient is an unsupported override. A/B-zero
   sites retain explicit assignments, provenance, charges, and Coulomb terms.
   `epsilon>0, sigma=0` remains invalid. Lorentz–Berthelot and geometric mixing
   yield zero epsilon without creating a pairwise interaction.
3. Amber snapshots previously exposed attribute names such as
   `bond_assignments` as family keys. Both snapshot constructors now expose
   `bond`, `angle`, `proper_torsion`, and the canonical `periodic_improper`.
   There are no aliases or duplicate interactions. Improper source tuple order
   and all periodic terms are preserved. The phenol source also revealed that
   its impropers do not all put the graph-central atom in position 3. The
   adapter validates a unique central site and records its 1-based source
   position per improper; one phenol term uses position 2.

## External Amber reference

`tests/fixtures/phenol_parmed_13239c2.prmtop` is copied from ParmEd's pinned
[`test/files/phenol.prmtop`](https://github.com/ParmEd/ParmEd/blob/13239c2516dd371317ffedeecb411864a5ccb365/test/files/phenol.prmtop).
It was retrieved through the pinned GitHub contents API after a raw download
timed out. SHA-256:
`4722fe1f53d89e9576841b74c2ad3a18494a466be8fb0a2995cd8d956d3d0c28`.
The upstream blob SHA is `daf004783938da29cf71d575a552922ab84fc9b9`.
The pinned [ParmEd README license section](https://github.com/ParmEd/ParmEd/blob/13239c2516dd371317ffedeecb411864a5ccb365/README.md#license)
states GNU LGPL terms and attributes the project to Jason Swails. The
[upstream fixture-addition commit](https://github.com/ParmEd/ParmEd/commit/eaf8db8f49044b3aa3749be8e20eeb3ba53ea415)
states that the phenol and biphenyl files came from David Mobley.

The reviewed independent phenol graph is `Oc1ccccc1` with explicit hydrogens:
six aromatic carbons, one hydroxyl oxygen, five aromatic hydrogens, and one
hydroxyl hydrogen. The stable mapping is `source_index i -> site_id 101+7*i`
for `i=0..12`; the checked source order is C1–C6, O1, H1–H6. The graph and
aromatic bond semantics are constructed from this chemical identity, not from
prmtop connectivity. The source has 13 atoms, 13 bonds, 19 angles, 26 distinct
proper paths, six ordered periodic impropers, and 23 effective 1–4 pairs.
All source charges sum to zero within the declared tolerance; the no-LJ
hydroxyl H retains +0.4186 e. Tests verify source exclusions and per-pair scales.

Representative source-side bond, angle, proper, and improper energies are
calculated with Amber coefficients at non-equilibrium geometry and compared
to converted ISLAND records. The improper angle is computed from four ordered,
noncoplanar Cartesian points, including a term whose central atom is in
position 2. LJ is checked from source A/B at a non-minimum separation;
electrostatics use decoded source charges and Amber's Coulomb constant.
These tests establish **conversion agreement**, not scientific suitability.
`production_validated=False` and `simulation_readiness="not_established"` remain.
Coordinates, if present, stay in Å; parameter lengths are nm.

The upstream fixture is an external Amber reference, but its exact GAFF/GAFF2
identity, version, and charge-method provenance are **not confirmed** by the
fixture file, its filename, or its atom-type labels. The same pinned repository
contains `test/files/parm/gaff.dat`, whose header identifies GAFF 1.8 and has
an `ho` zero-LJ record; this is corroborating chemistry, **not proof** of the
phenol file's preparation route. AmberTools was not installed, so the separate
GAFF/GAFF2-provenance acceptance condition remains unmet. A future Phase 4D2
fixture should be generated with recorded AmberTools versions, commands,
input chemistry, charge method, and checksums, then passed through this importer.
No automatic AmberTools runner is added here.
