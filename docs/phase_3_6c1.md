# Phase 3.6C1: stereochemistry and coordinate-provenance corrections

Phase 3.6C1 corrects three consistency gaps without changing the local-template
architecture or introducing whole-chain embedding.

## Final-graph stereochemistry is authoritative

All explicitly assigned final-graph tetrahedral centers are captured as stable
`site_id -> CIP` assignments before template extraction, whether or not tacticity
was requested. The tacticity sequence remains an additional constraint for its
controlled centers.

When a local fragment is created, every real neighbor maps back to a final stable
site. A temporary context cap maps explicitly to the neighboring repeat's stable
site. The fragment atom's tetrahedral tag is parity-adjusted for the permutation
between final-graph and fragment neighbor order. A final-polymer CIP label is not
copied onto the capped fragment because context caps can alter CIP priority.

After assembly, ISLAND creates a separate RDKit molecule with the generated
coordinates, removes stored stereochemical tags, assigns stereochemistry from 3D,
and compares every constrained site by stable ID. Errors report site ID, repeat
index, expected state, and observed state. If no assigned centers exist, metadata
records `status="not_applicable"`; it does not claim successful validation.

## Coordinate provenance follows coordinates

`ConformationResult.apply_to()` validates the complete coordinate mapping before
mutation, then applies both coordinates and normalized provenance. It remains
non-destructive by default. With `copy=False`, all fallible validation and
provenance construction occurs before modifying the target.

Both builder and direct-generator paths now use the same fields:

- `system.metadata["polymer"]["coordinates"]` identifies the actual source;
- `system.metadata["polymer"]["coordinate_generation"]` records method, seeds,
  attempts, rejections, rollbacks, steric diagnostics, and stereochemistry status.

Consequently, a 2D chemical build can receive local-template coordinates through
`apply_to()` and then be accepted by the Phase 3.6B1 random-walk generator. The
original 2D system remains unchanged.

## Explicit-hydrogen validation

The local-template generator itself enforces the explicit-hydrogen contract before
any template embedding. It reconstructs the authoritative chemical graph and
rejects heavy atoms with implicit or atom-property explicit hydrogens that are not
represented as sites. Merely finding one hydrogen somewhere in the system is not
sufficient.

The RDKit adapter now preserves `no_implicit_hydrogens` and
`explicit_hydrogen_count` metadata so bracket-atom semantics survive round trips.
Chemistry that genuinely requires no hydrogen is accepted; temporary hydrogens on
context caps do not become authoritative sites.

## Remaining scope

The Phase 3.6C restrictions remain: one finite linear atomistic chain, explicit
hydrogens when chemically required, single non-aromatic inter-repeat bonds, rigid
local templates, and no force-field optimization. Multicenter or cap environments
that cannot represent required tetrahedral geometry fail explicitly. Whole-chain
ETKDG, force fields, MD, LAMMPS, packing, and crosslinking are not fallbacks.
