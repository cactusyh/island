# Phase 4J3.1 — bidirectional RDKit isotope/radical identity

The reviewed J3 tip was `b5c908babb4ce56d776f4ab9b80c1bca73e9a861`.
Fetched main was still the J2 merge (`d7e30b0`), so the J4 branch starts from
that reviewed J3 tip and retains the unmerged dependency. Main is unchanged.

Reproduction: `[2H]O[2H] -> to_rdkit -> from_rdkit` lost isotope metadata,
changed reconstructed H masses to natural hydrogen and changed PCFF v3 `dw`
to `hw`. The new regression failed before implementation. The correction calls
`SetIsotope` and `SetNumRadicalElectrons` from the authoritative metadata.
No inference is made from mass, and no caller mass or metadata is modified.

Values must be nonboolean integers, nonnegative and representable in RDKit's
uint16 isotope / uint8 radical fields. Larger setter inputs were experimentally
observed to wrap silently (65536 -> 0, 256 -> 0); these are rejected with
`RDKitConversionError` identifying the site and field. Missing metadata means
zero; nonzero metadata is retained by `from_rdkit`, as introduced in J3.
Absent zero-valued fields remain absent in ordinary historical conversions.
The storage bounds are not a claim that every isotope/radical is chemically
supported; graph sanitization and backend-specific chemistry validation remain.

Before correction: 18 failures and 3 passes in the new tests. After correction:
82 focused RDKit/PCFF tests passed, including D2O, radicals, ordinary molecules,
malformed values and PCFF v3 typing. Real pinned-FRC D2O and ethanol round trips
preserve the native typing identities; methyl-radical identity is retained.
Logs and the real-source receipt are in `../island-validation/phase4j4-declared/`.
No source file, historical model, bundle or workflow record was rewritten.
Full ordinary-suite, lint and environment results are recorded with J4 delivery.

`production_validated=False`; `simulation_readiness="not_established"`.
