# Validation status after Phase 4J28

The J26 integrity corrections are software and numerical validation work. The
focused J26, J25 analysis and trajectory-store tests pass. The supported
environment is required for J17-J25 backend and OpenMM regressions; optional
dependency failures in the base environment are recorded separately.

The corrected PBC implementation supports linear bonded paths in orthorhombic
cells. It reports per-chain Re/Rg and uses equal chain weighting across frame
means. Branched, crosslinked and triclinic unwrapping returns an explicit
blocked diagnostic.

J28 adds graph-bound all-atom selections.  The retained J23 loader now carries
 authoritative final-graph bonds and molecule membership; explicit-H linear
 molecules use bonded placement for all-atom mass-weighted Rg while Re and
 orientation use explicit backbone endpoints.  Selection identities are
 persisted and trusted during record reconstruction.  This is bounded software
 integration and analytical control evidence, not production physical MD.

J29 exercises a retained seven-frame explicit-hydrogen PCFF workflow through
the real CLI and loader. All nine requested observables return `ok`; an
independent first-frame mass-weighted calculation agrees with graph-bound
selection and bonded unwrapping. Analysis-record relocation and separate
process reconstruction preserve identities and values. This is software and
numerical acceptance evidence; production trajectory validation, physical
equilibrium, and experimental agreement remain blocked.

No experimental predictive claim, equilibrium production claim, Tg result,
mechanical result, universal force-field claim or LAMMPS Ewald/PPPM claim is
released. The readiness flags remain `production_validated=False` and
`simulation_readiness="not_established"`.

TrajectoryStore v2 provides transactional persistence and opt-in export. The
J23 integration exports a validated completed workflow; it is not an online
writer attached to every MD step. No production-scale trajectory campaign,
physical equilibrium claim or experimental agreement is supported.
