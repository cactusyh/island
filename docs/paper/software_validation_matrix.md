# Software validation matrix

The machine-readable matrix separates contract, numerical, independent comparison
and physical evidence. Current `SOFTWARE_VERIFIED` entries describe tested code
paths; they do not establish production MD.

| Group | Evidence | Status | Boundary |
|---|---|---|---|
| Graph correctness | J18-J21 tests and identity/tamper records | SOFTWARE_VERIFIED | chemistry domains are bounded |
| Typing and assignment | J19/J22 tests, source and charge records | SOFTWARE_VERIFIED | no universal force-field coverage |
| Numerical backends | J23 Reference checks and receipts | NUMERICALLY_VERIFIED | LAMMPS KSPACE runtime blocked |
| Workflow reproducibility | J23 split/resume/relocation tests | SOFTWARE_VERIFIED | short Reference workflows |
| Analysis correctness | mapping, units, observable and store tests | SOFTWARE_VERIFIED | topology-dependent observables remain blocked |
| Physical validation | no independent converged cohort | BLOCKED | production readiness remains unestablished |

Synthetic and analytic controls are reported as software or numerical evidence.
No control is called polymer-property validation.
