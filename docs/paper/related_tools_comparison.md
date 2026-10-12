# Related tools comparison

This comparison uses verified repository capabilities and the supplied RadonPy
paper as a methodological reference. It does not claim universal superiority.

| Capability | ISLAND | RadonPy | msi2lmp | Foyer | AmberTools |
|---|---|---|---|---|---|
| Input representation | source-bound PSMILES/final graph | attachment SMILES workflow | molecular files | typed molecular topology | molecular files/Amber inputs |
| Polymer construction | linear/copolymer/stereo builders | automated chain growth | not evaluated as polymer builder | not evaluated | not evaluated |
| Crosslinking | deterministic graph plan before typing | workflow-specific polymerization | not evaluated | not evaluated | not evaluated |
| Packing | periodic multichain plan with stable remap | amorphous packing/compression | export conversion | not evaluated | not evaluated |
| Force-field assignment | unified PCFF/OPLS-AA/GAFF/GAFF2 records | GAFF2/LAMMPS workflow | LAMMPS coefficient conversion | SMARTS typing | GAFF/GAFF2 parameterization |
| Periodic MD | OpenMM Reference J23; LAMMPS export | LAMMPS workflow | LAMMPS input/data conversion | not evaluated | not evaluated |
| Analysis | provenance-bound observables and records | property post-processing | not evaluated | not evaluated | not evaluated |
| Checkpoint/resume | transactional, relocatable, separate-process | not evaluated here | not evaluated | not evaluated | not evaluated |
| Provenance | graph/source/charge/checkpoint hashes | workflow and source metadata | file conversion provenance | source XML metadata | artifact lineage |
| Validation strategy | contract, numerical, physical gates | convergence and experimental comparison in paper | conversion tests | typing tests | tool execution checks |

RadonPy is not used as a source of text, code, figures, data or conclusions.
Features marked `not evaluated` are intentionally not promoted to claims.
