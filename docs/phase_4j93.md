# J93 — LigParGen PE5 metadata gate

J93 is an additive precondition gate for the real LigParGen PE five-mer source.
It preserves the J92 receipt, result, and blocker report byte-for-byte. It
rehashes the two explicitly pinned LAMMPS files, inventories the external
directory, and checks for an authoritative input structure plus `request.json`,
`manifest.json`, `provider_metadata.json`, and `independent_audit.json`.

The current external directory contains only `LIG.lammps.lmp`,
`LIG.lammps.in`, and notebook checkpoint copies. No SMILES, MOL, SDF, PDB, or
required metadata file is present. J93 therefore stops before parsing source
coefficients, reconstructing a graph, assigning caps, assigning a formal
charge, or expanding `special_bonds`. The decision remains
`parameter_bundle_incomplete`.

The provider and version (`LigParGen 3.0.1`) and validated LAMMPS conventions
are inherited from the sealed J92 real-source report. This inheritance does not
establish the exact graph. Exact graph identity, terminal caps, formal charge,
pair14 count, and exclusion count remain unresolved (`null`, `null`, `null`,
zero, and zero respectively). J92's graph-level missing-term inventory remains
visible.

If an authoritative input structure and all metadata files become available,
they must be supplied as a complete hash-pinned manifest-bound bundle. The
existing J92 exact-graph adapter then validates graph construction, atom
mapping, every Class-I interaction instance, convention source lines,
deterministic pair14/exclusion expansion, and independent audit identities.
J93 never infers PSMILES, caps, charge, or graph identity from the LAMMPS atom
graph.

Run the current gate under the inherited execution fence:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j93-gate-fence -- \
  -m scripts.phase_4j93_metadata_gate emit \
  --directory /tmp/j93-gate \
  --source-dir /inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-opls-bundles/pe5 \
  --data-sha256 19aa11e34253a6838319cdbf183aacae35ace320563e3e9bf57c8f9b8e4a79f6 \
  --settings-sha256 89203559d4717fa50d1b1681c62db2fe0d2ddc264a7d18257584cf24e4cb3a63
```

Replay and `validate` use the same arguments plus the captured-input SHA-256.
Relocated replay may change only the source directory; both pinned files and
the reconstructed report must remain byte-identical.

All scientific and runtime gates remain closed. J93 executes no LAMMPS command,
MD, QM, energy, force, minimization, conversion, fitting, or physical
validation. It adds no GROMACS API.
