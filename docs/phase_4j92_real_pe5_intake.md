# J92 real LigParGen PE5 intake

This follow-up runs the J92 two-file adapter against the explicitly supplied
external directory `island-opls-bundles/pe5`. The source contains
`LIG.lammps.lmp` and `LIG.lammps.in`; it does not contain the required bundle
manifest, request, input structure, provider metadata, or independent audit.

The partial intake accepts the two files only with separate SHA-256 pins. It
parses source records and settings to report observed data, while always
returning `parameter_bundle_incomplete`. It cannot construct or qualify an
exact graph, enumerate graph obligations, or expand `special_bonds` without
the missing graph identity. The parser does not execute the inert `run 0`,
energy variable, or print commands found in the settings file.

Use the documented intake command under the inherited execution fence:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j92-real-fence -- \
  -m scripts.phase_4j92_intake ingest \
  --directory /tmp/j92-real-intake \
  --source-dir /inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-opls-bundles/pe5 \
  --data-sha256 19aa11e34253a6838319cdbf183aacae35ace320563e3e9bf57c8f9b8e4a79f6 \
  --settings-sha256 89203559d4717fa50d1b1681c62db2fe0d2ddc264a7d18257584cf24e4cb3a63
```

Replay and strict intake receipt validation use the same external originals
and pins plus the captured-input SHA-256:

```bash
PYTHONPATH=.:src ../island-validation/phase4g1-env/bin/python \
  -m scripts.run_phase_4j72_c1 --output /tmp/j92-real-replay-fence -- \
  -m scripts.phase_4j92_intake replay \
  --directory /tmp/j92-real-intake \
  --snapshot-sha256 SNAPSHOT_SHA256 \
  --source-dir /inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island-opls-bundles/pe5 \
  --data-sha256 19aa11e34253a6838319cdbf183aacae35ace320563e3e9bf57c8f9b8e4a79f6 \
  --settings-sha256 89203559d4717fa50d1b1681c62db2fe0d2ddc264a7d18257584cf24e4cb3a63
```

Add `-m scripts.validate_phase_4j92_receipt --intake-only` in place of the
intake module for strict receipt validation. Relocated replay may change only
`--source-dir`; all source bytes and pins must remain identical.
