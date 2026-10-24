# Phase 4J40-C3 — periodic preflight and pressure diagnostics

This additive phase audits the eleven C2 integration failures and adds a
fail-closed periodic cutoff preflight. It requires cutoff + declared safety
margin + neighbor skin to remain below half the smallest current or predicted
box length. It never changes a box, cutoff, density, force field, or target.

The C2 failures remain preserved. Their exact failure records, last checkpoint,
box lengths, volumes, densities and barostat configuration are audited by
`scripts/validate_phase_4j40_c3.py`. The retained OpenMM Reference histories do
not expose a validated full virial pressure. The target pressure and barostat
settings are therefore recorded separately as `pressure_unavailable`; no target
value is substituted.

A documented virial control is available for a backend that supplies complete
kinetic energy, configurational virial, volume and particle count:
`P=(2K+W)/(3V)`, converted from kJ mol^-1 nm^-3 to bar. The helper marks this
as measured virial pressure and explicitly records that it is independent of the
barostat target. It does not claim the retained C2 runs have such measurements.

No recovery dynamics were claimed: the C2 artifacts remain the actual failed
histories, and scientific acceptance remains blocked by pressure availability,
insufficient sampling, eleven cutoff failures and unresolved cross-engine
agreement. `production_validated=False`; `simulation_readiness="not_established"`.

Validation command:

```bash
PYTHONPATH=.:src ../island-validation/phase4e2-env/bin/python \
  -m scripts.validate_phase_4j40_c3 \
  --run /tmp/j40-c2-validation/run --output NEW/failure-audit.json
```
