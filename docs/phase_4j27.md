# Phase 4J27: PBC-correct radius of gyration

J27 corrects the remaining semantic defect in J26 trajectory analysis. For
periodic systems with explicit linear chain groups, radius of gyration is now
computed per chain and frame from mass-weighted unwrapped coordinates along the
bonded path. Chain values are equal-weight averaged per frame before block
statistics. Re, Re/Rg, chain orientation, nematic order and Rg therefore share
one periodic geometry policy.

Periodic analysis without chain groups returns an explicit blocked state. Branch,
crosslink and triclinic unwrapping remain blocked with diagnostics. The J26
TrajectoryStore wording is retained: v2 is transactional persistence with an
opt-in export of a validated completed J23 workflow, not an online writer
attached to every MD step.

No production-scale trajectory campaign, equilibrium physical claim, Tg scan,
mechanical test or experimental comparison is performed.

```text
production_validated=False
simulation_readiness="not_established"
```
