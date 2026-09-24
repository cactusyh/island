# Phase 4D1: validated resolved Amber topology import

Phase 4D1 imports **already parameterized** standard fixed-charge Amber `prmtop`
data. It does not run AmberTools, derive GAFF/GAFF2 types, or replace ISLAND's
chemical graph. Install the optional parser with `pip install '.[amber]'`.
`import island.forcefields` remains possible without ParmEd or RDKit; ParmEd is
loaded only when the adapter is called.

```python
from island.forcefields import import_amber_prmtop, ParameterizedSystem

result = import_amber_prmtop(
    system, "existing.prmtop",
    {0: 101, 1: 205, 2: 310},  # every source atom index -> stable site ID
    source="origin and preparation of this topology",
    force_field="unknown", charge_method="unknown",
)
snapshot = ParameterizedSystem.from_amber_import(system, result)
```

The map is exact, bijective, and **zero-based**. Array positions are not stable
site IDs. The adapter verifies atom counts, elements, mapped bond connectivity,
and per-source-atom numerical LJ type indices. It preserves ISLAND bond orders,
formal charges, stereochemistry, repeat-unit provenance, coordinates, and
coordinate-source metadata. A prmtop cannot verify those chemical attributes;
they remain authoritative in `MolecularSystem`. Coordinates are not required.

The result holds resolved per-site LJ, bond, angle, proper torsion, and periodic
improper selections, independently validated charges, a verified nonbonded
policy, source exclusions and effective 1–4 pairs. It does **not** fabricate an
`AtomTypingRuleSet` or SMARTS match history. It uses a separate
`island_amber_resolved_v1` result schema so Phase 4A–4C signatures stay unchanged.
Assignments retain stable IDs, source row references, source numerical LJ type
indices, and source atom-type names. The SHA-256 of the source bytes, a separate
deterministic mapping signature, and the signed mapping content,
resolved numerical content, charge result, policy, parser/adapter versions, and
declared provenance enter the result signature. `validate_integrity(system)`
checks content and graph compatibility; `is_compatible_with(system)` returns
`False` on malformed/stale results. `copy.deepcopy` and `dataclasses.replace`
are supported, and snapshot creation copies the system and assignments.
The checksum records which file was imported; reuse validation does not reread
the source file, so reimport after changing a prmtop.

## Supported source subset and conversion

| Source item | Import behavior |
| --- | --- |
| Standard `AmberParm`, atomistic fixed-charge sites | Supported with exact mapping, element and bond checks |
| Decoded ParmEd charges | Elementary charge, no second Amber raw-charge division; each component and total checked against ISLAND formal charge (default absolute tolerance `1e-4 e`) |
| LJ 12–6 | Source `Rmin/2` Å and epsilon kcal/mol converted to sigma nm and epsilon kJ/mol; numerical type indices and source A/B coefficients checked against Lorentz–Berthelot mixing |
| Harmonic bonds | Amber `k(r-r0)^2` to ISLAND `0.5*k'(r-r0)^2`: `k'=2*4.184*100*k`, Å to nm for `r0` |
| Harmonic angles | Amber `k(theta-theta0)^2` to ISLAND `0.5*k'(theta-theta0)^2`: `k'=2*4.184*k`; angle displacement is radians, equilibrium angle stored in degrees |
| Proper torsions | Each term `k[1+cos(n*phi-phase)]`; kcal/mol to kJ/mol, phase remains degrees; all source terms for one reversed-equivalent graph path form one record, including zero terms |
| Periodic impropers | Same term conversion, but source atom order is preserved and atom 3 is central; no proper reversal normalization |
| Exclusions and 1–4 | Exact source exclusion list and effective, non-`ignore_end` end pairs checked against shortest-bond-path 1–2/1–3/1–4 policy; duplicate effective pairs rejected; SCEE/SCNB map to independent Coulomb/LJ scale factors |

The LJ equation is `4*epsilon*((sigma/r)^12-(sigma/r)^6)`; `sigma` is the
zero-crossing distance, not the energy-minimum distance. `Rmin=2*(Rmin/2)` and
`sigma=Rmin/2^(1/6)`. Source A/B coefficients must equal
`4*epsilon*sigma^12` and `4*epsilon*sigma^6` in source units. Comparison is by
numerical LJ type pair, not textual atom-type label or an atom-by-atom matrix.
If source A/B are both zero, a prmtop cannot recover sigma from them; a zero LJ
record is accepted only when ParmEd supplies a positive, unambiguous `Rmin/2`.
Unrepresentable zero-length or zero bond/angle force constants are rejected.

### Explicitly unsupported

The importer rejects Chamber/Amoeba variants, CMAP, polarizability, virtual or
extra sites, 10–12 hydrogen-bond coefficients, extra/alternate 1–4 LJ tables,
pair-specific overrides, non-Lorentz–Berthelot source coefficients, nonuniform
SCEE/SCNB, duplicated effective 1–4 interactions, missing mandatory graph bond,
angle or proper-path terms, nonrepresentable exclusions, and other detected
nonstandard interaction families. It does not infer missing scientific terms.
Impropers are required only when explicitly present in the source (and must
match authoritative `Topology.impropers` if that list is populated); it does
not create one at every degree-three center. Complex Amber topologies whose
effective nonbonded behavior differs from ISLAND's shortest-path policy are
rejected rather than approximated.

`force_field` and `charge_method` are declarations only. Their default is
`unknown`; a filename never establishes a method or version. Even with complete
source-term coverage, `production_validated=False` and
`simulation_readiness="not_established"`. Import does not prove chemical
suitability, full scientific term coverage, or MD readiness. ISLAND coordinates
remain Å; imported parameter lengths use nm. No coordinates are rescaled.

## Reproducible fixtures and independent checks

`tests/test_amber_import.py` builds tiny **synthetic**, hand-checkable ParmEd
structures at test runtime. It specifies source atom types and LJ radius/epsilon,
harmonic constants, zero/multiple torsion terms, an ordered periodic improper,
and `ignore_end` behavior; `AmberParm.from_structure(...).save(...)` produces the
temporary prmtop. It checks source versus converted bond, angle, and improper
energies away from equilibrium with independently written formulas. Run
`pytest -q tests/test_amber_import.py` to regenerate and validate these fixtures.
`examples/import_amber_topology.py` generates a separate tiny synthetic chain.
Neither fixture is real GAFF/GAFF2 data.

A traceable real GAFF/GAFF2 topology containing an improper was **not obtained**
in this environment. AmberTools was absent; a candidate small upstream ParmEd
fixture was identified, but fetching its raw content timed out, and its GAFF
provenance/improper content could not be verified. Production-reference acceptance,
including independent real-source
proper/improper/LJ energy checks, remains unverified. A future fixture should
record AmberTools version, exact input structure, antechamber/parmchk2/tleap
commands, charge method, file checksum, expected source parameters, and license.
It must be checked against this importer without weakening unsupported-feature
diagnostics.

Primary format and parser references:

- [Amber file formats](https://ambermd.org/FileFormats.php)
- [ParmEd Amber topology documentation](https://parmed.github.io/ParmEd/html/amber.html)
- [ParmEd topology-object API](https://parmed.github.io/ParmEd/html/api/parmed/parmed.topologyobjects.html)

Phase 4D2 can automate AmberTools execution and pass its resulting prmtop plus
an explicit source-index map into this boundary. That automation and any
production force-field claim are outside Phase 4D1.
