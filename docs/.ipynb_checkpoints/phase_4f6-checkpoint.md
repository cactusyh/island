# Phase 4F6 — explicit generic fragment charges

## Base and preserved work

Branch `codex/phase-4f6-generic-fragment-charges` starts from fetched main
`a7b57da05dce90f485c887f42f27f9fe47a48e33`, containing the reviewed Phase 4F4.1
prerequisites. Phase 4F5 remains preserved separately on
`codex/phase-4f5-pe-charge-transfer`, commit
`36b522f272bfc7e49f3a80f768fb6b1173c75427`, with its frozen models, held-out
artifacts and research conclusions. It is not a dependency of this generic path.
Main, historical signed records and unrelated checkpoint directories were not
modified. No new PE fit or polymer-name dispatch is introduced.

## Supported chemical contract

The constructor uses existing `RepeatUnit` semantics and requires explicit
`[*:1]` head and `[*:2]` tail placeholders, each single-bonded to a heavy atom.
Each is replaced by an ordinary hydrogen, retaining its role, original dummy
index and bonded-parent identity in metadata. No isotopes label caps. Original
repeat heavy indices and native parent-H groups are recorded explicitly.
Both attachments may share one parent; both cap contributions are then handled
independently. Only the small computational fragment is embedded (ETKDGv3, fixed
seed). The original polymer is never rebuilt by charge assignment.

V1 supports connected, nonperiodic, explicit-H, closed-shell linear homopolymers
with hydrogen-terminated ends and unambiguous builder provenance. Rings and
branches inside a repeat are allowed. Support is graph based, not selected by
polymer name or literal species-specific PSMILES. The demonstrated PE, PEO,
polystyrene-type and PMMA-type repeats all use the same implementation.

Assigned tetrahedral/bond stereochemistry, isotopes, radicals, atomic formal
charges, disconnected repeats, non-single attachment chemistry, other end groups,
copolymers and networks are rejected. Potential but **unassigned** stereocenters
are not silently assigned a tacticity. The data contract conservatively supports
H/C/N/O/F/Cl/Br/S/P valence patterns; the selected backend can be narrower (the
existing Amber adapter limits elements and supported aromatic MOL2 types).
A chemically valid provided-charge fragment does not imply Amber support.
The target ceiling is 1,000 explicit sites; AM1-BCC remains limited to 100 sites.
Historical PE/PEO observation and charge-reference schemas remain unchanged.

## Separate generation, conservation, transfer and compatibility

Public APIs are in `island.fragment_charges`:

- `prepare_capped_fragment(psmiles, seed=2026)` -> `CappedFragment`.
- `FragmentChargeBackend.calculate(fragment)` is the small structural contract.
- `ProvidedFragmentChargeBackend(charges, source=...).calculate(fragment)` checks
  exact fragment-site coverage and finite nonboolean charges in elementary-charge
  units. It retains the supplied vector and labels compatibility **unverified**.
- `AM1BCCFragmentChargeBackend(force_field="gaff2", work_root=...,
  amberhome=..., timeout_seconds=600).calculate(fragment)` ->
  `FragmentChargeResult`.
- `create_fragment_template(calculation, conservation_policy=...,
  charge_tolerance=0.002, hydrogen_policy="verified_parent_group_equal_v1")` ->
  `FragmentChargeTemplate`.
- `assign_fragment_charges(template, system, force_field=None,
  charge_tolerance=1e-12)` -> `FragmentChargeAssignment`.
- `save_fragment_record(record, path)` / `load_fragment_record(path,
  expected_cache_key=...)` provide exclusive persistence and explicit cache reuse.

The charge-only Amber adapter uses the existing `prepare_input`, discovery,
process-group/timeout runner, generated-name lineage, full MOL2 chemistry/stereo
validation, fatal-diagnostic checks and SQM convergence/completion evidence.
It runs **antechamber only**; it does not create a prmtop, run an importer,
fabricate an `ImportedAmberResult`, or change the full parameterization engine.
The capped fragment is the complete QM molecule. All artifacts, source coordinates,
commands, settings, package/executable/selected force-field identities and raw
post-BCC charges remain in the charge record and retained directory. Failure
artifacts remain available; failed calculations cannot become templates. No retry,
seed change, charge fallback or automatic normalization is implemented.

The generic mapping, records and conservation do not select a force field.
AM1-BCC generation here requires an explicit **GAFF or GAFF2** selection.
`validate_force_field_compatibility(template, {"family": ..., "data_sha256": ...})`
checks that exact family/data identity. A matching identity documents the backend
route, **not scientific suitability**. If a target has `metadata["force_field"]`,
that identity must be supplied and checked rather than ignored. An unbound mapping
records compatibility as unassessed. Provided charges cannot satisfy this
compatibility check merely by existing. OPLS/PCFF requests are rejected; neither
backend is implemented and AM1-BCC is never an implicit substitute.

### Conservation is an explicit, separate policy

`strict` retains every fragment charge and rejects a molecular residual exceeding
the declared fragment tolerance. `uniform_fragment_l2_v1` explicitly projects the
**complete** fragment vector onto neutral total charge:

```
delta = -fsum(raw_fragment_charges) / N_explicit_fragment_atoms
q_work[i] = q_raw[i] + delta
```

This minimizes the unweighted squared per-atom correction under the charge-sum
constraint. All atoms, including both caps and native H, participate. Raw charges,
offset, corrected charges and residual are retained. There is no last-atom repair.
Projection verifies conservation at 1e-12 e; it cannot bypass missing atoms,
incompatible chemistry or failed QM evidence. A strict fragment accepted with a
small residual may subsequently fail the independently declared target charge
bound as its residual accumulates with DP. It is not silently projected then.

### Connection transfer and H groups

The transfer policy is `removed_cap_to_bonded_parent_v1`:

1. Resolve original heavy sites by repeat/source index and verify their atomic
   identity, mass, full internal bond graph and oriented inter-repeat connections.
2. For each formed connection, add each removed cap charge to its own parent.
3. At real hydrogen-terminated ends, retain the appropriate cap contribution in
   that parent's H group.
4. Add native-H and retained-cap charges; verify the target H count and parent
   membership, then divide that group total equally among those target hydrogens.

This named `verified_parent_group_equal_v1` transformation makes no arbitrary
individual-H correspondence claim. Removed caps never enter the retained H group.
Heavy sites with no H remain valid. Temporary cap IDs are fragment-local and are
never copied into the target topology. The algorithm checks exact final coverage,
finite assignments, `Q_chain = DP * Q_fragment` within 1e-12 e and the separately
declared target neutrality tolerance. It never repairs the complete target total.
No `AtomSite`, topology, coordinate frame or caller metadata is mutated.

## Records, cache and error boundaries

Data-only versioned records use schemas `island_capped_fragment_v1`,
`island_fragment_charges_v1`, `island_fragment_template_v1` and
`island_fragment_assignment_v1`. They use the existing strict tagged JSON/checksum
utilities, preserving integer IDs without executable deserialization. Payload
access is owned; validators recompute graph/mapping identities, charge totals,
policies, corrections and assignments. Rechecksummed contradictions are rejected
with `FragmentChargeError`. AM1 evidence checks bind commands, method/family,
required settings, zero stage exit status, logs/checksums, completion and exact
MOL2 name/element/connectivity/charge coverage. Live creation additionally uses
RDKit's existing full aromatic/bond-semantic check; offline loading checks aromatic
or neutral Kekule pi-valence consistency without importing RDKit.

A cache key covers the entire validated calculation and fragment (including
geometry/construction provenance), backend/version evidence, conservation tolerance
and hydrogen policy. It is not a PSMILES-only or DP-specific key. Reusing a saved
key does not launch tools; a changed calculation/policy cannot masquerade as that
cache entry. Historical artifact paths are provenance only during offline loading;
embedded evidence is checked without accessing original files. These checks
establish internal consistency, not computational authenticity.

Construction and live Amber chemistry validation require RDKit. Offline record
loading, template comparison and mapping require only core ISLAND/NumPy; no
RDKit, ParmEd, OpenMM, SciPy or QM execution. A future RESP backend needs its own
validated calculation evidence and compatibility rule through this same small
backend contract. Neighbor-containing fragments and alternative caps require an
explicit mapping/cap-policy version. There is no plugin registry or implicit
method substitution.

## Methodological reference

The cap-collapse operation is informed by RadonPy's published workflow and its
actual `connect_mols` implementation. We inspected
[RadonPy core/poly.py at 5d148935](https://github.com/RadonPy/RadonPy/blob/5d14893515376a4518e9f1373a1ebc4bb756db14/radonpy/core/poly.py#L127)
and its [BSD-3-Clause license](https://github.com/RadonPy/RadonPy/blob/5d14893515376a4518e9f1373a1ebc4bb756db14/LICENSE).
Source SHA256 is `1cb8ce58358a9a129435b7cf9af4d529941e7f176dad2eeb5aa355231e077cf8`;
license SHA256 is `0f660cb23fa6c593637c177cda5b08942a3e433c6a0a0cbb2860fcedf07c72f8`.
ISLAND's implementation is independent and uses external stable-ID records rather
than mutating charge properties during molecular connection.

[Hayashi et al., npj Computational Materials 8, 222 (2022)](https://doi.org/10.1038/s41524-022-00906-4)
describes RadonPy's RESP/HF/6-31G(d) charge protocol. This phase uses the existing
AM1-BCC backend and does **not** reproduce that RESP protocol or its scientific
validation. Isolated capped-repeat charges may differ substantially from charges
in a polymer environment. Conservation by construction is not transferability.

## Declared and executed acceptance

```
python -m scripts.validate_fragment_charges --declare-only --output /path/to/NEW
python -m scripts.validate_fragment_charges --execute-declared \
  --output /path/to/NEW --amberhome /path/to/ambertools
python examples/generic_fragment_charges.py /path/to/NEW/ps/template.json \
  --cache-key EXPECTED_SHA256 --dp 4
```

The declaration was published before calculations in
`island-validation/phase4f6-fragments/`. All four fragments use seed 2026,
GAFF2/AM1-BCC, 600 s/stage, no outer retry, `uniform_fragment_l2_v1`, explicit
parent-H symmetrization, fragment raw-strict comparison at 0.002 e, and target
charge checks at 1e-12 e. Targets are DP1/2/3/50. Actual AmberTools is 24.8,
conda build `cuda_None_nompi_py310h834fefc_101`; selected GAFF2 data is 2.2.20
(March 2021). Individual executable version banners remain unavailable; package,
executable and data hashes are recorded. No complete long chain was embedded in
3D or sent to QM. Graph-only targets use the builder's non-3D path.

| Repeat | Fragment atoms | Raw total (e) | Uniform offset per atom (e) | DP50 sites | DP50 total (e) |
| --- | ---: | ---: | ---: | ---: | ---: |
| `[*:1]CC[*:2]` | 8 | 0.002 | -0.00025 | 302 | -1.374e-15 |
| `[*:1]CCO[*:2]` | 9 | 0.000001 | -1.1111111e-7 | 352 | 2.061e-15 |
| `[*:1]CC(c1ccccc1)[*:2]` | 18 | 0.000999 | -0.0000555 | 802 | -6.245e-15 |
| `[*:1]CC(C)(C(=O)OC)[*:2]` | 17 | 0.001003 | -0.000059 | 752 | -4.184e-15 |

All four actual fragment calculations passed, and all 16 mappings had exact
coverage and met their numerical charge bounds. Each saved template was reused
across all DPs with no additional charge calculation. Source and template records
remain unchanged. These are software/command/conservation results, not independent
scientific accuracy evidence.

### Interoperability gate remains unmet

One PE DP3 GAFF2 **provided-charge** parameterization was attempted, recording the
fragment method, conservation/transfer/H policies and assignment identity as its
charge source. It ran antechamber/parmchk2/tleap without SQM, then failed import
under the separately predeclared **1e-12 e** molecular-charge tolerance. The
retained values explain the failure:

- Assigned total: `-6.938893903907228e-17 e`.
- Typed MOL2 total: `-1.3877787807814457e-17 e`.
- Prmtop total read by ParmEd: `-5.487781479240894e-10 e`.
- Maximum source-to-prmtop per-site difference: `2.743890670231508e-10 e`.

The written prmtop/readback precision exceeds that exceptionally strict total
bound, despite near-neutral assignments and much smaller per-site differences
than the existing serialization allowance. No charges or tolerances were changed,
no retry was made, and no successful preparation was fabricated. The fragment/
mapping gate passed; the requested aggregate gate did not (CLI exit 1). Outcome
fields explicitly separate `complete_fragment_mapping`, `integration_passed` and
`complete`. A later interoperability experiment must predeclare a justified
serialization-aware import tolerance; this run remains failed.

## Verification and remaining work

The ordinary suite passed **984 tests**, with **9 existing opt-in skips**
(80.06 s). The focused suite passed **27 tests**. It covers independent cap
algebra, shared parents, odd/even DP, rings/branches/heteroatoms, H bookkeeping,
ID/insertion-order independence, ownership, cache checks, incompatible graphs and
force fields, malformed vectors, rechecksummed records, optional dependency
isolation, synthetic stage failures and truthful CLI gates. Synthetic runner
injection is not counted as QM acceptance. Ruff and pip check passed, as did the
new cached-template example and the existing projection-inspection example.

[The acceptance evidence](phase_4f6_acceptance.json) records tool identities,
artifact checksums and numerical outcomes. Final read-only revalidation checked
109 new case files, 212 original Phase 4F3 source files, 19 Phase 4F4 outputs,
22 preserved Phase 4F5 outputs and 104 Phase 4F5 held-out source files. Historical
files and signatures remained unchanged. All four retained AM1-BCC templates and
their DP50 assignments also loaded in a process blocking RDKit, ParmEd, OpenMM
and SciPy imports.

The next charge-backend extension should be a separately validated RESP adapter
with explicit electronic-structure settings and force-field compatibility. The
next mapping study should test neighbor-containing fragments against independent
held-out polymer environments. Neither extension, automatic production charge
selection, OPLS/PCFF, packing, dynamics nor long-chain QM is implemented here.
All records retain `production_validated=False` and
`simulation_readiness="not_established"`.
