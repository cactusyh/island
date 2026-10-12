# Physical validation protocol

Physical validation is a future gate. The current J23 Reference fixtures remain
software and numerical checks with `production_validated=False` and
`simulation_readiness="not_established"`.

For each supported representative system, record PSMILES, DP, chain count,
crosslink state, box, family, source hash, charge method, parameter coverage,
protocol, replica seeds and trajectory/chunk hashes. Reject a system when any
declared force-field term is missing or a source identity is mixed.

The first physical gates require actual per-frame box volumes and traceable mass,
energy/temperature/pressure stability, density block convergence, RDF consistency,
Rg/Re and anisotropy convergence, topology and crosslink connectivity preservation,
and agreement across independent replicas. Target thermostat/barostat values are
stored separately from measured observables.

RDF requires an explicit stable-ID selection, orthorhombic or declared triclinic
minimum-image algorithm, frame count and block uncertainty. MSD requires validated
unwrapped coordinates or a declared image reconstruction. Unsupported periodic
geometry returns `BLOCKED`/`NOT_EVALUATED`.

Tg requires a predeclared multi-temperature protocol, equilibrated windows,
independent replicas, hysteresis checks and uncertainty. Mechanical properties
require a deformation and stress convention. Neither claim is released by this
software package.
