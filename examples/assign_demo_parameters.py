"""Assign synthetic Phase 4B parameters and show explicit limitations."""

from island.builders import build_linear_polymer
from island.chemistry import from_smiles
from island.forcefields import (
    ParameterAssignmentEngine,
    RDKitSmartsAtomTypingEngine,
    island_demo_parameters_v1,
    island_demo_v1_ruleset,
)

ruleset = island_demo_v1_ruleset()
library = island_demo_parameters_v1()
typing_engine = RDKitSmartsAtomTypingEngine()
parameter_engine = ParameterAssignmentEngine()

polymer = build_linear_polymer("[*:1]CC[*:2]", dp=5, generate_3d=False)
typing = typing_engine.type_system(polymer, ruleset)
assigned = parameter_engine.assign(polymer, typing, ruleset, library)

print(library.source)
for family, coverage in assigned.coverage.items():
    print(
        f"{family}: required={coverage.required} assigned={coverage.assigned} "
        f"missing={coverage.missing} ambiguous={coverage.ambiguous}"
    )
first_site = assigned.site_assignments[min(assigned.site_assignments)]
print(
    "example selection:",
    first_site.site_ids,
    first_site.parameter_id,
    first_site.source,
)
print("charges:", assigned.charges_status)
print("simulation readiness:", assigned.simulation_readiness)

# Methanol is typed by island_demo_v1, but the PE-only synthetic parameter
# library intentionally has no oxygen terms. Diagnostic mode reports the gaps.
unsupported = from_smiles("CO", random_seed=2026)
unsupported_typing = typing_engine.type_system(unsupported, ruleset)
partial = parameter_engine.assign(
    unsupported, unsupported_typing, ruleset, library, strict=False
)
print("intentional missing case complete:", partial.complete_supported_scope)
for diagnostic in partial.diagnostics[:5]:
    print(
        "missing diagnostic:",
        diagnostic.family,
        diagnostic.site_ids,
        diagnostic.atom_types,
        diagnostic.reason,
    )

print("This synthetic assignment is not scientifically validated or MD-ready.")
