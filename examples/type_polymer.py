"""Type a polymer graph with the illustrative Phase 4A SMARTS rules."""

from island.builders import build_linear_polymer
from island.forcefields import (
    RDKitSmartsAtomTypingEngine,
    island_demo_v1_ruleset,
)

system = build_linear_polymer("[*:1]CC[*:2]", dp=5, generate_3d=False)
result = RDKitSmartsAtomTypingEngine().type_system(
    system,
    island_demo_v1_ruleset(),
)

for site_id, assignment in sorted(result.assignments.items()):
    print(
        f"site={site_id:3d} type={assignment.atom_type:20s} "
        f"selected={assignment.selected_rule_ids} "
        f"matched={assignment.matched_rule_ids}"
    )
print(f"complete={result.complete}")
print(f"typing_signature={result.typing_signature}")
