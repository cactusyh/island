"""End-to-end capped PE preparation; requires actual AmberTools and ParmEd."""

from island.builders import build_linear_polymer
from island.exceptions import AmberToolsUnavailableError
from island.forcefields import (
    AmberToolsOptions,
    AmberToolsParameterizationEngine,
)

system = build_linear_polymer(
    "[*]CC[*]", dp=3, coordinate_method="local_templates",
    template_seed=2026, assembly_seed=2026,
)
charges = {site_id: 0.0 for site_id in system.topology.sites}
options = AmberToolsOptions(
    force_field="gaff2", charge_method="provided", provided_charges=charges,
    retain_success_artifacts=True,
)
try:
    result = AmberToolsParameterizationEngine().parameterize(system, options)
except AmberToolsUnavailableError as error:
    print("AmberTools integration not run:", error)
    print("No reference output or parameterized snapshot was generated.")
else:
    snapshot = result.to_parameterized_system(system)
    print("preparation signature:", result.record_signature)
    print("imported signature:", result.imported_result.result_signature)
    print("tool versions:", result.record["tool_versions"])
    print("artifact directory:", result.record["artifact_dir"])
    print("production validated:", snapshot.metadata["aggregate"]["production_validated"])
    print("simulation readiness:", snapshot.metadata["aggregate"]["simulation_readiness"])
    print("Zero example charges are software-test inputs, not a physical PE model.")
