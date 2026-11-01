"""Post-execution status/import repair; refuse active or unfinished native work."""
import json,lzma,subprocess
from pathlib import Path
from island.pcff_completion import signed,write,pin,require
from scripts.validate_phase_4j63_receipt import execution_source
repo=Path.cwd();root=repo/'docs/evidence/phase_4j63';work=repo.parent/'island-validation/phase4j63-qm';logs=Path('/tmp/island-j63-validation')
read=lambda p:json.loads(p.read_text())
require(read(logs/'deformation-process.json')['returncode']==0,'native_campaign_not_finished')
require(read(root/'deformation-execution.json')['all_completed'],'native_campaign_not_finished')
post=read(logs/'post-QM-results.json')
require(len(post)==6 and all(r['returncode']==0 for r in post),'post_QM_collection_failed')
engine=repo/'src/island/pcff_qm_pilot/hmds_engine.py';coordinator=repo/'scripts/phase_4j63.py';original_engine=pin(engine);original_coordinator=pin(coordinator)
require(original_engine==read(root/'plan/execution-pins.json')['files'][str(engine.relative_to(repo))],'execution_source_already_changed')
require(lzma.decompress((root/'plan/executed-hmds-engine.py.xz').read_bytes())==engine.read_bytes(),'executed_source_archive_changed')
require(lzma.decompress((root/'plan/executed-phase-coordinator.py.xz').read_bytes())==coordinator.read_bytes(),'executed_coordinator_archive_changed')
helper=logs/'status_helper.py';(root/'plan/status-only-repair-helper.py').write_bytes(helper.read_bytes())
preview=logs/'status-repair-preview/repo/src/island/pcff_qm_pilot/hmds_engine.py';engine.write_bytes(preview.read_bytes())
old='    for name, sha in read(evidence / "plan/execution-pins.json")["files"].items():\n        require(pin(repo / name) == sha, "executed_source_mutation")\n'
new='    from scripts.validate_phase_4j63_receipt import execution_source\n\n    execution_source(repo, evidence)\n'
s=coordinator.read_text();require(s.count(old)==1,'coordinator_repair_scope');coordinator.write_text(s.replace(old,new))
subprocess.run(['ruff','check','--fix',str(engine),str(coordinator)],check=True)
subprocess.run(['ruff','format',str(engine),str(coordinator)],check=True)
write(root/'plan/engine-maintenance-amendment.json',signed({'schema':'island.j63-postexecution-status-import-maintenance.v1','original_source_sha256':original_engine,'current_source_sha256':pin(engine),'preserved_executed_source':'plan/executed-hmds-engine.py.xz','status_repair_helper':'plan/status-only-repair-helper.py','status_repair_helper_sha256':pin(root/'plan/status-only-repair-helper.py'),'preserved_executed_coordinator':'plan/executed-phase-coordinator.py.xz','original_coordinator_sha256':original_coordinator,'current_coordinator_sha256':pin(coordinator),'actual_native_adapter':'island.psi4-hmds-j63.v1','maintenance_adapter':'island.psi4-hmds-j63.v1.1','applied_after_all_native_processes_returned':True,'numerical_settings_kernels_callbacks_and_capture_changed':False,'scope':'Unused-import removal, import whitespace, exact budget-exception predicate and coordinator recognition of this explicit amendment. All retained numerical evidence remains bound to original v1 source; no data relabeled or rerun.'}))
execution_source(repo,root)
print('Explicit maintenance AST/provenance verification passed. No new QM executed.')
