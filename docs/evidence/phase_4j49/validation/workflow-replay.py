import json,shutil
from pathlib import Path
from island.charge_references.records import unpack
from island.pdms import restore_pdms_system
from island.forcefields.opls_registry import OPLSAssets,use_opls_assets,prepare_polymer_forcefield,canonical_registry
from island.forcefields.selection.inventory import digest
from island.pdms_lammps import LAMMPSOptions,pin
from scripts.phase_4j49 import MODES
from scripts.phase_4j49_workflow import run_control
from scripts.phase_4j47 import write
source=Path('/tmp/island-j49-validation/audit-sealed')
out=Path('/tmp/island-j49-validation/final-evidence-audit')
shutil.copytree(source,out)
summary=json.loads((out/'summary.json').read_text())
assets=OPLSAssets('island_ffref/pcff/PDMS_benchmark-main.zip','8b0e8b74303b39b1b0650be3a15b8dd201b20027465f95ffc3f2989da123fea4')
engine=LAMMPSOptions('/tmp/island-j49-exact/build/lmp',pin('/tmp/island-j49-exact/build/lmp'))
for dp in (0,2,4,6):
 system=restore_pdms_system(unpack((out/f'dp{dp}/system.json').read_text()))
 with use_opls_assets(assets):p=prepare_polymer_forcefield(system,forcefield='opls')
 control=json.loads((out/f'components-dp{dp}.json').read_text())['exact-erfc-ewald']
 directory=out/f'workflow-dp{dp}'
 shutil.rmtree(directory)
 summary['diagnostic_workflows'][str(dp)]=run_control(system,p,directory,engine=engine,static_control=control,solver=MODES['matched-ewald'])
 print(dp,summary['diagnostic_workflows'][str(dp)]['bounded_diagnostic_workflow_passed'],flush=True)
summary.pop('identity')
write(out/'summary.json',{**summary,'identity':digest(summary)})
write(out/'diagnostic-workflows.json',summary['diagnostic_workflows'])
write(out/'canonical-registry.json',canonical_registry())
