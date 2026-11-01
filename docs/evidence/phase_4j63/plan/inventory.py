import os,sys,json,platform,importlib.metadata as md
from pathlib import Path
import psi4
from island.pcff_completion import signed,write,pin,digest
from island.pcff_qm_pilot.control_engine import molecule,request
repo=Path.cwd();e=repo/'docs/evidence/phase_4j63';scratch=repo.parent/'island-validation/phase4j63-qm/inventory-scratch';scratch.mkdir(parents=True,exist_ok=True);os.chdir(scratch)
c=json.loads((repo/'docs/evidence/phase_4j62/HMDS-center.json').read_text())['center'];s={**c,'basis':'def2-tzvp','grid':[150,974],'thresholds':[1e-10,1e-8]};m=molecule(psi4,request(s,c['coordinates']));bases={}
for name in ('def2-tzvp','def2-tzvpp','def2-qzvp','def2-universal-jkfit'):
 b=psi4.core.BasisSet.build(m,'BASIS',name,quiet=True);shells=[{'center':b.shell_to_center(i),'am':b.shell(i).am,'pure':b.shell(i).is_pure(),'exponents':[b.shell(i).exp(j) for j in range(b.shell(i).nprimitive)],'coefficients':[b.shell(i).coef(j) for j in range(b.shell(i).nprimitive)]} for i in range(b.nshell())];p=Path(psi4.core.get_datadir())/'basis'/(name+'.gbs');bases[name]={'nbf':b.nbf(),'nshell':b.nshell(),'shells':shells,'shell_identity':digest(shells),'file':str(p),'sha256':pin(p)}
write(e/'engine-environment.json',signed({'engine':'Psi4','version':psi4.__version__,'interpreter':sys.executable,'interpreter_sha256':pin(sys.executable),'library_sha256':pin(psi4.core.__file__),'dependencies':{n:md.version(n) for n in ('numpy','scipy','optking')},'bases':bases,'only_basis_construction_no_QM':True,'other_engines':'not_executed'}));print({k:(v['nbf'],v['nshell']) for k,v in bases.items()})
