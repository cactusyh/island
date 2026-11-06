from pathlib import Path
import json
from island.forcefields.assignment_core import normalize
from island.forcefields.bundle_adapters.schema import sha
from scripts.phase_4j76_source_resolution import scope
from scripts.phase_4j84_bundles import spec_from_data
scope()
root=Path('tests/fixtures/phase_4j87'); (root/'graphs').mkdir(parents=True,exist_ok=False)
rows=json.loads(Path('tests/fixtures/phase_4j85/organic-motifs.json').read_text())['motifs']
files=[]
for spec in rows:
    name=spec['id']; relative='graphs/'+name+'.json'; path=root/relative
    raw=(json.dumps({'schema':'island.j87.exact-graph-fixture.v1','id':name,'spec':spec},indent=2)+'\n').encode(); path.write_bytes(raw)
    files.append({'id':name,'path':relative,'sha256':sha(raw),'size':len(raw),'graph_identity':normalize(spec_from_data(spec)).identity})
request={'schema':'island.j87.opls-qualification-request.v1','family':'OPLS-AA','profile_id':'j82_contract_opls_aa','graph_identities':sorted(f['graph_identity'] for f in files),'fixtures':sorted(files,key=lambda f:f['id']),'expected_version':None,'expected_conventions':{},'expected_forms':[]}
p=root/'request.json';p.write_text(json.dumps(request,indent=2)+'\n')
s=Path('scripts/phase_4j87_opls.py');s.write_text(s.read_text().replace('"PENDING"','"'+sha(p.read_bytes())+'"'))
print('request_sha256',sha(p.read_bytes()))
