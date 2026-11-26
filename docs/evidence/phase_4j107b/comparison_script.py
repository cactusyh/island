import gzip, json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from hashlib import sha256
from scripts.benchmark_final_graphs import execute

import time
for _ in range(60):
    summary_path = Path('/tmp/island-j107b/full-panel/summary.json')
    if summary_path.exists() and json.loads(summary_path.read_text()).get('complete'):
        break
    time.sleep(2)
else:
    raise TimeoutError('Final paired summary is not complete')

old = Path('docs/evidence/phase_4j107/results.jsonl.gz')
new = Path('/tmp/island-j107b/full-panel/results.jsonl.gz')
def read(p):
    with gzip.open(p, 'rt') as h:
        return {(r['panel'],r['id'],r['field'],r['mode']):r for line in h if (r:=json.loads(line))}
a,b=read(old),read(new)
assert a.keys()==b.keys()
gains=[]; regressions=[]; type_changes=[]
for k,r in b.items():
    before=a[k]
    if before['assigned'] and not r['assigned']: regressions.append(r)
    if not before['assigned'] and r['assigned']: gains.append(r)
    if before['assigned'] and r['assigned'] and before['type_counts']!=r['type_counts']:
        type_changes.append({'key':k,'psmiles':r.get('psmiles'),'before':before['type_counts'],'after':r['type_counts']})
root=Path('/tmp/island-j107b/v1-attribution');root.mkdir(exist_ok=True)
fields=('id','panel','psmiles','smiles','end_groups','field','mode','dp','seed','backend_timeout')
jobs=[{**{k:r[k] for k in fields if k in r},'field_dir':'/tmp/island-j104-work/polymer-field'} for r in gains if r['field']=='opls']
with ThreadPoolExecutor(2) as pool:
    controls=list(pool.map(lambda j:execute(j,root,180),jobs))
from island.forcefields.emc import _field_sources
v1_sources = _field_sources('opls/2024/opls-aa', None, Path('/tmp/island-j104-work/polymer-field'))[3]
report={'v1_control_field_sha256': v1_sources, 'baseline_sha256':sha256(old.read_bytes()).hexdigest(),'current_sha256':sha256(new.read_bytes()).hexdigest(),
        'gains':[{'key':[r[k] for k in ('panel','id','field','mode')],'psmiles':r.get('psmiles')} for r in gains],
        'regressions':regressions,'type_changes':type_changes,'current_code_v1_controls':controls,
        'budget_caveat':'J107 outer budget was 120s. Current outer budget is 180s, graph backend command budget 30s. PCFF numerical/profile sources unchanged; OPLS v2 gains are attributed using current-code v1 controls.'}
Path('/tmp/island-j107b/comparison.json').write_text(json.dumps(report,indent=2)+'\n')
print('gains',len(gains),'regressions',len(regressions),'type changes',len(type_changes),'opls v1 controls',len(controls),flush=True)
