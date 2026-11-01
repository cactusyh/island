# Read-only progress report; launches no native calculation and changes no scientific record.
import json,time
from pathlib import Path
root=Path.cwd(); work=root.parent/'island-validation/phase4j65-qm'; out=Path('/tmp/island-j65-validation')
a=json.loads((work/'attempts.json').read_text())['calls']; print(time.strftime('%H:%M:%S UTC',time.gmtime()),len(a),'counted attempts;',sum(r['status']=='completed' for r in a.values()),'completed')
for k,r in a.items():
 if r['status']!='completed':
  p=work/'outputs'/k/'psi4.out';print('Active:',k,'elapsed seconds:',round(time.time()-r['started_unix']))
  if p.exists():print('\n'.join(p.read_text(errors='replace').splitlines()[-5:]))
for n in ['primary-numerical-gates.json','trial-fit.json']:
 p=root/'docs/evidence/phase_4j65'/n
 if p.exists():
  r=json.loads(p.read_text());print(n,{k:r[k] for k in ['observations','joints_authorized','fit_executed','validation_authorized','numerical_eligible'] if k in r})
for i in range(2):
 p=out/f'full-{i}.log';s=p.read_text(errors='replace').splitlines();print('Suite shard',i,s[-1][-180:] if s else 'starting')
