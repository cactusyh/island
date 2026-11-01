import json,time
from pathlib import Path
root=Path.cwd();work=root.parent/'island-validation/phase4j65-qm';out=Path('/tmp/island-j65-validation');previous=None
def read(p):return json.loads(p.read_text())
while True:
 a=read(work/'attempts.json')['calls'];done=sum(r['status']=='completed' for r in a.values());active=[k for k,r in a.items() if r['status']!='completed'];text=f"{time.strftime('%H:%M:%S UTC',time.gmtime())}: {done} completed / {len(a)} counted attempts (cap32)"
 if active:
  key=active[-1];text+=f"; active {key}, {round(time.time()-a[key]['started_unix'])}s"
 print(text,flush=True)
 p=out/'pipeline-results.json'
 if p.exists():
  rows=read(p);current=rows[-1]['name']
  if current!=previous:print('Completed stage:',current,'exit',rows[-1]['returncode'],flush=True);previous=current
 time.sleep(60)
