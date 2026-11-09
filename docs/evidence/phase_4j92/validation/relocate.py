import json,shutil,subprocess,tempfile,os,sys
from pathlib import Path
from scripts.phase_4j92_intake import ROOT,replay
from scripts.phase_4j73_source_intake import pin
repo=Path.cwd(); dst=Path(tempfile.mkdtemp(prefix='island-j92-relocated-'))
for n in ['src','scripts','docs/evidence/phase_4j92']:
 src=repo/n; target=dst/n; target.parent.mkdir(parents=True,exist_ok=True); shutil.copytree(src,target,dirs_exist_ok=True)
env={**os.environ,'PYTHONPATH':f'.:src'}
run=subprocess.run([sys.executable,'-m','scripts.phase_4j92_intake','replay','--directory','docs/evidence/phase_4j92'],cwd=dst,env=env,capture_output=True,text=True)
print(json.dumps({'destination':str(dst),'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr,'relocated_replay':run.returncode==0,'external_source_asset_root_required':True},sort_keys=True))
raise SystemExit(run.returncode)
