import json
import subprocess
from pathlib import Path


def test_numeric_characterization_from_retained_frame(tmp_path):
 src=Path('/tmp/j40-c2-validation/run/reference/gaff2-pe/cross-engine.json');forces=src.parent/'forces.dump'
 if not src.exists(): return
 out=tmp_path/'x.json';subprocess.run([__import__('sys').executable,str(Path('scripts/validate_phase_4j40_c11.py').resolve()),'--comparison',str(src),'--forces',str(forces),'--output',str(out)],check=False)
 d=json.loads(out.read_text());assert d['relative_energy_error']>0
