import json
from pathlib import Path

from scripts.validate_phase_4j40_c10 import main


def test_actual_polymer_record_is_numeric(tmp_path,monkeypatch):
 src=Path('/tmp/j40-c2-validation/run/reference/gaff2-pe/cross-engine.json')
 if not src.exists(): return
 out=tmp_path/'r.json';monkeypatch.setattr('sys.argv',['x','--cross-engine',str(src),'--output',str(out)])
 try: main()
 except SystemExit: pass
 d=json.loads(out.read_text()); assert d['first_failing_level']=='gaff2-pe_linear';assert d['actual_comparisons'][0]['energy_error_kj_mol']>0
