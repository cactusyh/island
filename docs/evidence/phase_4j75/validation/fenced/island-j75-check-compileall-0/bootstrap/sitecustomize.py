import os,sys,importlib.util,datetime
try:
 if not hasattr(datetime,'UTC'): datetime.UTC=datetime.timezone.utc
 spec=importlib.util.spec_from_file_location('scripts.phase_4j72_c1_fence','/inspire/ssd/project/sais-suiren-foundation-model-prediction/public/yh/island/scripts/phase_4j72_c1_fence.py')
 module=importlib.util.module_from_spec(spec)
 sys.modules[spec.name]=module
 spec.loader.exec_module(module)
 module.install()
except BaseException as exc:
 os.write(2,('C1 fence startup failed: '+repr(exc)+'\n').encode())
 os._exit(97)
