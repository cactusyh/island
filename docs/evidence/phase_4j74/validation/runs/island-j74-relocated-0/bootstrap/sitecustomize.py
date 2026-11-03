import os,sys,importlib.util,datetime
try:
 if not hasattr(datetime,'UTC'): datetime.UTC=datetime.timezone.utc
 spec=importlib.util.spec_from_file_location('scripts.phase_4j72_c1_fence','/tmp/island-j74-focused-1-tmp/test_relocation_and_offline_re0/repository/scripts/phase_4j72_c1_fence.py')
 module=importlib.util.module_from_spec(spec)
 sys.modules[spec.name]=module
 spec.loader.exec_module(module)
 module.install()
except BaseException as exc:
 os.write(2,('C1 fence startup failed: '+repr(exc)+'\n').encode())
 os._exit(97)
