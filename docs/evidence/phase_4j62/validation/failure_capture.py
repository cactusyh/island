"""Retain nonzero subprocess diagnostics before pytest fixture teardown."""
import json,os,subprocess,time,traceback
from pathlib import Path
import pytest
_original=subprocess.run

def _text(v):
 if isinstance(v,bytes):return v.decode('utf8',errors='replace')
 return v if isinstance(v,str) else repr(v)

def retained_run(*args,**kwargs):
 result=None;error=None
 try:
  result=_original(*args,**kwargs);return result
 except (subprocess.CalledProcessError,subprocess.TimeoutExpired) as exc:
  error=exc;raise
 finally:
  r=error or result
  if r is not None and (error or r.returncode):
   directory=Path(os.environ['J62_DIAGNOSTICS'])/f'{os.getpid()}-{time.time_ns()}'
   directory.mkdir(parents=True)
   (directory/'command.json').write_text(json.dumps({'argv':_text(args[0] if args else kwargs.get('args')),'cwd':str(kwargs.get('cwd',Path.cwd())),'returncode':getattr(r,'returncode',None),'error':type(error).__name__ if error else None},indent=2)+'\n')
   (directory/'stdout.log').write_text(_text(getattr(r,'stdout',None)))
   (directory/'stderr.log').write_text(_text(getattr(r,'stderr',None)))

def pytest_configure(config):subprocess.run=retained_run

@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item,call):
 outcome=yield;report=outcome.get_result()
 if report.failed:
  directory=Path(os.environ['J62_DIAGNOSTICS'])/f'pytest-{os.getpid()}-{time.time_ns()}';directory.mkdir(parents=True)
  (directory/'failure.txt').write_text(report.nodeid+'\n'+str(report.longrepr)+'\n'+str(report.sections))
  # --basetemp is persistent; record all artifact paths before teardown. No environment dump.
  base=item.config.option.basetemp
  if base:
   (directory/'native-artifact-paths.json').write_text(json.dumps([str(p) for p in Path(base).rglob('*') if p.is_file()],indent=2)+'\n')
