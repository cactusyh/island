"""Seal truthful J72 source evidence and terminal software coverage.

Never launches QM, fitting, MD, validation or conversion. Earlier failed and
interrupted runs remain separate from selected terminal software results.
"""
import importlib.metadata
from copy import deepcopy
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from xml.etree import ElementTree as ET

repo = Path.cwd()
out = Path('/tmp/island-j72-validation')
root = repo / 'docs/evidence/phase_4j72'
v = root / 'validation'
final = out / 'guarded-final'
results = json.loads((final / 'suite-results.json').read_text())
assert len(results) == 8 and all(r['returncode'] in (0, 1) for r in results), [(r['shard'],r['returncode']) for r in results]
# Preserve all prior scripts/logs, including inaccurate pre-correction plans,
# under an explicit verification history. No PDFs or native fixture directories.
for p in sorted(out.iterdir()):
    if p.is_file():
        shutil.copyfile(p, v / p.name)
for directory in ('terminal-checks','scope-corrected-checks','guarded','guarded-final'):
    dest=v/directory;dest.mkdir(exist_ok=True)
    for p in sorted((out/directory).iterdir()):
        if p.is_file():shutil.copyfile(p,dest/p.name)
base = json.loads((out / 'collection.json').read_text())
new = [n for n in (out / 'j72-collection-final.txt').read_text().splitlines() if n.startswith('tests/') and '::' in n]
assert len(new) == 173 and len(new) == len(set(new))
nodes = base['nodes'] + new
assert len(base['nodes'])==5660 and len(nodes) == len(set(nodes))
collection = {'nodes': nodes, 'base_commit': 'e0c82269164781e0aa12881ebb6b308e21735e84', 'historical_nodes': len(base['nodes']), 'J72_nodes': len(new), 'focused_reruns_add_unique_coverage': False, 'historical_node_shards': json.loads((final/'collection.json').read_text())['node_shards']}
(root / 'collection.json').write_text(json.dumps(collection, indent=2) + '\n')
sys.path[:0] = [str(repo), str(repo / 'src')]
from scripts.validate_phase_4j72_receipt import terminal_sources
selected=terminal_sources(root,nodes)
merged=ET.Element('testsuites');suite=ET.SubElement(merged,'testsuite',name='J72 unique terminal results')
for _,case in selected.values():suite.append(deepcopy(case))
ET.ElementTree(merged).write(v/'terminal-results.xml',encoding='utf-8',xml_declaration=True)
(v/'terminal-selection.json').write_text(json.dumps({n:p for n,(p,_) in selected.items()},indent=2)+'\n')
terminal=['validation/terminal-results.xml']
(v / 'junit-files.json').write_text(json.dumps(terminal, indent=2) + '\n')
actual = []; skips = 0
for name in terminal:
    for c in ET.parse(root / name).iter('testcase'):
        assert c.find('failure') is None and c.find('error') is None, (name,c.attrib)
        skips += c.find('skipped') is not None
        actual.append(c.attrib['classname'].replace('.', '/') + '.py::' + c.attrib['name'])
assert sorted(actual) == sorted(nodes) and len(actual) == len(set(actual))
counts = {'unique_tests': len(nodes), 'passed': len(nodes) - skips, 'skipped': skips, 'failure': 0, 'error': 0, 'historical_base_nodes': len(base['nodes']), 'J72_nodes': len(new), 'focused_overlap_coverage_added': 0, 'scope_compliance_is_not_software_test_success': True}
(v / 'terminal-test-counts.json').write_text(json.dumps(counts, indent=2) + '\n')
checks = json.loads((out / 'historical-and-environment-checks.json').read_text())
for c in json.loads((out/'scope-corrected-checks/checks-terminal.json').read_text()):
    c['log']='scope-corrected-checks/'+c['log'];checks.append(c)
assert all(c['returncode'] == 0 for c in checks), checks
histories = []
for p in sorted(set(out.glob('j72-*.xml')) | set(out.glob('guard-compatibility*.xml')) | set(out.glob('discovery-*.xml')) | set(out.glob('guard-corrections*.xml'))):
    cases = list(ET.parse(p).iter('testcase'))
    selected=p.name=='guard-corrections-terminal.xml'
    histories.append({'file': p.name, 'tests': len(cases), 'failure': sum(c.find('failure') is not None for c in cases), 'error': sum(c.find('error') is not None for c in cases), 'terminal':selected,'additional_unique_coverage':173 if selected else 0,'overlapping_historical_nodes':78 if selected else 0})
deviation=json.loads((root/'execution-scope-report.json').read_text())
(v / 'verification-history.json').write_text(json.dumps({'focused_runs': histories, 'initial_full_suite':{'results':'suite-results.json','logs':['full-0.log','full-1.log'],'returncodes':[-15,-15],'terminal_coverage_added':0,'unreported_failure_marker':'full-0.log contains F before termination; no terminal traceback/XML was emitted, so no specific cause is asserted','reason':'scope guard gap allowed J31 QM/MD; stopped, audited and quarantined'},'provisional_guarded_suite':{'results':'guarded/suite-results.json','returncodes':[-15]*8,'terminal_coverage_added':0,'reason':'stopped to extend fence to module/session fixtures'},'final_guarded_suite':{'raw_results':'guarded-final/suite-results.json','raw_failures_retained':True,'whole_file_correction':'guard-corrections-terminal.xml','correction_reason':'fence blocked read-only find_spec(psi4) and pinned J37 provenance reconstruction; corrected discovery and exact read-only entry point, preserving native/MD/refit blocks','selection_rule':'all 45 J35-C1 and 33 J37 nodes replaced unconditionally by whole-file corrections; all 173 J72 nodes use final rerun; other baseline nodes retain their original results','corrected_historical_nodes':78,'new_J72_unique_nodes':173,'terminal_result_source_mapping':'terminal-selection.json'},'generation_failure':{'file':'adjudication-first.log','reason':'archive metadata has no chemical_domain; explicit unknown fallback added'},'first_replay_failure':{'file':'replay-first.log','reason':'JSON hierarchy integer keys normalize to strings; canonical string keys added'},'isolated_pip_correction':{'file':'scope-corrected-checks/pip-PDF-scope-corrected.log','reason':'analysis PYTHONPATH exposed island egg-info without numpy to PDF-only environment; corrected isolated -I pip invocation passes without environment changes'},'earlier_lint_failure':{'file':'ruff-before-terminal.log','reason':'test import ordering; corrected'},'corrected_focused_failure':{'file':'j72-scope-terminal.log','reason':'stale assertion expected zero executions; replaced with independently reconstructed actual counts'},'first_sealing_failure':{'file':'diff-seal-first.log.xz','reason':'byte-exact pytest failure padding triggered diff whitespace checks; explicit J72-only attributes preserve these two raw artifacts without changing historical attributes or source checks'},'patch_context_failure':'one unapplied patch failed context matching; no changes executed; corrected patch applied','historical_failures_and_timing_disclosures':'all J54-J71 evidence unchanged','actual_execution_scope':deviation,'earlier_zero_count_plans_are_superseded_by':'../execution-scope-report.json'},indent=2)+'\n')
(v/'advisory-checks.json').write_text(json.dumps([{'name':'ambient-python-pip','returncode':1,'log':'ambient-pip.log','reason':'pre-existing torch missing fsspec and incompatible sympy; supported environments pass; no repair'}],indent=2)+'\n')
(v/'environment.json').write_text(json.dumps({'analysis_executable':sys.executable,'analysis_python':sys.version,'packages':{n:importlib.metadata.version(n) for n in ('numpy','scipy','pytest','rdkit','openmm')},'optional_PDF_reader':'/tmp/island-j71-validation/pdf-env; PyMuPDF1.26.7','thread_limits':{'OPENBLAS_NUM_THREADS':1,'OMP_NUM_THREADS':1,'OPENMM_CPU_THREADS':1},'full_suite_software_workers':8,'test_scope_fence':'scripts.phase_4j72_test_policy plus fail-closed descendant startup','phase_execution_scope':deviation,'reserved_validation':'untouched'},indent=2)+'\n')
doc=repo/'docs/phase_4j72.md';s=doc.read_text();s=s.replace('Verification totals and terminal checks are recorded in `validation/terminal-test-counts.json` and the receipt.',f"The unique terminal total is **{len(nodes):,} tests: {len(nodes)-skips:,} passed, {skips:,} skipped, zero failures/errors** ({len(base['nodes']):,} historical nodes plus {len(new)} J72 nodes). All **8,835 historical files** are byte-preserved. Focused reruns add zero unique coverage. Terminal details are recorded in `validation/terminal-test-counts.json` and the receipt. These are software results, not a claim of compliance with the violated allocation.")
s=s.replace('later tests extend mutation coverage.', 'the terminal J72 run passes all 173 tests, including complete outer-chain re-signing, strict count schemas, the child fence, quarantined execution reconstruction and terminal-result source closure.')
doc.write_text(s)
sys.path[:0]=[str(repo),str(repo/'src')]
from island.pcff_completion import pin
from scripts.validate_phase_4j72_receipt import FILES
(v/'implementation-files.json').write_text(json.dumps({n:pin(repo/n) for n in FILES},indent=2)+'\n')
subprocess.run(['git','add',*FILES,'docs/evidence/phase_4j72'],check=True)
with (v/'diff-seal-final.log').open('w') as log:
    p=subprocess.run(['git','diff','--cached','--check'],stdout=log,stderr=subprocess.STDOUT,check=False)
assert p.returncode==0
checks.append({'name':'diff-seal-final','command':['git','diff','--cached','--check'],'cwd':str(repo),'returncode':p.returncode,'log':'diff-seal-final.log'})
(v/'checks.json').write_text(json.dumps(checks,indent=2)+'\n')
env={**os.environ,'PYTHONPATH':'.:src','OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1','OPENMM_CPU_THREADS':'1'}
cmd=[sys.executable,'scripts/validate_phase_4j72_receipt.py','docs/evidence/phase_4j72']
for label,suffix in (('seal',['--seal']),('strict-receipt',[])):
    with (out/(label+'.log')).open('w') as log:
        p=subprocess.run(cmd+suffix,env=env,stdout=log,stderr=subprocess.STDOUT,check=False)
    assert p.returncode==0,(out/(label+'.log')).read_text()
print(counts,flush=True)
print((out/'strict-receipt.log').read_text(),flush=True)
