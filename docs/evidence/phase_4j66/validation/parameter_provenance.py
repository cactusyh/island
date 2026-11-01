"""Expose inherited graph/type/source bindings without changing any source row."""
from pathlib import Path
from island.pcff_completion import pin,signed,write
from island.pcff_qm_pilot.hmds_baseline import IDENTITY,corrections
from island.pcff_qm_pilot.hmds_hamiltonian import load_model,parts
from island.pcff_qm_pilot.hmds_precision import read
repo=Path.cwd();root=repo/'docs/evidence/phase_4j66';plan=read(root/'plan/scientific-plan.json');old=read(repo/'docs/evidence/phase_4j65/plan/scientific-plan.json');model=load_model(repo)
p=parts(model,old['center_coordinates']);c=corrections(model,old['center_coordinates'])
rows=[]
for spec in plan['parameter_specifications']:
 name=spec['name'];queries=[row for row in p['query_to_columns'] if name in row['columns']]
 rows.append({**spec,'origin':'explicit missing original-source coupling; trial value is research only' if queries else 'additive J64 quadratic research correction; not replacement source row','graph_applicable_source_queries':queries,'baseline_dependencies':[row for row in c['source_dependencies'] if row['column']==name]})
write(root/'parameter-provenance.json',signed({'Hamiltonian':IDENTITY,'model_identity':model['identity'],'graph_identity':model['graph_identity'],'atom_order':model['atom_order'],'source_hashes':model['source_hashes'],'source_record_sha256':model['source_record_sha256'],'source_record_identity':model['source_record_identity'],'source_typing_identity':model['source_typing_identity'],'charge_identity':model['charge_identity'],'charge_origin':model['charge_origin'],'conventions':model['conventions'],'equation_authority_sha256':model['J60_equation_authority_sha256'],'units':model['units'],'specifications':rows,'energy_intercept':{'units':'kcal/mol','scope':'single energy-only training nuisance; no physical source row or force derivative'},'fixed_MM_arrays':'physical-arrays.npz:fixed_energy,fixed_force; per22stable configurations','original_author_artifact_status':'not_yet_supplied','confirmation_authority':'28 user-attested confirmations; no author artifacts created','source_rows_changed':False,'producer_receipt_sha256':pin(repo/'docs/evidence/phase_4j65/receipt.json')}))
