"""Lossless repack: deflated NPZ containers defeat outer archive deduplication."""
from pathlib import Path
import tempfile,shutil
import numpy as np
from island.pcff_completion import inventory,write,signed,pin
from island.pcff_qm_pilot.bundle import pack,unpack
from island.pcff_qm_pilot.hmds_precision import read
root=Path('docs/evidence/phase_4j66')
with tempfile.TemporaryDirectory(prefix='j66-lossless-canonical-') as temp:
 target=Path(temp)/'records';target.mkdir()
 for file in (root/'bootstrap').iterdir():
  if file.suffix=='.npz':
   with np.load(file,allow_pickle=False) as stored:np.savez(target/file.name,**{k:stored[k] for k in stored.files})
  else:shutil.copyfile(file,target/file.name)
 # The exact numerical arrays and all JSON native solver states are unchanged.
 for file in (root/'bootstrap').iterdir():
  if file.suffix=='.npz':
   with np.load(file) as a,np.load(target/file.name) as b:assert a.files==b.files and all(np.array_equal(a[k],b[k]) for k in a.files)
  else:assert file.read_bytes()==(target/file.name).read_bytes()
 write(root/'bootstrap-file-inventory.json',inventory(target))
 pack(target,root/'bootstrap-solves-repacked.tar.xz')
 with tempfile.TemporaryDirectory(prefix='j66-archive-check-') as verification:
  unpack(root/'bootstrap-solves-repacked.tar.xz',Path(verification)/'records')
  assert inventory(Path(verification)/'records')==inventory(target)
shutil.move(root/'bootstrap-solves.tar.xz',Path('/tmp/island-j66-validation/deflated-container-bootstrap.tar.xz'))
shutil.move(root/'bootstrap-solves-repacked.tar.xz',root/'bootstrap-solves.tar.xz')
for name in ('block-uncertainty.json','compute-accounting.json'):
 data=read(root/name);data.pop('identity')
 key='full_records_archive_sha256' if name.startswith('block') else 'full64_records_archive_sha256'
 data[key]=pin(root/'bootstrap-solves.tar.xz');write(root/name,signed(data))
write(root/'archive-packaging.json',signed({'representation':'Canonical NPZ containers stored without inner deflate, then existing tar.xz bundle; every matrix array element and JSON solver record independently compared exactly before collection. No numerical rerun or change.','exact_arrays_and_solver_records_preserved':True,'raw_states':64,'canonical_problems':120,'archive_sha256':pin(root/'bootstrap-solves.tar.xz')}))
print((root/'bootstrap-solves.tar.xz').stat().st_size)
