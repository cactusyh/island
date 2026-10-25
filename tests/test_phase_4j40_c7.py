from scripts.validate_phase_4j40_c7 import inventory
from scripts.validate_phase_4j40_c7_receipt import validate


def test_real_kspace_inventory():
 r=inventory('../island-validation/phase4j22/lmp-kspace'); assert r['status']=='available'; assert len(r['sha256'])==64; assert 'Large-scale' in r['version']

def test_missing_kspace_is_asset_unavailable(tmp_path):
 assert inventory(str(tmp_path/'missing'))['status']=='asset_unavailable'

def test_receipt_blocks_polymer_closure():
 r=validate('docs/evidence/phase_4j40_c7/receipt.json'); assert r['scientific_acceptance'] is False
