import pytest

from island.analysis.long_range_reference import (
 closure_record,
 pairwise_coulomb_control,
)
from island.exceptions import AnalysisError


def test_pairwise_control():
 r=pairwise_coulomb_control(1,-1,1);assert r['energy_kj_mol']<0 and r['force_norm_kj_mol_nm']>0

def test_invalid_pairwise_control():
 with pytest.raises(AnalysisError): pairwise_coulomb_control(1,1,0)

def test_toy_agreement_does_not_validate_polymer():
 r=closure_record(controls=[{'status':'passed'}],polymer_cells=[{'status':'blocked'}],executable_manifest={'status':'available'},tolerance=1e-5)
 assert r['status']=='ewald_unvalidated' and not r['scientific_acceptance']

def test_polymer_validation_requires_all_cells():
 r=closure_record(controls=[{'status':'passed'}],polymer_cells=[{'status':'passed'},{'status':'blocked'}],executable_manifest={},tolerance=1e-5)
 assert not r['polymer_closure_passed']
