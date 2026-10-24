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

def test_control_matrix_requires_all_physical_controls():
 from island.analysis.long_range_reference import control_matrix
 assert {x['name'] for x in control_matrix()} == {'neutral_point_charge','counterion_pair','pairwise_reference','small_polymer_fragment','box_scaling','cutoff_tolerance_scaling'}

def test_control_comparison_rejects_metadata_mutation():
 from island.analysis.long_range_reference import compare_control
 a={k: 1 for k in ('coordinates','box','charges','masses','atom_order','exclusions','one_four_scaling','lj_mixing','cutoff','reciprocal_tolerance')}
 a.update({k: 0.0 for k in ('energy','coulomb','reciprocal','real_space','force_norm','virial','finite_volume_pressure')})
 b=dict(a,box=2)
 assert compare_control(openmm=a,lammps=b,tolerance=1e-5)['status']=='blocked'

def test_control_comparison_passes_only_complete_equal_record():
 from island.analysis.long_range_reference import compare_control
 a={k: 1 for k in ('coordinates','box','charges','masses','atom_order','exclusions','one_four_scaling','lj_mixing','cutoff','reciprocal_tolerance')}
 a.update({k: 0.0 for k in ('energy','coulomb','reciprocal','real_space','force_norm','virial','finite_volume_pressure')})
 assert compare_control(openmm=a,lammps=dict(a),tolerance=1e-5)['status']=='passed'
