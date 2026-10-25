import pytest

from island.analysis.ewald_diagnosis import compare_frame, diagnose
from island.exceptions import AnalysisError


def record():
 fields={k:k for k in ('graph','coordinates','box','atom_order','masses','charges','charge_origin','exclusions','one_four_scaling','lj_mixing','cutoff','switching','reciprocal_tolerance','units','precision')}
 return dict(fields,**{k:0.0 for k in ('total','coulomb','real_space_coulomb','reciprocal','lj','bonded','angle','torsion','improper','force_norm','virial','pressure_proxy')})
def test_equal_components_close(): assert compare_frame(openmm=record(),lammps=record(),tolerance=1e-5)['status']=='passed'
def test_box_mutation_blocks():
 b=record();c=record();c['box']='changed';assert compare_frame(openmm=b,lammps=c,tolerance=1e-5)['status']=='blocked'
def test_charge_mutation_blocks():
 b=record();c=record();c['charges']='changed';assert compare_frame(openmm=b,lammps=c,tolerance=1e-5)['status']=='blocked'
def test_bad_tolerance_rejected():
 with pytest.raises(AnalysisError): compare_frame(openmm=record(),lammps=record(),tolerance=0)
def test_incomplete_matrix_stays_unvalidated():
 r=diagnose({'pairwise':{'status':'passed'}});assert r['status']=='blocked' and r['ewald_unvalidated']
