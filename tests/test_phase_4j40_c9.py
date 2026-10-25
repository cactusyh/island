import pytest

from island.analysis.ewald_root_cause import isolate, one_variable_experiment
from island.exceptions import AnalysisError


def test_first_failing_system_is_explicit():
 r=isolate({'point_charges':{'status':'passed'},'counterion_system':{'status':'blocked'}})
 assert r['first_failing_system']=='counterion_system' and r['ewald_unvalidated']
def test_matrix_order_is_preserved():
 r=isolate({x:{'status':'passed'} for x in ('point_charges','counterion_system','pairwise_control','small_fragment','linear_polymer','aromatic_polymer','packed_network')})
 assert r['status']=='closed' and r['first_failing_system'] is None
def test_variable_mutation_is_recorded_not_repaired():
 a={'identity':'a','status':'blocked','box':'a'};b={'identity':'b','status':'blocked','box':'b'}
 r=one_variable_experiment(baseline=a,variable='box_convention',changed=b)
 assert r['status']=='blocked' and not r['repair_justified'] and r['metadata_changed']==['box']
def test_unknown_variable_rejected():
 with pytest.raises(AnalysisError): one_variable_experiment(baseline={},variable='charge',changed={})
