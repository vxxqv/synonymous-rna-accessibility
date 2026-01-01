"""Boundary and information-limit tests for graph barrier certificates."""
from pathlib import Path
import sys
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from barrier_identification import barrier, sharp_barrier_bounds, threshold_status, additive_background_barrier


def test_unmeasured_bypass_prevents_false_separation():
    graph=[[1,2],[0,3],[0,3],[1,2]]
    lower=[1,0,0,1];upper=[1,0,1,1]
    assert sharp_barrier_bounds(graph,lower,upper,0,3)==(0,1)
    assert threshold_status(graph,lower,upper,0,3,.5)=='unresolved'
    assert barrier(graph,[1,0,0,1],0,3)==1
    assert barrier(graph,[1,0,1,1],0,3)==0


def test_endpoints_are_not_implicitly_retained():
    graph=[[1],[0,2],[1]]
    assert threshold_status(graph,[0,1,1],[1,1,1],0,2,.5)=='endpoints_not_certified'


def test_adjacent_and_coincident_endpoints_have_zero_barrier():
    assert sharp_barrier_bounds([[1],[0]],[-4,-8],[10,3],0,1)==(0,0)
    assert sharp_barrier_bounds([[]],[-2],[3],0,0)==(0,0)


def test_endpoint_separation_makes_bounds_sharp():
    graph=[[1],[0,2],[1]]
    assert sharp_barrier_bounds(graph,[2,-1,3],[4,1,5],0,2)==(1,5)
    assert barrier(graph,[2,1,3],0,2)==1
    assert barrier(graph,[4,-1,5],0,2)==5


def test_invalid_input_cannot_issue_certificate():
    with pytest.raises(ValueError): barrier([[]],[0],1,0)
    with pytest.raises(ValueError): sharp_barrier_bounds([[1],[]],[0,0],[1,1],0,1)
    with pytest.raises(ValueError): sharp_barrier_bounds([[],[]],[0,0],[1,1],0,1)
    with pytest.raises(ValueError): sharp_barrier_bounds([[1],[0]],[2,0],[1,1],0,1)
    with pytest.raises(ValueError): threshold_status([[1],[0]],[0,0],[1,1],0,1,float('nan'))
    with pytest.raises(ValueError): additive_background_barrier(float('nan'),0,1)


def test_additive_buffer_identity_and_perturbation_sharpness():
    assert additive_background_barrier(2,0,3)==0
    assert additive_background_barrier(2,0,1)==1
    graph=[[1],[0,2],[1]]
    assert barrier(graph,[2,-1,2],0,2)-barrier(graph,[1,0,1],0,2)==2
