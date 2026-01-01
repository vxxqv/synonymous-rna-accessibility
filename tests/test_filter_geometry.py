"""Counterexample-style checks of the exploratory filter geometry methods."""
import sys,itertools,math
from pathlib import Path
import pytest
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'analysis/expansion'));sys.path.insert(0,str(ROOT/'src'))
from retention_filtration import filtration
from conditional_vertex_null import conditional_null
from accessibility import components

def test_path_valley_bar_and_integrated_excess():
    r=filtration([[1],[0,2],[1]],[.9,.2,.8],full=True)
    assert r['longest_lifetime']==pytest.approx(.6)
    assert r['positive_bar_count']==1
    assert r['integrated_excess']==pytest.approx(.3)
    assert r['maximum_excess']==.5

def test_baseline_disconnection_is_not_extra_persistence():
    r=filtration([[1],[0],[3],[2]],[.9,.8,.7,.6],full=True)
    assert r['baseline_components']==2
    assert r['positive_bar_count']==0
    assert r['integrated_excess']==0

def test_equal_level_activation_has_no_positive_tie_bar():
    r=filtration([[1],[0,2],[1]],[.5,.5,.5],full=True)
    assert r['positive_bar_count']==0
    assert r['longest_lifetime']==0

def test_conditional_null_exact_known_path():
    g=[[1],[0,2],[1]];labels,_=components(g)
    r=conditional_null(g,labels,{0,2},False,1729)
    assert r['exact'] and r['population_size']=='3'
    assert r['null_mean']==pytest.approx(1/6)
    assert r['observed_excess']==.5
    assert r['statistic_varied']

def test_component_conditioning_retains_serine_split_denominator():
    g=[[1],[0],[3],[2]];labels,_=components(g)
    r=conditional_null(g,labels,{0,2},False,1729)
    assert r['connection_denominator']=='2'
    assert r['observed_excess']==0
    assert r['null_max']==0

def test_degree_null_exposes_an_invariant_control():
    g=[[1],[0,2],[1]];labels,_=components(g)
    r=conditional_null(g,labels,{0,2},True,1729)
    assert r['population_size']=='1'
    assert not r['statistic_varied']
