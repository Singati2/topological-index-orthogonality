"""Wuzi family reduces to classical indices at four parameter settings.

Expected values are pre-computed literal expressions; see EXPECTED below.
"""
import math

import networkx as nx
import pytest

from src.wuzi_index import wuzi


# Pre-computed expected values for each (graph, classical-index) pair.
# P_5: edges (1,2),(2,2),(2,2),(2,1).
# C_6: edges all (2,2), 6 edges.
# K_{2,3}: edges all (3,2), 6 edges.
EXPECTED = {
    "P5": {
        "graph": nx.path_graph(5),
        "M2":  12.0,                                # 2+4+4+2
        "R":   math.sqrt(2.0) + 1.0,                # 2/sqrt(2) + 2/2
        "SCI": 2.0 / math.sqrt(3.0) + 1.0,          # 2/sqrt(3) + 2*0.5
        "H":   7.0 / 3.0,                           # 2*(2/3) + 2*(2/4)
    },
    "C6": {
        "graph": nx.cycle_graph(6),
        "M2":  24.0,
        "R":   3.0,
        "SCI": 3.0,                                 # 6/sqrt(4)
        "H":   3.0,                                 # 6 * 2/4
    },
    "K23": {
        "graph": nx.complete_bipartite_graph(2, 3),
        "M2":  36.0,                                # 6 edges * 6
        "R":   6.0 / math.sqrt(6.0),
        "SCI": 6.0 / math.sqrt(5.0),
        "H":   12.0 / 5.0,                          # 6 * 2/5
    },
}


@pytest.mark.parametrize("name", list(EXPECTED))
def test_wuzi_1_0_0_equals_M2(name):
    G = EXPECTED[name]["graph"]
    expected = EXPECTED[name]["M2"]
    assert wuzi(G, 1.0, 0.0, 0.0) == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("name", list(EXPECTED))
def test_wuzi_neg_half_0_0_equals_R(name):
    G = EXPECTED[name]["graph"]
    expected = EXPECTED[name]["R"]
    assert wuzi(G, -0.5, 0.0, 0.0) == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("name", list(EXPECTED))
def test_wuzi_0_neg_half_0_equals_SCI(name):
    G = EXPECTED[name]["graph"]
    expected = EXPECTED[name]["SCI"]
    assert wuzi(G, 0.0, -0.5, 0.0) == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("name", list(EXPECTED))
def test_2_times_wuzi_0_neg_1_0_equals_H(name):
    # H(G) = 2 * W(G; 0, -1, 0) by definition.
    G = EXPECTED[name]["graph"]
    expected = EXPECTED[name]["H"]
    assert 2.0 * wuzi(G, 0.0, -1.0, 0.0) == pytest.approx(expected, rel=1e-10)
