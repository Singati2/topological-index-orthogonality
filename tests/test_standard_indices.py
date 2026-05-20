"""Regression tests for the 30 baseline indices on P_5 and C_6.

Reference values are pre-computed by hand; see comments for each
edge-by-edge breakdown.
"""
import math

import networkx as nx
import pytest

from src.standard_indices import compute_all


def test_compute_all_returns_30_keys_and_finite_values_on_P5():
    values = compute_all(nx.path_graph(5))
    assert len(values) == 30
    for name, v in values.items():
        assert math.isfinite(v) or math.isnan(v), (
            f"{name} returned non-finite, non-NaN value {v!r}"
        )


def test_M1_on_P5():
    # P_5 edge degree pairs: (1,2), (2,2), (2,2), (2,1).
    # M_1 = sum (d_u + d_v) = 3 + 4 + 4 + 3 = 14
    values = compute_all(nx.path_graph(5))
    assert values["M1"] == pytest.approx(14.0)


def test_M2_on_P5():
    # M_2 = sum (d_u * d_v) = 2 + 4 + 4 + 2 = 12
    values = compute_all(nx.path_graph(5))
    assert values["M2"] == pytest.approx(12.0)


def test_F_on_P5():
    # F = sum (d_u^2 + d_v^2) = (1+4)+(4+4)+(4+4)+(4+1) = 26
    values = compute_all(nx.path_graph(5))
    assert values["F"] == pytest.approx(26.0)


def test_R_on_P5():
    # R = sum 1/sqrt(d_u*d_v)
    #   = 1/sqrt(2) + 1/sqrt(4) + 1/sqrt(4) + 1/sqrt(2)
    #   = sqrt(2) + 1
    values = compute_all(nx.path_graph(5))
    assert values["R"] == pytest.approx(math.sqrt(2.0) + 1.0, rel=1e-10)


def test_W_on_P5():
    # W(P_n) = n*(n^2 - 1)/6; W(P_5) = 5*24/6 = 20
    values = compute_all(nx.path_graph(5))
    assert values["W"] == pytest.approx(20.0)


def test_M1_M2_on_C6():
    # C_6 is 2-regular with 6 edges, all (2,2). M_1 = M_2 = 24.
    values = compute_all(nx.cycle_graph(6))
    assert values["M1"] == pytest.approx(24.0)
    assert values["M2"] == pytest.approx(24.0)


def test_R_on_C6():
    # R(C_n) = n/2 for n >= 3.
    values = compute_all(nx.cycle_graph(6))
    assert values["R"] == pytest.approx(3.0, rel=1e-10)


def test_W_on_C6():
    # W(C_n) for even n = n^3/8 = 27 for n = 6.
    values = compute_all(nx.cycle_graph(6))
    assert values["W"] == pytest.approx(27.0)
