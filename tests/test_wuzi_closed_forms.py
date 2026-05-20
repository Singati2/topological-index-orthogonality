"""Wuzi closed forms on standard graph families match direct computation.

Closed forms are pre-derived in src/wuzi_analytical.py and mirrored here
as small Python functions so the test file is self-contained.
"""
import math

import networkx as nx
import pytest

from src.wuzi_index import wuzi


PARAMETER_TRIPLES = [
    (1.0, 0.0, 0.0),
    (-0.5, 0.0, 0.0),
    (0.5, -0.5, 1.0),
    (0.0, 0.0, 2.0),
    (-1.0, -1.0, 1.0),
]


def closed_form_Kn(n: int, a: float, b: float, _g: float) -> float:
    """K_n is (n-1)-regular with m = n(n-1)/2 edges; gamma irrelevant on regular graphs."""
    m = n * (n - 1) // 2
    return m * (2.0 ** b) * ((n - 1) ** (2 * a + b))


def closed_form_Cn(n: int, a: float, b: float, _g: float) -> float:
    """C_n is 2-regular; gamma irrelevant on regular graphs."""
    return n * (4.0 ** (a + b))


def closed_form_Sn(n: int, a: float, b: float, g: float) -> float:
    """S_n: all n-1 edges have degree pair (1, n-1)."""
    return (n - 1) * ((n - 1) ** a) * (n ** b) * math.exp(g * (n - 2) / n)


def closed_form_Kpq(p: int, q: int, a: float, b: float, g: float) -> float:
    """K_{p,q}: all p*q edges have degree pair (p, q)."""
    return p * q * ((p * q) ** a) * ((p + q) ** b) * math.exp(
        g * abs(p - q) / (p + q)
    )


@pytest.mark.parametrize("a,b,g", PARAMETER_TRIPLES)
@pytest.mark.parametrize("n", [4, 5, 6])
def test_closed_form_Kn(n, a, b, g):
    actual = wuzi(nx.complete_graph(n), a, b, g)
    expected = closed_form_Kn(n, a, b, g)
    assert actual == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("a,b,g", PARAMETER_TRIPLES)
@pytest.mark.parametrize("n", [4, 5, 6])
def test_closed_form_Cn(n, a, b, g):
    actual = wuzi(nx.cycle_graph(n), a, b, g)
    expected = closed_form_Cn(n, a, b, g)
    assert actual == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("a,b,g", PARAMETER_TRIPLES)
@pytest.mark.parametrize("n", [4, 5, 6])
def test_closed_form_Sn(n, a, b, g):
    # nx.star_graph(k) has k+1 vertices: one centre + k leaves.
    G = nx.star_graph(n - 1)
    actual = wuzi(G, a, b, g)
    expected = closed_form_Sn(n, a, b, g)
    assert actual == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("a,b,g", PARAMETER_TRIPLES)
@pytest.mark.parametrize("p,q", [(2, 3), (3, 4), (2, 4)])
def test_closed_form_Kpq(p, q, a, b, g):
    G = nx.complete_bipartite_graph(p, q)
    actual = wuzi(G, a, b, g)
    expected = closed_form_Kpq(p, q, a, b, g)
    assert actual == pytest.approx(expected, rel=1e-10)
