"""Wuzi ratio-bound constants collapse to a single value at exact reductions.

At Wuzi-classical reductions (e.g. (alpha, beta, gamma) = (1, 0, 0)
gives W = M_2), the ratio psi(i, j) / h(i, j) over the 10 admissible
(i, j) pairs on Delta <= 4 should be a single value (the identity
constant relating Wuzi to the classical index).
"""
import math

import pytest


def psi(i: int, j: int, alpha: float, beta: float, gamma: float) -> float:
    return ((i * j) ** alpha) * ((i + j) ** beta) * math.exp(
        gamma * abs(i - j) / (i + j)
    )


ADMISSIBLE = [(i, j) for i in range(1, 5) for j in range(i, 5)]


def ratio_extremes(alpha: float, beta: float, gamma: float, h_fn):
    ratios = [psi(i, j, alpha, beta, gamma) / h_fn(i, j) for i, j in ADMISSIBLE]
    return min(ratios), max(ratios)


def test_M2_reduction_at_1_0_0():
    """At (alpha, beta, gamma) = (1, 0, 0), W = M_2; ratio is 1.0 everywhere."""
    c_min, c_max = ratio_extremes(
        1.0, 0.0, 0.0, lambda i, j: i * j
    )
    assert c_min == pytest.approx(1.0, rel=1e-10)
    assert c_max == pytest.approx(1.0, rel=1e-10)


def test_Randic_reduction_at_neg_half_0_0():
    """At (-0.5, 0, 0), W = R; ratio is 1.0 everywhere."""
    c_min, c_max = ratio_extremes(
        -0.5, 0.0, 0.0, lambda i, j: 1.0 / math.sqrt(i * j)
    )
    assert c_min == pytest.approx(1.0, rel=1e-10)
    assert c_max == pytest.approx(1.0, rel=1e-10)


def test_H_reduction_at_0_neg_1_0():
    """At (0, -1, 0), W = H/2; ratio is 0.5 everywhere."""
    c_min, c_max = ratio_extremes(
        0.0, -1.0, 0.0, lambda i, j: 2.0 / (i + j)
    )
    assert c_min == pytest.approx(0.5, rel=1e-10)
    assert c_max == pytest.approx(0.5, rel=1e-10)


def test_SCI_reduction_at_0_neg_half_0():
    """At (0, -0.5, 0), W = SCI; ratio is 1.0 everywhere."""
    c_min, c_max = ratio_extremes(
        0.0, -0.5, 0.0, lambda i, j: 1.0 / math.sqrt(i + j)
    )
    assert c_min == pytest.approx(1.0, rel=1e-10)
    assert c_max == pytest.approx(1.0, rel=1e-10)
