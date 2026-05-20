"""is_caterpillar predicate, loaded from scripts/10_wuzi_extremal_search.py.

The scripts/ directory is not a Python package, so we load the function
via importlib.util.spec_from_file_location.
"""
import importlib.util
import os

import networkx as nx
import pytest


def _load_is_caterpillar():
    here = os.path.dirname(__file__)
    script_path = os.path.abspath(os.path.join(
        here, "..", "scripts", "10_wuzi_extremal_search.py"
    ))
    spec = importlib.util.spec_from_file_location(
        "_ext_search_for_tests", script_path
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.is_caterpillar


is_caterpillar = _load_is_caterpillar()


@pytest.mark.parametrize("n", [3, 5, 7, 10])
def test_path_is_caterpillar(n):
    assert is_caterpillar(nx.path_graph(n)), f"P_{n} should be a caterpillar"


@pytest.mark.parametrize("n", [3, 5, 7])
def test_star_is_caterpillar(n):
    # nx.star_graph(k) has k+1 vertices total: 1 centre + k leaves.
    assert is_caterpillar(nx.star_graph(n - 1)), f"S_{n} should be a caterpillar"


def test_handmade_caterpillar_8_vertices():
    # Internal path: 0 - 1 - 2 - 3, with leaves attached at 1, 2, 2, 3.
    G = nx.Graph()
    G.add_edges_from([
        (0, 1), (1, 2), (2, 3),
        (1, 4),
        (2, 5), (2, 6),
        (3, 7),
    ])
    assert is_caterpillar(G)


def test_order_7_spider_is_not_caterpillar():
    # Centre vertex 0, three length-2 limbs: 0-1-2, 0-3-4, 0-5-6.
    # Removing leaves leaves the internal subgraph {0, 1, 3, 5} with edges
    # (0, 1), (0, 3), (0, 5) -- a star K_{1,3}, not a path.
    G = nx.Graph()
    G.add_edges_from([
        (0, 1), (1, 2),
        (0, 3), (3, 4),
        (0, 5), (5, 6),
    ])
    assert not is_caterpillar(G), (
        "order-7 spider with three length-2 limbs is NOT a caterpillar"
    )
