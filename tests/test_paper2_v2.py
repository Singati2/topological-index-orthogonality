"""Paper 2 (v2) checks: span certificate, exact identities, numerically safe
partial correlation, and consistency of the committed result files and
generated tables with the manuscript."""
from __future__ import annotations

import filecmp
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from src import orthogonality as ortho
from src.standard_indices import ALL_INDICES, compute_all

PROJECT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
RESULTS = os.path.join(PROJECT, "results")
BID18 = ["M1", "M2", "mM1", "mM2", "F", "R", "SCI", "H", "GA", "AG",
         "ABC", "ABS", "AZI", "SO", "SO_red", "Alb", "Sigma", "redM1"]


def _load_script(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(PROJECT, "scripts", name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _random_g4_graphs(k=40, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    while len(out) < k:
        n = int(rng.integers(4, 14))
        G = nx.gnm_random_graph(n, int(rng.integers(n - 1, 2 * n)), seed=int(rng.integers(1e9)))
        if G.number_of_edges() == 0:
            continue
        G = nx.convert_node_labels_to_integers(G.subgraph(max(nx.connected_components(G), key=len)).copy())
        if G.number_of_edges() and max(d for _, d in G.degree()) <= 4:
            out.append(G)
    return out


# ---------------------------------------------------------------- certificate

def test_data_free_certificate_rank_10():
    cert = _load_script("70_bid_span_certificate")
    F = cert.f_matrix(BID18)
    assert F.shape == (18, 10)
    assert np.linalg.matrix_rank(F, tol=1e-9) == 10


def test_universal_identity_holds_on_random_graphs():
    cert = _load_script("70_bid_span_certificate")
    F10 = cert.f_matrix(cert.BASIS10)
    graphs = _random_g4_graphs()
    for lab, f in cert.VARIANTS.items():
        fv = np.array([f(i, j) for (i, j) in cert.P4])
        c = np.linalg.solve(F10.T, fv)
        for G in graphs:
            d = dict(G.degree())
            z = sum(f(*sorted((d[u], d[v]))) for u, v in G.edges())
            base = compute_all(G)
            approx = sum(c[k] * base[b] for k, b in enumerate(cert.BASIS10))
            assert approx == pytest.approx(z, rel=1e-8, abs=1e-8), lab


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_exact_baseline_identities(seed):
    for G in _random_g4_graphs(15, seed):
        b = compute_all(G)
        n, m = G.number_of_nodes(), G.number_of_edges()
        assert b["Sigma"] == pytest.approx(b["F"] - 2 * b["M2"])
        assert b["redM1"] == pytest.approx(b["M1"] - 4 * m + n)


def test_span_certified_pcor_exact_zero_and_noise():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(400, 6)) * np.array([1, 10, 100, 1e3, 1e4, 1e5])
    y = X[:, 0] + rng.normal(size=400)
    z = X @ np.arange(1, 7) + 7.0
    r = ortho.span_certified_pcor(z, X, y)
    assert r["in_span"] and r["pcor"] == 0.0
    w = rng.normal(size=400)
    r2 = ortho.span_certified_pcor(w, X, y)
    assert not r2["in_span"] and abs(r2["pcor"]) < 0.2
    lo, hi = ortho.fisher_ci(r2["pcor"], r2["df"])
    assert lo < r2["pcor"] < hi


# --------------------------------------------------------------- result files

def _csv(name):
    p = os.path.join(RESULTS, name)
    if not os.path.exists(p):
        pytest.skip(f"missing {p}")
    return pd.read_csv(p)


def test_certificate_holds_on_all_datasets():
    c = _csv("bid_span_certificate.csv")
    assert set(c.dataset) == {"esol", "freesolv", "lipophilicity", "bbbp"}
    assert c.certificate_holds.all()
    assert (c.rank_BID18 == c.rank_mij).all()
    assert (c.max_degree <= 4).all()


def test_census_in_span_and_controls_not():
    v = _csv("bid_variant_screen.csv")
    bid = v[v.family != "non-BID control"]
    assert bid.in_span.all()
    assert (bid.pcor_certified == 0).all()
    assert bid.rel_resid_on_baseline30.max() < 1e-10
    ctrl = v[v.family == "non-BID control"]
    assert (~ctrl.in_span).all() and ctrl.rel_resid_on_baseline30.min() > 1e-4
    assert v[v.family == "published BID"]["index"].nunique() == 15


def test_loyola_grid_pcor_exactly_zero():
    for ds in ["esol", "freesolv", "lipophilicity"]:
        g = _csv(os.path.join(ds, "wuzi_grid.csv"))
        assert len(g) == 100
        assert g.in_span.all()
        assert (g.partial_corr_target == 0).all()


def test_generated_tables_match_committed():
    """docs/tables_v2 must be exactly what scripts/75 produces from results/."""
    committed = os.path.join(PROJECT, "docs", "tables_v2")
    if not os.path.isdir(committed):
        pytest.skip("no committed tables")
    tmp = tempfile.mkdtemp()
    try:
        shutil.copytree(os.path.join(PROJECT, "results"), os.path.join(tmp, "results"))
        os.makedirs(os.path.join(tmp, "scripts"))
        shutil.copytree(os.path.join(PROJECT, "src"), os.path.join(tmp, "src"))
        for f in ("70_bid_span_certificate.py", "75_make_v2_tables.py"):
            shutil.copy(os.path.join(PROJECT, "scripts", f), os.path.join(tmp, "scripts"))
        subprocess.run([sys.executable, os.path.join(tmp, "scripts", "75_make_v2_tables.py")],
                       check=True, capture_output=True)
        for f in os.listdir(committed):
            assert filecmp.cmp(os.path.join(committed, f), os.path.join(tmp, "docs", "tables_v2", f),
                               shallow=False), f
    finally:
        shutil.rmtree(tmp)


def test_dataset_checksums():
    """Raw MoleculeNet files must match data/SHA256SUMS (verified copies)."""
    import hashlib
    sums = os.path.join(PROJECT, "data", "SHA256SUMS")
    if not os.path.exists(sums):
        pytest.skip("no SHA256SUMS")
    for line in open(sums).read().split("\n"):
        if not line.strip():
            continue
        digest, name = line.split()
        with open(os.path.join(PROJECT, "data", name), "rb") as fh:
            assert hashlib.sha256(fh.read()).hexdigest() == digest, name


def test_curated_counts():
    exp = {"esol": (1115, 628), "freesolv": (639, 236), "lipophilicity": (4200, 3472), "bbbp": (1955, 1693)}
    for name, (n, g) in exp.items():
        p = os.path.join(PROJECT, "data", "curated", f"{name}_curated.csv")
        if not os.path.exists(p):
            pytest.skip("no curated data")
        df = pd.read_csv(p)
        assert len(df) == n and df.graph_wl_hash.nunique() == g, name


# ------------------------------------------------------------- v2.1 additions

def test_certify_tool_data_free_and_identity():
    import math
    from src.certify import certify, phi_matrix, BID18 as B18
    assert np.linalg.matrix_rank(phi_matrix(B18), tol=1e-9) == 10
    graphs = _random_g4_graphs(25, 7)
    rep = certify(graphs, f=lambda a, b: math.sqrt(a * a + b * b + a * b))   # Euler Sombor
    assert rep["data_free_certified"] and rep["certified"] and rep["in_span"]
    assert rep["identity_max_error"] < 1e-9
    # hyper-Zagreb (a+b)^2 = F + 2 M2 exactly; coefficients are unique only on a
    # 10-index basis (on the 18-index baseline they are unique up to ker Phi^T)
    basis10 = ["M1", "M2", "mM1", "mM2", "F", "R", "SCI", "H", "GA", "AG"]
    rep2 = certify(graphs, f=lambda a, b: (a + b) ** 2, baseline=basis10)
    c = rep2["coefficients"]
    assert c["F"] == pytest.approx(1.0, abs=1e-8) and c["M2"] == pytest.approx(2.0, abs=1e-8)
    assert all(abs(v) < 1e-8 for k, v in c.items() if k not in ("F", "M2"))


def test_certificate_generality_results():
    g = _csv("certificate_generality.csv")
    assert len(g) >= 13
    ok = g[g.dataset != "HIV"]
    assert ok.certificate_holds.all() and (ok.max_rel_resid_mij_on_BID18 < 1e-8).all()
    hiv = g[g.dataset == "HIV"].iloc[0]
    assert (not hiv.certificate_holds) and hiv.n_maxdeg_gt4 > 0 and hiv.rank_BID18 < hiv.rank_mij
    if "restricted_certificate_holds" in g.columns:
        assert bool(hiv.restricted_certificate_holds)


def test_chemistry_baseline_screen_results():
    s = _csv("screen_vs_chemistry_summary.csv")
    assert set(s.dataset) == {"esol", "freesolv", "lipophilicity", "bbbp"}
    assert int(s.cand_pass_both.sum()) == 0
    c = _csv("screen_vs_chemistry_candidates.csv")
    assert (c.verdict_both != "PASS").all()


def test_graph_oracle_floor_consistent_with_curation():
    o = _csv("graph_oracle_floor.csv").set_index("dataset")
    cur = _csv("benchmark_v2_curation.csv").set_index("dataset")
    for d in ["esol", "freesolv", "lipophilicity"]:
        # in-sample group-mean RMSE equals sqrt(SS_within / N) reported by scripts/72
        assert o.loc[d, "oracle_in_sample_rmse"] == pytest.approx(
            cur.loc[d, "graph_only_rmse_floor(sqrt SSwithin/N)"], rel=1e-3)
