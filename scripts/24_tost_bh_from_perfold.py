"""Deterministic recomputation of Table 11/12 inference from per-fold scores.

Paper 2, matched-model-class benchmark (Section 6.3.2). This script closes a
reproducibility gap: the paired-fold bootstrap CIs, the TOST equivalence
verdicts, and the Benjamini--Hochberg multiple-comparison correction that
appear in the manuscript had NO committed producer. ``scripts/18_pipeline_
then_lasso.py`` only emits per-fold primary-metric scores and a two-sided
paired t-test; the TOST/bootstrap CSV (``results/paper2_table11_tost_
bootstrap.csv``) was generated outside the committed pipeline and the BH
correction was done by hand.

This script recomputes ALL of that deterministically from the existing
per-fold output, with a fixed bootstrap seed, and does NOT re-run the
benchmark.

Input
-----
results/paper2_pipeline_then_lasso_perfold.csv
    Long format: dataset, task, config, fold, metric, value
    (RMSE for regression datasets, ROC-AUC for BBBP; one row per fold.)

Outputs (NEW files; nothing existing is overwritten)
----------------------------------------------------
results/paper2_table11_tost_bootstrap_recomputed.csv
results/paper2_bh_correction_recomputed.csv
results/paper2_recompute_verification.txt   (console report, also saved)

Method (matches the conventions already used in the manuscript)
---------------------------------------------------------------
* Paired difference per fold:  diff = metric(config) - metric(full).
* Paired t-test: scipy.stats.ttest_rel on the finite paired folds,
  two-sided, df = (#finite folds) - 1.   For the matched 5-fold design df=4.
* Equivalence margin: +/-5% of the full-baseline mean primary metric
  (RMSE for regression, AUC for BBBP), i.e. equiv_bound = 0.05 * mean(full).
* t-based 90% CI: mean +/- t_{0.95, df} * SE,  SE = std(diff, ddof=1)/sqrt(n).
  The TOST equivalence verdict is "the 90% CI lies entirely within
  (-bound, +bound)", which is exactly the two-one-sided-tests procedure at
  alpha = 0.05 (Schuirmann).  The two one-sided p-values are also reported.
* Paired-fold bootstrap: resample the finite paired differences with
  replacement B = 10000 times (fixed seed), percentile 95% CI of the mean.
  With only 5 folds the bootstrap is a sensitivity summary, not independent
  confirmation of the t-based verdict (stated in the manuscript).
* Benjamini--Hochberg step-up at q = 0.05 over the family of tests with a
  FINITE p-value.  The degenerate Lipophilicity ``combined_pruned`` row
  (zero features in 4/5 folds -> only 1 finite fold -> no t-test) is
  EXCLUDED from the family, giving m = 7.  Bonferroni is reported alongside.
"""
from __future__ import annotations

import csv
import math
import os

import numpy as np
from scipy import stats

# ----------------------------------------------------------------------------
THIS_DIR = os.path.abspath(os.path.dirname(__file__))
PROJECT = os.path.abspath(os.path.join(THIS_DIR, ".."))
RESULTS = os.path.join(PROJECT, "results")

PERFOLD_CSV = os.path.join(RESULTS, "paper2_pipeline_then_lasso_perfold.csv")
EXISTING_TOST_CSV = os.path.join(RESULTS, "paper2_table11_tost_bootstrap.csv")

OUT_TOST = os.path.join(RESULTS, "paper2_table11_tost_bootstrap_recomputed.csv")
OUT_BH = os.path.join(RESULTS, "paper2_bh_correction_recomputed.csv")
OUT_REPORT = os.path.join(RESULTS, "paper2_recompute_verification.txt")

BOOTSTRAP_SEED = 20260531   # canonical seed (original artifact had no committed seed)
N_BOOTSTRAP = 10000
EQUIV_PCT = 0.05            # +/-5% equivalence margin
ALPHA = 0.05               # TOST one-sided alpha AND BH/Bonferroni q
CONFIGS = ["pairwise_pruned", "combined_pruned"]   # tested vs "full"
DATASET_ORDER = ["esol", "freesolv", "lipophilicity", "bbbp"]

# Manuscript Table 11 reported rounded p-values, for the verification report.
MANUSCRIPT_P = {
    ("esol", "pairwise_pruned"): 0.020,
    ("esol", "combined_pruned"): 0.022,
    ("freesolv", "pairwise_pruned"): 0.748,
    ("freesolv", "combined_pruned"): 0.672,
    ("lipophilicity", "pairwise_pruned"): 0.026,
    ("bbbp", "pairwise_pruned"): 0.163,
    ("bbbp", "combined_pruned"): 0.194,
}


def load_perfold(path):
    """Return {(dataset, config): np.array of per-fold values ordered by fold},
    plus {dataset: task}."""
    by_key = {}
    task_of = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            ds, task, cfg = row["dataset"], row["task"], row["config"]
            fold = int(row["fold"])
            try:
                val = float(row["value"])
            except ValueError:
                val = float("nan")
            by_key.setdefault((ds, cfg), {})[fold] = val
            task_of[ds] = task
    arrays = {}
    for (ds, cfg), folds in by_key.items():
        arrays[(ds, cfg)] = np.array([folds[k] for k in sorted(folds)], dtype=float)
    return arrays, task_of


def bootstrap_ci(diffs, seed, B, lo=2.5, hi=97.5):
    """Percentile bootstrap CI of the mean of paired differences."""
    rng = np.random.default_rng(seed)
    n = len(diffs)
    idx = rng.integers(0, n, size=(B, n))
    means = diffs[idx].mean(axis=1)
    return float(np.percentile(means, lo)), float(np.percentile(means, hi))


def fmt(x):
    return "nan" if (x is None or (isinstance(x, float) and not math.isfinite(x))) else f"{x:+.4f}"


def main():
    lines = []   # verification report

    def out(s=""):
        print(s)
        lines.append(s)

    arrays, task_of = load_perfold(PERFOLD_CSV)

    out("=" * 78)
    out("RECOMPUTE: TOST / bootstrap / BH from per-fold scores")
    out(f"  input : {os.path.relpath(PERFOLD_CSV, PROJECT)}")
    out(f"  bootstrap seed = {BOOTSTRAP_SEED}, B = {N_BOOTSTRAP}, equiv = +/-{EQUIV_PCT:.0%}")
    out("=" * 78)

    tost_rows = []
    pvalue_family = []   # (dataset, config, p) for finite-p tests only

    for ds in DATASET_ORDER:
        full = arrays.get((ds, "full"))
        if full is None:
            continue
        full_finite = full[np.isfinite(full)]
        equiv_bound = EQUIV_PCT * float(np.mean(full_finite))
        for cfg in CONFIGS:
            arr = arrays.get((ds, cfg))
            if arr is None:
                continue
            diffs_all = arr - full
            mask = np.isfinite(diffs_all)
            diffs = diffs_all[mask]
            n = int(mask.sum())

            rec = {
                "dataset": ds, "config": cfg, "n_folds": n,
                "mean_paired_diff": "", "t_stat": "", "p_two_sided": "",
                "df": "", "t_ci90_lo": "", "t_ci90_hi": "",
                "boot_ci95_lo": "", "boot_ci95_hi": "",
                "equiv_bound_pm5pct": f"+/-{equiv_bound:.4f}",
                "tost_p": "", "tost_t_equivalent": "n/a",
                "boot_within_bound": "n/a", "note": "",
            }

            if n < 2:
                rec["note"] = ("combined screen degenerates to zero features in "
                               f"{int((~mask).sum())}/{len(diffs_all)} {ds} folds")
                tost_rows.append(rec)
                continue

            mean_diff = float(np.mean(diffs))
            sd = float(np.std(diffs, ddof=1))
            se = sd / math.sqrt(n)
            df = n - 1

            # paired t-test vs full (recomputed from per-fold; matches script 18)
            t_stat, p_two = stats.ttest_rel(arr[mask], full[mask])
            t_stat = float(t_stat)
            p_two = float(p_two)

            # t-based 90% CI
            tcrit = stats.t.ppf(1 - ALPHA, df)   # 0.95 quantile -> 90% two-sided
            ci_lo = mean_diff - tcrit * se
            ci_hi = mean_diff + tcrit * se

            # TOST: two one-sided tests at alpha; equivalent iff 90% CI within bound
            t_low = (mean_diff - (-equiv_bound)) / se   # H0: diff <= -bound
            p_low = float(stats.t.sf(t_low, df))
            t_high = ((equiv_bound) - mean_diff) / se    # H0: diff >= +bound
            p_high = float(stats.t.sf(t_high, df))
            tost_p = max(p_low, p_high)
            tost_equiv = (ci_lo > -equiv_bound) and (ci_hi < equiv_bound)

            # paired-fold bootstrap (fixed seed per (dataset,config) for determinism)
            row_seed = BOOTSTRAP_SEED + abs(hash((ds, cfg))) % 100000
            b_lo, b_hi = bootstrap_ci(diffs, row_seed, N_BOOTSTRAP)
            boot_within = (b_lo > -equiv_bound) and (b_hi < equiv_bound)

            rec.update({
                "mean_paired_diff": f"{mean_diff:+.4f}",
                "t_stat": f"{t_stat:+.3f}", "p_two_sided": f"{p_two:.4f}",
                "df": df,
                "t_ci90_lo": f"{ci_lo:+.4f}", "t_ci90_hi": f"{ci_hi:+.4f}",
                "boot_ci95_lo": f"{b_lo:+.4f}", "boot_ci95_hi": f"{b_hi:+.4f}",
                "tost_p": f"{tost_p:.4f}",
                "tost_t_equivalent": "yes" if tost_equiv else "no",
                "boot_within_bound": "yes" if boot_within else "no",
            })
            tost_rows.append(rec)
            pvalue_family.append((ds, cfg, p_two))

    # ---- write recomputed TOST/bootstrap CSV ----
    fields = ["dataset", "config", "n_folds", "mean_paired_diff", "t_stat",
              "p_two_sided", "df", "t_ci90_lo", "t_ci90_hi",
              "boot_ci95_lo", "boot_ci95_hi", "equiv_bound_pm5pct",
              "tost_p", "tost_t_equivalent", "boot_within_bound", "note"]
    with open(OUT_TOST, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in tost_rows:
            w.writerow(r)
    out(f"\nwrote {os.path.relpath(OUT_TOST, PROJECT)}  ({len(tost_rows)} rows)")

    # ---- Benjamini-Hochberg + Bonferroni over the finite-p family ----
    m = len(pvalue_family)
    out("\n" + "-" * 78)
    out(f"MULTIPLE-COMPARISON CORRECTION  (family size m = {m}, q = {ALPHA})")
    out("  Degenerate Lipophilicity combined_pruned EXCLUDED (no finite p-value).")
    out("-" * 78)

    ordered = sorted(pvalue_family, key=lambda t: t[2])   # ascending p
    # BH step-up: largest k with p_(k) <= k/m*q ; reject all ranks <= k
    bh_k = 0
    for i, (_, _, p) in enumerate(ordered, start=1):
        if p <= i / m * ALPHA:
            bh_k = i
    bonf_thresh = ALPHA / m

    bh_rows = []
    out(f"  {'rank':>4} {'dataset':>14} {'config':>16} {'raw p':>9} "
        f"{'BH thr':>9} {'BH rej':>7} {'Bonf rej':>9}")
    for i, (ds, cfg, p) in enumerate(ordered, start=1):
        bh_thr = i / m * ALPHA
        bh_reject = i <= bh_k
        bonf_reject = p <= bonf_thresh
        bh_rows.append({
            "dataset": ds, "config": cfg, "raw_p_two_sided": f"{p:.4f}",
            "rank": i, "family_size_m": m, "q": ALPHA,
            "bh_critical_value": f"{bh_thr:.4f}", "bh_reject": bh_reject,
            "bonferroni_critical_value": f"{bonf_thresh:.4f}",
            "bonferroni_reject": bonf_reject,
        })
        out(f"  {i:>4} {ds:>14} {cfg:>16} {p:>9.4f} "
            f"{bh_thr:>9.4f} {str(bh_reject):>7} {str(bonf_reject):>9}")

    n_bh = sum(r["bh_reject"] for r in bh_rows)
    n_bonf = sum(r["bonferroni_reject"] for r in bh_rows)
    out(f"\n  BH rejects {n_bh}/{m} at q={ALPHA};  Bonferroni rejects {n_bonf}/{m}.")

    with open(OUT_BH, "w", newline="") as f:
        fields_bh = ["dataset", "config", "raw_p_two_sided", "rank",
                     "family_size_m", "q", "bh_critical_value", "bh_reject",
                     "bonferroni_critical_value", "bonferroni_reject"]
        w = csv.DictWriter(f, fieldnames=fields_bh)
        w.writeheader()
        for r in bh_rows:
            w.writerow(r)
    out(f"  wrote {os.path.relpath(OUT_BH, PROJECT)}")

    # ---- VERIFICATION against existing artifacts & manuscript ----
    out("\n" + "=" * 78)
    out("VERIFICATION")
    out("=" * 78)

    # (1) recomputed p vs manuscript Table 11 rounded p
    out("\n[1] paired-t p-values: recomputed vs manuscript Table 11 (rounded)")
    pmap = {(ds, cfg): p for ds, cfg, p in pvalue_family}
    for key, mp in MANUSCRIPT_P.items():
        rp = pmap.get(key)
        if rp is None:
            out(f"    {key[0]:>14} {key[1]:>16}: recomputed=NA  manuscript={mp}")
            continue
        ok = abs(round(rp, 3) - mp) <= 0.001
        out(f"    {key[0]:>14} {key[1]:>16}: recomputed={rp:.4f}  "
            f"manuscript={mp:.3f}  {'OK' if ok else 'DIFF'}")

    # (2) recomputed TOST/CI/bound vs existing artifact CSV
    out("\n[2] recomputed vs existing paper2_table11_tost_bootstrap.csv")
    if os.path.exists(EXISTING_TOST_CSV):
        existing = {}
        with open(EXISTING_TOST_CSV, newline="") as f:
            for row in csv.DictReader(f):
                existing[(row["dataset"], row["config"])] = row

        def pf(s):
            if s is None:
                return None
            s = s.strip().replace("±", "").replace("+/-", "")
            try:
                return float(s)
            except ValueError:
                return None

        rec_by_key = {(r["dataset"], r["config"]): r for r in tost_rows}
        out(f"    {'dataset/config':>30} {'field':>12} {'existing':>10} {'recomp':>10} {'match':>6}")
        for (ds, cfg), er in existing.items():
            rr = rec_by_key.get((ds, cfg))
            if rr is None:
                continue
            checks = [
                ("mean_diff", pf(er.get("mean_paired_diff")), pf(rr["mean_paired_diff"]), 0.0005),
                ("t_ci90_lo", pf(er.get("t_ci90_lo")), pf(rr["t_ci90_lo"]), 0.0005),
                ("t_ci90_hi", pf(er.get("t_ci90_hi")), pf(rr["t_ci90_hi"]), 0.0005),
                ("equiv_bnd", pf(er.get("equiv_bound_pm5pct")), pf(rr["equiv_bound_pm5pct"]), 0.0005),
                ("boot_lo", pf(er.get("boot_ci95_lo")), pf(rr["boot_ci95_lo"]), 0.005),
                ("boot_hi", pf(er.get("boot_ci95_hi")), pf(rr["boot_ci95_hi"]), 0.005),
            ]
            for name, ev, rv, tol in checks:
                if ev is None and rv is None:
                    continue
                if ev is None or rv is None:
                    flag = "n/a"
                else:
                    flag = "OK" if abs(ev - rv) <= tol else "DIFF"
                out(f"    {ds + '/' + cfg:>30} {name:>12} {fmt(ev):>10} {fmt(rv):>10} {flag:>6}")
            # TOST verdict string match
            ev_t = (er.get("tost_t_equivalent") or "").strip()
            rv_t = rr["tost_t_equivalent"]
            out(f"    {ds + '/' + cfg:>30} {'TOST':>12} {ev_t:>10} {rv_t:>10} "
                f"{'OK' if ev_t == rv_t else 'DIFF':>6}")
        out("\n    NOTE: bootstrap CIs use a NEW fixed seed (the original artifact had no")
        out("    committed seed); agreement within ~0.005 Monte-Carlo tolerance is expected.")
    else:
        out("    existing artifact not found; skipping comparison.")

    # (3) BH conclusion vs corrected manuscript claim
    out("\n[3] BH conclusion vs corrected manuscript claim")
    out(f"    manuscript claim: 'none of the nominally-significant entries survives at q=0.05'")
    out(f"    recomputed      : BH rejects {n_bh}/{m} at q=0.05  ->  "
        f"{'MATCH' if n_bh == 0 else 'MISMATCH'}")
    nominal = [k for k in pmap if pmap[k] <= ALPHA]
    out(f"    nominally significant (raw p<=0.05): {len(nominal)} entries "
        f"{[f'{k[0]}/{k[1]}' for k in nominal]}")

    with open(OUT_REPORT, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"\n[report written to {os.path.relpath(OUT_REPORT, PROJECT)}]")


if __name__ == "__main__":
    main()
