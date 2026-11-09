"""Reproduce every result and figure: ``python experiments/run_all.py`` (add ``--quick`` for a smoke run)."""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from ddselect import plotting as P  # noqa: E402
from ddselect.bath import BathConfig  # noqa: E402
from ddselect.report import write_tables  # noqa: E402
from ddselect.sequences import Sequence  # noqa: E402
from ddselect.simulate import CentralSpinBath, pick_device  # noqa: E402
from ddselect.study import (default_grid, evaluate_condition, run_condition,  # noqa: E402
                            select_spectroscopy)

torch.set_num_threads(6)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="tiny smoke run (< 1 minute)")
    args = ap.parse_args()
    quick = args.quick
    out_res, out_fig = ROOT / "results", ROOT / "figures"
    if quick:
        out_res, out_fig = ROOT / "results_quick", ROOT / "figures_quick"
    out_res.mkdir(exist_ok=True); out_fig.mkdir(exist_ok=True)

    t0 = time.time()
    seeds = [10, 11, 12, 13, 14] if not quick else [10]  # disjoint from the pilot seeds 0-2
    epss = [0.0, 0.02, 0.05, 0.1] if not quick else [0.05]
    budgets = [100, 300, 1000, 3000, 10_000, 100_000] if not quick else [300, 3000]
    reps = 30 if not quick else 3
    cfg = BathConfig(n_nuclei=8 if not quick else 4)
    grid = default_grid(quick)
    rows, truth, pred_rows = [], [], []
    for seed in seeds:
        for eps in epss:
            tab = run_condition(cfg, seed, eps, grid)
            rows += evaluate_condition(tab, grid, budgets, reps, seed)
            for i, c in enumerate(grid.cands):
                for j, t in enumerate(grid.storage_times):
                    truth.append({"seed": seed, "eps": eps, "sequence": c.name, "family": c.family,
                                  "storage_time_us": float(t), "fidelity": float(tab.fid[i, j])})
            if eps == 0.0 or quick:
                _, pred = select_spectroscopy(tab, grid, None, np.random.default_rng(0))
                for i, c in enumerate(grid.cands):
                    for j, t in enumerate(grid.storage_times):
                        pred_rows.append({"seed": seed, "family": c.family, "i": i, "j": j,
                                          "exact": float(tab.fid[i, j]), "pred": float(pred[i, j])})
            print(f"seed {seed} eps {eps} done ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(out_res / "regret_all.csv", index=False)
    pd.DataFrame(truth).to_csv(out_res / "truth_fidelity.csv", index=False)

    # --- surrogate calibration -----------------------------------------------------------
    pr = pd.DataFrame(pred_rows)
    pr.to_csv(out_res / "surrogate_predictions.csv", index=False)
    sp = []
    for (_, _), g in pr.groupby(["seed", "j"]):
        if g.exact.std() > 1e-6:
            sp.append(spearmanr(g.exact, g.pred)[0])
    cal = {"rmse": float(np.sqrt(((pr.pred - pr.exact) ** 2).mean())),
           "spearman": float(np.nanmean(sp)),
           "rmse_by_family": {f: float(np.sqrt(((g.pred - g.exact) ** 2).mean())) for f, g in pr.groupby("family")}}
    P.fig_pred_vs_exact(pr.pred.values, pr.exact.values, pr.family.values, out_fig / "fig4_surrogate_calibration.png", cal)

    # --- figures 1-3 ---------------------------------------------------------------------
    eps_focus = 0.05 if 0.05 in epss else epss[0]
    shots_focus = 1000 if not quick else 300
    P.fig_regret_bars(df, shots_focus, out_fig / "fig1_regret_by_method.png")
    slopes = P.fig_shot_scaling(df, eps_focus, out_fig / "fig2_regret_vs_shots.png")
    P.fig_pareto(df, eps_focus, out_fig / "fig3_cost_quality_pareto.png")

    # --- bath-size scaling: simulation cost and selector regret --------------------------
    dev = pick_device()
    size_rows, size_regret = [], []
    sizes = [4, 5, 6, 7, 8, 9, 10] if not quick else [3, 4]
    for n in sizes:
        for seed in seeds[:3]:
            bath = CentralSpinBath(BathConfig(n_nuclei=n), seed, dev)
            seq = Sequence("cpmg", 8)
            bath.memory_fidelity(seq, 5.0, 0.0)  # warm-up
            if dev.type == "cuda":
                torch.cuda.synchronize()
            ts = time.time()
            for _ in range(3):
                bath.memory_fidelity(seq, 5.0, 0.0)
            if dev.type == "cuda":
                torch.cuda.synchronize()
            size_rows.append({"n_nuclei": n, "dim": 2**n, "seed": seed,
                              "seconds_per_sequence": (time.time() - ts) / 3})
    for n in ([4, 6] if not quick else [3]):
        for seed in seeds[:3]:
            tab = run_condition(BathConfig(n_nuclei=n), seed, eps_focus, grid)
            for r in evaluate_condition(tab, grid, [300, 3000] if not quick else [300], reps, seed):
                if r["method"] in ("direct successive halving", "spectroscopy (filter function)",
                                   "spectroscopy (infinite shots)", "best single sequence"):
                    size_regret.append({**r, "n_nuclei": n})
    sz = pd.DataFrame(size_rows)
    sz.to_csv(out_res / "size_scaling_time.csv", index=False)
    sr = pd.DataFrame(size_regret)
    main_sel = df[(df.eps == eps_focus) & df.method.isin(
        ["direct successive halving", "spectroscopy (filter function)", "spectroscopy (infinite shots)",
         "best single sequence"]) & df.shots.isin([0, 300, 3000, float("inf")])].copy()
    main_sel["n_nuclei"] = cfg.n_nuclei
    sr = pd.concat([sr, main_sel[sr.columns]]) if len(sr) else main_sel
    sr.to_csv(out_res / "size_scaling_regret.csv", index=False)
    size_slope = P.fig_size_scaling(sz, out_fig / "fig5_simulation_cost_scaling.png")

    # --- summary ------------------------------------------------------------------------
    def agg(method: str, eps: float, shots: float) -> dict:
        v = df[(df.method == method) & (df.eps == eps) & (df.shots == shots)].groupby("seed").regret.mean()
        f = df[(df.method == method) & (df.eps == eps) & (df.shots == shots)].groupby("seed").fidelity.mean()
        return {"regret_mean": float(v.mean()), "regret_std": float(v.std(ddof=1)) if len(v) > 1 else 0.0,
                "fidelity_mean": float(f.mean()), "fidelity_std": float(f.std(ddof=1)) if len(f) > 1 else 0.0}

    methods = sorted(df.method.unique())
    table = {m: {str(e): {str(s): agg(m, e, s) for s in sorted(df[df.method == m].shots.unique())} for e in epss}
             for m in methods}
    cross = {}
    for e in epss:
        cross[str(e)] = None
        for s in budgets:
            a = agg("spectroscopy (filter function)", e, s)["regret_mean"]
            b = agg("direct successive halving", e, s)["regret_mean"]
            if b <= a:
                cross[str(e)] = s
                break
    summary = {"seeds": seeds, "pulse_errors": epss, "budgets": budgets, "reps_per_budget": reps,
               "n_nuclei": cfg.n_nuclei, "n_candidates": len(grid.cands),
               "storage_times_us": [float(t) for t in grid.storage_times],
               "n_probes": len(grid.probes()), "loglog_slopes": slopes, "crossing_budget_halving_beats_spectroscopy": cross,
               "surrogate_calibration": cal, "simulation_cost_slope": size_slope,
               "device": str(dev), "table": table, "runtime_seconds": time.time() - t0}
    with open(out_res / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=1)
    write_tables(out_res)
    print(f"done in {time.time() - t0:.0f} s; results in {out_res}")


if __name__ == "__main__":
    main()
