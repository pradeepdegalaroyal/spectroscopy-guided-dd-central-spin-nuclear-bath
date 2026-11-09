"""Figures: one accent colour plus greys; hatched bars, Pareto scatter, log-log fits."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ACCENT = "#c2410c"
GREYS = ["#222222", "#555555", "#888888", "#bbbbbb"]
HATCHES = ["", "///", "xxx", "..."]
DPI = 200


def _style(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=8)


def fig_regret_bars(df: pd.DataFrame, shots: int, path: Path) -> None:
    """Grouped hatched bars: regret per method, one bar per pulse error."""
    methods = ["fixed Hahn echo", "fixed CPMG-8", "fixed XY4-8", "best single sequence",
               "direct successive halving", "spectroscopy, no Larmor lines",
               "spectroscopy, no pulse calibration", "spectroscopy (filter function)"]
    sub = df[(df.shots == shots) | (df.shots == 0)]
    epss = sorted(sub.eps.unique())
    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    w = 0.8 / len(epss)
    for k, e in enumerate(epss):
        means, stds = [], []
        for m in methods:
            v = sub[(sub.method == m) & (sub.eps == e)].groupby("seed").regret.mean()
            means.append(v.mean()); stds.append(v.std(ddof=1) if len(v) > 1 else 0.0)
        colors = [ACCENT if m.startswith("spectroscopy") else GREYS[1] for m in methods]
        ax.bar(np.arange(len(methods)) + (k - (len(epss) - 1) / 2) * w, means, w, yerr=stds, capsize=2,
               color=colors, edgecolor="black", hatch=HATCHES[k % 4], linewidth=0.6,
               alpha=1.0 - 0.12 * k, label=f"pulse error {e:g}")
    ax.set_xticks(range(len(methods)))
    short = {"fixed Hahn echo": "fixed\nHahn echo", "fixed CPMG-8": "fixed\nCPMG-8", "fixed XY4-8": "fixed\nXY4-8",
             "best single sequence": "best single\nsequence", "direct successive halving": "direct\nsuccessive\nhalving",
             "spectroscopy, no Larmor lines": "spectroscopy\nno Larmor\nlines",
             "spectroscopy, no pulse calibration": "spectroscopy\nno pulse\ncalibration",
             "spectroscopy (filter function)": "spectroscopy\nfilter\nfunction"}
    ax.set_xticklabels([short[m] for m in methods], fontsize=7.5)
    ax.set_ylabel("mean regret in memory fidelity (dimensionless)", fontsize=9)
    ax.set_title(f"Regret of each selection rule at {shots} total shots (mean +- std over bath draws)", fontsize=10)
    ax.legend(fontsize=8, frameon=False, title="hatch = flip-angle error", title_fontsize=8)
    _style(ax)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)


def fig_shot_scaling(df: pd.DataFrame, eps: float, path: Path) -> dict[str, float]:
    """Log-log regret versus shot budget with fitted slopes; returns the fitted slopes."""
    sub = df[(df.eps == eps) & np.isfinite(df.shots) & (df.shots > 0)]
    spec = {"direct successive halving": (GREYS[0], "o"), "direct uniform": (GREYS[2], "s"),
            "spectroscopy (filter function)": (ACCENT, "^")}
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    slopes: dict[str, float] = {}
    for m, (c, mk) in spec.items():
        g = sub[sub.method == m].groupby("shots").regret.agg(["mean", "std"]).reset_index()
        ax.errorbar(g.shots, g["mean"].clip(lower=1e-4), yerr=g["std"], color=c, marker=mk, capsize=2, lw=1.2,
                    label=m)
        if m.startswith("direct"):
            fit = g[(g.shots >= 300) & (g.shots <= 1e4) & (g["mean"] > 0)]
            a, b = np.polyfit(np.log10(fit.shots), np.log10(fit["mean"]), 1)
            slopes[m] = float(a)
            xs = np.array([300, 1e4])
            ax.plot(xs, 10 ** (b + a * np.log10(xs)), "--", color=c, lw=0.9)
            ax.text(3e3, 10 ** (b + a * np.log10(3e3)) * (2.4 if m.endswith("uniform") else 0.35),
                    f"slope {a:.2f}", color=c, fontsize=8)
    floor = df[(df.eps == eps) & (df.method == "spectroscopy (infinite shots)")].groupby("seed").regret.mean().mean()
    ax.axhline(floor, color=ACCENT, ls=":", lw=1)
    ax.text(110, floor * 1.12, "surrogate floor (infinite shots)", color=ACCENT, fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("total measurement shots per decision table", fontsize=9)
    ax.set_ylabel("mean regret in memory fidelity", fontsize=9)
    ax.set_title(f"Regret versus shot budget, flip-angle error {eps:g}", fontsize=10)
    ax.legend(fontsize=8, frameon=False, loc="lower left")
    _style(ax)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)
    return slopes


def pareto_mask(cost: np.ndarray, qual: np.ndarray) -> np.ndarray:
    """True for points not dominated in (lower cost, lower regret)."""
    keep = np.ones(len(cost), dtype=bool)
    for i in range(len(cost)):
        keep[i] = not np.any((cost <= cost[i]) & (qual <= qual[i]) & ((cost < cost[i]) | (qual < qual[i])))
    return keep


def fig_pareto(df: pd.DataFrame, eps: float, path: Path) -> None:
    """Cost-quality scatter of (method, budget) points with the Pareto front highlighted."""
    sub = df[(df.eps == eps) & np.isfinite(df.shots) & (df.shots > 0)]
    g = sub.groupby(["method", "shots"]).regret.mean().reset_index()
    g = g[g.method.isin(["direct successive halving", "direct uniform", "spectroscopy (filter function)",
                         "spectroscopy, no Larmor lines", "spectroscopy, no pulse calibration"])]
    mask = pareto_mask(g.shots.values.astype(float), g.regret.values)
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for m, mk in [("direct successive halving", "o"), ("direct uniform", "s"), ("spectroscopy (filter function)", "^"),
                  ("spectroscopy, no Larmor lines", "v"), ("spectroscopy, no pulse calibration", "D")]:
        h = g[g.method == m]
        ax.scatter(h.shots, h.regret.clip(lower=1e-4), marker=mk, s=34,
                   color=ACCENT if m.startswith("spectroscopy") else GREYS[1], alpha=0.75, label=m)
    front = g[mask].sort_values("shots")
    ax.step(front.shots, front.regret.clip(lower=1e-4), where="post", color="black", lw=1.4)
    ax.scatter(front.shots, front.regret.clip(lower=1e-4), s=110, facecolors="none", edgecolors="black", lw=1.4,
               label="Pareto front")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("measurement cost: total shots", fontsize=9)
    ax.set_ylabel("quality: mean regret in fidelity (lower is better)", fontsize=9)
    ax.set_title(f"Cost versus quality of selection rules, flip-angle error {eps:g}", fontsize=10)
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    _style(ax)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)


def fig_pred_vs_exact(pred: np.ndarray, exact: np.ndarray, fam: np.ndarray, path: Path, stats: dict) -> None:
    """Scatter of surrogate-predicted versus exact fidelity, coloured by sequence family."""
    fig, ax = plt.subplots(figsize=(5.2, 4.8))
    for f, c, mk in [("free", ACCENT, "*"), ("cpmg", GREYS[0], "o"), ("udd", GREYS[2], "s"), ("xy4", GREYS[3], "^")]:
        m = fam == f
        ax.scatter(exact[m], pred[m], s=14 if f != "free" else 40, color=c, marker=mk, alpha=0.7,
                   edgecolors="black" if f == "xy4" else "none", linewidths=0.3, label=f)
    ax.plot([0.6, 1.02], [0.6, 1.02], color="black", lw=0.8)
    ax.set_xlabel("exact memory fidelity (quantum nuclear bath)", fontsize=9)
    ax.set_ylabel("filter-function surrogate prediction", fontsize=9)
    ax.set_title("Surrogate calibration at infinite shots, no pulse error", fontsize=10)
    ax.text(0.97, 0.04, f"RMSE {stats['rmse']:.3f}\nmean Spearman {stats['spearman']:.2f}", fontsize=8,
            transform=ax.transAxes, ha="right", va="bottom")
    ax.legend(fontsize=8, frameon=False, loc="upper left", title="sequence family", title_fontsize=8)
    _style(ax)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)


def fig_size_scaling(sz: pd.DataFrame, path: Path) -> float:
    """Log-log exact-simulation time per sequence versus Hilbert-space dimension."""
    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    g = sz.groupby("dim").seconds_per_sequence.agg(["mean", "std"]).reset_index()
    ax.errorbar(g.dim, g["mean"], yerr=g["std"], color=ACCENT, marker="o", capsize=2, label="exact simulation")
    big = g[g.dim >= 256]  # small sizes are launch-overhead limited
    if len(big) < 2:
        big = g
    a, b = np.polyfit(np.log10(big.dim), np.log10(big["mean"]), 1)
    xs = np.array([big.dim.min(), big.dim.max()])
    ax.plot(xs, 10 ** (b + a * np.log10(xs)), "--", color="black", lw=0.9)
    ax.text(0.03, 0.9, f"fitted slope {a:.2f} (dimension 256 and above)", fontsize=8, transform=ax.transAxes)
    ax.set_xscale("log", base=2); ax.set_yscale("log")
    ax.set_xticks(list(g.dim)); ax.set_xticklabels([str(int(d)) for d in g.dim], fontsize=8)
    ax.minorticks_off()
    ax.set_xlabel("nuclear Hilbert-space dimension 2^N", fontsize=9)
    ax.set_ylabel("wall-clock time per six-state sequence (s)", fontsize=9)
    ax.set_title("Cost of exact central-spin simulation (GPU if available)", fontsize=10)
    _style(ax)
    fig.tight_layout(); fig.savefig(path, dpi=DPI); plt.close(fig)
    return float(a)
