"""Derive markdown tables and the win map from ``results/summary.json``."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

SPEC = "spectroscopy (filter function)"
HALV = "direct successive halving"


def write_tables(results_dir: Path) -> None:
    """Write ``win_map.csv`` and ``tables.md`` next to the summary."""
    s = json.loads((results_dir / "summary.json").read_text())
    t, budgets, epss = s["table"], s["budgets"], s["pulse_errors"]
    e_focus = 0.05 if 0.05 in epss else epss[0]
    b_focus = min(budgets, key=lambda b: abs(b - 1000))

    def cell(m: str, e: float, sh: float) -> str:
        d = t[m][str(e)][str(float(sh))]
        return f"{d['regret_mean']:.3f} ± {d['regret_std']:.3f}"

    wins = [{"eps": e, "shots": b, "spectroscopy_regret": t[SPEC][str(e)][str(float(b))]["regret_mean"],
             "halving_regret": t[HALV][str(e)][str(float(b))]["regret_mean"],
             "spectroscopy_wins": t[SPEC][str(e)][str(float(b))]["regret_mean"] < t[HALV][str(e)][str(float(b))]["regret_mean"]}
            for e in epss for b in budgets]
    pd.DataFrame(wins).to_csv(results_dir / "win_map.csv", index=False)

    ms = ["direct uniform", HALV, SPEC, "spectroscopy, no Larmor lines", "spectroscopy, no pulse calibration"]
    lines = [f"Table A: regret versus shot budget at flip-angle error {e_focus:g} (mean ± std over bath draws)", "",
             "| Selection rule | " + " | ".join(f"{b:,}" for b in budgets) + " |", "|---|" + "---|" * len(budgets)]
    lines += [f"| {m} | " + " | ".join(cell(m, e_focus, b) for b in budgets) + " |" for m in ms]
    lines += ["", f"Table B: regret at {b_focus:,} shots versus flip-angle error", "",
              "| Selection rule | " + " | ".join(f"ε = {e:g}" for e in epss) + " |", "|---|" + "---|" * len(epss)]
    for m in ["fixed Hahn echo", "fixed CPMG-8", "fixed XY4-8", "best single sequence"]:
        lines.append(f"| {m} | " + " | ".join(cell(m, e, 0) for e in epss) + " |")
    lines += [f"| {m} | " + " | ".join(cell(m, e, b_focus) for e in epss) + " |" for m in ms]
    lines += ["", f"Table C: model-error floor (unlimited shots) versus model-free selection at {budgets[-1]:,} shots", "",
              "| Selection rule | " + " | ".join(f"ε = {e:g}" for e in epss) + " |", "|---|" + "---|" * len(epss)]
    for m in ["spectroscopy (infinite shots)", "spectroscopy, no Larmor lines (infinite shots)"]:
        lines.append(f"| {m} | " + " | ".join(cell(m, e, float("inf")) for e in epss) + " |")
    lines.append(f"| {HALV}, {budgets[-1]:,} shots | " + " | ".join(cell(HALV, e, budgets[-1]) for e in epss) + " |")
    lines += ["", "Table D: does the surrogate have lower mean regret than successive halving? (yes / no)", "",
              "| Budget | " + " | ".join(f"ε = {e:g}" for e in epss) + " |", "|---|" + "---|" * len(epss)]
    for b in budgets:
        lines.append(f"| {b:,} | " + " | ".join(
            "yes" if t[SPEC][str(e)][str(float(b))]["regret_mean"] < t[HALV][str(e)][str(float(b))]["regret_mean"] else "no"
            for e in epss) + " |")
    (results_dir / "tables.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
