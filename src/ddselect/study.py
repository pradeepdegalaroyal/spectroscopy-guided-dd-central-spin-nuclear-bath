"""Study logic: exact truth tables, simulated experiments, selection rules and regret."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .bath import BathConfig
from .filterfunction import FilterFunctionModel, dephasing_fidelity, make_basis
from .sequences import Sequence, candidate_set, standard_baselines
from .simulate import CentralSpinBath

N_CAL = 16  # pulses in the pulse-error calibration trains
T_CAL = 0.05  # us, short enough that bath dephasing is negligible


@dataclass
class StudyGrid:
    """Candidate library, storage times and probe grid shared by all selectors."""

    cands: list[Sequence]
    storage_times: np.ndarray
    probe_pulses: tuple[int, ...] = (2, 4, 8, 16, 32)
    probe_times: np.ndarray = field(default_factory=lambda: np.geomspace(1.5, 64.0, 10))

    def probes(self) -> list[tuple[Sequence, float]]:
        """CPMG probe sequences (pulse number x total time) used for spectroscopy."""
        return [(Sequence("cpmg", n), float(t)) for n in self.probe_pulses for t in self.probe_times]


def default_grid(quick: bool = False) -> StudyGrid:
    """Full study grid, or a tiny one for smoke tests."""
    if quick:
        return StudyGrid(candidate_set(8), np.array([4.0, 16.0]), (2, 4, 8), np.geomspace(2.0, 32.0, 4))
    return StudyGrid(candidate_set(32), np.array([3.0, 4.5, 7.0, 10.0, 15.0, 22.0, 33.0, 50.0]))


@dataclass
class Tables:
    """Exact quantities for one (bath, pulse error) condition: everything a 'lab' can measure."""

    seed: int
    eps: float
    fid: np.ndarray  # (n_cand, n_T) true memory fidelity
    probe_c: np.ndarray  # (n_probe,) true sigma_x coherence of the CPMG probes
    cal_f: dict[str, float]  # six-state fidelity of pulse-error calibration trains
    larmor: dict[str, float]


def build_tables(bath: CentralSpinBath, eps: float, grid: StudyGrid) -> Tables:
    """Run the exact simulations that define truth and probe outcomes."""
    fid = np.array([[bath.memory_fidelity(c, float(t), eps) for t in grid.storage_times] for c in grid.cands])
    probe_c = np.array([bath.coherence_x(s, t, eps) for s, t in grid.probes()])
    cal = {
        "x": bath.memory_fidelity(Sequence("cpmg", N_CAL), T_CAL, eps),
        "xy4": bath.memory_fidelity(Sequence("xy4", N_CAL), T_CAL, eps),
    }
    return Tables(bath.seed, eps, fid, probe_c, cal, bath.cfg.larmor())


def run_condition(cfg: BathConfig, seed: int, eps: float, grid: StudyGrid) -> Tables:
    """Build the exact tables for one bath draw and pulse error."""
    return build_tables(CentralSpinBath(cfg, seed), eps, grid)


def _binom(rng: np.random.Generator, shots: int, p: float) -> float:
    return rng.binomial(shots, float(np.clip(p, 0.0, 1.0))) / shots


# --- selection rules --------------------------------------------------------------------
def select_uniform(tab: Tables, shots: int, rng: np.random.Generator) -> np.ndarray:
    """Model-free: equal shots on every (candidate, storage time) cell, pick best per time."""
    n_c, n_t = tab.fid.shape
    per = max(shots // (n_c * n_t), 1)
    est = np.array([[_binom(rng, per, tab.fid[i, j]) for j in range(n_t)] for i in range(n_c)])
    return est.argmax(axis=0)


def select_halving(tab: Tables, shots: int, rng: np.random.Generator) -> np.ndarray:
    """Model-free successive halving, run independently for each storage time."""
    n_c, n_t = tab.fid.shape
    rounds = int(np.ceil(np.log2(n_c)))
    out = np.zeros(n_t, dtype=int)
    for j in range(n_t):
        alive = np.arange(n_c)
        succ, tot = np.zeros(n_c), np.zeros(n_c)
        budget = shots / n_t
        for _ in range(rounds):
            if len(alive) == 1:
                break
            per = max(int(budget / rounds / len(alive)), 1)
            for i in alive:
                succ[i] += rng.binomial(per, float(np.clip(tab.fid[i, j], 0, 1)))
                tot[i] += per
            alive = alive[np.argsort(-succ[alive] / tot[alive])[: max(len(alive) // 2, 1)]]
        out[j] = alive[0]
    return out


def select_spectroscopy(
    tab: Tables,
    grid: StudyGrid,
    shots: int | None,
    rng: np.random.Generator,
    use_lines: bool = True,
    use_calibration: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Fit a filter-function surrogate to CPMG probes and pick the predicted-best candidate.

    ``shots=None`` means infinite statistics (exact probe values): the model-error floor.
    """
    probes = grid.probes()
    n_meas = len(probes) + (2 if use_calibration else 0)
    per = None if shots is None else max(shots // n_meas, 1)
    if per is None:
        c_hat = tab.probe_c.copy()
    else:
        c_hat = np.array([2 * _binom(rng, per, (1 + c) / 2) - 1 for c in tab.probe_c])
    model = FilterFunctionModel(make_basis(tab.larmor, use_lines))
    model.fit(probes, c_hat, per if per is not None else 10**9)
    g = {"x": 1.0, "xy4": 1.0}
    if use_calibration:
        for k, f in tab.cal_f.items():
            f_hat = f if per is None else _binom(rng, per, f)
            g[k] = float(np.clip(f_hat, 0.05, 1.0)) ** (1.0 / N_CAL)
    n_c, n_t = len(grid.cands), len(grid.storage_times)
    pred = np.zeros((n_c, n_t))
    for i, c in enumerate(grid.cands):
        pen = 1.0 if c.family == "free" else g[c.pulse_family()] ** c.n
        for j, t in enumerate(grid.storage_times):
            pred[i, j] = dephasing_fidelity(model.predict_coherence(c, float(t))) * pen
    return pred.argmax(axis=0), pred


def fixed_choice(grid: StudyGrid, seq: Sequence) -> np.ndarray:
    """Choice vector that always uses ``seq``."""
    idx = [c.name for c in grid.cands].index(seq.name)
    return np.full(len(grid.storage_times), idx)


def regret(tab: Tables, choice: np.ndarray) -> tuple[float, float]:
    """(mean selected fidelity, mean regret to the per-time oracle) over storage times."""
    cols = np.arange(tab.fid.shape[1])
    sel = tab.fid[choice, cols]
    return float(sel.mean()), float((tab.fid.max(axis=0) - sel).mean())


def evaluate_condition(tab: Tables, grid: StudyGrid, budgets: list[int], reps: int, seed: int) -> list[dict]:
    """Average regret of every method over ``reps`` shot-noise realisations per budget."""
    rows: list[dict] = []
    rng = np.random.default_rng(1000 * seed + int(1e4 * tab.eps) + 7)
    names = [c.name for c in grid.cands]

    def add(method: str, shots: float, fs: list[float], rs: list[float]) -> None:
        rows.append({"seed": tab.seed, "eps": tab.eps, "shots": shots, "method": method,
                     "fidelity": float(np.mean(fs)), "regret": float(np.mean(rs))})

    n_t = tab.fid.shape[1]
    f, r = regret(tab, tab.fid.argmax(axis=0))
    add("oracle", 0, [f], [r])
    f, r = regret(tab, np.full(n_t, int(tab.fid.mean(axis=1).argmax())))
    add("best single sequence", 0, [f], [r])
    for k, s in standard_baselines().items():
        if s.name in names:
            f, r = regret(tab, fixed_choice(grid, s))
            add(k, 0, [f], [r])
    for label, kw in (("spectroscopy (infinite shots)", {}), ("spectroscopy, no Larmor lines (infinite shots)", {"use_lines": False})):
        f, r = regret(tab, select_spectroscopy(tab, grid, None, rng, **kw)[0])
        add(label, float("inf"), [f], [r])
    variants = {
        "spectroscopy (filter function)": {},
        "spectroscopy, no Larmor lines": {"use_lines": False},
        "spectroscopy, no pulse calibration": {"use_calibration": False},
    }
    for s in budgets:
        for name, kw in variants.items():
            fs, rs = zip(*[regret(tab, select_spectroscopy(tab, grid, s, rng, **kw)[0]) for _ in range(reps)])
            add(name, s, list(fs), list(rs))
        fs, rs = zip(*[regret(tab, select_uniform(tab, s, rng)) for _ in range(reps)])
        add("direct uniform", s, list(fs), list(rs))
        fs, rs = zip(*[regret(tab, select_halving(tab, s, rng)) for _ in range(reps)])
        add("direct successive halving", s, list(fs), list(rs))
    return rows
