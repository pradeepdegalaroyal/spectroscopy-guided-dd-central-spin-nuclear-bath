"""Fast tests: sequence definitions, exact simulation limits, surrogate and selectors."""
import numpy as np
import pytest
import torch

from ddselect.bath import BathConfig, draw_bath
from ddselect.filterfunction import FilterFunctionModel, make_basis, switching_ft
from ddselect.sequences import Sequence, candidate_set
from ddselect.simulate import CentralSpinBath
from ddselect.study import default_grid, evaluate_condition, regret, run_condition, select_spectroscopy

CPU = torch.device("cpu")


def _bath(n: int = 3, seed: int = 0) -> CentralSpinBath:
    return CentralSpinBath(BathConfig(n_nuclei=n), seed, CPU)


def test_pulse_positions_are_sorted_and_inside():
    for s in candidate_set(16):
        f = s.fractions()
        assert len(f) == s.n and np.all(np.diff(f) > 0) and (len(f) == 0 or (f[0] > 0 and f[-1] < 1))


def test_udd2_equals_cpmg2():
    assert np.allclose(Sequence("udd", 2).fractions(), Sequence("cpmg", 2).fractions())


def test_zero_time_gives_unit_fidelity():
    b = _bath()
    assert b.memory_fidelity(Sequence("free", 0), 1e-9) == pytest.approx(1.0, abs=1e-5)
    assert b.memory_fidelity(Sequence("cpmg", 4), 1e-9) == pytest.approx(1.0, abs=1e-4)


def test_uncoupled_bath_is_lossless():
    cfg = BathConfig(n_nuclei=3, hyperfine_scale=0.0, flipflop_scale=0.0)
    b = CentralSpinBath(cfg, 1, CPU)
    assert b.memory_fidelity(Sequence("cpmg", 4), 7.0, 0.0) == pytest.approx(1.0, abs=1e-4)


def test_free_evolution_bounded_by_dephasing_floor():
    f = _bath().memory_fidelity(Sequence("free", 0), 50.0, 0.0)
    assert 2 / 3 - 0.02 <= f <= 1.0 + 1e-4


def test_pulse_error_hurts_cpmg_more_than_xy4():
    b = _bath()
    t = 0.05
    assert b.memory_fidelity(Sequence("xy4", 16), t, 0.1) > b.memory_fidelity(Sequence("cpmg", 16), t, 0.1)


def test_exact_ideal_dephasing_for_ising_bath():
    """With purely parallel couplings and no Zeeman term, the echo coherence is exactly one."""
    cfg = BathConfig(n_nuclei=3, transverse_range=(0.0, 0.0), flipflop_scale=0.0)
    bb = CentralSpinBath(cfg, 2, CPU)
    assert bb.coherence_x(Sequence("cpmg", 2), 5.0, 0.0) == pytest.approx(1.0, abs=1e-4)
    assert bb.coherence_x(Sequence("free", 0), 5.0, 0.0) < 0.999


def test_switching_ft_matches_free_formula():
    w = np.array([0.7, 2.0])
    y = switching_ft(Sequence("free", 0), 3.0, w)
    assert np.allclose(y, (np.exp(1j * w * 3.0) - 1) / (1j * w))


def test_filter_function_fit_recovers_synthetic_spectrum():
    basis = make_basis(BathConfig().larmor())
    model = FilterFunctionModel(basis)
    model.weights = np.zeros(len(basis.labels)); model.weights[2] = 0.05
    probes = [(Sequence("cpmg", n), t) for n in (2, 4, 8) for t in (2.0, 6.0, 20.0)]
    c = np.array([model.predict_coherence(s, t) for s, t in probes])
    fit = FilterFunctionModel(basis)
    fit.fit(probes, c, 10**9)
    assert fit.predict_coherence(Sequence("cpmg", 16), 10.0) == pytest.approx(
        model.predict_coherence(Sequence("cpmg", 16), 10.0), abs=0.03)


def test_bath_draw_reproducible():
    a, b = draw_bath(BathConfig(), 5), draw_bath(BathConfig(), 5)
    assert a.species == b.species and np.allclose(a.a_par, b.a_par)


def test_selectors_and_regret_quick():
    grid = default_grid(quick=True)
    tab = run_condition(BathConfig(n_nuclei=3), 0, 0.02, grid)
    rows = evaluate_condition(tab, grid, [300], 2, 0)
    methods = {r["method"] for r in rows}
    assert {"oracle", "direct successive halving", "spectroscopy (filter function)"} <= methods
    assert all(r["regret"] >= -1e-9 for r in rows)
    ch, _ = select_spectroscopy(tab, grid, None, np.random.default_rng(0))
    assert regret(tab, tab.fid.argmax(axis=0))[1] == pytest.approx(0.0)
    assert ch.shape == (len(grid.storage_times),)
