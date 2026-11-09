"""Central-spin model of an electron spin in a GaAs-like dot with a mesoscopic nuclear bath.

Hamiltonian (rad/us, rotating frame of the electron, electron S_z = +-1/2)::

    H_pm = sum_k w_k Iz_k  +-  1/2 sum_k (A_k Iz_k + B_k Ix_k)  +  sum_{k<l} J_kl (Ix_k Ix_l + Iy_k Iy_l)

with nuclear Larmor frequencies ``w_k`` of 75As, 69Ga, 71Ga (spin-1/2 approximation),
parallel / transverse hyperfine couplings ``A_k``, ``B_k`` and flip-flop couplings ``J_kl``.
Each simulated nucleus stands for a group of nuclei, so couplings are effective values.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

GAMMA_MHZ_PER_T = {"As75": 7.31, "Ga69": 10.22, "Ga71": 12.98}
SPECIES_PROB = {"As75": 0.5, "Ga69": 0.3, "Ga71": 0.2}
TWO_PI = 2.0 * np.pi

SX = np.array([[0, 1], [1, 0]], dtype=complex) / 2
SY = np.array([[0, -1j], [1j, 0]], dtype=complex) / 2
SZ = np.array([[1, 0], [0, -1]], dtype=complex) / 2


@dataclass(frozen=True)
class BathConfig:
    """Parameters of the nuclear bath (frequencies in rad/us)."""

    n_nuclei: int = 8
    field_tesla: float = 0.1
    hyperfine_scale: float = TWO_PI * 0.4
    transverse_range: tuple[float, float] = (0.3, 1.5)
    flipflop_scale: float = TWO_PI * 0.03

    def larmor(self) -> dict[str, float]:
        """Nuclear Larmor frequency of each species in rad/us."""
        return {k: TWO_PI * g * self.field_tesla for k, g in GAMMA_MHZ_PER_T.items()}


@dataclass(frozen=True)
class BathDraw:
    """One random realisation of the couplings."""

    species: tuple[str, ...]
    omega: np.ndarray
    a_par: np.ndarray
    b_perp: np.ndarray
    j_ff: np.ndarray


def draw_bath(cfg: BathConfig, seed: int) -> BathDraw:
    """Sample species, hyperfine and flip-flop couplings for a seed."""
    rng = np.random.default_rng(seed)
    names = list(SPECIES_PROB)
    species = tuple(rng.choice(names, size=cfg.n_nuclei, p=[SPECIES_PROB[s] for s in names]))
    lar = cfg.larmor()
    omega = np.array([lar[s] for s in species])
    a_par = cfg.hyperfine_scale * rng.uniform(0.2, 1.0, cfg.n_nuclei) * rng.choice([-1, 1], cfg.n_nuclei)
    lo, hi = cfg.transverse_range
    b_perp = np.abs(a_par) * rng.uniform(lo, hi, cfg.n_nuclei)
    j = cfg.flipflop_scale * rng.normal(size=(cfg.n_nuclei, cfg.n_nuclei))
    j = np.triu(j, 1)
    return BathDraw(species, omega, a_par, b_perp, j)


def _site_op(op: np.ndarray, k: int, n: int) -> np.ndarray:
    out = np.array([[1.0 + 0j]])
    for s in range(n):
        out = np.kron(out, op if s == k else np.eye(2))
    return out


def conditional_hamiltonians(draw: BathDraw) -> tuple[np.ndarray, np.ndarray]:
    """Return (H_up, H_down): nuclear Hamiltonians for electron S_z = +1/2 and -1/2."""
    n = len(draw.species)
    dim = 2**n
    ix = [_site_op(SX, k, n) for k in range(n)]
    iy = [_site_op(SY, k, n) for k in range(n)]
    iz = [_site_op(SZ, k, n) for k in range(n)]
    base = np.zeros((dim, dim), dtype=complex)
    hf = np.zeros((dim, dim), dtype=complex)
    for k in range(n):
        base += draw.omega[k] * iz[k]
        hf += 0.5 * (draw.a_par[k] * iz[k] + draw.b_perp[k] * ix[k])
        for m in range(k + 1, n):
            if draw.j_ff[k, m] != 0.0:
                base += draw.j_ff[k, m] * (ix[k] @ ix[m] + iy[k] @ iy[m])
    return base + hf, base - hf
