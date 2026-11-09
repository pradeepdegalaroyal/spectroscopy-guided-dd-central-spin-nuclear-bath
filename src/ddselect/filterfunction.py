"""Classical Gaussian-noise filter-function surrogate fitted from DD coherence probes.

For classical dephasing noise with one-sided spectrum ``S(w)`` and ideal pi pulses the
electron coherence is ``exp(-chi)`` with ``chi = int S(w) |Y(w, T)|^2 dw`` and ``Y`` the
Fourier transform of the sign-switching function. ``S`` is expanded in a non-negative basis
of Lorentzian lines (broad continuum at w = 0 plus narrow lines at known nuclear Larmor
frequencies) and the weights are fitted by weighted non-negative least squares.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from scipy.optimize import nnls

from .sequences import Sequence

_OMEGA = np.geomspace(0.02, 600.0, 3000)
_DOMEGA = np.gradient(_OMEGA)


def switching_ft(seq: Sequence, total_time: float, omega: np.ndarray = _OMEGA) -> np.ndarray:
    """Fourier transform Y(w) of the +-1 switching function of ``seq`` over ``total_time``."""
    edges = np.concatenate([[0.0], seq.fractions() * total_time, [total_time]])
    signs = (-1.0) ** np.arange(len(edges) - 1)
    ph = np.exp(1j * np.outer(omega, edges))
    return ((ph[:, 1:] - ph[:, :-1]) * signs).sum(axis=1) / (1j * omega)


def lorentzian(center: float, width: float, omega: np.ndarray = _OMEGA) -> np.ndarray:
    """Normalised Lorentzian line of half-width ``width`` at ``center``."""
    return (width / np.pi) / ((omega - center) ** 2 + width**2)


@dataclass
class Basis:
    """Spectral basis functions sampled on the frequency grid."""

    spectra: np.ndarray  # (n_basis, n_omega)
    labels: list[str]
    rows: dict = field(default_factory=dict)  # cache of design rows


def make_basis(larmor: dict[str, float], use_lines: bool = True) -> Basis:
    return _make_basis(tuple(sorted(larmor.items())), use_lines)


@lru_cache(maxsize=None)
def _make_basis(larmor_items: tuple, use_lines: bool) -> Basis:
    """Continuum Lorentzians at w = 0 plus (optionally) lines at the nuclear Larmor frequencies."""
    spectra, labels = [], []
    for g in np.geomspace(0.1, 30.0, 8):
        spectra.append(lorentzian(0.0, g))
        labels.append(f"cont{g:.2g}")
    if use_lines:
        for name, w in larmor_items:
            for g in (0.3, 1.0, 3.0):
                spectra.append(lorentzian(w, g))
                labels.append(f"{name}:{g}")
    return Basis(np.array(spectra), labels)


def design_row(seq: Sequence, total_time: float, basis: Basis) -> np.ndarray:
    """chi contribution of every basis function for ``seq`` at ``total_time``."""
    f = np.abs(switching_ft(seq, total_time)) ** 2
    return (basis.spectra * f * _DOMEGA).sum(axis=1)


class FilterFunctionModel:
    """Fit the basis weights to coherence probes and predict coherence of new sequences."""

    def __init__(self, basis: Basis) -> None:
        self.basis = basis
        self.weights = np.zeros(len(basis.labels))

    def fit(self, probes: list[tuple[Sequence, float]], c_hat: np.ndarray, shots: int) -> None:
        """Weighted NNLS fit of chi = -ln(c) to measured coherences."""
        c = np.clip(np.asarray(c_hat, dtype=float), 0.03, 0.999)
        chi = -np.log(c)
        sigma_c = np.sqrt(np.clip(1 - c**2, 1e-3, None) / max(shots, 1) + 1e-6)
        sigma_chi = sigma_c / c
        g = np.array([self._row(s, t) for s, t in probes])
        scale = g.max()
        w, _ = nnls(g / scale / sigma_chi[:, None], chi / sigma_chi)
        self.weights = w / scale

    def _row(self, seq: Sequence, total_time: float) -> np.ndarray:
        key = (seq.name, float(total_time))
        if key not in self.basis.rows:
            self.basis.rows[key] = design_row(seq, total_time, self.basis)
        return self.basis.rows[key]

    def predict_coherence(self, seq: Sequence, total_time: float) -> float:
        """Predicted exp(-chi) for a candidate sequence."""
        return float(np.exp(-self._row(seq, total_time) @ self.weights))


def dephasing_fidelity(c: float) -> float:
    """Six-state mean fidelity of a pure-dephasing qubit channel with x-y coherence ``c``."""
    return (2.0 + c) / 3.0
