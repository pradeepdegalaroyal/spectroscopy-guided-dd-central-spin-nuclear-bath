"""Exact statevector simulation of DD sequences on the electron + nuclear-bath system.

The nuclear bath starts in the maximally mixed state (high temperature, as for unpolarised
GaAs nuclei); it is sampled exactly by propagating every nuclear basis state.
"""
from __future__ import annotations

import numpy as np
import torch

from .bath import BathConfig, conditional_hamiltonians, draw_bath
from .sequences import Sequence

_S2 = 2**-0.5
CARDINAL_STATES = np.array(
    [[_S2, _S2], [_S2, -_S2], [_S2, 1j * _S2], [_S2, -1j * _S2], [1, 0], [0, 1]], dtype=complex
)


def pick_device() -> torch.device:
    """CUDA if available, otherwise CPU."""
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


class CentralSpinBath:
    """Electron spin coupled to ``n_nuclei`` nuclear spins-1/2, propagated exactly."""

    def __init__(self, cfg: BathConfig, seed: int, device: torch.device | None = None) -> None:
        self.cfg, self.seed = cfg, seed
        self.device = device or pick_device()
        self.draw = draw_bath(cfg, seed)
        self.dim = 2**cfg.n_nuclei
        h_up, h_dn = conditional_hamiltonians(self.draw)
        self._eig = []
        for h in (h_up, h_dn):
            e, v = np.linalg.eigh(h)
            vt = torch.tensor(v, dtype=torch.complex64, device=self.device)
            self._eig.append((e, vt, vt.conj().T.contiguous()))

    # -- propagation -------------------------------------------------------------------
    def _free(self, psi: torch.Tensor, tau: float) -> torch.Tensor:
        """Propagate psi[(electron), nuclear, columns] for time tau under H_up / H_down."""
        if tau <= 0:
            return psi
        out = torch.empty_like(psi)
        for s, (e, v, vh) in enumerate(self._eig):
            ph = torch.tensor(np.exp(-1j * e * tau), dtype=torch.complex64, device=self.device)
            out[s] = v @ (ph[:, None] * (vh @ psi[s]))
        return out

    @staticmethod
    def _pulse_matrix(axis: int, eps: float, device: torch.device) -> torch.Tensor:
        theta = np.pi * (1.0 + eps)
        sig = np.array([[0, 1], [1, 0]], dtype=complex) if axis == 0 else np.array([[0, -1j], [1j, 0]])
        r = np.cos(theta / 2) * np.eye(2) - 1j * np.sin(theta / 2) * sig
        return torch.tensor(r, dtype=torch.complex64, device=device)

    def _run(self, seq: Sequence, total_time: float, eps: float, states: np.ndarray) -> torch.Tensor:
        """Final state tensor of shape (2, dim, n_states * dim)."""
        d, ns = self.dim, len(states)
        eye = torch.eye(d, dtype=torch.complex64, device=self.device)
        st = torch.tensor(states, dtype=torch.complex64, device=self.device)  # (ns, 2)
        psi = torch.einsum("se,mc->emsc", st, eye).reshape(2, d, ns * d)
        times = np.concatenate([[0.0], seq.fractions() * total_time, [total_time]])
        mats = [self._pulse_matrix(a, eps, self.device) for a in seq.axes()]
        for j in range(len(times) - 1):
            psi = self._free(psi, times[j + 1] - times[j])
            if j < len(mats):
                psi = torch.einsum("ab,bmc->amc", mats[j], psi)
        return psi

    # -- observables -------------------------------------------------------------------
    def memory_fidelity(self, seq: Sequence, total_time: float, eps: float = 0.0) -> float:
        """Mean fidelity of the electron state over six cardinal states (even-n sequences)."""
        d = self.dim
        psi = self._run(seq, total_time, eps, CARDINAL_STATES)
        st = torch.tensor(CARDINAL_STATES, dtype=torch.complex64, device=self.device)
        psi = psi.reshape(2, d, len(CARDINAL_STATES), d)
        amp = torch.einsum("se,emsc->msc", st.conj(), psi)  # overlap with ideal state
        return float((amp.abs() ** 2).sum(dim=(0, 2)).mean() / d)

    def coherence_x(self, seq: Sequence, total_time: float, eps: float = 0.0) -> float:
        """Expectation of sigma_x for the initial state |+x> (a dephasing probe)."""
        d = self.dim
        psi = self._run(seq, total_time, eps, CARDINAL_STATES[:1]).reshape(2, d, d)
        sx = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex64, device=self.device)
        val = torch.einsum("emc,ef,fmc->", psi.conj(), sx, psi).real
        return float(val / d)
