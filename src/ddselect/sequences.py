"""Dynamical-decoupling (DD) sequence families used as candidates and probes.

All sequences act on a central electron spin with instantaneous pulses. A sequence is
defined by a family, a pulse number ``n`` and (after scaling by the storage time ``T``) the
pulse times and rotation axes.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

AXIS_X, AXIS_Y = 0, 1


@dataclass(frozen=True)
class Sequence:
    """A DD sequence: ``family`` in {free, cpmg, udd, xy4} with ``n`` pi pulses (n even)."""

    family: str
    n: int

    @property
    def name(self) -> str:
        return "free" if self.family == "free" else f"{self.family}-{self.n}"

    def fractions(self) -> np.ndarray:
        """Pulse positions as fractions of the storage time, strictly inside (0, 1)."""
        j = np.arange(1, self.n + 1)
        if self.family == "free":
            return np.zeros(0)
        if self.family in ("cpmg", "xy4"):
            return (j - 0.5) / self.n
        if self.family == "udd":
            return np.sin(np.pi * j / (2 * self.n + 2)) ** 2
        raise ValueError(f"unknown family {self.family}")

    def axes(self) -> list[int]:
        """Rotation axis (0 = x, 1 = y) of every pulse."""
        if self.family == "xy4":
            return [AXIS_X if k % 2 == 0 else AXIS_Y for k in range(self.n)]
        return [AXIS_X] * self.n

    def pulse_family(self) -> str:
        """Group used for pulse-error calibration: sequences with the same axes share it."""
        return "xy4" if self.family == "xy4" else "x"


def candidate_set(max_pulses: int = 32) -> list[Sequence]:
    """Candidate library: free evolution, CPMG, UDD and XY4 with n up to ``max_pulses``."""
    cands = [Sequence("free", 0)]
    n = 2
    while n <= max_pulses:
        cands.append(Sequence("cpmg", n))
        cands.append(Sequence("udd", n))
        if n >= 4:
            cands.append(Sequence("xy4", n))
        n *= 2
    return cands


def standard_baselines() -> dict[str, Sequence]:
    """Fixed sequences a practitioner would use without characterising the bath."""
    return {
        "fixed Hahn echo": Sequence("cpmg", 2),
        "fixed CPMG-8": Sequence("cpmg", 8),
        "fixed XY4-8": Sequence("xy4", 8),
    }
