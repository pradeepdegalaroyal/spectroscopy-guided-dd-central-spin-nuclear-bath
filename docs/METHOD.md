# Method

This document specifies the dynamical-decoupling (DD) selection study for an electron spin in a gallium arsenide (GaAs) quantum dot (QD).

## System (the "truth")

Electron spin (S_z = +-1/2) in the rotating frame, coupled to N = 8 nuclear spins-1/2 whose species (75As, 69Ga, 71Ga; probabilities 0.5, 0.3, 0.2) set the Larmor frequencies at 0.1 T. For electron state s = +-1/2 the nuclear Hamiltonian is

H_s = sum_k w_k Iz_k + s * sum_k (A_k Iz_k + B_k Ix_k) + sum_{k<l} J_kl (Ix_k Ix_l + Iy_k Iy_l)

with parallel couplings A_k, transverse couplings B_k = |A_k| r_k (r_k uniform in [0.3, 1.5]) and flip-flop couplings J_kl. Frequencies are in rad per microsecond. The nuclei start maximally mixed and are sampled exactly by propagating every nuclear basis state. Pulses are instantaneous rotations by pi(1 + eps) about x or y, with flip-angle error eps in {0, 0.02, 0.05, 0.1}. Each propagation segment uses the eigendecomposition of H_+ and H_- (computed once per bath).

**Memory fidelity.** Mean over the six cardinal electron states of the overlap of the final reduced electron state with the ideal one, for even pulse number. This is the decision metric.

## Candidates, storage times, probes

- Candidates (15): free evolution; Carr-Purcell-Meiboom-Gill (CPMG) and Uhrig dynamical decoupling (UDD) sequences with n = 2, 4, 8, 16, 32 pulses; XY4-type four-pulse phase-alternating trains (x, y, x, y) repeated to n = 4, 8, 16, 32 pulses.
- Storage times (8): 3, 4.5, 7, 10, 15, 22, 33, 50 microseconds, standing for different heralding waits. A selector must return one sequence per storage time.
- Spectroscopy probes (50): CPMG-n, n in {2, 4, 8, 16, 32}, at ten times from 1.5 to 64 microseconds. The probe state is |+x>; the measured quantity is the sigma_x expectation. The probe times do not coincide with the storage times.
- Calibration probes (2): six-state fidelity of a 16-pulse x-train and a 16-pulse XY4-train at 0.05 microseconds, which isolates pulse error from bath dephasing.

## Selection rules

1. **Spectroscopy surrogate (proposed).** Fit chi = -ln(c) of the probe coherences to chi = sum_b beta_b * integral of S_b(w) |Y(w,T)|^2 dw, where Y is the Fourier transform of the pulse sign-switching function and S_b are Lorentzians: eight broad lines at zero frequency plus three widths at each of the three nuclear Larmor frequencies (known from the applied field). Weighted non-negative least squares with shot-noise weights. Predicted memory fidelity of a candidate = (2 + exp(-chi)) / 3 times g^n, where g is the per-pulse factor from the calibration probe of that pulse family. The best predicted candidate is chosen for each storage time. Shots are split equally over the 52 probes.
2. **Model-free direct trial (classical baseline).** Spend the shots measuring candidate fidelities. Two variants: uniform allocation over all 15 x 8 cells, and successive halving per storage time (the stronger baseline).
3. **Standard fixed sequences (literature practice).** Hahn echo (CPMG-2), CPMG-8, XY4-8, and the best single sequence over all storage times (an oracle that fixes one sequence).
4. **Oracle.** The best candidate from the exact table.

The "standard quantum baseline" of this study is therefore the fixed standard DD sequences plus the exactly simulated quantum bath as truth; the "classical baseline" is the model-free direct trial. The surrogate itself is the classical-noise model under test.

## Metrics

Regret in memory fidelity averaged over storage times; mean selected fidelity; surrogate root-mean-square error (RMSE) and mean Spearman rank correlation against exact fidelities at unlimited shots; crossing budget (smallest budget at which successive halving is at least as good as the surrogate).

## Ablations and robustness sweeps

- Ablation A: spectral basis without nuclear-Larmor lines (continuum only).
- Ablation B: no pulse-error calibration (g = 1).
- Noise/robustness sweep: flip-angle error eps in {0, 0.02, 0.05, 0.1}; shot budget in {100, 300, 1000, 3000, 10^4, 10^5}; unlimited shots as the model-error floor.
- Size check: bath size N = 4, 6, 8 for the regret, N = 4 to 9 for the exact-simulation cost.

## Seeds and statistics

Five bath draws (seeds 10-14) set species, hyperfine and flip-flop couplings. For each draw and each budget the shot-noise realisation is repeated 30 times; reported spreads are the standard deviation of the 30-repeat mean across the five draws.
