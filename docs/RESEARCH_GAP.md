# Research gap

## State of the art (from the supplied literature)

An electron spin in a gallium arsenide (GaAs) quantum dot (QD) is a natural memory node for a quantum network, but its coherence is limited by the surrounding nuclear spins. Dynamical decoupling (DD), a train of control pulses that refocuses slow environmental fields, is the standard tool for extending it, and the literature offers many pulse sequences (Hahn echo, Carr-Purcell-Meiboom-Gill (CPMG), Uhrig dynamical decoupling (UDD), XY-type phase-alternating trains). The supplied papers do not study DD of QD spins directly; they frame the problem from neighbouring angles:

- Degen et al., 2017 (doi:10.1103/RevModPhys.89.035002) present quantum sensing, with spin qubits as a leading platform, as the experimentalist's toolbox for measuring noise and fields. Spin-qubit measurements are the natural way to learn about the bath.
- Streltsov et al., 2017 (doi:10.1103/RevModPhys.89.041003) treat coherence as a quantifiable resource whose dynamical evolution must be controlled.
- Aharonovich et al., 2016 (doi:10.1038/nphoton.2016.186) and Wang et al., 2020 (doi:10.1038/s41566-019-0532-1) place QD emitters and integrated photonics at the centre of scalable quantum technologies; Nusyura et al., 2026 (doi:10.1007/s42484-026-00402-8) list quantum memories among the open bottlenecks.
- Yousif et al., 2027 (doi:10.1007/978-3-032-39708-9_27) show that device performance is ultimately priced in shots per second and tail latency, and McClean et al., 2016 (doi:10.1088/1367-2630/18/2/023023) show that a measurement-hungry hybrid loop can be made far cheaper with the right classical optimiser. Measurement budget is a first-class cost.
- Abdiwi et al., 2027 (doi:10.1007/978-3-032-39708-9_30) improve spin coherence by device design, whereas Gnatowski et al., 2027 (doi:10.1007/978-3-032-37936-8_15) show that an exact optimiser can beat a standard heuristic for compilation; neither addresses choosing control for a given bath.

## Related work

| Paper | Year | Role | What it does | What it does NOT cover |
|---|---|---|---|---|
| Degen et al. (10.1103/RevModPhys.89.035002) | 2017 | landmark | Principles and methods of quantum sensing with spin qubits | A shot-budgeted rule for picking a decoupling sequence |
| Streltsov et al. (10.1103/RevModPhys.89.041003) | 2017 | landmark | Coherence as a resource: quantification, manipulation, dynamics | Measurement-limited protection of a concrete memory |
| Aharonovich et al. (10.1038/nphoton.2016.186) | 2016 | landmark | Solid-state single-photon emitters including QDs | Spin-memory control and its calibration |
| McClean et al. (10.1088/1367-2630/18/2/023023) | 2016 | landmark | Hybrid quantum-classical loops, measurement-cheap optimisers | Control selection for spin memories |
| Wang et al. (10.1038/s41566-019-0532-1) | 2020 | review | Integrated photonic quantum technologies | Matter-spin memories and their decoupling |
| Kjaergaard et al. (10.1146/annurev-conmatphys-031119-050605) | 2020 | review | Superconducting-qubit hardware, gates and readout | Nuclear-spin baths of semiconductor spins |
| Ashida et al. (10.1080/00018732.2021.1876991) | 2020 | review | Effective non-Hermitian descriptions of open quantum systems | When a classical-noise description of a finite bath fails |
| Abdiwi et al. (10.1007/978-3-032-39708-9_30) | 2027 | recent | Spin-orbit and topological device design for longer coherence | Control-sequence choice on a fixed device |
| Gnatowski et al. (10.1007/978-3-032-37936-8_15) | 2027 | recent | Exact routing that beats a standard heuristic compiler | Control of spin memories |
| Nusyura et al. (10.1007/s42484-026-00402-8) | 2026 | adjacent | Review of photonic quantum information, memories, loss | Quantified memory-control costs |
| Yousif et al. (10.1007/978-3-032-39708-9_27) | 2027 | adjacent | Maps device metrics to shots per second and latency | Spin coherence control |
| Das and Roy Choudhury (10.1016/j.rse.2026.115665) | 2026 | adjacent | Statistics-limited sensing at very low photon numbers | Spin control under finite statistics |
| Sallah et al. (10.1016/j.jpcs.2026.114199) | 2027 | adjacent | Spectroscopic characterisation of defect emitters | Coherent spin control |

## The gap

Recent device-level papers price quantum hardware in shots, landmark and review papers establish spin-qubit measurement as a sensor of the environment, and the open problem named in the reviews is how to protect a memory when characterisation time is scarce. Experimentally, a group running a GaAs QD spin node faces exactly this: every candidate sequence (and every storage time demanded by a network heralding window) costs shots to test, while a noise-spectroscopy model costs a fixed number of probe shots but assumes classical Gaussian noise, which a small, strongly coupled nuclear bath violates. No paper in the brief, and to the author's knowledge no simple benchmark, tells an experimenter at which shot budget the spectral model is the better investment than direct trial of sequences, or how large the model-misspecification floor is for a genuinely quantum bath.

**Testable question.** For an electron spin coupled to an exactly simulated quantum nuclear bath (8 spins-1/2, three GaAs isotopes, hyperfine and flip-flop couplings), with finite flip-angle error, at what total shot budget does a filter-function surrogate fitted to CPMG probes choose decoupling sequences with lower regret than model-free selection of the same candidates?

**Hypotheses (fixed before the final run; the budget grid was set from a pilot on bath seeds 0-2, the final run uses disjoint seeds 10-14).**

- H1. At 1000 shots or fewer per decision table, the surrogate has lower mean regret than the best model-free selector (successive halving) at every flip-angle error studied.
- H2. With unlimited shots the surrogate retains a positive regret floor (greater than 0.01 in fidelity) caused by classical-model misspecification, whereas the model-free selector's regret goes to zero; hence a crossing budget exists below 10^5 shots.
- H3. Removing the nuclear-Larmor lines from the spectral basis raises the floor, and removing the pulse-error calibration raises regret when the flip-angle error is 0.05 or more.

**Success metric.** Regret = oracle memory fidelity minus fidelity of the selected sequence, averaged over eight storage times; fidelity is the mean over six cardinal electron states of the exact state fidelity. Reported as mean +- standard deviation over five bath draws.

**Novelty in one sentence.** A shot-budgeted, exactly simulated comparison of a classical noise-spectroscopy surrogate against model-free selection of decoupling sequences for a quantum nuclear-spin bath, locating the crossing budget and the misspecification floor.

**Scale statement.** All results are exact statevector simulations of 8 nuclear spins-1/2 (Hilbert-space dimension 256 times 2 for the electron) with instantaneous pulses; couplings are effective values chosen to give sub-microsecond free dephasing and visible nuclear-resonance structure, not a fit to a particular device.
