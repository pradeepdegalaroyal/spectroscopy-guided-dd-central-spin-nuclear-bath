Table A: regret versus shot budget at flip-angle error 0.05 (mean ± std over bath draws)

| Selection rule | 100 | 300 | 1,000 | 3,000 | 10,000 | 100,000 |
|---|---|---|---|---|---|---|
| direct uniform | 0.206 ± 0.021 | 0.187 ± 0.023 | 0.098 ± 0.009 | 0.030 ± 0.003 | 0.009 ± 0.003 | 0.001 ± 0.000 |
| direct successive halving | 0.169 ± 0.016 | 0.144 ± 0.018 | 0.065 ± 0.005 | 0.018 ± 0.005 | 0.005 ± 0.002 | 0.001 ± 0.000 |
| spectroscopy (filter function) | 0.153 ± 0.018 | 0.076 ± 0.019 | 0.039 ± 0.011 | 0.022 ± 0.005 | 0.023 ± 0.011 | 0.024 ± 0.012 |
| spectroscopy, no Larmor lines | 0.158 ± 0.019 | 0.069 ± 0.023 | 0.064 ± 0.019 | 0.073 ± 0.020 | 0.071 ± 0.019 | 0.063 ± 0.019 |
| spectroscopy, no pulse calibration | 0.215 ± 0.017 | 0.276 ± 0.013 | 0.233 ± 0.007 | 0.212 ± 0.010 | 0.211 ± 0.007 | 0.206 ± 0.018 |

Table B: regret at 1,000 shots versus flip-angle error

| Selection rule | ε = 0 | ε = 0.02 | ε = 0.05 | ε = 0.1 |
|---|---|---|---|---|
| fixed Hahn echo | 0.159 ± 0.021 | 0.158 ± 0.021 | 0.150 ± 0.020 | 0.111 ± 0.022 |
| fixed CPMG-8 | 0.180 ± 0.024 | 0.185 ± 0.023 | 0.206 ± 0.019 | 0.238 ± 0.020 |
| fixed XY4-8 | 0.180 ± 0.024 | 0.180 ± 0.024 | 0.176 ± 0.024 | 0.150 ± 0.024 |
| best single sequence | 0.064 ± 0.020 | 0.064 ± 0.020 | 0.063 ± 0.019 | 0.064 ± 0.017 |
| direct uniform | 0.099 ± 0.009 | 0.105 ± 0.010 | 0.098 ± 0.009 | 0.086 ± 0.011 |
| direct successive halving | 0.062 ± 0.007 | 0.071 ± 0.009 | 0.065 ± 0.005 | 0.064 ± 0.011 |
| spectroscopy (filter function) | 0.049 ± 0.015 | 0.044 ± 0.014 | 0.039 ± 0.011 | 0.081 ± 0.025 |
| spectroscopy, no Larmor lines | 0.065 ± 0.020 | 0.077 ± 0.018 | 0.064 ± 0.019 | 0.091 ± 0.008 |
| spectroscopy, no pulse calibration | 0.047 ± 0.010 | 0.253 ± 0.015 | 0.233 ± 0.007 | 0.299 ± 0.023 |

Table C: model-error floor (unlimited shots) versus model-free selection at 100,000 shots

| Selection rule | ε = 0 | ε = 0.02 | ε = 0.05 | ε = 0.1 |
|---|---|---|---|---|
| spectroscopy (infinite shots) | 0.043 ± 0.042 | 0.026 ± 0.019 | 0.022 ± 0.013 | 0.077 ± 0.064 |
| spectroscopy, no Larmor lines (infinite shots) | 0.064 ± 0.020 | 0.064 ± 0.020 | 0.063 ± 0.019 | 0.134 ± 0.017 |
| direct successive halving, 100,000 shots | 0.000 ± 0.000 | 0.001 ± 0.000 | 0.001 ± 0.000 | 0.001 ± 0.001 |

Table D: does the surrogate have lower mean regret than successive halving? (yes / no)

| Budget | ε = 0 | ε = 0.02 | ε = 0.05 | ε = 0.1 |
|---|---|---|---|---|
| 100 | yes | no | yes | no |
| 300 | yes | no | yes | no |
| 1,000 | yes | yes | yes | no |
| 3,000 | no | no | no | no |
| 10,000 | no | no | no | no |
| 100,000 | no | no | no | no |
