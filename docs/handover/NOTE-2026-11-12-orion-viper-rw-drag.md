# NOTE 2026-11-12 Orion -> viper: rw-drag (localise the drag deficit at P 3.33 tau 10)

FINAL 84cd20888 (rw-amp1e6 build: Release, MPI, -ffp-contract=off). Host orion004, 36 cores. 36/36 cases "ok", 0 FATAL.
Wall: cases_drag 30 s, cases_scale 30 s (concurrent). The pre-check passed: arm overrides present, dt_retry = false, nx2 = 4. RW_C=20 changes only
c_light (10 -> 20) and therefore dt / implicit_cfl (2.3 -> 9.4-10.2).

## Reading (eps(G) from the 4-mode N->inf limit)
| arm | edd | vetf | vs base |
|---|---|---|---|
| base vimp=true (control) | -3.061e-4 (resid 0.09) | -3.015e-4 (0.08) | reproduces viper's -3.15e-4 |
| (i) gas_src=conservative | -1.787e-1 (0.10) | -1.824e-1 (0.11) | 600x worse; slow-mode amp 0.14 (not the fix) |
| (ii) vimp/cs2_vimp/vimp_bc false | +2.335e-3 (0.00) | +2.397e-3 (0.00) | opposite sign, 8x larger |
| (iii) time_scheme=be | -9.19e-3 (0.52) | -1.04e-2 (0.50) | dt-biased, as expected |
| P 1.667 | -1.79e-5 (0.30) | -1.81e-5 (0.18) | predicted x1/2; observed x1/17 |
| c = 20 (P 3.33) | -3.1e-6 (0.69) | -5.7e-6 (0.63) | predicted x1/4; observed x1/100 / x1/53 (resid large: eps(G) is not a good one-term fit there) |

- **conservative2 + vimp is the best of the variants, not the culprit.** Removing the dA sharing (i) or the implicit velocity (ii) makes the drag
  error far larger. (ii) also flips its sign. The drag term itself, on both sides, depends strongly on how the velocity enters the coupling.
- **The scaling is not linear in E0/(rho c^2).** Halving P cuts eps by about 17 (~P^4); doubling c cuts the acoustic offset |dw4| from 1.6e-3 to
  8e-5 / 1.6e-4. The fit residual rises to 0.6-0.7 at c = 20, so there the offset is near the fit floor and no longer drag-shaped.
  So the deficit is a higher-order effect of radiation-dominated, optically-thick coupling at c ~ 2-4 c_ph, not a fixed O(E0/rho c^2) factor.
  The c = 20 arm also runs at 4x larger implicit_cfl: part of the c-scaling may be the dt/t_diff ratio, not c itself.
- The base edd and vetf ladders are the rw-edd / rw-p333 runs again: 1-mode omegas are bitwise those NOTEs' values.
  The 4-mode limits are edd 13.700935 - 4.087871i, vetf 13.346034 - 3.661931i.

## python3 ana_drag.py runs cases_drag.txt
```
c = 10.0
== edd P 1.667 vimp=true
  N   64  1-mode 10.200661-3.708785j  4-mode 10.204083-3.709615j  slow amp 1.22e-03  rms 3.8e-05
  N  128  1-mode 10.206693-3.715340j  4-mode 10.207422-3.715445j  slow amp 3.50e-04  rms 4.4e-06
  N  256  1-mode 10.208316-3.717142j  4-mode 10.208532-3.717159j  slow amp 1.04e-04  rms 1.2e-06
  inf: 4-mode 10.208902-3.717730j  1-mode 10.208857-3.717742j  ref 10.208996-3.717690j  dw4 = -9.38e-05-4.01e-05j  eps(G) = -1.794e-05 (resid 0.30)
== edd P 3.333 rm_implicit_gas_src=conservative vimp=true
  N   64  1-mode 13.384753-6.090710j  4-mode 13.467320-4.974092j  slow amp 1.47e-01  rms 1.3e-05
  N  128  1-mode 13.417692-6.071745j  4-mode 13.482403-4.972493j  slow amp 1.45e-01  rms 4.5e-06
  N  256  1-mode 13.428311-6.062067j  4-mode 13.486914-4.970435j  slow amp 1.44e-01  rms 1.9e-06
  inf: 4-mode 13.488418-4.969750j  1-mode 13.431851-6.058841j  ref 13.701320-4.086365j  dw4 = -2.13e-01-8.83e-01j  eps(G) = -1.787e-01 (resid 0.10)
== edd P 3.333 rm_implicit_vimp=false
  N   64  1-mode 13.734822-3.959022j  4-mode 13.724589-4.015798j  slow amp 1.13e-02  rms 4.2e-06
  N  128  1-mode 13.718597-4.022860j  4-mode 13.713002-4.051935j  slow amp 5.68e-03  rms 4.3e-07
  N  256  1-mode 13.710132-4.054691j  4-mode 13.707193-4.069410j  slow amp 2.86e-03  rms 2.7e-07
  inf: 4-mode 13.705257-4.075235j  1-mode 13.707310-4.065302j  ref 13.701320-4.086365j  dw4 = +3.94e-03+1.11e-02j  eps(G) = +2.335e-03 (resid 0.00)
== edd P 3.333 rm_time_scheme=be vimp=true
  N   64  1-mode 13.330583-4.484628j  4-mode 13.429498-4.276092j  slow amp 3.85e-02  rms 3.1e-06
  N  128  1-mode 13.530196-4.287534j  4-mode 13.569652-4.185152j  slow amp 1.92e-02  rms 1.1e-06
  N  256  1-mode 13.619227-4.187235j  4-mode 13.636628-4.136753j  slow amp 9.54e-03  rms 5.0e-07
  inf: 4-mode 13.658953-4.120619j  1-mode 13.648903-4.153802j  ref 13.701320-4.086365j  dw4 = -4.24e-02-3.43e-02j  eps(G) = -9.188e-03 (resid 0.52)
== edd P 3.333 vimp=true
  N   64  1-mode 13.679731-4.091637j  4-mode 13.684286-4.083043j  slow amp 1.82e-03  rms 5.1e-06
  N  128  1-mode 13.695168-4.089714j  4-mode 13.696520-4.086434j  slow amp 6.07e-04  rms 1.0e-05
  N  256  1-mode 13.699142-4.089324j  4-mode 13.699831-4.087512j  slow amp 3.82e-04  rms 3.7e-06
  inf: 4-mode 13.700935-4.087871j  1-mode 13.700466-4.089194j  ref 13.701320-4.086365j  dw4 = -3.85e-04-1.51e-03j  eps(G) = -3.061e-04 (resid 0.09)
== vetf P 1.667 vimp=true
  N   64  1-mode 10.152481-3.312568j  4-mode 10.155258-3.313149j  slow amp 1.07e-03  rms 3.7e-05
  N  128  1-mode 10.157775-3.318737j  4-mode 10.158369-3.318829j  slow amp 3.18e-04  rms 4.7e-06
  N  256  1-mode 10.159167-3.320392j  4-mode 10.159343-3.320408j  slow amp 9.33e-05  rms 1.2e-06
  inf: 4-mode 10.159668-3.320934j  1-mode 10.159631-3.320944j  ref 10.159739-3.320887j  dw4 = -7.05e-05-4.72e-05j  eps(G) = -1.814e-05 (resid 0.18)
== vetf P 3.333 rm_implicit_gas_src=conservative vimp=true
  N   64  1-mode 13.041761-5.324009j  4-mode 13.144682-4.456989j  slow amp 1.38e-01  rms 1.3e-05
  N  128  1-mode 13.070936-5.310923j  4-mode 13.158868-4.456193j  slow amp 1.36e-01  rms 4.4e-06
  N  256  1-mode 13.080269-5.303593j  4-mode 13.162996-4.454419j  slow amp 1.35e-01  rms 1.8e-06
  inf: 4-mode 13.164372-4.453828j  1-mode 13.083380-5.301150j  ref 13.346366-3.660631j  dw4 = -1.82e-01-7.93e-01j  eps(G) = -1.824e-01 (resid 0.11)
== vetf P 3.333 rm_implicit_vimp=false
  N   64  1-mode 13.376956-3.550225j  4-mode 13.366393-3.596152j  slow amp 1.06e-02  rms 4.4e-06
  N  128  1-mode 13.362215-3.605876j  4-mode 13.356625-3.629379j  slow amp 5.33e-03  rms 2.9e-07
  N  256  1-mode 13.354438-3.633389j  4-mode 13.351542-3.645285j  slow amp 2.68e-03  rms 2.5e-07
  inf: 4-mode 13.349848-3.650587j  1-mode 13.351845-3.642560j  ref 13.346366-3.660631j  dw4 = +3.48e-03+1.00e-02j  eps(G) = +2.397e-03 (resid 0.00)
== vetf P 3.333 rm_time_scheme=be vimp=true
  N   64  1-mode 12.984015-4.028192j  4-mode 13.079259-3.858083j  slow amp 3.66e-02  rms 3.0e-06
  N  128  1-mode 13.178740-3.846899j  4-mode 13.217656-3.762854j  slow amp 1.82e-02  rms 1.0e-06
  N  256  1-mode 13.266493-3.753988j  4-mode 13.283808-3.712390j  slow amp 9.09e-03  rms 4.4e-07
  inf: 4-mode 13.305859-3.695568j  1-mode 13.295744-3.723018j  ref 13.346366-3.660631j  dw4 = -4.05e-02-3.49e-02j  eps(G) = -1.043e-02 (resid 0.50)
== vetf P 3.333 vimp=true
  N   64  1-mode 13.325986-3.663332j  4-mode 13.330350-3.656704j  slow amp 1.59e-03  rms 3.0e-05
  N  128  1-mode 13.340693-3.663039j  4-mode 13.341982-3.660514j  slow amp 5.74e-04  rms 1.3e-05
  N  256  1-mode 13.344398-3.663014j  4-mode 13.345021-3.661577j  slow amp 3.43e-04  rms 3.6e-06
  inf: 4-mode 13.346034-3.661931j  1-mode 13.345633-3.663006j  ref 13.346366-3.660631j  dw4 = -3.33e-04-1.30e-03j  eps(G) = -3.015e-04 (resid 0.08)
```

## RW_C=20 python3 ana_drag.py runs_c20 cases_scale.txt
```
c = 20.0
== edd P 3.333 vimp=true
  N   64  1-mode 6.186121-6.031833j  4-mode 6.185293-6.031945j  slow amp 9.68e-04  rms 1.9e-07
  N  128  1-mode 6.169865-6.046440j  4-mode 6.169650-6.046441j  slow amp 2.30e-04  rms 4.8e-08
  N  256  1-mode 6.165594-6.050383j  4-mode 6.165534-6.050370j  slow amp 5.50e-05  rms 1.9e-08
  inf: 4-mode 6.164162-6.051680j  1-mode 6.164170-6.051697j  ref 6.164119-6.051746j  dw4 = +4.34e-05+6.64e-05j  eps(G) = -3.125e-06 (resid 0.69)
== vetf P 3.333 vimp=true
  N   64  1-mode 6.733351-6.372800j  4-mode 6.731628-6.373566j  slow amp 9.27e-04  rms 3.2e-07
  N  128  1-mode 6.715739-6.393838j  4-mode 6.715147-6.393920j  slow amp 2.04e-04  rms 9.1e-08
  N  256  1-mode 6.711232-6.399240j  4-mode 6.711076-6.399294j  slow amp 4.31e-05  rms 9.9e-09
  inf: 4-mode 6.709720-6.401085j  1-mode 6.709730-6.401040j  ref 6.709868-6.401137j  dw4 = -1.49e-04+5.15e-05j  eps(G) = -5.740e-06 (resid 0.63)
```

Run trees: Orion /orion/ptmp/jinma/rwdrag_1112/rw_drag (runs, runs_c20).
