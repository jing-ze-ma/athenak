# NOTE 2026-11-01 Caltech -> viper (USER): tag-only floor arm results; the 1.2e-5 bar shell difference is (A), the intended correction, and it sits in floor-supported cells with local p < 1e-6 bar

Answers TASK 0d9a039e.

## Runs
- Binary: tracers-defaults-1101 @ 8ccae663, athena_gpu_dhj_trcdef_8ccae663_fmadoff, md5 de6bf3f204de89d84ef0ebce296d707b.
- Both arms: fresh from t = 0 to 0.5 rot (55077 s), r1m_trc deck (T1 deck), xe_tracer = false, 1 rank. They shared one H200 (02-03 GPU 1).
  - **tagonly**: floor_tracer = true, floor_tracer_fmax = 1e-12 (deck gates_1101/tagonly.athinput).
  - **fl**: floor_tracer = true, defaults (margin 2, fmax 0.01).
- Both ran clean: 8700 cycles, mean dt 6.24 s, 14.8 RKG stages, 0.20 s/cycle (shared GPU).
- hst means over 0.4-0.5 rot:

| arm | Mfloor net [g/s] | Mdrain [g/s] | gross [g/s] | Mdot_top [g/s] |
|---|---|---|---|---|
| fl | 9.83e13 | 5.53e13 | 1.54e14 | 4.6e10 |
| tagonly | 1.56e14 | 1.3e7 (i.e. none) | 1.56e14 | 3.4e11 |

- Bins: /resnick/groups/carnegie_poc/jingze/trc_1101/runs/tg_{fl,tagonly}/bin/dhj.mhd_w_bcc.00004.bin (t = 0.5 rot).
- fl on 8ccae663 reproduces the T1 fl arm (89213fcb) to all printed digits, as expected: the x_e fix does nothing with xe off.

## tag_where.py output (as requested: `tag_where.py <fl bin> <tagonly bin>`)
```
file runs/tg_fl/bin/dhj.mhd_w_bcc.00004.bin: total tag mass fraction 5.261e-09
| log p [bar] | tag mass share | tag/mass in band | share at rho<2rho_fl | share 2..4 | share >4 |
|---|---|---|---|---|---|
| -9..-8 | 0.059 | 9.89e-01 | 0.991 | 0.009 | 0.000 |
| -8..-7 | 0.270 | 5.62e-01 | 0.934 | 0.066 | 0.000 |
| -7..-6 | 0.622 | 6.06e-01 | 0.984 | 0.016 | 0.000 |
| -6..-5 | 0.048 | 2.25e-02 | 1.000 | 0.000 | 0.000 |
| -5..-4 | 0.000 | 1.43e-05 | 1.000 | 0.000 | 0.000 |
| -4..-3 | 0.000 | 3.85e-14 | 0.000 | 0.000 | 1.000 |

vs runs/tg_tagonly/bin/dhj.mhd_w_bcc.00004.bin: per band, mass ratio (A/B - 1) and tag/mass of B
| -9..-8 | +1.796e-01 | 9.907e-01 |
| -8..-7 | -7.740e-02 | 6.587e-01 |
| -7..-6 | -1.043e-01 | 6.564e-01 |
| -6..-5 | -6.391e-03 | 3.551e-02 |
| -5..-4 | -8.861e-04 | 2.998e-04 |
| -4..-3 | +1.766e-03 | 6.522e-07 |
| -3..-2 | +1.334e-04 | 1.701e-09 |
```

tagonly alone (tag location):
```
file runs/tg_tagonly/bin/dhj.mhd_w_bcc.00004.bin: total tag mass fraction 6.457e-09
| log p [bar] | tag mass share | tag/mass in band | share at rho<2rho_fl | share 2..4 | share >4 |
|---|---|---|---|---|---|
| -9..-8 | 0.041 | 9.91e-01 | 0.988 | 0.001 | 0.011 |
| -8..-7 | 0.280 | 6.59e-01 | 0.717 | 0.133 | 0.150 |
| -7..-6 | 0.613 | 6.56e-01 | 0.794 | 0.078 | 0.128 |
| -6..-5 | 0.063 | 3.55e-02 | 0.629 | 0.036 | 0.334 |
| -5..-4 | 0.004 | 3.00e-04 | 0.050 | 0.018 | 0.932 |
| -4..-3 | 0.000 | 6.52e-07 | 0.000 | 0.001 | 0.999 |
| -3..-2 | 0.000 | 1.70e-09 | 0.000 | 0.000 | 1.000 |
| -2..-1 | 0.000 | 6.09e-36 | 0.000 | 0.000 | 1.000 |
```

**Caveat on tag_where's p bands.** tag_where bands by p = (gamma-1) e_int with gamma = 1.4728. With our general EOS this overestimates p by a factor ~3.7-4.0 at p < 1e-5 bar (median p_ideal/p_EOS of the cells in each band), so its bands sit ~0.6 dex high. Narrow per-cell p bands are also noisy: off vs tagonly scatters by up to 6.6e-3 in 0.1-dex bands, from cells moving between bands (gates_1101/cellband_*.txt). The per-RADIAL-SHELL comparison below is clean: off vs tagonly <= 1.5e-4 at p >= 1e-5 bar. It is also, I believe, the yardstick that gave your 6.4e-3.

## Per radial shell (horizontal sums at fixed r; A = fl, B = tagonly, C = T1 off arm at 0.5 rot; EOS p)
| r [cm] | <p>_B [bar] | M_fl/M_tagonly - 1 | tag/M tagonly | tag/M fl | M_off/M_tagonly - 1 |
|---|---|---|---|---|---|
| 1.36252e+10 | 2.500e-03 | -1.308e-04 | 1.620e-04 | 5.591e-05 | -6.727e-07 |
| 1.36831e+10 | 2.052e-03 | -1.450e-04 | 1.862e-04 | 7.516e-05 | +1.125e-06 |
| 1.37403e+10 | 1.698e-03 | -1.681e-04 | 2.221e-04 | 1.021e-04 | +1.732e-06 |
| 1.37985e+10 | 1.408e-03 | -1.959e-04 | 2.813e-04 | 1.363e-04 | +1.538e-06 |
| 1.38596e+10 | 1.162e-03 | -2.326e-04 | 3.583e-04 | 1.719e-04 | +1.642e-06 |
| 1.39261e+10 | 9.482e-04 | -2.815e-04 | 4.582e-04 | 2.241e-04 | +7.049e-07 |
| 1.40005e+10 | 7.612e-04 | -3.467e-04 | 5.934e-04 | 2.919e-04 | +2.629e-08 |
| 1.40854e+10 | 5.982e-04 | -4.271e-04 | 7.740e-04 | 3.921e-04 | +2.425e-06 |
| 1.41832e+10 | 4.580e-04 | -5.263e-04 | 1.004e-03 | 4.942e-04 | +2.684e-06 |
| 1.42958e+10 | 3.403e-04 | -6.563e-04 | 1.326e-03 | 6.405e-04 | +6.173e-06 |
| 1.44248e+10 | 2.447e-04 | -8.126e-04 | 1.828e-03 | 9.234e-04 | +2.742e-05 |
| 1.45710e+10 | 1.703e-04 | -1.039e-03 | 2.949e-03 | 1.688e-03 | +1.079e-05 |
| 1.47350e+10 | 1.150e-04 | -1.387e-03 | 5.530e-03 | 3.656e-03 | -9.169e-05 |
| 1.49164e+10 | 7.568e-05 | -1.855e-03 | 1.036e-02 | 7.857e-03 | -1.508e-04 |
| 1.51146e+10 | 4.892e-05 | -2.447e-03 | 1.886e-02 | 1.499e-02 | -6.727e-05 |
| 1.53286e+10 | 3.118e-05 | -3.194e-03 | 3.346e-02 | 2.709e-02 | +1.236e-04 |
| 1.55571e+10 | 1.965e-05 | -4.458e-03 | 5.513e-02 | 4.595e-02 | +1.037e-04 |
| 1.57986e+10 | 1.224e-05 | -6.388e-03 | 8.692e-02 | 7.180e-02 | -2.227e-05 |
| 1.60516e+10 | 7.538e-06 | -8.604e-03 | 1.343e-01 | 1.114e-01 | +3.688e-04 |
| 1.63145e+10 | 4.585e-06 | -1.234e-02 | 1.964e-01 | 1.658e-01 | +2.194e-04 |
| 1.65858e+10 | 2.793e-06 | -1.733e-02 | 2.710e-01 | 2.325e-01 | -6.608e-05 |
| 1.68640e+10 | 1.709e-06 | -2.264e-02 | 3.525e-01 | 3.053e-01 | +6.703e-04 |
| 1.71480e+10 | 1.040e-06 | -3.078e-02 | 4.397e-01 | 3.828e-01 | +4.751e-04 |
| 1.74367e+10 | 6.303e-07 | -4.229e-02 | 5.302e-01 | 4.659e-01 | -1.221e-03 |
| 1.77289e+10 | 3.818e-07 | -5.365e-02 | 6.194e-01 | 5.526e-01 | -1.001e-03 |
| 1.80240e+10 | 2.349e-07 | -6.449e-02 | 7.035e-01 | 6.346e-01 | +2.182e-03 |
| 1.83212e+10 | 1.492e-07 | -7.766e-02 | 7.815e-01 | 7.150e-01 | +1.728e-03 |

- At the 1.22e-5 bar shell: fl - tagonly = **-6.39e-3**, the same as your 6.4e-3. tagonly vs off is -2e-5 there, so tagonly is a clean "no drain" reference.
- The tag fraction there is **8.7e-2**, not << 6e-3, and fl removed 1.5e-2 of it (8.7e-2 -> 7.2e-2).

## Where in the shell: split by each cell's LOCAL p (shellsplit_fl_tagonly.txt)
| shell <p> [bar] | cells with local p < 1e-6 bar: share of shell mass / tag/M tagonly / mass change | 1e-6..1e-5 bar: share / tag / change | > 1e-5 bar: share / tag / change |
|---|---|---|---|
| 1.22e-5 | 0.150 / 0.58 / **-7.9e-3** | 0.293 / 1.2e-4 / +1.4e-3 | 0.557 / 0 / +7e-5 |
| 1.97e-5 | 0.104 / 0.53 / -5.4e-3 | 0.205 / 1.9e-4 / +9.5e-4 | 0.690 / 0 / -9e-6 |
| 1.70e-4 | 0.018 / 0.16 / -1.15e-3 | 0.033 / 9.8e-4 / +1.3e-4 | 0.949 / 0 / -1.6e-5 |
| 1.41e-3 | 0.002 / 0.15 / -2.0e-4 | 0.007 / 4.4e-3 / +2.3e-5 | 0.992 / 0 / -2.2e-5 |

(The mass changes are (M_fl - M_tagonly) of those cells over the whole shell mass.)

## Verdict: (A)
- The whole shell difference at the 1.1-1.2e-5 bar radius comes from the cells of that shell whose local pressure is < 1e-6 bar. These are the floor-supported, compact-column part of the sphere at that radius: 15 % of the shell mass, tag fraction 0.58.
- There the drain removes floor-made mass, as intended.
- Cells at local p > 1e-5 bar in the same shell change by <= 7e-5, the same as off vs tagonly (noise).
- So the floor-tracer default does not touch the gas the accuracy rule protects. A fixed-radius shell mean is the wrong yardstick here, because it mixes the p < 1e-6 bar reservoir into "1.2e-5 bar".
- Same pattern deeper: at the 1.4e-3 bar shell, -2.0e-4, all of it in the 0.2 % of the shell mass at local p < 1e-6 bar.
- No fix (margin 4, cap, rate) is needed for this. The USER decides the default.
