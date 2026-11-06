# NOTE 2026-11-06 viper -> Caltech: early KE in the real Athena++ 1-D BSG runs (reply to TASK-2026-11-06-caltech-viper-athpp-early-ke)

**Answer: no t=0 kick in real Athena++.** In Linhao's actual code (his `bsg_star.cpp` pgen, his inner heat-bath BC), KE rises smoothly and
monotonically from 0. There is no sharp peak in the first ~1000 s and no coherent outward pulse. At ~1000 s the KE is 9.5e39 (af1) /
8.7e40 (tc4) erg in a pi/3 sr wedge, 70-2000x below our `vstar_n` arms (6.7e42 / 1.9e43). Interior velocities are INWARD (infall at the
inner wall, r = 20-24 Rsun) and below 1.3 km/s at 1000 s.

## Arms (Athena++ e47fabb7 + his bsg_star.cpp; only output code was patched; hydro + radiation are his code unmodified)
| arm | nx1 | nmu | angle_flag | taucell | dt | notes |
|---|---|---|---|---|---|---|
| A_tc4 | 384 uniform 20-80 Rsun | 4 | 0 (Athena++ default) | 4 | ~212 s (cfl 0.3) | his unmodified code + BC; Jacobi error_limit 1e-6, nlimit 100 |
| A_af1 | 384 | 4 (nzeta 4) | 1 (spherical angular flux) | 4 | ~212 s | same, angle_flag 1 = closer to a multi-D spherical transport |

1-D spherical polar mesh, theta pi/3..2pi/3, phi 0..pi/3. The .hst volume IS the pi/3 sr wedge, so the numbers need no scaling.
No density seed (1-D). Units: rho 1e-9, T 1e5, l Rsun, mu 0.62 (t_unit 1.90e4 s, E_unit 4.51e36 erg).
These are the Raven audit decks (`/raven/ptmp/jinma/athpp_audit_1015/h2h/runs/A_*`). I re-ran them for 3 d on viper with output every
100 s (hst) / ~500 s (dumps), from a viper rebuild of the same source. The rebuild's .hst matches the Raven runs to all printed digits
(checked at 1055/1892/2727/3560/4391 s).
Data: `/viper/ptmp2/jinma/athpp_ke_1106/` (ana.txt, ana2.txt, scripts, outputs).

## KE (whole column; KE at r < 40 Rsun is 98-100 % of it), max |v_r| inside r < 45 Rsun, first-cell luminosity
L_* = 5.58e38 erg/s (IC). L_wall = 4 pi r^2 F_r1(lab) in the first active cell, r = 20.08 Rsun. Cell-centred: the face flux is not
in the dumps. Negative = into the bath.

| t | A_tc4 KE [erg] | max\|v_r\| [km/s] @ r | L_wall/L_* | A_af1 KE [erg] | max\|v_r\| @ r | L_wall/L_* |
|---|---|---|---|---|---|---|
| 0 | 0 | 0 | +1.0 | 0 | 0 | +1.0 |
| 212 s | 7.3e39 | | | 4.2e39 | | |
| 634 s | 3.8e40 | 0.95 (in) @ 20.4 | -48 | 9.1e39 | 0.52 (in) @ 24.3 | -38 |
| 1055 s | 8.7e40 | 1.22 (in) @ 22.1 | -67 | 9.5e39 | 0.59 (in) @ 22.1 | -49 |
| 2100 s | 2.9e41 | 2.0 (in) @ 22.0 | -94 | 1.6e40 | 0.72 (in) @ 22.0 | -63 |
| 3560 s | 7.5e41 | 3.2 (in) @ 21.8 | -106 | 3.1e40 | 0.92 (in) @ 21.8 | -59 |
| 2 h | 2.3e42 | 5.3 (in) @ 22.0 | -158 | 6.6e40 | 1.02 (in) @ 22.4 | -63 |
| 3 h | 4.0e42 | 7.6 (in) @ 22.7 | -222 | 1.08e41 | 1.41 (in) @ 23.2 | -71 |
| 6 h | 8.4e42 | 14.3 (in) @ 24.5 | -312 | 2.2e41 | 2.6 (in) @ 24.6 | -67 |
| 12 h | 9.1e42 | 28.8 (in) @ 27.9 | -298 | 2.3e41 | 4.3 (in) @ 28.5 | -77 |
| 1 d | 4.3e42 | 53 (in) @ 33.4 | -308 | 7.4e40 | 5.5 (in) @ 33.7 | -106 |
| 2 d | 2.5e42 | 72 (out) @ 43.7 | -161 | 7.0e41 | 4.5 (in) @ 20.7 | -264 |
| 3 d | 1.0e42 | 24 (in) @ 44.0 | -589 | 1.3e42 | 6.2 (in) @ 21.3 | -148 |

Every 500 s to 3 h: `ana.txt` (both arms grow smoothly at every sample). The global max |v| is ~35-42 km/s at the outer edge
(78-80 Rsun, the near-vacuum top), and it carries negligible KE.

## Reading
1. **No kick.** A_af1 rises to ~9e39 in the first ~600 s (3 steps), stays flat to ~1300 s, then grows smoothly (~linear, 1e41 at 3 h).
   A_tc4 grows smoothly from the first step. Neither has a peak-then-decay like our `vstar_n` arms.
2. What KE there is comes from the **bath draining** through his inner BC. L_wall drops from +1 L_* to -40..-70 L_* within the first
   ~600 s (the known Athena++ wall drain; see athpp_audit_1015 NOTE H2 / athpp_fix_1016). The gas falls INWARD onto the wall, and that
   infall front moves out (r 20 -> 34 Rsun by 1 d). In AthenaK the same BC feeds +2.5 L_* instead.
3. So the coherent radial kick right after t = 0 in our `vstar_n` arms is not in his 1-D code. A_af1 (closer to multi-D transport) at
   7e40 by ~2 h (af0: 2.3e42) fits his 3-D Fig. 2 level (1e40 -> 7e40 at 0.5 d). Our KE peak at ~1000 s (6.7e42 / 1.9e43) is
   70-2000x above his, and its sign is opposite (outward kick vs inward infall). That points to our v* mimic, not to the real v*
   ordering.
4. Caveat: this is 1-D and unseeded. His 3-D run may use a different nmu / angle setup. Neither changes the absence of a t=0 transient
   here.
