# NOTE DeltaAI -> raven: AG Car A he_ic_rad = sph1d / re_atm (TASK-2026-10-11-deltaai-icgen) -- DONE

Picked up 2026-10-10 23:38 CDT; result 23:58 CDT. Job **3357616** (ghx4-interactive gh095, 1 node x 4 GH200, 1 rank/GPU,
MPICH_GPU_SUPPORT_ENABLED=1, plain srun), wall 5 min 54 s. Binary icgen-1011 3d8861c7 (`athena_hes_gpu_3d8861c794a9`, md5 df2e42fa04b237143d0813f4115abc92).
Inputs md5 as listed; every key in a per-run deck copy (nlim 120, ndiag 1, implicit_cfl 1251.3, output2/3/7 dcycle 10,
output6 dcycle 100000, no rst, rad_m1/ic_rad_diag 2 + the arm keys). **All 6 runs rc 0, FATAL 0, NON-CONVERGED 0.**

## Answer to the main question
**gha (re_atm + hse) is at the g0 level, not the gh level**: Picard 3.53 vs g0 3.82 vs gh 5.33 (end-of-run means); per window,
gha 4.20 / 3.20 / 3.10 / 3.15 / 3.55 / 4.00 vs g0 4.90 / 3.30 / 3.15 / 3.65 / 3.90 / 4.00, while gh goes to 6.00 from cycle 61.
The floor-gas infall confirms the cause: shell-mean v_r at 2.0 R_ph at cycle 120 is -17.4 km/s for gha (g0 -17.3) vs **-51.0 km/s for gh**.
gha keeps the gh photosphere (static to 0.47 km/s, rho drift 9 % at cycle 120; g0 2.9 km/s, 96 %) and the gh L_top (+1.0 % step 1).

## Per run
s/cycle = median of per-cycle wall excluding rotation cycles and the last; rotation = mean of rotation cycles; all = mean.

| arm | Picard mean / max (solves) | s/cycle median (rotation, all) | sph1d solve | ic_rad run_c0 max\|<E>/<J>-1\| (top) |
|---|---|---|---|---|
| g0 | 3.82 / 9 (120) | 0.207 (0.49, 0.275) | - | 1.061 (1.0796) |
| gha | 3.53 / 7 (120) | 0.201 (0.47, 0.267) | 5.21 s (331 VET passes, 29 outer) | 0.0177 (1.0177) |
| gha2 | 3.53 / 7 (120) | 0.200 (0.46, 0.266) | 5.20 s | 0.0177 (1.0177) |
| ghafr | 3.37 / 7 (123) | 0.206 (0.88, 0.300) | 5.22 s + ic_rad 2.19 s | 0.0289 (1.0289) |
| gra | 3.52 / 6 (120) | 0.208 (0.94, 0.316) | 0.29 s (18 passes, 1 outer) | 0.0326 (1.0326) |
| gh | 5.33 / 6 (120) | 0.229 (1.00, 0.339) | 6.10 s (340 passes) | 0.0183 (1.0183) |

- ghafr `ic_rad: done, 2.186 s (fs 0.422 s, relax 3 steps 0.689 s), L_top = 2.8430856815e+38 (L/L_IC = 0.999997)`.
- The rotation-cycle means vary 0.46-1.0 s between arms with the same code path; single rotations, so treat as noise (median steady cost is 0.200-0.208 for all but gh 0.229).
- Setup lines, e.g. gha: `he_ic_rad = sph1d (gas re_atm, hse on): r_b = 2.43206e+13 (tau 3002.59), blend tau 30 .. 3000, 1343 solve points + 64 core rays,
  331 VET passes, 29 outer, g_top = H/J 0.97864, f_top = K/J 0.957965, 5.21367 s`; R_ph (tau_R 2/3, run kappa) 2.70169e13 -> 2.70237e13.

## Picard (hst) mean per 20-cycle window
| arm | 1-20 | 21-40 | 41-60 | 61-80 | 81-100 | 101-120 |
|---|---|---|---|---|---|---|
| g0 | 4.90 | 3.30 | 3.15 | 3.65 | 3.90 | 4.00 |
| gha | 4.20 | 3.20 | 3.10 | 3.15 | 3.55 | 4.00 |
| gha2 | 4.20 | 3.20 | 3.10 | 3.15 | 3.55 | 4.00 |
| ghafr | 3.55 | 3.30 | 3.15 | 3.05 | 3.10 | 3.60 |
| gra | 3.95 | 3.25 | 3.10 | 3.25 | 3.60 | 4.00 |
| gh | 4.70 | 4.60 | 4.65 | **6.00** | **6.00** | **6.00** |

## L_top / L_top(g0, cycle 0)
| cycle | g0 | gha | gha2 | ghafr | gra | gh |
|---|---|---|---|---|---|---|
| 0 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 | 1.00000 |
| 1 | **1.53558** | 1.01036 | 1.01036 | **0.99997** | 1.00995 | 1.01034 |
| 2 | 1.44485 | 1.00550 | 1.00550 | 0.99995 | 1.00481 | 1.00485 |
| 3 | 1.24912 | 1.00221 | 1.00221 | 0.99990 | 1.00158 | 1.00187 |
| 4 | 1.10827 | 1.00078 | 1.00078 | 0.99986 | 1.00037 | 1.00064 |
| 5 | 1.04033 | 1.00020 | 1.00020 | 0.99983 | 1.00000 | 1.00016 |
| 6 | 1.01441 | 0.99998 | 0.99998 | 0.99982 | 0.99990 | 0.99996 |
| 7 | 1.00573 | 0.99988 | 0.99988 | 0.99981 | 0.99987 | 0.99988 |
| 8 | 1.00289 | 0.99984 | 0.99984 | 0.99981 | 0.99987 | 0.99984 |
| 9 | 1.00185 | 0.99982 | 0.99982 | 0.99981 | 0.99988 | 0.99982 |
| 10 | 1.00138 | 0.99981 | 0.99981 | 0.99981 | 0.99988 | 0.99981 |
| 20 | 1.00051 | 0.99983 | 0.99983 | 0.99985 | 0.99999 | 0.99983 |
| 60 | 1.00035 | 0.99996 | 0.99996 | 0.99997 | 1.00020 | 0.99996 |
| 120 | 1.00046 | 1.00004 | 1.00004 | 1.00004 | 1.00041 | 1.00003 |

ghafr (sph1d + fs_relax) is the only flat one from cycle 0 (min 0.99981).

## Photosphere hydrostatics (0.97-1.028 R_ph, 50 shells) and floor-gas infall
max |v_r| (cellwise, km/s) / max |<rho>/<rho>(cycle 0) - 1| (shell means) / shell-mean v_r at 2.0 R_ph (km/s, shell at 1.987 R_ph):

| arm | cycle 10 | cycle 30 | cycle 60 | cycle 120 |
|---|---|---|---|---|
| g0 | 0.987 / 0.052 / -1.19 | 1.995 / 0.164 / -4.10 | 2.416 / 0.390 / -8.47 | 2.922 / 0.955 / -17.30 |
| gha | 0.031 / 0.0009 / -1.34 | 0.090 / 0.0042 / -4.24 | 0.158 / 0.0116 / -8.62 | 0.472 / 0.089 / -17.45 |
| gha2 | (= gha) | | | |
| ghafr | 0.036 / 0.0009 / -1.45 | 0.094 / 0.0043 / -4.35 | 0.160 / 0.0120 / -8.73 | 0.459 / 0.087 / -17.56 |
| gra | 0.922 / 0.051 / -1.34 | 1.937 / 0.160 / -4.24 | 2.405 / 0.387 / -8.62 | 2.837 / 0.929 / -17.45 |
| gh | 0.031 / 0.0009 / **-4.01** | 0.089 / 0.0043 / **-12.16** | 0.156 / 0.0115 / **-24.67** | 0.478 / 0.090 / **-51.01** |

- The static photosphere comes from the HSE re-integration: gra (re_atm, no hse) behaves like g0 there.
- re_atm removes the 3x faster floor infall of gh (gas = re): at 2 R_ph gha/ghafr/gra follow g0.
- The photosphere is static over 60 cycles (Raven's window) but drifts by cycle 120 for all sph1d+hse arms (rho 9 %, v_r 0.47 km/s), 10x below g0.

## gha vs gha2
**Bitwise**: 44/44 bin files md5-identical (and array-identical), user.hst and hydro.hst identical.

## GPU-h
5.9 min x 4 GH200 = **0.39 GPU-h**.

## Files
In git: docs/handover/deltaai-icgen-1011/ (job.sub, analyse.py, analysis.txt). On DeltaAI: /work/nvme/bivj/jma20/icgen_1011/runs/3357616/.
