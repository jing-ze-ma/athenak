# NOTE 2026-11-09 ORION -> viper: raycheap13 -- per-stage frozen operator pieces (1-D / 2-D pulse, cs transport-only)

Answer to `TASK-2026-11-09-orion-raycheap13.md`. No code changes. Binary /orion/ptmp/jinma/raycheap11_1109/athena_addb642d0.

**Summary**
- **Part 3 (cs transport-only, T1 recipe): the lag keys are inert; the answer stays 1st order or worse.** S1 (recon
  picard), S2 (closure pass) and S3 (both) are bitwise identical to S0 in m1_e at c 0.15 and 0.01875: d 8.979e-04
  8.028e-04 5.848e-04, p 0.16 / 0.46, same Picard counts. (The raycheap14 rerun of S0/S3 on dea3e3732 gives the same
  numbers, S0 bitwise across the two builds.) S4 (+ thin_switch none, thin_relax 0) changes the solution (max|dE|/max E
  0.5-0.6% vs S0) and is worse: d 2.9e-03 2.4e-03 1.8e-03, p 0.30 / 0.41. S5 (+ vet_col lag, frozen) moves it 15% and
  does not converge (d grows with refinement, p -0.62 / -0.29), at ~0.66-0.82x wall.
- **2-D Cartesian diffusion pulse: the multi-D defaults are 1st order, and lifting BOTH freezes restores 2nd order.**
  D0 (defaults, recon_lag step + closure_lag step) p 1.22 / 1.12 / 1.06; **D3 (recon_lag picard + closure_lag pass) p 1.96 /
  1.98 / 1.99**. Either toggle ALONE does nothing: D1 (recon picard) and D2 (closure pass) are bitwise identical to D0 at the
  coarsest and finest levels. In D0-D2 the Picard loop takes exactly 3 passes every solve (the frozen operator is
  linear: nothing to iterate); D3 iterates 7-10 passes (mean 9.7 -> 7.0), never non-converged (cap 200 not hit). D3 costs
  ~2x the wall of D0.
- **1-D pulse: P0-P5 are all bitwise identical (p 1.95 / 1.98 / 1.99).** In the 1-D Cartesian path (implicit_flux central,
  faces not `all`) neither implicit_recon_lag nor implicit_closure_lag changes anything, also with implicit_recon = plm_dc
  (extra arms P4/P5). So 1-D cannot show the freeze; 2-D can.
- **Refused:** 2-D with `implicit_flux_faces = all` needs `implicit_flux = berthon | blend` (set berthon) and
  `takes implicit_recon = dc only` (the planned plm_dc 2-D arms D4/D5 were refused at start and dropped).

## Decks
- Part 1: the raycheap12 1-D pulse decks (dt_fixed 0.24/2^L, tlim 2.4, dt_retry false, transport implicit, implicit_tol
  1e-12, maxit 200, hesdirk2, closure m1; 1 rank). The deck itself names `implicit_recon_lag = picard`.
- Part 2: same deck made 2-D: nx2 = nx1 = 128 (x2 in [-0.5, 0.5], periodic), one meshblock, 1 rank;
  implicit_allow_multid = true, implicit_flux_faces = all, implicit_flux = berthon, implicit_solver = bicgstab, and the
  deck's `implicit_recon_lag` line REMOVED so D0 gets the multi-D default (log: closure_lag=step in D0/D1). The pgen pulse
  on a 2-D mesh is the 2-D Gaussian exp(-r^2/w^2), periodic, centre (0.5, 0).
- d = |E_k - E_k+1|_1/|E_finest|_1 between successive levels; all runs rc 0, FATAL 0, NOT ADMISSIBLE 0, DT_RETRY 0, and
  every run ends at tlim (cycles 10/20/40/81/161, one round-off step at L0/L3/L4 as in raycheap12).

## Parts 1-2 (hesdirk2, closure m1, dt 0.24 ... 0.015)

| arm | d (4 pairs) | p (3) | Picard mean/max/non-converged per level | wall (s) |
|---|---|---|---|---|
| P0 deck (recon_lag picard, closure_lag default) | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 4.89/5/0 4.54/5/0 4.22/5/0 3.98/5/0 3.99/5/0 | 1/1/1/0/1 |
| P1 recon_lag=step | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 4.89/5/0 4.54/5/0 4.22/5/0 3.98/5/0 3.99/5/0 | 1/0/1/1/1 |
| P2 closure_lag=step | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 4.89/5/0 4.54/5/0 4.22/5/0 3.98/5/0 3.99/5/0 | 0/1/1/1/1 |
| P3 both step | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 4.89/5/0 4.54/5/0 4.22/5/0 3.98/5/0 3.99/5/0 | 1/1/1/1/1 |
| P4 (extra) implicit_recon=plm_dc, recon_lag picard | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 4.89/5/0 4.54/5/0 4.22/5/0 3.98/5/0 3.99/5/0 | 1/0/1/1/1 |
| P5 (extra) plm_dc, recon_lag step | 5.573e-04 1.438e-04 3.654e-05 9.216e-06 | 1.95 1.98 1.99 | 4.89/5/0 4.54/5/0 4.22/5/0 3.98/5/0 3.99/5/0 | 0/1/1/1/1 |
| **D0** multi-D defaults (step / step) | 2.541e-03 1.094e-03 5.039e-04 2.413e-04 | 1.22 1.12 1.06 | 3.00/3/0 3.00/3/0 3.00/3/0 2.99/3/0 2.99/3/0 | 4/5/8/13/26 |
| D1 recon_lag=picard | 2.541e-03 1.094e-03 5.039e-04 2.413e-04 | 1.22 1.12 1.06 | 3.00/3/0 3.00/3/0 3.00/3/0 2.99/3/0 2.99/3/0 | 3/5/9/15/27 |
| D2 closure_lag=pass | 2.541e-03 1.094e-03 5.039e-04 2.413e-04 | 1.22 1.12 1.06 | 3.00/3/0 3.00/3/0 3.00/3/0 2.99/3/0 2.99/3/0 | 4/5/7/14/24 |
| **D3** recon_lag=picard + closure_lag=pass | 8.844e-04 2.279e-04 5.789e-05 1.459e-05 | 1.96 1.98 1.99 | 9.74/12/0 9.05/11/0 8.37/10/0 7.84/10/0 6.99/9/0 | 8/12/19/32/56 |

## Part 3 (cs Ar4 transport-only, S0-S5; job rc13cs, 3 p.exclusive nodes, 24 ranks per run, all concurrent-ish)

| arm | d (3 pairs) | p (2) | Picard mean/max/non-conv (c 0.15 ... 0.01875) | wall (s) | vs S0 (m1_e) |
|---|---|---|---|---|---|
| S0 T1 recipe (step / step) | 8.979e-04 8.028e-04 5.848e-04 | 0.16 0.46 | 13.96/24/0 10.76/19/0 7.93/16/0 6.04/12/0 | 308/470/694/1125 | - |
| S1 recon_lag=picard | 8.979e-04 8.028e-04 5.848e-04 | 0.16 0.46 | same | 296/444/673/1095 | bitwise |
| S2 closure_lag=pass | 8.979e-04 8.028e-04 5.848e-04 | 0.16 0.46 | same | 304/454/707/1128 | bitwise |
| S3 both | 8.979e-04 8.028e-04 5.848e-04 | 0.16 0.46 | same | 297/481/732/1178 | bitwise |
| S4 S3 + thin_switch none + closure_thin_relax 0 | 2.928e-03 2.378e-03 1.791e-03 | 0.30 0.41 | 13.96/24/0 10.77/19/0 7.93/16/0 6.04/12/0 | 311/491/734/1183 | 4.7e-03 / 6.2e-03 |
| S5 S4 + time2_vet_col lag, vet_col_every 1e6 | 5.393e-04 8.278e-04 1.011e-03 | -0.62 -0.29 | 12.45/24/0 9.25/19/0 6.53/16/0 4.73/12/0 | 252/360/517/743 | 1.5e-01 / 1.5e-01 |

All rc 0, FATAL 0, NOT ADMISSIBLE 0, DT_RETRY 0; all end at tlim 20 (cycles 68/136/271/542). Picard never hit the cap
(max 24 of 200), no limit cycles, implicit_maxit not raised. Logs: closure_lag=step in S0/S1, pass in S2-S5.
On this deck (unlike the 2-D pulse) the per-stage freezes are not what holds the order down; the Picard loop already
iterates 6-14 passes in S0 and the frozen pieces do not move the converged answer.
