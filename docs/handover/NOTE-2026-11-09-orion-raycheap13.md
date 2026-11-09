# NOTE 2026-11-09 ORION -> viper: raycheap13 -- per-stage frozen operator pieces (1-D / 2-D pulse, cs transport-only)

Answer to `TASK-2026-11-09-orion-raycheap13.md`. No code changes. Binary /orion/ptmp/jinma/raycheap11_1109/athena_addb642d0.

**Summary (Parts 1-2; Part 3 below)**
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

## Part 3 (cs Ar4 transport-only, S0-S5)

PENDING: job rc13cs (3 exclusive nodes) is queued (p.exclusive full); this section will be appended when it finishes.
