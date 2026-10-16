# NOTE 2026-10-16 raven -> caltech: accretor 3-D (+-20 deg) half-resolution benchmark on 4 A100

Reply to TASK-2026-10-16-raven-accretor-3d.md (08593683). Prerequisite (b), the 2-D v3 test: see
NOTE-2026-10-16-raven-accretor-bigdom.md. Up to t 0.025 it passes everything except one flag: dt drops to 2.4-3.3e-6 before
impact, driven by the thin atmosphere and not by the nozzle or the mask. Its impact part is still queued. The USER
authorised the 3-D benchmark meanwhile.

## Build
- accretor-bigdom-1015 @ 08625ed7 (code 049c1f39), full fresh nofma build (job 31213702), same recipe as the 2-D one.
  Binary `/raven/ptmp/jinma/accbig_1016/bin/athena_a100_ry_per_accretor_08625ed7_nofma`,
  **md5 5fcf9c42a05bc0e66b42160b307126f8**.
- 2-D bitwise gate on Raven (job 31213837, part 1): 08625ed7 vs efe33e6a, v3 2-D deck (+dt_report), fresh, 50 cycles,
  4 ranks. The hst (hydro and user), the bins (hydro_w and m1) and the cycle/dt/dt_report lines are all **IDENTICAL**.

## Benchmark run (job 31213837, gpudev, 1 node, 4 A100-SXM4-40GB, 4 ranks x 1 MeshBlock 420x114x256)
Deck: b3d_half.athinput with @HDIR@ replaced (dt_report = true is in the deck), plus `time/nlim=250 time/ndiag=10`,
and the run used `-t 00:10:00`. It ended on the wall limit at **cycle 245, t 3.451e-3** (stream on from t 0). The rst
files (7.1 + 10.8 GB) are deleted. Dir: /raven/ptmp/jinma/accbig_1016/b3d/run1.

| quantity | Raven 4 A100 | Caltech 4 H200 |
|---|---|---|
| steady s/cycle (cycles 140-230) | **1.642** (cycles 20-110: 2.9 incl. switch-on; init ~60 s) | 0.80 |
| dt and setting cell | 1.4096e-5 -> 1.4066e-5 flat. The rank-0 dt_report cell is an envelope cell at rest, r 8.786, phi 0.003 rad (as Caltech). | 1.408e-5, r 8.786 |
| GPU memory (nvidia-smi peak) | **24.98 GB per A100** (of 40; fits, 15 GB margin) | 24.6 GB/H200 |
| GPU-h per orbit = (0.473169/dt) x s/cycle x 4/3600 | (0.473169/1.4067e-5 = 33,637 cycles) x 1.642 x 4/3600 = **61.4 A100-h = 15.3 node-h per orbit** (pre-impact dt) | ~30 H200-h |
A100/H200 per cycle: 2.05x.

## Checks
| check | result | verdict |
|---|---|---|
| No FATAL / NON-CONVERGED / dt retries | None. Picard mean 2.67, max 5, NC 0, 0 retries. | PASS |
| Mnzf rate (col 25) = 254.843 within 0.5 % after the fill | 246.4 in the first interval (fill), then **254.934 -> 254.886** (t 0.0004-0.0034): +0.04 % -> +0.017 %. Min - Mnzf constant at 0.41955-0.41976. | PASS |
| Nozzle z-Gaussian | First nozzle cell r 14.5107, phi 0.03375 rad, t 3.45e-3: rho / (0.2493 exp(-z_f^2/(2 x 1.0054^2))) = **0.99998** at every theta with \|z_f\| < 3.02 (z_f = 14.47214 cos theta). Beyond: 2.2e-4 at \|z_f\| 3.04 (the edge cell, 8 % of the Gaussian), then 2.8e-6, 1.4e-7, 3e-8, 1e-8 = floor. | PASS |
| Mask quiet at all theta | The 3-D lobe (exact Phi(r,theta,phi) < Phi_L1, x > d_L1; 345,268 cells) spans theta 76.8-103.5 deg. Rho = 1e-8 (floor) at every theta, \|v\| < 0.034 km/s, e/rho <= 1490. Outside the lobe at r > 16.3: rho <= 1.09e-8, \|v\| < 3.9 km/s. Mmsk (col 9) 6.8e-6 -> 9.0e-5, linear and smooth (2.6e-2 per time unit). | PASS |
| Donor irradiation at all theta | Eirr (col 26) grows linearly, 4.9e4 -> 1.23e6 (3.7e8 per time unit). Ltop 1.44e5 -> 1.90e5, as on Caltech. The m1 bin lobe E at output time has median 0.68 (theta 80/90/100 deg: 0.62/0.67/0.62), so it is present and symmetric. That is only ~10 % of a T_d^4 = 6.787, because the bin is written after the implicit radiation step that follows the reset (2-D: 16-50 %). The held value itself is right (log: E = 6.78662). | PASS (present at all theta; output E is the post-transport value) |
| Back-sink row artefact (also in 2-D) | 1876 cells at r 14.665 (the first back-sink row behind the nozzle), theta 83.8-96.2, phi -9..+13 deg: floor rho with \|v\| up to 1.05e6 km/s in the hydro_w bin. Same as in 2-D: reset every stage, not seen by dt. | flag |

## Cost estimates for a longer 3-D test on Raven (gpu partition, 4 A100/node, 24 h links)
- **Half resolution (49.0 M cells, 1 node)**: 15.3 node-h (61 A100-h) per orbit at the pre-impact dt, i.e. one orbit
  per ~15 h on one node: **1 orbit fits in one 24-h link**. If dt drops by 2-3x as in the 2-D v3 run (thin atmosphere,
  impact, retries), expect 30-45 node-h per orbit (2 links).
- **Full resolution (840 x 228 x 2048 = 392 M cells)**: 8x cells and 2x cycles (dt ~ dx) = 16x, about **245 node-h (980
  A100-h) per orbit** at the pre-impact dt, assuming ideal weak scaling. Memory ~2 GB per M cells (25 GB per 12.3 M), so
  >= 32 A100 = 8 nodes, i.e. ~31 h wall per orbit (2 links). BiCGSTAB iteration growth with resolution and MPI halos
  will add to this. With the 2-3x dt loss: 500-750 node-h per orbit.
- For comparison, the user's planning numbers were 2500 H200-h per orbit (full, old domain, x1.3 big domain) and 300
  (half). Our pre-impact half-resolution figure is 61 A100-h ~ 30 H200-h, matching Caltech.

## Recommendations before a longer 3-D run
1. Make the physics 3-D (049c1f39 is plumbing only, as the design doc 14.1 says):
   - **gravity** (no vertical gravity now: the star is a theta-invariant cylinder and stream/atmosphere gas is not
     bound in z);
   - **Coriolis/centrifugal** with the full 3-D Omega x v and the cylindrical-radius centrifugal term;
   - **IC** (column replicated in theta -> needs a spherical/Roche-equipotential envelope);
   - **WB reference and ambient** (theta-invariant now).
2. **Theta walls**: reflecting at +-20 deg reflects any vertically spreading stream and atmospheric outflow back into
   the domain. Consider outflow (or a sponge layer) at the walls once the stream spreads in z.
3. Fix the r 14.6x **back-sink row** momentum (zero v in sink cells after the radiation/floor step).
4. Resolve the **2-D v3 dt issue first**: thin atmosphere gas at r 9.1-9.6 is accelerated to ~940 km/s on both sides of
   the star (phi ~90/270 deg), plausibly by the donor light. That halves dt before impact, and in 3-D the same
   atmosphere covers all theta.
5. **dt_report**: reduce over ranks (MPI MINLOC); at present it reports rank 0's local minimum only
   (ry_per_accretor.cpp:1367).
6. Decide whether the back-sink sector should radiate at T_d (it doubles Eirr in 2-D).
