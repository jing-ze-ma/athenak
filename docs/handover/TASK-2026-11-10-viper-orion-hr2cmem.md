# TASK 2026-11-10 viper -> ORION: hr2cmem -- Cartesian (vet_sc) order ladders of the half-range blend with implicit_hr_cmem

User GO 11:40 PST. CPU work. Rotation off. Answer as NOTE hr2cmem on bsg-files-1009.

## Why
On a beam-free plane-parallel Cartesian pulse with the production Cartesian keys (vet_sc + halfrange + fse thin
switch + hc), the TIME order collapses whenever the face weight sits at a fixed 0 < w < 1: the E error between dt
levels is flat (~5.5e-4). w = 0 (central), w = 1 everywhere (`w1`), and DTC (w -> 0 as dt -> 0) are second order.
The thin-switch corrections (hc, js, thin_freeze) are NOT the cause. Hypothesis: the face-flux memory f0x stores the
BLENDED applied flux, so the central branch relaxes toward the half-range flux at a rate ~ 1/dt. `implicit_hr_cmem`
(F_c in f0x, applied flux in fax) is the candidate fix; we need its C and T ladders.

Raven A100 numbers that motivate it (binary 2aeefd3d, CART_1D=1, pk12, T 1..16 at n = 128, p of the 3 finest pairs):

| arm | E p | F1 p |
|---|---|---|
| cen | 2.02 2.03 2.03 | 1.82 1.86 1.88 |
| base | 0.87 0.05 0.16 | 0.14 0.02 0.07 |
| w1 | 1.83 1.69 1.52 | 1.88 1.70 1.37 |
| dtc | 2.06 2.07 2.04 | 1.81 1.86 1.89 |
| nofse / nocj / nocorr / norecon / frz0 / plm / plmh / gjump | 0.6-1.0, then 0.02-0.16 | ~0 |

Space (C, 32..256/512) of the same plane-parallel pulse: cen E 1.91 1.77 1.47 (512); base 1.76 0.92; plm 1.69 0.99
(earlier decks; numbers with the beam-free decks are being re-measured).

## Code
Branch **fork/hr2prod-1110, commit 2aeefd3df** (fork/hr2nd-1110 b0cbe4fd5 + fork/vetgdup-1110 a19a3d835 + two
opt-in keys: `implicit_hr_wscale_x1/_lat` (diagnostic) and `implicit_thin_gface = max | mean | jump`). Keys absent:
bitwise (CPU, cs deck); cs np 1/3/4 bitwise with the new keys on. `implicit_hr_cmem` comes from raycheap18 (bb58d80c2).
CPU build: `cmake -B build -D Athena_ENABLE_MPI=ON -D CMAKE_BUILD_TYPE=Release && make -C build -j` (PROBLEM none;
the pgen is the built-in rad_m1_beam / thick_pulse).

## Inputs (no /viper or /raven paths)
`docs/handover/scripts/hr2cmem/` on this branch: `cart.py` (deck maker + evaluator), `cb_FIN3_r4.athinput` (the
template: only its `<rad_m1>` block is used), `bin_convert.py` (reader), `run.sh` (driver).

Deck (cart.py `deck`, with `CART_1D=1`): plane-parallel Gaussian E pulse (`m1_test = thick_pulse`, `pulse_1d = true`,
width 0.25, bg 1e-3) in a static pure scatterer, box x1 in [-1, 1] (outflow; the vet_sc sweep axis), nx2 = 4 periodic
(x2 in [-1, 1]), one MeshBlock, c = 1e3, tlim = 0.3/c. `<rad_m1>` = the cb_FIN3_r4 block with kappa_p = kappa_f = 0,
kappa_s = 12.8 (pk12; pk1: 1.28, pk128: 128), time2_k1_redist 0, implicit_tol 1e-10, tol_rel 0, rotation OFF
(vet_gd_rotate_every/dt 0, vet_rotate_avg 1, _ema 0), vet_nmu 8 x vet_nphi 32 (vet_ring removed), and **no
implicit_beam\*** keys (the cb deck's boundary beams are deleted). C: nx1 = n, dt = 0.5 dx/c. T: n = 128, dt =
4 dx/(c l). Output: one m1 bin at tlim.

Arms (`~`-joined MOD keys in cart.py): `cen`, `base`, `w1`, `dtc`, `cm`, `cm~gjump`, `cm~plmh`, `cm~plmh~gjump`,
`nofse~cm`, `cm~dtc`.
- cm = implicit_hr_cmem true; plmh = thin_corr none + thin_recon dc + implicit_hr_recon plm, q = h,
  lim smooth, implicit_hr_face cell (q = h is refused with the default sf face ratios of vet_sc); gjump =
  implicit_thin_gface jump; w1 = thin_h0 -2, h1 -1 (w = 1 on every face); dtc = implicit_blend_dtc 1e-6.

## Runs
```
export BIN=<athena> CART_H=$PWD/cart_runs CART_1D=1
for a in cen base w1 dtc cm cm~gjump cm~plmh cm~plmh~gjump nofse~cm cm~dtc; do
  bash run.sh pk12 $a C "32 64 128 256 512"
  bash run.sh pk12 $a T "1 2 4 8 16"
done
```
Then pk1 and pk128 for `cen base cm cm~plmh~gjump` (C and T). Smoke each arm first (nlim 2) if you prefer.

## Analysis
`cart.py eval CASE ARM MODE LEVELS` prints, for m1_e, m1_f1, m1_f2, d_l = mean|q_l - R q_{l+1}| / mean|q_l -
<q_l>| (R = 2:1 restriction along x1 for C, identity for T) and p = log2(d_l/d_{l+1}). Report per arm one row:
`| arm | C E p | C F1 p | T E p | T F1 p |` (all pairs), plus rc/FATAL/NON-CONVERGED counts from the logs.

## Questions
1. Does cm restore T order on its own (base -> cm), and with gjump / plmh?
2. Is cm second order in SPACE where plmh is on (cm~plmh~gjump C)?
3. Does pk128 (diffusive) or pk1 (streaming) behave differently?
