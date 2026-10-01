# ck_spherical face-closure variants vs the exact reference (cksph_test_0930)

Branch `ck-sph-closure-fix` (worktree `/viper/ptmp2/jinma/cksph_test_0930/wt`), not pushed, not
merged.  All numbers below come from serial CPU runs on the login node.  They use the binaries
`/viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_4882f92d` (modes 0-4, ck_sph_top) and `_4b1da214`
(mode 5); `_4b1da214` has the same source tree as the final HEAD `847b74ff`.  The run dirs are
`runs/<case>_<arm>/`.  The tables are made by `ana/table.py` / `ana/summ.py`, and `runv.sh` is
the matrix driver.  The cases are those of RESULTS.md (static, `rt_test_freeze`, no restart, so
no rotation applies), plus:
- `grad_*`: the W121 grid with T = 2500 (1 - 0.3 sin(pi xi/2)).  It measures the deep diffusive
  flux, at faces where every (b,g) has tau-from-top > 50.
- `nw10_*`: W121 at 2000 K, `rt_test_dt = 10` s, 10 cycles, the full ck_implicit Newton.

## Recommendation

**Accuracy: `ck_sph_face = 5` (variable Eddington factor, VEF) + `ck_sph_top = 1`.**  Reason: it
is the only variant with L within 0.4 % on both grids, the thin-shell deposit within 2.4 %, and
the deep flux equal to the Eddington face (1e-8).  Its one weakness is an **unresolved** sharp
photosphere (the st slab, tau per cell 5e6), where it gives +6.9 %.  In the production grids the
photosphere of every g is resolved, so this case does not arise there.

**If the formal-solution cost is not acceptable: `ck_sph_face = 4` (blend face + cone path) +
`ck_sph_top = 1`.**  This is the cheap geometric flavour.  It is exact on every sharp case and
matches the deep flux, but L is +2.6-3.0 % high on W121 (+4.2 % on the old grid).  The deposit
is within +0.2..+5.0 %, except the top cell of the old 3000 K grid (1.17).

## The variants (one key, `problem/ck_sph_face`; 0 stays bitwise the current code)

All modes keep A x (u - d) (x the chain's flux weight) single-valued at every face, so the
deposit telescopes.  They are implemented in the tm (FRM = 1) chain kernel, explicit and JAC
(ck_implicit).  They are refused with ck_impl_lin / ck_impl_jac_lin.

| mode | name | what it does |
|---|---|---|
| 0 | Eddington (default) | S and A D continuous (K = J/3) |
| 1 | dilute | = ck_sph_dilute prototype: d continuous, A_a u_a = A_b u_b + (A_a - A_b) d_a |
| 2 | V2 thick/thin blend | family e (u_a - u_b) + (d_a - d_b) = 0, e = 1 - exp(-dtau_f/mu), dtau_f = the two half cells at the face; e = 1 Eddington, e = 0 dilute |
| 3 | V3a dilute + cone path | up crossing of each half cell absorbs along dtau/(mu (1 + cos th_c)); emission and down ray unchanged; cos^2 th_c = T^4 R (R = relation reflection above the face, T = plain half-cell transmission: 1 - R = A_ph/A through a transparent shell) |
| 4 | V3a blend + cone path | 2 + 3 |
| 5 | V3b VEF (Auer 1971 / Mihalas ch. 7) | f = K/J per (b,g) from a p-z formal solution through the column (8 wall rays + the wall-edge ray + one tangent ray per node, S = B, vacuum top), Newton-lagged (`ck_vef_every` cycles).  Both chains of a pair at mu_eff = sqrt f in the cell (then the pair is exactly the moment system with K = f J).  Faces: Phi f J and A mu_eff (u - d) continuous, Phi = q r^2, ln Phi = int (3f-1)/(r f) dr.  Top: Auer's H = h J as a reflection d = zeta u.  Wall: Marshak-type J + beta H = B (1/2 + beta/4), beta = J-/H- of the formal incoming field |

`problem/ck_sph_top` (V4, independent of the face): 0 = the ghost cell (historical, the hot IC
ghost), 1 = d_top = 0, 2 = kappa and B at the top active cell (column depth still from the ghost
pressure).  With mode 2 the tridiagonal gets the d_top(B_ie) entry.

## Results

Columns:
- **sharp L**: L_code/L_exact for tr x = 1.2 / 1.7 / 3 and st x = 1.2 / 1.7.
- **W121 / old L**: total L_code/L_exact, then the worst band 0-9 and band 10.
- **deposit**: Q_code/Q_exact over the thin shell (tau_v(max g) < 40), shown as the min-max
  without the top cell, with [the top cell] in brackets.
- **deep flux**: mean (max) of F/F_Eddington - 1 over the 22 deep faces of `grad`.
- **Newton**: max passes / total nonconv / max budget residual over 10 cycles at 10 s.

| arm | sharp L | W121 2000 K L | W121 3000 K L | old 2000 K L | old 3000 K L | dep W121 2000 | dep W121 3000 | dep old 2000 | dep old 3000 | deep flux | Newton |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 Eddington | 1.0909 / 1.2590 / 1.5000 / 1.0893 / 1.2541 | 1.2243 (b7 1.239; b10 1.161) | 1.2481 (b6 1.253; b10 1.214) | 1.2360 (b7 1.287; b10 -41.9) | 1.3381 (b6 1.366; b10 1.135) | 0.849-1.078 [0.80] | 0.900-1.049 [0.83] | 0.054-0.948 [-1.27] | 0.619-0.848 [-0.75] | 0 | 3/0/2.8e-04 |
| 1 dilute | 1.0000 / 0.9998 / 1.0000 / 0.9939 / 0.9928 | 0.9660 (b8 0.952; b10 0.914) | 0.9809 (b0 0.979; b10 0.967) | 0.8978 (b9 0.863; b10 -41.6) | 0.9363 (b0 0.928; b10 0.830) | 0.865-0.924 [0.85] | 0.899-0.944 [0.87] | 0.103-0.883 [-1.22] | 0.647-0.830 [-0.72] | -2.2e-03 (3.4e-03) | 3/0/2.8e-04 |
| 2 blend | 1.0000 / 0.9998 / 1.0000 / 1.0000 / 1.0000 | 0.9683 (b8 0.955; b10 0.916) | 0.9830 (b7 0.982; b10 0.971) | 0.9019 (b9 0.867; b10 -41.8) | 0.9403 (b0 0.932; b10 0.835) | 0.860-0.923 [0.85] | 0.893-0.948 [0.86] | 0.097-0.882 [-1.23] | 0.643-0.831 [-0.72] | -1e-12 (1.2e-10) | 3/0/2.7e-04 |
| 3 dilute+cone | 1.0000 / 0.9998 / 1.0000 / 0.9939 / 0.9928 | 1.0333 (b7 1.035; b10 0.979) | 1.0302 (b7 1.032; b10 1.028) | 1.0423 (b4 1.046; b10 -41.5) | 1.0459 (b6 1.048; b10 0.962) | 0.994-1.033 [0.98] | 1.004-1.044 [1.01] | 0.196-0.998 [-1.14] | 0.723-0.948 [-0.66] | -2.2e-03 (3.4e-03) | 3/0/2.7e-04 |
| 4 blend+cone | 1.0000 / 0.9998 / 1.0000 / 1.0000 / 1.0000 | 1.0295 (b6 1.031; b10 0.974) | 1.0261 (b7 1.028; b10 1.025) | 1.0399 (b4 1.044; b10 -41.6) | 1.0417 (b6 1.043; b10 0.959) | 0.991-1.033 [0.98] | 1.008-1.047 [1.00] | 0.190-0.999 [-1.14] | 0.720-0.949 [-0.67] | -1e-12 (1.2e-10) | 3/0/2.7e-04 |
| 0 + top 1 | 1.0909 / 1.2592 / 1.5000 / 1.0893 / 1.2541 | 1.2243 (b7 1.239; b10 1.216) | 1.2481 (b6 1.253; b10 1.214) | 1.2386 (b7 1.287; b10 1.218) | 1.3387 (b6 1.366; b10 1.219) | 0.861-1.078 [0.81] | 0.904-1.049 [0.83] | 0.881-0.958 [0.89] | 0.904-0.950 [1.09] | 0 | 3/0/2.8e-04 |
| **4 + top 1** | 1.0000 / 1.0000 / 1.0000 / 1.0000 / 1.0000 | 1.0295 (b6 1.031; b10 1.026) | 1.0261 (b7 1.028; b10 1.025) | 1.0425 (b4 1.044; b10 1.040) | 1.0423 (b6 1.043; b10 1.042) | 1.002-1.033 [0.99] | 1.008-1.047 [1.00] | 1.005-1.025 [1.02] | 1.005-1.050 [1.17] | -1e-12 (1.2e-10) | 3/0/2.7e-04 |
| 4 + top 2 | 1.0000 / 0.9999 / 1.0000 / 1.0000 / 1.0000 | 1.0295 (b6 1.031; b10 1.026) | 1.0261 (b7 1.028; b10 1.025) | 1.0410 (b4 1.043; b10 1.016) | 1.0418 (b6 1.043; b10 1.021) | 1.001-1.033 [0.99] | 1.008-1.046 [1.00] | 0.823-0.996 [0.79] | 0.677-0.898 [0.50] | -1e-12 (1.2e-10) | 3/0/2.6e-04 |
| 5 VEF | 0.9999 / 0.9995 / 1.0007 / 1.0692 / 1.0688 | 1.0036 (b0 1.004; b10 0.884) | 1.0040 (b1 1.005; b10 1.004) | 0.9965 (b9 0.945; b10 -57.4) | 1.0000 (b2 0.991; b10 0.888) | 0.976-1.003 [0.98] | 0.986-1.023 [0.98] | -0.069-0.977 [-1.86] | 0.679-0.835 [-1.18] | +7e-09 (1.1e-08) | 3/0/2.4e-04 |
| **5 VEF + top 1** | 0.9999 / 1.0001 / 1.0007 / 1.0692 / 1.0688 | 1.0036 (b0 1.004; b10 1.004) | 1.0040 (b1 1.005; b10 1.004) | 1.0008 (b0 1.001; b10 1.001) | 1.0010 (b1 1.001; b10 1.001) | 0.993-1.003 [0.99] | 0.992-1.024 [0.99] | 1.000-1.000 [1.00] | 1.000-1.023 [1.14] | +7e-09 (1.1e-08) | 3/0/2.4e-04 |

**Eddington face vs diffusion.**  In `grad` the Eddington face itself matches the non-grey
discrete diffusion flux to 9.5e-5 at the deep faces.

**Bitwise checks.**
- Mode 0 is bitwise the Eddington face.  `cmp` of the dumps against the prototype-era `_eo`
  runs (binary e7c50db5, key off): identical for all 11 cases and nw10 rt.txt with
  `_7e70fe86`; for w_T2000, o_T2000, tr_x1.7, grad and nw10 with `_4882f92d`; and for w_T2000,
  o_T3000, grad and nw10 (dump and rt.txt) with the final `_4b1da214`.
- Mode 1 is bitwise the prototype ck_sph_dilute (all 11 cases + rt.txt).
- 4 + top 1 is bitwise identical between `_4882f92d` and `_4b1da214`.
- Every mode changes the result where it should.  Examples: st 0.994 -> 1.000 (blend), the
  deposit 0.92 -> 1.00 (cone), band 10 on the old grid -41.9 -> 1.04 (top 1), L 0.97 -> 1.004
  (VEF).

## Findings per variant

- **V2 blend (mode 2).**
  - It recovers the Eddington deep flux exactly (1e-10; dilute: -0.22 % mean, -0.34 % max).
  - It removes the half-cell error of the sharp st cases (0.994 -> 1.000).
  - It keeps the dilute emergent L (0.968 / 0.983; old grid 0.90 / 0.94).
  - It does nothing for the deposit (still 8-15 % under-cooled).
  - The dilute deep-flux error is a local truncation error of order (A_a - A_b)/A per face
    (0.3 %), not a thin-region effect.  It is present even where cells have tau 1e4.
- **V3a cone path (modes 3, 4).**
  - It fixes the W121 deposit (0.85-0.92 -> 0.99-1.05).
  - It moves L from -3 % to +3 % (old grid -10/-6 % -> +4 %).
  - It is exact on the sharp cases.
  - It is cheap: one extra exp per (cell, chain) in each pass, no new storage.
  - What remains is geometric: the up field of an extended isothermal atmosphere is broader than
    the cone of a sharp photosphere, which a single R-based cone cannot represent.
- **V4 top datum.**
  - `top 1` removes the hot-ghost corruption.  On the old grid, band 10 goes from -41.9 to 1.04
    and the deposit of the top 5 cells from 0.05-0.9 to 1.00-1.05 (top cell 1.02 / 1.17).
  - On W121, where the ghost is not hot, it changes only band 10 at 2000 K (0.97 -> 1.03).
  - `top 2` (the top active cell's kappa and B) radiates back down at the full top T through
    the ghost-pressure column.  That is a model of the unresolved column above, not vacuum, so
    it cannot be judged against this vacuum reference.  On the old grid it moves the top-cell
    deposit to 0.79 / 0.50.
- **V3b VEF (mode 5).**
  - With the formal-solution f it reproduces the formal solution's own L: +0.4 % (W121), +0.1 %
    (old grid).
  - It recovers most of the "grazing-chord" deficit that RESULTS.md called unrepresentable by a
    radial column.  The formal solution contains the chords, and f/q pass them to the column.
  - What made it right (all needed; each was measured wrong without it):
    1. Both chains at mu_eff = sqrt f.  With the Gauss pair scaled by sqrt(3f), the pair's K/J
       is not f in semi-thick cells: L was 0.940, deposit ~0.98-1.02.
    2. The top reflection zeta from h.
    3. The Marshak-type wall.  With u = B at the wall, the tr cases over-emitted 7 % (2 pi mu_eff
       B = 1.155 pi B).
  - Remaining weakness: an unresolved photosphere inside a single ultra-thick cell (st: +6.9 %).
    The cell's centre f is 1/3 and a single mu_eff = 1/sqrt 3 emits 1.155 pi B.  A possible
    remedy, not done: fall back to the Gauss pair in optically thick cells.

## VEF lag and cost

- **Lag.**
  - nw10 run (T falls by up to 297 K in 10 cycles at the top): refreshing f every 5 / 10 cycles
    instead of every cycle changes the final T by at most 0.37 / 0.83 K.
  - The Newton pass counts are identical (3/0), and so are the budget residuals (2.2e-4 max).
  - In production T changes far less per cycle, so a refresh every O(10-100) cycles looks
    adequate.  This is not tested on a relaxed production state.
- **Cost (operation count, not a timing).**  GPU timing was not allowed in this task.
  - Per column and refresh, the formal solution does ~ 88 (b,g) x [(ncore + 1 + n + 2) rays x
    2 (n+1) segments x <= 2 exp].  For W121 (n = 76) that is ~1.3e6 exp + 88 x 76 k-table
    lookups.
  - One tm pass is ~176 x 76 = 1.3e4 half-layer expm1 + 88 x 76 lookups.
  - So a refresh costs about 100 tm passes of exponentials and one pass of table look-ups.
  - With ~3 passes per cycle, an overhead of <= 20 % needs a refresh every ~150 cycles at the
    exp-bound estimate.  Two cheap cuts are possible, neither implemented:
    - Restrict the formal solve to nodes above the level where every (b,g) has tau-from-top > 50
      (f = 1/3 below).  On W121 that is ~22 of 78 nodes, ~10x cheaper.
    - Use fewer core rays.
  - Must be timed on apudev before any production use.  The formal kernel is one thread per
    (column, b, g) with O(n^2) serial work.
  - Memory: 2 x (nmb, 88, n1+3, n3, n2) Reals for s / gamma, plus a 14x work array of the same
    shape (prototype layout; the work array could be thread-private).
- **Error of the cheap geometric flavour vs the formal-solution version** (mode 4 + top 1 minus
  mode 5 + top 1):
  - L: +2.6 / +2.2 points on W121, +4.2 / +4.1 points on the old grid.
  - Deposit: up to +4.7 % vs up to +2.4 % in the thin shell.
  - Sharp cases: 4 is exact, 5 is +6.9 % on the unresolved slab.
  - Deep flux: both equal the Eddington face.

## Newton at larger dt

- At `rt_test_dt = 30` s, modes 0, 4 + top 1 and 5 + top 1 all converge: 6-8 passes, 0 nonconv,
  per-cycle counts identical to within one pass.
- At 100 s every arm hits maxit 16 in most cycles, including mode 0.  Out of 10 cycles, the
  nonconv counts are:
  - Modes 0 / 1 / 2 / 3 / 4: 3 / 9 / 4 / 7 / 7.
  - With top 1: mode 0 9, mode 4 10, mode 5 8 (5 alone: 8).
- Mode 0 itself goes from 3 to 9 with top 1, so the extra failures come from that datum.  With
  vacuum above, the top cells cool fast: T drops by up to 297 K over the 100 s of nw10 with
  5 + top 1.  The failures are not from the face.
- Pre-existing limit: the prototype already failed identically at 300 s.

## Porting the winner to ck_impl_lin / ck_impl_jac_lin (not done)

The stored-factorisation path is `rt_chain_ck_lin` / `_lin1` and `ck_lin_build` in
two_stream_rt.hpp (~7000-7600 now, the old ~6901-7229 / ~7842 / ~8058).
- **Mode 4.**
  - `ck_lin_build` would form R per face with the blend map; e needs dtau_f from `ckkro_g` and
    dz, which are already stored.
  - It would store per (face, chain) the up-crossing transmission tu, or R_f and e.  That is
    +1-2 Views of the `lP_g` size (today 5 per chain).
  - The residual kernel's pass 1 needs `rn/sv` replaced by (r1, c1) and tu, and pass 2 needs
    `ub/db` from (al, be) and the up ray with tu.
  - The jac_lin window recursion (`rt_chain_ck_jlin`) needs the same al/be/c1 and wr = tu^2
    tr^2 c1 Q, as in the tm JAC block.
- **Mode 5.**
  - Everything for mode 4, plus the per-(cell, b, g) mu_eff and per-face gamma from the formal
    kernel, which the lin kernels only read.
  - The stored triple (`ckc0/ci/co`) must be formed with dtau sqrt(3) mu/s.  It already is when
    the storing pass runs the tm kernel.
  - The top zeta and the wall beta go into the lin kernels' first and last rows.
- **ck_sph_top.**  The lin kernels read `cktpf_g * Bb(ie+1)`.  Mode 1 is a zero there; mode 2
  needs B(ie) and the extra diagonal term.
- In all cases the lin path must be re-gated bitwise at key 0 and re-checked for the Jacobian
  sign drop (negative off-diagonals are dropped).

## Open questions for the user

1. Is the formal-solution cost acceptable, at a refresh every ~50-150 cycles, after a GPU
   timing?  If not, mode 4 + top 1 is the fallback, at the +3-4 % L cost.
2. The top datum.  Should the thermal down-stream at the top be vacuum (top 1, which matches
   the exact reference and removes the hot-ghost artefact) or the column-above model (top 0/2)?
   This is a modelling choice for WASP-121b, where the domain top is at ~1e-8 bar.
3. Should mode 5 fall back to the Gauss pair in optically thick cells, for unresolved
   photospheres?  The production grids resolve them.
4. ck_dif_dtau handover columns are refused for modes 3-5 (R = 1 datum).  Is that fine, given
   that the production keys have it off (default 0)?
