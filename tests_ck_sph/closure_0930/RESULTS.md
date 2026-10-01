# ck_spherical thermal closure vs an EXACT reference (cksph_test_0930)

**Verdict: CONFIRMED.**  The spherical face closure of the correlated-k thermal two-stream
(`problem/ck_spherical = true`, the tm form `ck_sweep_form = 1` that production runs) over-emits.
Below a transparent shell of area ratio x = A_top/A_ph, a sharp photosphere emits
**2x/(1+x) sigma T^4**, reproduced to 1e-4 by the real kernel.  On static isothermal atmospheres
it over-emits by **22-25 % on the WASP-121b 1x grid** and by **24-34 % on the old dhj grid**.
The plane-parallel form (`ck_spherical = false`) emits A_top F = x sigma T^4 A_ph, so it over-emits
by 67-72 % (W121 grid) and 104-136 % (old grid).

**Mechanism.**  The face conditions "S = (u+d)/2 continuous, A (u-d) continuous" are the spherical
moment equations with the Eddington closure K = J/3, whose geometric term (3K - J)/r is then zero.
J therefore does not dilute through a transparent shell.  With d = 0 at the top, S = const and
A u = const give d_ph = B (1-x)/(1+x) < 0 at the photosphere, and the photosphere loses
u - d = 2x/(1+x) B.  Every sweep form (4-pass probe, tm, sd) and the Jacobian carry the same face.

## Setup (the AthenaK kernel itself, not a re-implementation)

- Binaries: `/viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_7314aad3` (test hooks only) and
  `athena_dhj_cpu_e7c50db5` (+ prototype; key off it is BITWISE the first binary: `cmp` of the
  dumps for w_T2000, w_T3000 and tr_x1.7 are identical).  Branch `ck-sph-closure-fix`, off
  rt-integration fef0ea83.  Both are serial CPU runs on the login node.
- Input: `base_w121.athinput` is the WASP-121b 1x production input (`w121prod_1x.athinput`: ck
  hiT2 11-band table, ck_nquad 2, ck_implicit, tm, ck_beam_sph, pcut 1e8) moved to a
  spherical-polar 76 x 4 x 8 mesh.  It uses the same stretched radial grid, x1 1.126575e10 ..
  1.631037e10.  Test keys: `rt_test_freeze`, `rt_test_dt = 1e-3` (static), `rt_test_T` (uniform T
  on the IC density), `ck_int_at_cut = false` (no T_int), and `ck_impl_lin = ck_impl_jac_lin =
  ck_beam_par = false`.  The lin keys make no difference to the answer: with the production lin
  keys on, w_T2000/w_T3000 give the same ratios, 1.2243 / 1.2481.  No restart was used: this is
  a static RT test on the pgen IC density (ic_w121_1x.txt), so no rotation applies.
- Old dhj grid: x1 9.44e9 .. 2.0556e10, nx1 128, the 4-coefficient stretch, and the
  picket-fence IC (ap 9.44e9, grav 942).
- New test hooks (default off): `problem/ck_dump_kap` appends kappa(b,g) of every cell and the
  g weights to the `ck_dump_file` column dump.  `problem/rt_test_rstep` / `rt_test_rho_top` set
  a density step.
- Exact reference (`ana/exact.py`): L = 4 pi^2 sum_b B_b sum_g w_g int (1 - e^-tau_bg(p)) 2p dp.
  tau is taken along straight chords through the dumped kappa_bg * rho of the same cells.  Rays
  hitting the inner wall get I = B (isothermal), which is the code's bottom datum.  L_code is
  4 pi r_top^2 F_b(top face).
- Check: sum pi B_b / sigma T^4 = 1.000005.

## 1. Sharp photosphere (single T = 2000 K, 11 bands)

The shell is rho = 1e-14 (tau_v <= 0.07 in every b, g).  Below it the photosphere is either the
inner wall itself (tr) or a rho = 0.1 slab (st, tau per cell >= 5.6e6).  Exact L/(A_ph sigma T^4)
= 1.0000 in every case.

| case | x | Eddington face (sph=true) | 2x/(1+x) | plane-parallel (sph=false) | ck_sph_dilute |
|---|---|---|---|---|---|
| tr (wall), 1.2 | 1.2000 | 1.09091 | 1.09091 | 1.20001 | **1.00000** |
| tr (wall), 1.7 | 1.7000 | 1.25905 | 1.25925 | 1.69970 | **0.99984** |
| tr (wall), 3   | 3.0000 | 1.50002 | 1.50000 | 3.00001 | **1.00002** |
| st (slab), 1.2 | 1.2049 | 1.08927 | 1.09292 | 1.20488 | **0.99391** |
| st (slab), 1.7 | 1.7024 | 1.25415 | 1.25991 | 1.70238 | **0.99277** |

In the st cases the remaining -0.6/-0.7 % is the half-cell at the photosphere face.

## 2. Static isothermal atmospheres: L_code / L_exact per band and total

| run | b0 | b1 | b2 | b3 | b4 | b5 | b6 | b7 | b8 | b9 | b10 | TOTAL | L_exact/(A_in pi B) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| W121 1x grid, 2000 K, sph=true (Eddington face) | 1.2188 | 1.2093 | 1.2083 | 1.2324 | 1.2079 | 1.2321 | 1.2311 | 1.2387 | 1.1849 | 1.1880 | 1.1611 | **1.2243** | 1.2525 |
| W121 1x grid, 2000 K, sph=false (plane-par.) | 1.6462 | 1.6290 | 1.6282 | 1.6940 | 1.6296 | 1.6940 | 1.6918 | 1.7139 | 1.5729 | 1.5824 | 1.5885 | **1.6736** | 1.2525 |
| W121 1x grid, 2000 K, sph=true + ck_sph_dilute | 0.9687 | 0.9626 | 0.9618 | 0.9689 | 0.9605 | 0.9686 | 0.9680 | 0.9703 | 0.9520 | 0.9524 | 0.9142 | **0.9660** | 1.2525 |
| W121 1x grid, 3000 K, sph=true (Eddington face) | 1.2377 | 1.2451 | 1.2413 | 1.2499 | 1.2490 | 1.2501 | 1.2531 | 1.2476 | 1.2424 | 1.2406 | 1.2138 | **1.2481** | 1.2216 |
| W121 1x grid, 3000 K, sph=false (plane-par.) | 1.6812 | 1.6976 | 1.6936 | 1.7204 | 1.7125 | 1.7226 | 1.7307 | 1.7178 | 1.6946 | 1.6918 | 1.6327 | **1.7159** | 1.2216 |
| W121 1x grid, 3000 K, sph=true + ck_sph_dilute | 0.9794 | 0.9831 | 0.9803 | 0.9816 | 0.9831 | 0.9811 | 0.9821 | 0.9795 | 0.9807 | 0.9795 | 0.9674 | **0.9809** | 1.2216 |
| old dhj grid, 2000 K, sph=true (Eddington face) | 1.2247 | 1.1957 | 1.1857 | 1.2629 | 1.1934 | 1.2626 | 1.2594 | 1.2867 | 1.1328 | 1.1135 | -41.9462 | **1.2360** | 2.3158 |
| old dhj grid, 2000 K, sph=false (plane-par.) | 1.9412 | 1.8438 | 1.8368 | 2.1510 | 1.8381 | 2.1542 | 2.1372 | 2.2772 | 1.6218 | 1.6251 | -41.9398 | **2.0448** | 2.3158 |
| old dhj grid, 2000 K, sph=true + ck_sph_dilute | 0.9100 | 0.8989 | 0.8902 | 0.9022 | 0.8956 | 0.9023 | 0.9014 | 0.9054 | 0.8875 | 0.8626 | -41.6492 | **0.8978** | 2.3158 |
| old dhj grid, 3000 K, sph=true (Eddington face) | 1.2855 | 1.3212 | 1.3058 | 1.3489 | 1.3435 | 1.3522 | 1.3658 | 1.3339 | 1.3088 | 1.3010 | 1.1351 | **1.3381** | 2.0079 |
| old dhj grid, 3000 K, sph=false (plane-par.) | 2.1026 | 2.2356 | 2.2278 | 2.4066 | 2.3498 | 2.4263 | 2.4844 | 2.3624 | 2.1998 | 2.1818 | 1.8138 | **2.3609** | 2.0079 |
| old dhj grid, 3000 K, sph=true + ck_sph_dilute | 0.9280 | 0.9397 | 0.9311 | 0.9400 | 0.9424 | 0.9406 | 0.9430 | 0.9311 | 0.9337 | 0.9305 | 0.8304 | **0.9363** | 2.0079 |

In the last column the photosphere sits above x1min, so L_exact/(A_in pi B) > 1.

- Band 10 on the old grid at 2000 K is not a closure effect.  B_10(2000 K) is tiny, and the
  top ghost datum dominates it: I_down(top) = (1 - e^-dtau) B(T_ghost), with T_ghost = 3415 K and
  p_ghost = 7.8e-7 bar left over from the IC ghost cells.  It carries a negligible share of L.
- The same ghost datum explains the top 2-3 cells of the old grid, where the code heats and
  the exact solution cools (see section 3).

**What remains with the fix (-2 to -10 %) is not the face closure.**  A radial column cannot
represent the long grazing chords through an extended semi-transparent atmosphere; these
chords raise the exact limb photosphere.  The deficit grows with H/r: it is larger on the old
grid and at 2000 K.  The sharp tests (section 1) have no such chords, and there the fix is exact.

## 3. Deposit (net radiative heating) in the thin shell vs exact 4 pi kappa rho (J - B)

The exact J comes from 48-point Gauss rays per (b, g) at the cell centres.  The ratio is
Q_code / Q_exact; both are cooling (Q < 0) in an isothermal atmosphere.  The cells are those
with tau_v(max g) < 40.

| case | Eddington face | plane-parallel | ck_sph_dilute |
|---|---|---|---|
| W121, 2000 K (r/r_in 1.18 -> 1.42) | 1.08 -> 0.80 | 0.87 -> 0.69 | 0.92 -> 0.85 |
| W121, 3000 K (1.22 -> 1.42) | 1.05 -> 0.83 | 0.87 -> 0.73 | 0.92 -> 0.87 |
| old grid, 2000 K (1.86 -> 2.17) | 0.95 -> **-1.27** | 0.76 -> -1.33 | 0.88 -> -1.22 |
| old grid, 3000 K (2.00 -> 2.17) | 0.85 -> -0.75 | 0.73 -> -0.78 | 0.83 -> -0.72 |

The negative values in the top cells are the hot-ghost top datum, identical in all three
closures.  Under the Eddington face the deposit swaps sign across the shell: it over-cools the
lower thin shell by up to 8 % and under-cools the top by 20 %.  Under the dilute face it is a
uniform 8-15 % under-cooling.  The cause is that the diluted up stream is absorbed as if it were
hemispherically isotropic (path 1/mu), whereas the real field is forward-peaked.  That would
need a variable Eddington (path) factor, which is not attempted here.

## 4. Prototype fix: `problem/ck_sph_dilute` (branch ck-sph-closure-fix, commit e7c50db5)

The face condition for hemispherically isotropic streams: the down intensity is invariant
(d_below = d_above), and the down rays that miss the sphere below turn round into the up ray.
The up ray satisfies A_a u_a = A_b u_b + (A_a - A_b) d_a.
- A (u - d) is still single-valued at every face, so the deposit telescopes as before.
- An isotropic field (u = d) passes unchanged.
- A u is conserved through a transparent shell.
- tm maps: R <- 1 - rho (1 - R), Sc <- rho Sc, with rho = A_b/A_a.  Solve: u_b = R d_a + Sc,
  d_b = d_a.  There is no divide at the face.
- The JAC tridiagonal uses the same maps (al = 1, be = 0, c = rho).

Status:
- Default false, and bitwise off.
- It is implemented in the tm (FRM = 1) chain kernel only, explicit and ck_implicit, and is
  REFUSED at startup with ck_impl_lin / ck_impl_jac_lin.  The stored-factorisation kernels still
  carry the old face: two_stream_column_ck.hpp ~6901-7229 and rt_chain_ck ~7842 / 8058.
- Newton behaviour (w121 grid, 2000 K, rt_test_dt = 10 s, 10 cycles): 3 passes and 0 nonconv
  per cycle, both keys off and on.  Budget residual 5e-5 in both.  At dt = 300 s both fail
  identically (maxit 16), which is pre-existing.
- Thick-limit cost: in a column with a T gradient (runs/grad_*), the deep diffusive flux is
  0.2-0.3 % below the Eddington face.  The Eddington face matches the non-grey diffusion flux to
  1e-4 there.  The emergent flux is 0.788x, which is (1+x)/2x.

## Files

- Runs: `/viper/ptmp2/jinma/cksph_test_0930/runs/<case>/`, with dump.txt (column + kappa) and
  exact.txt.  Suffixes: `_eo`/`_sph`/`_sphS` = Eddington face, `_pp` = plane-parallel,
  `_dil` = prototype.
- Scripts: `run.sh`, `runall.sh`, `keys.sh`, `base_w121.athinput`, `ana/exact.py`.  They are
  also committed on the branch under `tests_ck_sph/closure_0930/`.
