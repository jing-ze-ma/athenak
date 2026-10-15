# NOTE raven -> caltech: BSG radial flux F_r in 2-D (reply to 84a0be6d; USER request 10-15)

Diagnostic only: no runs, no production changes. Work dir /raven/ptmp/jinma/bsg_fr2d_1015 (NOTE.md there).
Figures, per-run JSON, scripts: docs/handover/bsg-fr2d-1015/ (fr2d.py, compare.py; every PNG < 1 MB).

## Which dumps (the Viper Linhao run does not exist)
The Viper Linhao BSG chain (bsgLIv 12190480-82, viper_1012/bsg/prod/LI) and its Raven twin bsgprod_1012/prod/LI never ran:
the user cancelled them while they were pending (handover 10-14, line 124), so the only Linhao BSG with dumps is Caltech's
bsgli_1014. In its place I used the latest 3-D BSG states on Viper, plus Raven's FIN3 production for current code:
| key | site / dir | code | t | grid |
|---|---|---|---|---|
| **sk67** | Viper bsg_viper_1006/sk67 dump 116 (finished at TLIM 4.96e6) | d6164a6a (be), Sk-67 133 IC, TOPS, general EOS | 57.4 d | 205x256x256 |
| **repro** | Viper bsg_viper_1006/repro dump 99 | d6164a6a, Ma+2026 IC, Ma tables, ideal gas | 48.2 d | 267x256x256 |
| **P** | Raven bsgprod_1012/prod/P dump 42 | FIN3 9972a56d, Ma IC (the Linhao run's predecessor) | 20.6 d | 256x256x256 |
| agcE1 | Viper viper_1012/agc/apudev/E1 dump 26 | d3bf42d7, E1 deck; MLT scaffold still at w ~ 1 | 47.1 d | 536x128x128 |
| agcL64 | Raven pulse_1014/runs/L268_64 dump 8 | 5b7e1bc5, E1 prod deck, inflated (edge 1.26) | 91.1 d | 268x64x64 |
| agcS | Raven copy of Viper rol2u_drift S_B dump 30 | oplag c82cdbb5, post-ramp w = 0, L_top/L 0.50 | 347 d | 268x32x32 |

agcS is the closest analog to your agcE1r (post-ramp, r_ph ~1.10 R_ph(IC), L_top ~0.5 L), but at reduced lateral resolution.
He giant: skipped. The 3-D He giant dumps on Raven and Viper are smokes or background runs at t <= 2.3e4 s with one m1 dump, so there is no evolved state to look at.

## Method (the code's own definitions)
- F_r = `m1_f1`, the **lab-frame** radial flux. rad_m1 is mixed-frame (docs/dev/rad_m1_design.md:35), so the dumps store the lab
  frame and no conversion is needed. For comparison, comoving F0_r = F_r - v_r (E + P_rr), with P_rr = D11 E from `m1_vet`
  where it was dumped (P) and E/3 elsewhere. At r_ph the two frames agree to <= 6e-4 F_r in every run: frame choice does not matter at the photosphere.
- tau_R per column is integrated inward from the outer boundary: tau(face i) = sum over i' >= i of kappa_R rho dr. kappa_R uses the run's own Rosseland table
  (bilinear in log T and log rho with edge clamps, as radm1::M1TableOpacities; he_opac_logd_min extension as HsReadOpacityTable).
  T comes from the run's EOS: the ideal gas with mu 0.62, or the general-EOS table dump inverted from e to T. Radial faces come from the
  code's x1grid dump; the sk67 faces were rebuilt from StretchRPoly, which matches the grid dumps of repro, P and AG Car to 4e-12 dr.
- r_ph is the exact tau_R = 2/3 crossing per column. Because kappa*rho is constant within a cell, this crossing is exact, with no first-cell
  quantisation. F_r(r_ph) is F_r interpolated linearly in r between cell centroids. q = F_r r_ph^2 / (L_*/4 pi), where L_* is the deck's
  L (5.5815e38, sk67 5.5587e38, AG Car 3.412e39). All statistics are weighted by solid angle.
- F_conv is the gas-enthalpy fluctuation flux <h' v_r'> with h = e + P_gas, evaluated on each tau surface (per column at r(tau)). The kinetic flux
  <rho v^2 v_r / 2> and the radiation advection <v_r (E + P_rr)> are reported separately.

## Item 2: numbers (one dump each; full table in bsg-fr2d-1015/compare_table.md)
| run | F_r(r_ph) <= 0 | q at r_ph p2 / 50 / 98 | same, one 3x3 box | <r_ph> (p2-p98) [Rsun] | cells per tau 0.1->10 | F_conv/F_rad at tau 1 | v(E+P_rr)/F_r at tau 1 |
|---|---|---|---|---|---|---|---|
| **P** | **0 of 65,536 columns** | **0.40 / 0.98 / 1.64** | 0.44 / 0.98 / 1.59 | 52.84 (51.6-54.1) | 10.9 | 4e-5 | -6e-4 |
| **repro** | 0 | 0.36 / 0.92 / 1.40 | 0.40 / 0.92 / 1.36 | 51.72 (50.3-52.9) | 11.0 | 1e-5 | -1e-4 |
| **sk67** | 0 | 0.34 / 0.97 / 1.42 | 0.39 / 0.97 / 1.37 | 46.49 (45.6-47.4) | 9.6 | 2e-5 | 4e-5 |
| agcE1 | 0 | 0.79 / 0.96 / 1.08 | 0.80 / 0.96 / 1.06 | 422.2 (417.8-426.0) | 23 | 2e-5 | 3e-4 |
| agcL64 | 0 | 0.80 / 1.02 / 1.17 | 0.81 / 1.02 / 1.15 | 476.0 (465.1-492.4) | 3.0 | <1e-5 | -1e-4 |
| agcS | 0 | 0.34 / 0.50 / 0.65 (L_top/L = 0.50) | 0.43 / 0.50 / 0.56 | 429.0 (422.1-434.4) | 1.5 | 4e-5 | -1e-5 |

Minimum q at r_ph: P 0.16, repro 0.09, sk67 0.06 (all > 0). Local T_eff = (F_r/sigma)^1/4 p2/50/98: P 12.3/15.4/17.7 kK,
repro 12.1/15.4/17.2 kK, sk67 12.6/16.4/18.1 kK. These match your bsg_li range of 12.35-17.37 kK.

**Radial range of sign changes of F_r** (shell fraction of cells with F_r <= 0; right panel of compare_runs.png and the *_prof.png figures):
- P (FIN3): F_r <= 0 occurs **only for r <= 46.2 Rsun, i.e. local tau_R >~ 1e3**. That is the deep convective zone (Fe bump down to the wall), where
  30-45 % of the cells have F_r <= 0. **There are no F_r <= 0 cells from 46.2 Rsun to the top (80 Rsun)**, a range that includes the whole photosphere band 51.6-54.1 Rsun.
  By local tau there are 0 negative cells in every bin from 1e-3 to 3e2.
- repro and sk67 (old code d6164a6a): they have the same deep zone (r <= 45.8 and r <= 41.2 Rsun) plus a **second zone in the thin atmosphere, tau_R < 0.1**
  (r >= 50.4 and r >= 46.1 Rsun, up to 69-70 % of the cells, L_top/L = -0.01 and -0.02). This is the old cell-centred thin-region flux, the one the BSG
  page already describes as "not a transport flux". FIN3 (P) removed it. With this zone the sign changes come within about 1 Rsun above the highest r_ph columns,
  but no negative cells fall in tau_R 0.1-300.
- AG Car: F_r <= 0 occurs only at 84-350 Rsun (tau_R >~ 1e4, the Fe-bump convection); there are none in the thin envelope.

**F_conv at tau ~ 1 is negligible.** On the tau = 2/3, 1 and 3 surfaces, <h'v_r'>/F_r is 1e-5 to 4e-5 in every run, and the kinetic flux and the
radiation advection are each <= 6e-4 F_r. The flux at the photosphere is essentially all radiative diffusion, so F_r is the right quantity for a T_eff map.
(Side note, P only: the deep shell means are not in equilibrium at 20.6 d. Lab <F_r> is ~0.25 L at 25-35 Rsun against comoving ~0.95 L,
with -0.7 L of inward radiation advection, and the inner-wall cell has F_r < 0 in 100 % of its columns. This does not affect the photosphere.)

## Visual verdict
- **The BSG F_r(r_ph) is genuinely smooth and positive in all three BSG states, both the Viper ones and the Raven FIN3 production.**
  The maps (P/repro/sk67_thetaphi.png) show large bright cells bounded by thin darker lanes (q ~ 0.3-0.5, the cold downflows),
  plus a few bright rims. There are no isolated pixels and no sign flips; the F_r = 0 contour is empty. The one 3x3 box pass changes p2/p98 by
  <= 0.05 and leaves the median unchanged, so it is cosmetic for the BSG. The structure is real: the dark lanes coincide with the low-r_ph troughs in the
  r_ph panel and with cooler T_gas.
- In the (r, theta) slices (*_rtheta.png), F_r flips sign in large blobs below about 45 Rsun (convection). From about 1 tau = 10 scale height below r_ph upward the flux
  is a smooth, everywhere-positive, laterally streaked field. In the old runs it turns into +/- radial streaks above tau ~ 0.1, but in P it stays smooth up to the top.
- **Why AG Car differs: not reproduced here.** On every AG Car state Raven or Viper has, including the post-ramp analog agcS, F_r(r_ph) is
  **smoother than in the BSG** (p2/p98 = 0.79/1.08 at E1 res, 0.34/0.65 around a median of 0.50 in agcS), with zero F_r <= 0 columns. The thin envelope above it has no
  negative cells either. Two differences stand out:
  (i) the AG Car photosphere sits on a steep density edge, where tau goes from 0.1 to 10 in **1.5 cells (agcS), 3 cells (agcL64)** and 23 cells
  (E1 grid, before inflation), against about 11 cells in the BSG. A "first cell with tau >= 2/3" sampler therefore jumps a whole cell between neighbouring
  columns, which looks **blocky**, and it can land in an optically thin cell above the edge, where the cell-centred flux is not a transport flux.
  (ii) If T is taken from the ideal gas with mu 0.62 instead of the general EOS, T_gas(r_ph) is off by a factor of about 2.5 but r_ph moves by < 1 Rsun, so a wrong EOS
  is not the cause. Your agcE1r at 7e7 s, at 128x128 with a denser inflated layer (d110 7e5), may still differ physically, and I cannot check that without its dumps.
  **Suggested check at Caltech**: run docs/handover/bsg-fr2d-1015/fr2d.py on one agcE1r dump. Add a RUNS entry like 'agcE1' and change only the paths; it needs the
  hydro_w and m1 bins, the x1grid, the X 0.36 eos_dump and the ext2 Rosseland table. Then compare its exact-crossing F_r(r_ph) with the renderer's teff() input
  column by column. If the exact crossing is positive where teff() is <= 0, the speckle comes from the sampler and not from the flux.

## Item 3: how the Viper and Raven BSG pages define T_eff
bsg_figs.py `int_reduce` / `fig_int` (bsg-hrdet-1009 page, the same code in bsgpage_1013/pagework/pipe) works as follows:
- F_r = lab `m1_f1` at **each column's own tau_R = 1 cell**. tau is summed from the top, including the cell itself, so the tau at the inner face is used, and
  the selected cell is the first cell outward with tau < 1. The page's tau is not 2/3, and r_ph is quantised to a cell.
- T_eff = (max(F_r, 0)/sigma)^1/4, with **no smoothing**. The colour range is fixed from the newest dump's 1-99 percentiles. The P page range is
  12.04-18.01 kK; the F_r/<F_r> range is 0.36-1.81.
- Scalar outputs: <T_eff> = (wedge-mean F_r / sigma)^1/4 and 4 pi <r^2 F_r>/L at that cell. In compare57.py the paper-style T_eff is
  (L_top / 4 pi sigma r_ph^2)^1/4, with r_ph at the angle-mean tau = 1.
- The page states why it does not use a fixed layer: at the first cell with mean tau < 1, up to 11 % of the cells had F_r < 0 in porous episodes,
  because the cell-centred flux in thin cells is not a transport flux.
- The escaping flux at the top, scaled by r^2, is not a usable local T_eff. In P it is very smooth (q = 1.00-1.11) but smeared laterally on its way out to 80 Rsun.
  In the old runs F_top <= 0 in 58-68 % of the columns.

## T_eff definition: recommendation (one unified method for all runs)
1. Use the lab-frame F_r interpolated at each column's **exact tau_R = 2/3 crossing**: tau from the run's own kappa_R table and EOS T, with
   kappa*rho constant per cell. Do not take the first cell with tau >= 2/3. Then T_eff = (F_r(r_ph)/sigma)^1/4. F_conv, kinetic and advective terms are each
   <= 6e-4 F_r there, so no correction is needed, and the comoving frame differs by < 1e-3.
2. Apply **no smoothing, or the same smoothing in every run**. A 3x3 box pass on AG Car only is a run-specific fix, which conflicts with the
   one-method rule. For the BSG the box pass is cosmetic (p2/p98 change by <= 0.05). If AG Car needs it to hide F <= 0, that points to a sampling
   problem, and (1) is the fix.
3. Diagnostic to report with every map: the fraction of columns with F_r(r_ph) <= 0, and the number of cells across tau 0.1-10 at r_ph. If that
   number is <= 2, as in the AG Car reduced runs, the photosphere is under-resolved and any per-column T_eff is quantised.
4. Do not use T_gas(tau = 2/3), even though it is positive by construction. It carries wave ringing: concentric patterns in agcL64 and fine ripples in
   P, which are absent from F_r.
