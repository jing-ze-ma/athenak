---
name: exofms-cross-validation
description: "2026-08-21: Exo-FMS_column_ck built and run on identical columns. C++ kernel longwave agrees to 1.7% at TOA. THE COMPARISON FOUND A REAL BUG: stellar heating was losing 26% of the insolation (fixed, 7dbfa1fd) and THE GREY PICKET FENCE STILL LOSES ~59% -- production runs are affected."
metadata:
  node_type: memory
  type: project
  modified: 2026-08-21T21:00:00.000Z
---

The physics rung of the validation ladder in [[correlated-k-design]]. Everything lives in
`/viper/u2/jinma/ATHENAK/bench/exofms_compare/` -- README.md has the full recipe, the
matched namelist, Exo-FMS's converged output, and `ourrt.py`.

## Build traps

- `git clone` works; the codeload tarball URL 404s.
- gfortran from `gcc/14` plus MKL for LAPACK/BLAS (swap `-llapack -lblas` in `src/Makefile`
  for `-L${MKLROOT}/lib/intel64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm
  -ldl`), and `export LD_LIBRARY_PATH=$MKLROOT/lib/intel64:$LD_LIBRARY_PATH` to run.
- **Build SERIALLY.** `make -j` races on the `.mod` files and dies in IC_mod.
- DISORT is already commented out of the Makefile's object list; leave it that way.

## Setup traps

- **Insolation does NOT come from Tirr.** `Finc = (Rs*Rsun/(sm_ax*au))^2 * swflux`; Tirr
  only sets the initial condition. Leaving Rs and sm_ax at their HD189 defaults while
  changing Tirr gives a silently wrong run.
- **This also resolves the sw_flux normalisation question**: the files are stellar SURFACE
  flux in W/m2. The W121 file sums to 9.788e7 = a 6446 K star. Our own reader sidesteps
  this by using only the shape and renormalising to sigma T_irr^4.
- To match our insolation: Rs = 1.458, sm_ax = 0.02254 -> dilution 0.09052 -> Finc total
  = sigma T_irr^4 = 8.860e6 W/m2.
- **nstep = 10000 is NOT converged** for this case: OLR off from ASR + F_int by 23 %.
  200000 steps gets it to 0.1 %. Always check that balance before using a profile.
- Exo-FMS is MKS, our pgen cgs.

## Result

Our longwave net flux on their converged profile: **OLR agrees to 2.4 %** (5.244e6 vs
5.120e6 W/m2), 2-4 % everywhere the flux is large, and the base net flux is sigma T_int^4
to **0.04 %**. Ratios wander where the net flux has nearly vanished, which is a difference
of large numbers, not a disagreement. Sanity: their OLR gives Teff 3082.5 K against the
analytic (T_int^4 + mu_0 T_irr^4)^(1/4) = 3082.1 K.

Residual is most likely (a) different longwave closure -- theirs `lw_AA_L`, absorption
approximation; ours short characteristics with a diffusivity factor, (b) we lack He- and
H2- free-free that their n_cia = 7 carries, (c) their k-table interpolation is Bezier.

## DIRECT C++ COMPARISON (the caveat below is now closed)

`problem/ck_dump_file` + `ck_dump_m/j/k` dumps one column from the production kernel.
Exo-FMS was patched with `iIC = 7` (reads layer T from `IC_Tl.txt`) and given OUR grid via
`a_sh` = our face pressures in Pa, `b_sh` = zeros, so both codes run the identical column.
Use `adj_scheme = 'none'` (lowercase!), `nstep = 1`, `corr = .False.`.
Meshblock 0 is all NIGHTSIDE in this setup; m = 16, j = 9, k = 9 gives mu0 = +0.539.

**Longwave agrees to 1.7 % at TOA** (7.545e6 vs 7.678e6 W/m2) and ~1.7-2.4 % down to
1e-2 bar. Ratios fall to 0.6-0.75 between 0.08 and 3 bar where the net flux is small.

## THE BUG IT FOUND -- and one that is still live

Stellar heating was deposited as `kappa rho F exp(-tau)` with tau at the cell's LOWER face,
correct only for thin layers. Deposited/absorbed is `u e^-u/(1-e^-u)`, u = dtau/mu: 0.95 at
u = 0.1, 0.58 at u = 1, 0.31 at u = 2.

- Correlated-k was absorbing 3.513e6 W/m2 where (1-A) F* mu0 = 4.778e6: **26 % of the
  insolation gone**. FIXED in 7dbfa1fd by depositing the flux DIFFERENCE across the cell;
  now 4.722e6, 98.8 % of analytic, 99.2 % of Exo-FMS. No extra transcendental -- the upper
  face transmission is carried from the previous cell.
- **The grey picket fence had the same defect: it deposited ~76 % of F* mu0**, losing about
  a quarter of the insolation, near enough independent of the visible-band ratio from 0.1
  to 3 because the ~0.46 scale-heights-per-cell grid rather than the opacity sets
  dtau/mu ~ 1. **FIXED in b4e0953c at the user's instruction. THIS CHANGED THE ANSWER OF
  EVERY PRODUCTION RUN.** RT went 124.8 -> 132.9 ms/100cyc (+6.5 %). Split and monolithic
  still agree to 1.3e-6.
  - **`bench/polar_ab/fixed200.a` IS NO LONGER A VALID REFERENCE.** Use
    `bench/polar_ab/fixed200_swfix.a`. Every A/B script in that directory still points at
    the old one.
  - Every deep-hot-Jupiter run before b4e0953c absorbed about three quarters of the stellar
    flux it should have.

**A number I got wrong first time:** I initially reported the grey scheme as depositing
41 % / losing 59 %. That was a labelling error in my own analysis script -- the 0.41 was in
units of F*, not F* mu0, and with mu0 = 0.539 the exact answer IS 0.539. The real figure is
76 % deposited, 24 % lost, which is consistent with the 26 % measured for correlated-k on
the same grid, as it should be.

## The caveat that mattered, now closed

`ourrt.py` transcribes our scheme; it is NOT the C++ kernel. The chain is: C++ opacity and
continuum match this Python digit for digit on device, the C++ sweeps pass the isothermal
and transparent-slab limits exactly, and the transcribed scheme matches Exo-FMS to 2.4 %.
**A direct C++-vs-Exo-FMS comparison still needs a column-dump diagnostic in the pgen** --
that is the remaining gap, and it is a small one.
