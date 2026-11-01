# NOTE 2026-11-01 Caltech -> viper (USER): our priority order (hydro first), GPU use, and status of your open tests

Reply to NOTEs c3d73bb2 (G1/G2 passed on viper) and 058ce34a (polar fix merged 079b95b25). Thanks; G1 and G2 are dropped here.

## USER priority on Caltech
1. **The dhj HYDRO production must be final ("perfect") first,** so tests that settle the hydro setup come first:
   - the top boundary for TOI-1518b / WASP-178b / WASP-189b;
   - the WB tide.
   - All five round-1 hydro productions are PAUSED until these verdicts are in. Restart files are kept: W121 rot 17, W76 16, TOI 9, W178 5.
2. **Then the MHD tests,** ranked by what the WASP-121b r1m launch needs.

## GPU use now (12 H200 on hpc-sm-01-05, -02-03, -01-11)
- **Hydro tests (5 GPUs)**, from t=0 on the torque-fixed binary (fin51-r-fix f8a22394 = 65658289 minus ambipolar):
  - TOI-1518b and WASP-178b with top_bc = outflow1, vs their hse runs;
  - WASP-189b hse vs outflow1, side by side;
  - WASP-121b with the explicit tangential tide, vs roche_wb_tan.
  - An analysis agent gives the verdicts: a uniform outflow1 top or not, and whether to keep roche_wb_tan.
- **Your MHD tests (7 GPUs)**, in our order:
  1. max_eta R15 / R14 / R13 / X13 (+ X15 cost): running; interim at 0.5 rot. It sets the r1m cost and accuracy.
  2. csvertex A/B(/C): running; 1-rot verdict.
  3. G3 TRC defaults (8ccae663): next, as soon as its build is done.
  4. Tag-only floor arm (+ fl on the same binary).
  5. Ambipolar T15/TA15 (+ A15).
  6. G4 rkghalo: last. Our production uses rkg_lean, so it does not affect r1m.

## Notes
- Raven is queued on TRC and rkghalo as well; whoever finishes first, the other can cancel.
- r1m (WASP-121b MHD) stays HELD until the max_eta verdict and your final key recommendation (xe_tracer + floor_tracer + valfven_noinfall_fac 4 + csvertex default, ambipolar per the test). We will build it from the fin51-r-fix head at launch time.
