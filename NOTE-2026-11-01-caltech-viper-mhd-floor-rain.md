# NOTE 2026-11-01 Caltech -> viper (USER): MHD floor rain, STS phase 0, one-way top, round-1 status

## 1. MHD floor rain (please check your s10 run)
- In 10 G resistive MHD, the Alfven-speed density floor (valfven_max_cgs 2e6) props up gas in the non-hydrostatic top. That gas then rains down at v_r ~ -3 km/s into p > 1e-6 bar.
- Mfloor in the hst is a RATE (g/s). In our r1m smoke it was ~2e15 g/s and rising with B, which is ~15 % of the mass above 1 mbar per rotation. Your s10 most likely has the same, ~1e15 g/s, and it is not flagged anywhere.
- Fix: new key <mhd>/valfven_noinfall_fac (branch mhdfloor-1101, being merged into fin51-r-fix). It removes the infall momentum in floored cells.
  - fac 2 cuts Mfloor 33x, to a flat ~5e13 g/s, and raises dt 11 %.
  - B growth and everything at p > 1e-5 bar are unchanged (1-2 %).
  - fac 4 and 10 are being tested now.
- MHD winds are weaker than hydro because of magnetic drag, not because of the floor mass.

## 2. STS phase 0 (merged fin51-r-fix 899153cb)
- New keys in <mhd>: rkg_dt_safety = 0.5 (per-cell eigenvalue bound, 1.07x the true value) and rkg_lean = true (skips C2P, the U halo exchange and the dt computation in the intermediate RKG stages; uniform cs only).
- r1m: 175 -> 64 ms/cycle (STS 32 -> 14 stages). 300 rot: 18.6 -> 6.8 d on 1 H200.
- Differences are at the cfl-twin noise level.

## 3. One-way top (merged fin51-r-fix 2c1da681)
- Key problem/top_bc = outflow1. Ghosts are copied with v_r = max(v_r, 0) and a diode at the top face; no hydrostatic component. MHD ghost B: tangential copied, normal from div B = 0.
- Used where the top is on open equipotentials: WASP-121b and WASP-76b.
- Kept the hydrostatic top (full potential including the tide) where the whole top is inside the polar Roche lobe: TOI-1518b, WASP-178b and WASP-189b.
- The top sponge is kept.
- Correction to our earlier r(1e-6 bar) numbers: they came from a 3x anchor that wrongly used the 10x p_ref 1.8e-4. The correct 3x anchor is ap 1.141384e10, grav 1137.7664.

## 4. Round 1 (1x) on Caltech, 1 H200 per run, all to 300 rot
Status ~08:00 PST 11-01 (tlim 1e9; we stop each run by hand at 300 rot):

| run | settings | rot so far |
|---|---|---|
| WASP-121b hydro r1h2 | roche_wb_tan + outflow1 | ~91 |
| WASP-76b | outflow1 | ~91 |
| TOI-1518b | hse top | ~57 |
| WASP-178b | hse top | ~29 |
| WASP-189b | hse top | ~26 |

- WASP-121b MHD r1m launches after the fac smokes, with outflow1 + phase 0 + noinfall.
- Deck changes: cs_seam_flux = shared is now pinned explicitly (seam A/B: positive_donor == shared). BVALS_CC_SPLIT=1 is exported (3.6 % faster, bitwise).
- Open: L_z does not close in the first 0.2 rot (~-1e-4 Iw*Omega). We are adding an hst column for the applied tangential torque.
