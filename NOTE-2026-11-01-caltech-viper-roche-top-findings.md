# NOTE 2026-11-01 Caltech -> viper (USER): tide + top-boundary findings that affect EVERY planet package; p_ref catch acknowledged

Thanks for the p_ref catch (6243b948). Our 3x Roche reference used the 10x NRS1 pressure 1.8e-4 bar; the correct 3x anchor is yours (ap 1.141384e10, grav 1137.7664). Round 1 uses your 1x deck and your 1x Roche anchor (ap = x1min 1.144050e10, grav 1132.4682), re-derived here.

## 1. "Inside the Roche lobe" must mean inside the POLAR lobe radius, not below L1
- The lobe is narrowest at the poles: about 0.64 r_L1 for WASP-121b (polar 1.5596e10 cm vs L1 2.42e10).
- **A top at x1max/r_L1 = 0.80-0.92 is still on OPEN equipotentials over most of the top shell:** 98 % of the shell at 2.1e10 for WASP-121b.
- **Please re-check WASP-76b at 0.80 r_L1** and every package's x1max against the polar lobe radius. roche_anchor.py prints L1; the polar radius is where Phi_tot on the z axis equals Phi(L1).

## 2. But a top inside the polar lobe does not work either (WASP-121b, tide on)
- **The dayside 1e-6 bar isobar rises above the polar lobe radius:**
  - 1.564e10 at t=0 over the substellar point;
  - 1.73e10 after 1-3 rot;
  - estimated >~1.8e10 evolved, since day heating clipped even the Sept 3x no-tide run at its 1.595e10 top.
- **The user requires the 1e-6 bar level inside the domain everywhere.** So the top must stay high (~2.1-2.2e10) and sit on open equipotentials at the poles.
- **Tested x1max 1.49e10 (inside the polar lobe):** a drain of +3e13 to +5e13 g/s in both tide arms. That removes more than the whole atmosphere above 1e-6 bar per rotation. Rejected.

## 3. Hydrostatic top ghosts are wrong on open equipotentials (user point: there is no HSE beyond the lobe)
- The current top ghosts are an isothermal hydrostatic extrapolation in the full potential, with v_r copied, so they allow inflow.
- At 2.1e10 with the tide, the closed substellar cap (tiny effective g) acts as a dense reservoir and feeds -3e10 g/s of inflow.
- **Being built here (branch topbc-1101, key problem/top_bc = hse | outflow1, default hse = bitwise):**
  - pure zero-gradient copy from the last active cell, with v_r = max(v_r, 0);
  - no hydrostatic component anywhere;
  - energy recomputed with the ghost's own potential (etotgrav).
- **The top sponge is checked in the same pass:**
  - whether its band is radius-based, which would damp gas at p > 1e-6 bar on the tilted dayside;
  - whether it damps v_r against the outflow.
- Results follow in a NOTE; the round-1 decks will use the new top.

## 4. The explicit tangential tide has a truncation residual; a well-balanced version is being built
- **The on-minus-off excess is truncation error:** C32 vs C64, tide on vs off, from t=0 to 2 rot (w121hroche_1101/smoke/c64cmp_out.txt). The excess in horizontal velocity shrinks 4-5.5x at C64 at every time: deep shell 6, 6.5 -> 1.4 m/s; near-wall, 9 -> 1.7 m/s.
- **Being built here:** branch roche-wbtan-1101, an x2/x3 tidal (+ centrifugal) source taken as the discrete face difference of the equilibrium p(Phi), behind a key, default off.
- **Zonal winds at 1e-4..1e-1 bar** are nearly the same with and without the tide (e.g. 115 vs 134 m/s at 1e-2 bar, 3 rot). The tide and the top change the flow only above ~1e-4 bar.

## 5. Deep v_r growth = convective onset (FYI)
- It appears in both arms, at 30-190 bar, where the IC is marginally adiabatic (30-70 % of cells slightly superadiabatic, N^2 ~ -1e-9..-3.5e-9).
- The cells are grid scale (m ~ 20 at C32) and grow faster at C64. This is consistent with your "convective below ~8 bar". It is physical.

## 6. Tide cost
- +52 % per rotation on WASP-121b (dt 8.7 vs 13.6 s on 3x). It comes from the radial CFL where the bulge puts ~360 bar / 4600 K gas into the finest radial cells near r ~ 1.19e10.
- Per-cycle speed is unchanged.

## Asks
- In each package, report the polar Roche-lobe radius, x1max/r_polar, and the expected dayside r(1e-6 bar) with the tide.
- Keep the deck's top so that the dayside 1e-6 bar isobar is inside with margin; the outflow1 top handles the open poles.
