# NOTE 2026-11-01 Caltech -> viper (USER): hydro setup settled (uniform outflow1 + WB tide), the L_z leftover = start-up drag, MHD WB tide merged, max_eta final, HG and drag tests

## 1. Hydro round-1 setup: decided (USER)
**Top boundary:** problem/top_bc = outflow1 for ALL planets, including those whose top is inside the polar lobe.
- A/B hse vs outflow1 from t=0: TOI-1518b to 7.6 rot, WASP-178b to 4.3, WASP-189b to ~5, plus the pre-fix TOI/W189 pairs.
- dt identical (TOI 9.36/9.36 s, W178 9.28/9.28); no dips or failures.
- 1e-6 bar stays inside by 4-5 dex; mass, energy and L_z agree to 3 digits.
- Winds and T at p >= 1e-2 bar differ by far less than the run's own 1-rot variability.
- For these planets outflow1 acts as a closed lid: Mdot_top = 0 after rot ~1.5, where hse lets ~1e9 g/s fall in. dfloor events drop 2x.
- The W189 hot streak at lon +150 (0.1 bar) is in both arms, so the top does not cause it.

**WB tide:** keep problem/roche_wb_tan = true (with our torque fix 65658289). W121 WB vs explicit, matched at 8.7 rot:
- the explicit tide drives spurious deep flows (1 bar rms 413 vs 202 m/s at rot 8; 58 vs 10 m/s at rot 0.5);
- the deep state is systematically off, by ~3x the 1-rot variability at 1-10 bar;
- the tide torque has the opposite sign (explicit spins down, WB spins up);
- the upper atmosphere is the same.

**The rot 0-5 L_z offset (-5.1e-4 Iw0·Ω) is the hard-coded start-up Rayleigh drag** (pgen ~line 5594, "initial relaxation"):
- it damps all momentum components, tau = P/10 for t < 2 rot, a log ramp to P at 5 rot, then off;
- -Lrel/tau accounts for the offset to < 1 % in W121 (WB and explicit), TOI and W189;
- after rot 5, L_z closes against Ttide alone to 0.2-0.3 %.
- Branch relaxdrag-1101 (9d7db7fe, not merged) adds problem/relax_drag (default true = bitwise), relax_drag_rot and a Tdrag hst column.
- **USER is testing a start WITHOUT the drag** on all five planets (in progress). Early W121 at 3 rot: stable, dt unchanged, jet 3.4 km/s vs 0.2 km/s with the drag, L_z closes from t=0; but e/T-floor events are ~10x higher early. Verdict to follow.

**Hydro productions stay PAUSED until the drag verdict.**

## 2. MHD WB tide merged: fin51-r-fix dd781f22
- problem/roche_wb_tan for MHD (the probe uses B = 0). Shared helper src/pgen/wbtan_common.hpp: one probe and one force path for sp and cs, hydro and MHD.
- **Gates:** key off is bitwise (GPU r1m C16 + CPU suite); key on, B = 0, one stage: |v_h| 7.5e-14 cm/s vs 65 explicit.
- **r1m 0.5-rot smoke** (WB on vs explicit):
  - deep flows 2-4x smaller (1-30 bar: 8.5 vs 36 m/s);
  - Ttide varies (not frozen);
  - dt, RKG stages and div B unchanged.
- **L_z closes in both MHD arms** to <= 2e-7 Iw0·Ω once the start-up drag is included. There is no Maxwell-torque column; the residual bounds the boundary Maxwell torque at ~4e-7 Iw0·Ω/rot.
- Our r1m deck now has roche_wb_tan = true.

## 3. max_eta final (NOTE 8ae2b6ff)
- R14 is borderline: 0.61x cost, 5 RKG stages, T < 1 %, B_rms < 1.2 %, but E_B drifts +0.46 -> +0.85 % over 0.5 -> 1.5 rot.
- R13/X13 fail (the cap, not STS). Explicit loses at every cap.
- **Our recommendation to the USER: keep 1e15**, because of the growing E_B drift. USER to decide.

## 4. Status of your other open tests
- **Done:** G3 TRC (PASS, 6bbfb530), tag-only (A, 7437a9d4), G4 HALO v1 (bitwise but no saving, 96b69569), csvertex (PASS, 5e1ee1c4).
- **Running:** ambipolar A15 / T15 / TA15 (interim at 0.5 rot soon); the halo v2 vs rkg_lean timing (TASK bf6ae9e3) on a dedicated H200.

## 5. FYI
- **HG (He giant):** the 96-160 d big expansion and fallback is a start-up readjustment after the MLT -> resolved-convection switch. Infall does not quench convection; it re-compresses the envelope and damps the pulsation. A physicality check of the late phase (structure vs 1-D, thermal time, the ~40 d oscillation) is running.
- **Accretor moves to the cubed sphere** (USER: sp pole artifacts). csacc-1031 is rebased on fin51-r-fix; the impact point sits in panel 0, ~33 cells from the nearest seam. The tip tests and the cs accretor WB-tide port are in progress.
