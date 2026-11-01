# NOTE 2026-11-01 viper -> Caltech: the four dhj planet packages (WASP-76b, WASP-178b, TOI-1518b, WASP-189b) + Roche-top table

Completes the reply to TASK-2026-11-01-caltech-viper-more-dhj-planets (bed9c542e). WASP-121b 1x is in
NOTE-2026-11-01-viper-caltech-w121-1x-deck (6243b948).
- Files: docs/handover/dhj_planets_1101/<planet>/. Each README.md has every parameter, reference and choice, and each package has MD5SUMS.
- Pushed: decks, ICs, the new sw_flux files, gate summaries, scripts and grid keys.
- Not pushed: run outputs and restarts (viper /viper/ptmp2/jinma/dhj_planets_1101/<planet>/).
- The ck/CKnX dirs here hold only the new sw_flux files. Use your existing CK1X/CK3X sets (cia, ray, ck, CE_tables are unchanged) and add these.
- Also updated: seidel2026.md (arXiv:2606.22455, submitted 21 June 2026; the first search missed it), roche_top_table.md.

## Plan (user)
- **Round 1: everything at 1x.** **Round 2: 3x** for WASP-121b and WASP-178b (supersolar in the literature) and WASP-189b (check).
- Hydro first (lhllc, cb_sensor default, am_vector = true, problem/am_transport = false). Tide on in your decks.

| planet | Teq (ours / Seidel) | nominal Z (round 1 = 1x) | M_p rel. err | nx1, x1max [cm] | cost/rot vs W121 hfresh | rot recommended | tide/g top | x1max / r_polar |
|---|---|---|---|---|---|---|---|---|
| WASP-76b | 2216 / 2229 | 1x | 1.5 % | 112, 2.40e10 (tide-on anchor) | 1.31x (tide on, measured) | 300 (~24 h, 1 viper node) | 0.57 | 1.33 (open poles) |
| WASP-178b | 2452 / 2398 (\*) | 3x nominal + 1x arm | 7.2 % | 120, 1.8985e10 (3x) / 1.9325e10 (1x) | 4.44x (3x, measured) | 120 (~33 h) | 0.045 | 0.57 |
| TOI-1518b | 2546 / 2491 (\*\*) | 1x (+3x deck) | 25.7 % | 120, 2.0078e10 | 2.41x (measured) | 200 (~30 h) | 0.10-0.14 | 0.85 |
| WASP-189b | 2638 / 2638 | 1x (+3x deck) | 7.5 % | 124, 1.5953e10 (1x) | ~4.2x (scaled, not timed) | 150 (~40 h) | 0.03 | 0.51 |

- \* No published parameter set reproduces Seidel's 2398 K; Teff is the main uncertainty (the KELT SED value gives 2303 K). The user keeps 2452 K.
- \*\* All of the difference is a/R*: 4.109 (Simonnin+2025) vs Cabot's 4.291 used by Seidel.
- Cost is measured on 1 MI300A node, without the tide unless noted. With the tide, add your +52 %-type radial-CFL cost where the bulge is strong.

## Gates (viper)
- `-n` parse; CPU 2-cycle fresh start (r(p_ref) vs the anchor tool to 1e-4 or better; T vs table to <= 0.4 % in active cells).
- A short GPU smoke on apudev: no FATAL, no eos_fail/c2p/fofc/vceil. Floor counts are in each README.
- **WASP-76b is gated with the tide ON with c8aac483 binaries.** The 182477553 binaries ignore roche/roche_a, so **round 1 needs fin51-r-fix >= c8aac483c**. The others were gated tide-off.
- Bug found in the shared gate_ic.py: cubed-sphere corner ghosts were counted as active. The old 1-5 % T maxima came from them; fixed.

## Roche-top table (roche_top_table.md; your ask)

| planet | r_polar [cm] | x1max/r_polar | r(1e-6 bar) substellar t=0 | evolved est. | x1max/evolved |
|---|---|---|---|---|---|
| WASP-121b 1x (deck / roche grid) | 1.5595e10 | 1.379 / 1.431 | 1.637e10 | 1.873e10 | 1.148 / 1.192 |
| WASP-76b 1x (2.40e10) | 1.8034e10 | 1.331 | 1.738e10 | 2.016e10 | 1.191 |
| WASP-178b 3x / 1x | 3.3340e10 | 0.569 / 0.580 | 1.413e10 / 1.427e10 | 1.532e10 / 1.558e10 | 1.239 / 1.240 |
| TOI-1518b 1x (deck / roche grid) | 2.3626e10 | 0.850 / 0.860 | 1.562e10 | 1.720e10 | 1.168 / 1.181 |
| WASP-189b 1x / 3x | 3.0879e10 | 0.517 / 0.508 | 1.286e10 / 1.271e10 | 1.390e10 / 1.361e10 | 1.148 / 1.151 |

- r_polar / r_L1 = 0.64 for every planet; our W121 value matches yours. The evolved estimate carries your W121 rise (+0.166e10 after 1-3 rot, +0.236e10 evolved), scaled by each planet's height.
- **Question:** for W121 1x we get r(1e-6 bar) at t=0 = 1.637e10 substellar; you quote 1.564e10. The 3x anchor and IC give 1.597e10. Could yours come from the earlier 3x / 10x-p_ref setup? Our margins use our higher value, so they are conservative.
- WASP-121b and WASP-76b have open polar tops and need your outflow1 top. The other three have their whole top inside the polar lobe.

## Kepler consistency of roche_a (affects the tide strength)
roche_anchor derives G M_* from omega and roche_a. With literature a, P and M_* from different papers the implied M_* differs:
- WASP-121b +2.6 %;
- WASP-76b +3.3 %;
- WASP-178b +1.6 %;
- TOI-1518b +0.09 %;
- **WASP-189b +14 %** (2.32 vs 2.03 M_sun; the tide is small there anyway).

## Caveats (user: noted, not acted on now)
1. **Short UV is folded into band 10, and band 10 deposits it ~2 dex too deep.**
   - Stellar flux below 0.26 um: WASP-178b 16.5 %, KELT-20 14.4 %, WASP-189b 6.2 %, TOI-1518b 4.3 %, WASP-76b 0.5 %, WASP-121 ~0.5 % (ATLAS9).
   - 1-D RCE tests (viper /viper/ptmp2/jinma/uvband_1101/RESULTS.md): a grey UV band at 1e-6..1e-5 bar heats p < 1e-4 bar by +500..+2200 K and cools the 1-100 mbar Fe layer by 70-86 K for WASP-178b (~50 K for WASP-189b).
   - An LTE photoionization opacity instead deposits it at 1e-4..1e-2 bar and heats; the sign depends on the opacity.
   - So treat WASP-178b, and WASP-189b above 1e-4 bar, with care. No UV band is implemented (CK_NB = 11 is hard-coded; a grey UV chain would be ~2-3 days).
   - Also: Wang+2026 put 2.3 % of WASP-121's energy below 0.258 um (blackbody-like), vs our ATLAS9 fold of ~0.5 %. Not yet checked.
2. **The ck tables stop at 6100 K**, with the hiT2 extension above that. A UV-heated top could exceed it.
3. **Hydro has no drag,** so expect winds faster than the observed Fe winds, especially at the hot end. The trend comparison needs the MHD arms.
4. **Non-ideal MHD estimate** (W121 10 G, offline): at the Fe level and on the dayside / terminator the flow is near-ideal (Rm ~1e4-1e5), so Hall and ambipolar are negligible there (<1e-3 of advection). They matter only on the nightside above ~0.1 mbar (Hall/adv 3-7, AD/adv ~1) and at the top. UV photoionization changes eta by ~1 % at the Fe level, so it is not needed. An offline nightside advected-ionization estimate is running; an ambipolar-diffusion port (Cartesian -> sp -> cs) is in progress on viper.
