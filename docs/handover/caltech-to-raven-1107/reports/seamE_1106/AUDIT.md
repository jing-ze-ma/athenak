# Audit of the HG seam recipe (seamE-1106), independent, 11-06/07

Scope: r1 (sd10) = cs_cross_mono + Newton off + guard 4 + tcap 1.1 + time2_k1_redist 0.1; l1 (sd12) = r1 + implicit_fb_local rounds 4.
Control hg/off, same rst 139. Compared per sim time, bin_double dumps 1741-1743 (t0 + 0.33, 0.67, 1 d; 1742/1743 differ by 100 s between arms).
Scripts and raw output: /resnick/groups/carnegie_poc/jingze/audit_r1_1106.

## Noise
- r1n (implicit_tol x0.999) is NOT a noise member: bitwise with r1 over all shared rows. tol_eff = max(tol, min(0.1 r0, 1e-4)), so tol is inert. Same check shows implicit_seam_dump is state-inert.
- kp6n (he_gm +1 ulp) and r1p (rst_pert_erad 1e-10) are valid. r1p has dump 1741. L_top noise <= 3.5e-4 relative (r1 vs r1p), 7e-4 (kp6 vs kp6n).

## "Mode 7"
No arm has cs_cross_mono + mode 7 + fb_local + redist. r1, l1, r2, sdk are Newton OFF and do not set the key; the guard mode acts only on the Newton term, so it is moot there. Mode 7 was measured only without guard 4 and redist (fixg7): vmax 1394 km/s, NOTPOS 87 in t0+39000. The Newton-on combination is UNMEASURED.

## Better than control (measured, >> noise)
| | off | r1 / l1 |
|---|---|---|
| NOTPOS | 321 to t0+88000 | l1 12 to t0+130000 |
| vmax i100-600 at 1743 (km/s) | 543 (760 at 1742) | l1 135 (r1 131 at 1742) |
| v > 300 km/s cells | 2-4 | 0 |
| \|F\| > cE cells | ~1000 | ~450 |
| dfloor, vceil first log window | 458, 234 | 0, 1 |
| m1_pos_flr, m1_acl_e | 51, 245 (first window) | 0, 0 |
| energy residual / L | +4e-4 to +1.3e-3 | r1 up to +1.1e-4; l1 <= 7e-6 |
| wall per 1000 sim s | 26-32 s | 20-23 s |
| Picard mean | 8-10.7 | 5.5-6.1 |

Mass: guard 4 arms close dM = -int Mdot_top to ~1 %. Guard 0 arms (off, fixg7, s2) gain ~4.7e19 g in 15000 s (density floor, inferred; 7e-13 of M). Seam imprint (median |ln x/shell median|, seam vs interior; odd-even at seams): unchanged at 1741, no new seam pattern.

## Changed versus control (measured, >> noise; good or bad unknown)
1. L_top +1.0 % (1.354 vs 1.341e38, steady to 1 d; noise 3.5e-4). s1, s2, fixg7 give +0.7 %, so most is cs_cross_mono (or a key common to them); redist adds +0.3 %. fb_local adds nothing above noise. Net heating dEtot/dt 6.4 -> 4.9e36 erg/s (-23 %); part of that is the control's own energy leak, 1.3-1.8e35.
2. Near-photosphere structure GROWS: shell E at i 280-320 (113-117 Rsun) +1.45 %, +2.1 %, +2.5 % at 0.33, 0.67, 1 d; eint +0.46 %, +0.9 %, +1.1 %. Layer i 160-200: E -0.03 %, -0.06 %, -0.08 %; eint -0.08 %, -0.13 %, -0.17 %. Not saturated in 1 d. r1p shows the same, so it is not chaos.
3. Thin-layer grid-scale E noise: third difference of ln E at i >= 360 (r > 123 Rsun), interior, +50 to +100 % (0.024 vs 0.016) at all three dumps, in r1, r1p, l1; NOT in fixg7. So it comes from guard 4, Newton off, tcap, redist or the larger dt, not cs_cross_mono. Above the accuracy region by the user's rule; amplitude ~0.5 %.
4. Deep KE (i 0-100, 7-38 Rsun): -0.04 %, -0.21 %, -0.41 % at 0.33, 0.67, 1 d (vrms -0.09 %), growing and unexplained. fixg7 -0.13 % at 0.33 d, so not from guard or redist alone. Deep E, eint, rho within 7e-5, smaller than control's own breathing. Only one member beyond 0.33 d.
5. dt 100-111 s vs 77-87 s (1.3x): time error is larger per step. Mdot_top about -10 %.

## Code review (wt_seamE, line numbers as read at 8aa1a94be + neighbours)
OK by reading: x1 is one block on HG (6 x 1 x 2 x 2). redist columns are race-free on GPU; K1_E sum V conserved, E channel only, gas/F/momentum untouched; K1 is re-formed each step, and FSAL storage and restart hold no redist state. tcap default 2 is bitwise. cs_cross_mono has only 2 call sites (M1CsTanDer), all in ImplicitTransverseTerms, one face value, so telescoping is kept and the flux is single-valued; it is lagged RHS only. Guard 4 acts in the fluid halo only: no mass or energy added, and hst confirms. fb_local: mask grown through face ghosts, VlatTT() used by assembly, stencil and face update; tau_ten untouched; reset each solve.

Problems:
- rad_m1_time2.cpp:825 `sg = -1` first: donor is always tried deeper first. Replay: 88 % of moved energy comes from the deeper side, so a systematic outward transport. Prefer a symmetric split or say so in the comment. Small (below).
- rad_m1_time2.cpp:827 silently skips other blocks: a radial split is layout dependent with no FATAL. Add one.
- rad_m1_implicit.hpp (M1CsTanMono comment) "A smooth field never reaches the clamp" is FALSE. ord.py: on sin(3x), 0.4-3 % of seam positions clamp; Linf error drops to FIRST order, L1 stays 2nd (constant 1.5x). On real data, 9-17 % of seam-adjacent cells clamp at ALL depths including i < 100 (E used as proxy for D E; inferred): median change 4e-4 E deep.
- rad_m1_implicit.cpp ~13381 comment "the fixed point does not move" is FALSE as stated: tcap 1.1 vs 2, stop tol 1e-8: E median 9e-8 but max 1.1e-2 (i 240-360), C16 CPU. Path dependent somewhere; local.
- Newton off: "converged state is the same" not found on C16. Off vs on, both tol 1e-8: E median 3.6e-4, max 6.4e-3, at i >= 240; deep 1e-13. Stop error itself is NOT worse with Newton off (off 1e-9, on 3e-8 at i >= 240; deep 1e-14), so fewer passes do not mean less converged. Cause of the different fixed point not found (opacity lag passes suspected, not checked). Real state s2 vs fixg7 L_top differs 4e-4, about noise.
- Redist: stress th=0.999 is 97 % "unpaid" (rank 0: 212924 cells, 207051 unpaid), so that case exercises few paid moves. th=0.1 never fires on C16 (a vs g bitwise), so C16 cannot test the real path.

## Size of redist (INFERRED, offline replay of K1 in rst 139, dt 100, control-lineage state, first step)
Cell count depends on the k1_stiff move, whose k dt test was not evaluated: 19 cells (none moved) to 43291 (all moved; 6e-4 of cells), i 162-268. Upper case: gross 3.6 % of L dt per step (x0.7 for the final weight b1), median 0.44 E_cell, p95 3.4 E_cell. 72 % over 1 cell, 18 % over 2, tail to > 10; 0.8 % unpaid. Net radial transport peaks at 1.1e-3 L at i 171. Real counts in r1 were not obtained: the k1_redist report prints only at a clean end and arms are killed by timeout. Print it periodically.

## Tests
- Restart bitwise, CPU, C16, all r1 keys, th=0.999, 4 vs 2+2 cycles: IDENTICAL (redist fired). fb_local not covered. No GPU restart test.
- Layout: he_hg_xmono b2x2 no worse than he_hg by numbers (dens 3.0e-6 vs 3.5e-6, m1_e 1.7e-8 vs 1.6e-8), but the suite verdict is ERROR because no he_old is in bins, so NOWORSE was never computed.
- Layout OPEN: he_hg_k1rd b2x2 dens 3.9e-5, ener 2.2e-4, m1_e 9.0e-7 at cycle 1; baseline 3.5e-6, 1.9e-5, 1.6e-8; baseline 1-ulp spread 1.6e-5, 8.7e-5, 2.9e-8. Cycle 1 is backward Euler, so redist and the predictor cannot cause it (inferred): guard 4 or Newton off. Split running in li3. No layout case for fb_local exists.
- Keys-absent GPU gate: GATE.txt is sd4 only; sd10 onward not re-gated.

## Unchecked
Runs beyond 1.5 d; any resolution or dt convergence test saying which L_top is right; fb_local restart and layout; dt_retry with redist (no retry occurred in these arms; checkpoint order relative to FormStage not read); r2 (redist 0.01), sdk (sdirk2), time2_k1_fct arms; manufactured-solution order in the full code (only the formula tested).

## Verdict
Safe to keep testing; do NOT use for production yet. No conservation or seam harm found; stability, holes, floors and cost are clearly better. Fix or settle first: (a) layout excess of guard 4 / Newton off, (b) a fb_local layout and restart gate, (c) sd-current keys-absent GPU gate, (d) a longer r1 vs r1p vs off run for the growing near-photosphere E and deep KE, (e) false comments and the x1 FATAL.
