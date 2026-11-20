# NOTE 2026-11-20b viper -> orion (USER): clean-cc handover; viper stops clean-cc work, the test nodes are yours

USER decision 11-19 ~17:30 PST: viper stops the clean-cc (cell layout) line and consolidates on your staggered line.
The pHG test nodes vipa1073/1081/1186 (and the successor's, see NOTE 2026-11-20c) are yours. Nothing was deleted,
nothing pushed to clean-cc. All RSG arms below: rsg.01658, 3 nodes x 2 MI300A ranks; P = the live RSG production
(old prod binary). Runs under /viper/ptmp2/jinma/ccfull_1117/runs/<arm>.

## (a) Base clean-cc on the lbot base is stable and correct, but about 50x production cost

lx_baselong (2eda97438 = 65f075a38 + our keys, all off; prod keys), 1000 cycles, no collapse:
dt mean 94.5 (P 931), 116 retries (112 NON-CONVERGED, 4 eint<=0, ~11 per 100 cycles), 72 linear its/stage,
17.5 s/cycle (P 3.40) -> about 50x P per simulated second. Ledger -9.7e-5 (P 4.2e-5), L_top 6.50e38 (P 6.45e38),
L_bot 6.36e38 (P 6.85e38). The r9 base (245ea4418) lacks 65f075a38 and gives L_bot -1e40: do not use it.

## (b) Frozen pass-1 opacity: the converged state depends on the first iterate

With clean_opac_passes = 2 (our prod keys) the opacities and everything CleanCellFaceSpeeds derives from them
(al, RC/TH, AF) are evaluated from the first-iterate T and the pass-1 T, then frozen: the converged state is not the
fixed point of its own opacities. With clean_opac_passes = 1000 (every pass), xt_cons and base give the same states
(identical dt history, ledger -1.1e-4, L_top 6.553e38 over 21 cycles), but the first cycle needs 7 cuts and the
21-cycle dt mean is 20.2 (base 37.6). Ladder 3 / 4 / 6: dt mean 12.5 / 30.4 / 12.4, not monotone (still frozen
states). A 1000-cycle run of the every-pass solve was stopped at its first cycle by this decision: untested long.

## (c) First backward-Euler step: E < 0

- Prod keys: 510583 E<0 cells at dt 1000, deep interior (rho 1e-7..1e-5, inner ~30 % of x1, block interior).
  clean_xt_cons (first iterate with the right-hand side's exchange: G = 1 - beta^2/2 and the beta terms, a = v + Dv)
  removes them: 425 cells left, all in a thin zone i ~ 407-450 (rho 1e-10..1e-9).
- At the worst cell (gid 7 i 408 j 46 k 47): chi d 0.08 radially, h c chi 44 at dt 1000, a 170x opacity jump across
  its k-1 face; E-E row Z-signed and dominant (+100), but the E row's central y columns 0.5 c n are ~415 per unit y per
  radial face against a diagonal of 708, and the y rows couple back with sum s|dfY/dE| = 247 vs 1 + hc chi = 45; the
  all-low fallback row is Z but not row-dominant (-23.7). Work columns negligible (5e-11).
- clean_thin_low = 1 (faces with chi d < 1 and hc chi > 1 on both sides start low-order) removes these too:
  0 E<0 cells at dt 1000 (with xt_cons). clean_thin_ydonor (your lagged donor + d/(1 - l), restart-safe) leaves them
  (421) but fixes the stage ledger. clean_low_face_d (one D_nn per low-order face) changes nothing.

## (d) xt_cons long runs collapse through the write-back eint <= 0

lx_long (lbot + xt_cons + rst_k1=2) and lx_xtthin (lbot + xt_cons + thin_low=1) both fall to dt 4-15 s after about
100 cycles: from then on almost every retry is eint<=0 (the exact exchange of the write-back has no root with
eint > 0: rhs = E + eint - E_tr/(1 + h c kE) <= 0, i.e. transport drains more than the cell holds), about 10 per
100 cycles; ledger +0.06..+0.10, L_top 20-28 % low. thin_low does not prevent it. Over 21 cycles xt_cons alone has a
per-step ledger residual of +-1-5e-3 (base 1e-4) and E_old<0 at stage starts; both disappear with thin_ydonor /
thin_low or with every-pass opacities.

## (e) Code facts for your line

- rad_m1_clean.cpp:861: a write-back with cln_neint > 0 set dtr_fail without dtr_why (logged as reason "?").
  Fixed locally (dtr_why = 4, eint<=0).
- Deadlock: CleanCellVelocityUpdate set cln_fe / cln_fv from per-rank maxima (since b44bd47cd), so the Picard exit
  was rank-local; fixed with one MPI_MAX of {dvmax, vmax, femax, fvmax} + non-finite -> inf (f64af1747, 2f8f57445,
  06c9ee3cc). Also in clean-cc head.
- time2_k1_redist acts under scheme = clean (Time2FormStage applies it to the clean stages) although the clean
  preflight listed it as without effect; time2_face_fct is forced off under clean.
- A cubed-sphere energy balance must take the boundary fluxes through CleanCellRows (the solve's own face areas);
  pcoord area x f0x1 is off by a constant ~0.6 on cs (fixed in our diagnostic, fdc30cf54).
- Old prod restarts carry 37 slope channels, the clean run keeps 136, so clean starts with BE at full dt;
  clean_rst_k1 = 2 imports K1 without the face channels (=1 with faces gave ledger +0.2).

## (f) Local branches (not pushed) and run dirs

| branch (worktree) | tip | content |
|---|---|---|
| clean-cc-lbot-pos-viper (/viper/ptmp2/jinma/ccfull_1117/wt_pos) | fdc30cf54 | 65f075a38 + xt_cons, rst_k1, lin_continue, deadlock fix, diag_neg/diag_cell/ledger, bcg_fallback, low_face_d, thin_low, thin_ydonor, k1_* preflight fix |
| clean-cc-lbot-xtcons-viper (wt_lbx) | 2eda97438 | 65f075a38 + the r9 work |
| clean-cc-r9-xtcons-viper (wt_r9x) | e08f39ee3 | r9 base + the same |
| clean-cc-coldstart-diag-viper (/viper/ptmp2/jinma/ccint_1118/wt) | 2f8f57445 | int base + cold-start diagnostics, deadlock fix |
| clean-cc-int-viper (ccint_1118/wt_int) | 06c9ee3cc | int + deadlock fix |

Binaries /viper/ptmp2/jinma/ccint_1118/bin/he_star_m1_<hash>; run dirs /viper/ptmp2/jinma/ccfull_1117/runs
(lx_baselong, lx_long, lx_xtthin, rp_*, st_*, op*, led_*, det_*, ab_*); key files ccfull_1117/keys; restarts kept in
runs/*/rst (lx_long, lx_xtthin, lx_baselong, rsg_long_k2, det_D1, smk_long_k2: about 2 TB, for the user to decide).
The MI300A UCX result (default UCX bitwise over 3 x 100 cycles + a restart pair, -21 %) is in NOTE 2026-11-19q
and the follow-up runs det_*.

-- viper
