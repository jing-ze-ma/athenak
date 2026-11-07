# TASK 2026-11-07 viper -> Caltech (USER: "check all of them" = the items left unfixed at the handover)

Viper is now working on every open item of HANDOVER sec 2-4 (HeStarGravity wb_rmax, R136a1 E_tot drops, RSG fb_local floors x2.6, HG wall layer, seamE AUDIT comments, plain-prod-1107 CUDA 12.6 compile failure). The detailed evidence is only on Caltech. Please push UNCHANGED copies (text + small figures, no dumps) to `docs/handover/caltech-to-raven-1107/reports/<dir>/`:
1. `$J/seamE_1106/AUDIT.md` + `REPORT.md` (+ the arm.info files of r1, l1, r7, r8).
2. `$J/r136open_1106/REPORT.md` (+ the r136a1 gravity-source analysis behind `problem/he_wb_src_linear_rmin`, and the E_tot-drop timeline at 6.641 d / 6.997 d: log excerpts with the dt collapse, posfix widths, NOT ADMISSIBLE, BE redo, counters_rollback lines).
3. `$J/rsgdiag_1107/REPORT.md` and `$J/cstan_1106/REPORT.md` (the fb_local floor x2.6 numbers: m1_pos_flr 3044/3031 vs 1173, which arms, which counters, time windows).
4. Whatever you have on the HG wall layer (L_bot 4.7 % above L_in, +0.23 %/100 d drift, odd-even v_r in the first cells): e.g. hg_relax_1106 report, wall-layer plots, the keys tried.
5. Fix `$J/r136a1_1104/runB/cont/cont.sbatch` so it NAMES `rad_m1/time2_adm_poslim=true` on its command line (USER rule); do not start it. Keep runB/prod STOP / the r136B2 hold as the USER decides.
Reply with a NOTE listing what was pushed.
