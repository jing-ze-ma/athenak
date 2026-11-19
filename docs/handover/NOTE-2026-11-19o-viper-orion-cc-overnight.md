# NOTE 2026-11-19o viper -> orion (USER "tell orion"): Raven pHG262 is running HG prod; clean-cc RSG overnight results

Times PST.

## 1. GPUs
- **Raven pHG262 32307710 (6 nodes) is RUNNING HG production** (old prod run.cfg, from hegiant.00186, since 00:41).
  It is healthy: A100s at 75-79 %, hst t = 6.530e7 s, dump 02267, rst 00197. squeue shows only step .1
  `gpumon_wrap_raven.sh` because the run is inside that wrapper. Please do not start steps in 32307710/32307711.
  (run.log in raven_gap_1105/hg stopped updating at 04:08 PST although the run continues; viper is checking.)
- Viper pHG262 13199032 is no longer in the queue; successor 13255126 PENDING (yours).
- Viper (all busy): pRSG = RSG prod, pr136f = R136 prod, pHG = clean-cc tests (3 nodes) + T19L6.4 (2 nodes, new)
  + dhj W178/W189 30G (vipa1208), pool = dhj W76/TOI. T9L6.2 stopped (user).

## 2. clean-cc on the full RSG (rsg.01658), viper pHG
1. **All r9-based arms lacked 65f075a38** (the lbot boundary-stress fix; r9 = its parent 245ea4418). That was the
   L_bot -1e40 of every r9 arm. Port onto 65f075a38: 2eda97438 (local), bitwise = gb_lbot.
2. **xt_cons + rst_k1=2 on the lbot base, 21 cycles:** dt 103 vs 38 s mean, but 3x cost per cycle -> same throughput
   (7.9 vs 8.5 sim s per wall s); ledger -4.7e-3 vs -4.7e-6 (88 % of it in the first staged step on the imported K1);
   nreal 682 vs 692166.
3. **Long arms collapse:** r9 + xt_cons (1000 cycles) and lbot + xt_cons (stopped at 349 cycles): dt 4-13 s after
   ~cycle 100 vs prod 1000; ledger +0.10; L_top 18-28 % low. Retries are mostly the write-back eint<=0 rejection
   (rad_m1_clean.cpp:861, now labelled reason 4) plus NON-CONVERGED.
   The reference long arm with prod keys on the lbot base (no xt_cons) holds dt 92-116 s mean per 100 cycles, no collapse (running).
4. **Where E<0 comes from (dg_negloc, first BE at dt 1000):** 424 cells, ALL block interior (no face/edge cells),
   rho 1e-10..1e-9, r 1.30-1.36e14 cm (i 407-450). **Not near-vacuum block/panel edges**, so no E floor.
   The solution's h div F is 60-1000x E_old there. |F|>cE (nreal) is a separate outer near-vacuum family (i 466-476).
   clean_bcg_fallback=4 has no effect.
5. **Stencil at the worst cell (gid 7, i 408):** chi d 0.08 per cell but h c chi 44 (thin per cell, thick per step).
   E-E block is a Z-matrix and row-dominant (+100). The central y columns 0.5 c n give ~415 per unit y per radial face
   vs diagonal 708 (sign-indefinite). In the all-low pass the row is Z but not dominant (-24 at dt 1000), because the
   per-cell h+/- and D_nn differ across the face. A 170x chi jump sits across one k face.
6. **Prototypes** (local branch clean-cc-lbot-pos-viper, off by default): `clean_thin_low=theta` (faces with chi d < theta
   and h c chi > 1 go low-order: no central y columns) and `clean_low_face_d` (face-valued D_nn in the low-order
   Jacobian). Test = 21-cycle RSG with E<0 = 0 and eint<=0 = 0, then the IDORT/linear-wave order gate.

**Question:** did stagger (or your CPU runs) see the thin-per-cell / thick-per-step regime break positivity,
and how did you limit the y coupling there?
