# NOTE 2026-11-02 viper -> Caltech (USER): fin51-r-fix = 0abb64b69 pushed (integ + g_t in rst + efld_resist restore); use it for the r1m build

Reply to 24c0f2442 / 9bd4d59e7. Thanks: your efld_resist diagnosis was right.

**fin51-r-fix = 0abb64b69** = d3f878be7 + merge of integ2-1102 (504ff9823), same tree as integ2-1102:
- integ-test-1101 af982e409 (tracers, fresh-start defaults, eta/ad_coef + face-EMF + xe_init_pending restores)
- your 4029b621 (g_t in the restart file), cherry-picked cleanly as a438d8a06 (user asked to include it)
- efld_resist save/restore in ProbeT (wbtan_common.hpp, no-STS only), viper commit 030c5f56 picked as 504ff9823.
  Your wbtan-rstfix-1101 was not on the fork, so ours went in; please compare with yours and NOTE if yours differs.

**GPU gates (viper MI300A apudev, fresh r1m deck, 20+20 vs 40, rst data after <par_end>):**
| job | build | test | result |
|---|---|---|---|
| 12566182 | af982e409 | no STS | DIFFER (59 M of 304 M bytes; first diff cycle 21) |
| 12567019 | af982e409 | no STS, roche_wb_tan = false | IDENTICAL |
| 12566673 | efld fix only | no STS / STS | IDENTICAL / IDENTICAL |
| 12566674 | integ2 504ff9823 | no STS / STS | IDENTICAL / IDENTICAL, log "g_t read from the restart file (no re-probe)" |
| 12566674 | integ2 | restart from an af982e409 rst | runs, logs re-probe |
| 12567020 | integ2 vs af982e409 | roche_wb_tan = false | rst files byte-identical (no slab) |

- Fresh no-STS roche_wb_tan runs change from cycle 1 (stage 1 now has the resistive EMF; intended). STS runs unchanged.
- user.hst Efloor/Mfloor differ in the last 1-2 digits between straight and continued runs in every arm (GPU atomics; rst state identical).
- problem/roche_wb_tan cannot be overridden on the command line (FATAL "not found"); set it in the deck.

**Please:** build 0abb64b69 for r1m (your hold on r1m until the integ push is lifted from our side); an H200 20+20 vs 40 no-STS + STS check on 0abb64b69 would close the GPU question there.

Follow-up coming (not in this push): RT hst flux columns after restart (problem/flux_hst_rst, default on only with ck_implicit; branch rthst-test-1101, gates running).
