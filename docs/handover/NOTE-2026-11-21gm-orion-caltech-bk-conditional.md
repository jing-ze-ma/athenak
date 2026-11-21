# NOTE 2026-11-21gm orion -> caltech (coordinator): BI_bsgtau1000 -- keep the ray-cut branch tag; CONDITIONAL arm BK_bsgtau1000thin0

1. **For BI_bsgtau1000 (and BJ):** raycut.txt must keep the **full** "CLEAN ray cut cycle C shell A -> B (want N, tau_top|R)" lines, including the "(want N, tau_top)" or "(want N, R)" tag. The NOTE fn sed strips it; please use plain `grep 'ray cut cycle'` for these arms.
2. **Conditional, no need to ask us:** if BI's log shows **any** ray-cut line with the **", R)"** tag after cycle 0, the thin-cell / gradient branch is binding instead of the tau cut. In that case start **BK_bsgtau1000thin0** on the next free node right away, with the same settings as BI (B_v203, fresh, tlim 1.12e6, same outputs, m1prof 0.5 d + end, full raycut tags):
   - deck **bsg_diag_tau1000thin0.athinput** (committed here, md5 19e58e11) = bsg_diag_tau1000 + `vet_col_lat_cut_thin = 0.0`, so the thin-cell walk is off and the cut is pure tau.
   If BI shows only "tau_top" tags, BK is not needed.
   Tell us in your arm line which case happened (first cycle with ", R)").

-- orion (coordinator)
