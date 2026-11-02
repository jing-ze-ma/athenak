# NOTE 2026-11-02 Caltech -> viper (USER): prod-1102 = the single production branch; please run the HIP bitwise gate

Branch `prod-1102` @ afa17af29 on the fork. It is fin51-r-fix c38bd1981 (your head, merged in) plus the Caltech 11-02 branches below.
Every new key is read only when named and is bitwise when absent; the one exception is `<hydro>/split_bval_fast`, default ON (hydro only, gated bitwise).
Caltech did not push to fin51-r-fix.

## Merged
- cbcont-1102 34ec05e93: `<hydro>/lhllc_cb_cont`, `<mhd>/lhlld_cb_cont`.
- ckhesdirk / cbck-test-1102 09f3a03c7 and ck2nd-1102 b028ee690: `problem/ck_time_scheme = hesdirk2`, `ck_t2_start = probe`, `ck_t2_seed`, `ck_impl_xpass`, `ck_cad_drift`, `ck_vef_dt`, ...
- rkgfuse-fin51-1102 6f7aba385: `<mhd>/rkg_strang`, `<mhd>/rkg_strang_fuse`, `problem/ck_strang`, `<mhd>/tracer_strang`.
- ckbval-1102 c674cead5 (incl. f5d86e6de): `<hydro>/split_bval_fast` default ON, `<mhd>/split_bval_fast`, `<hydro>/split_c2p_once`, `<hydro>/c2p_idempotent`.
- lhllc2nd-1102 68ebbe44d (incl. lwave-nghbr-1102 49e677563): test keys `<time>/dt_fixed`, `<outputN>/bin_double`, balsara hydro mode; `lhllc_c1` / `plm_c1` (not recommended, default off).
- csaccmerge-test-1102 34c7b76e8 (= your 24d044fc1) and fin51-r-fix c38bd1981.

## Caltech changes in the merge
1. FATAL for `ck_time_scheme = hesdirk2` together with `problem/ck_strang` (deep_hot_jupiter_rt.cpp).
2. 54e74ee0b: RyPerFinal frees wbt_ga2_/wbt_ga3_ (you have it as 4c5f540f4).
3. c5f1d72f2, **ck_t2_seed restart fix**. The hesdirk2 stage seed was capped in CkWarmSeed against ck_ei, the previous call's last iterate.
   ck_ei is not in the rst, and the restart rebuild passes overwrite it. So a restart whose first full call is seeded was not bitwise:
   fresh r1m 300 straight != 150 + 150; r1m restart split at 8967 also differed. Bisect: identical without ck_t2_seed.
   The cap (ck_impl_dtmax x e) now uses the cell's current internal energy and is applied where the seed is formed (DhjCkT2Stage).
   This changes ck_t2_seed results vs ck2nd-1102. Key-absent paths are unchanged.
4. dhj WB tangential tide restart fixes (wbtan-rstfix 0f9dc36ec, emffix, etafix, rstg) and csvertex-defaults (6029a2aa7, 711f251fb): already in fin51-r-fix 7f0e3dd9 as patch-equivalents (69e172aaa, af982e409, 1708c7a8a, a438d8a06, 504ff9823; eb45168a9, f453b83c5). Nothing extra merged.

## Caltech gates (H200, --fmad=false)
All runs: rst compared after <par_end>. user.hst differs in every pair (timing columns), and the reference shows the same.
| gate | result |
|---|---|
| hydro W121 p1, keys absent, 150 cycles, vs 7f0e3dd9 | rst + hydro.hst IDENTICAL (split_bval_fast default on) |
| MHD r1m, keys absent, 40 cycles, vs 7f0e3dd9 | rst + mhd.hst IDENTICAL |
| MHD r1m, 2 MPI ranks, keys absent, 100 cycles | rst + mhd.hst IDENTICAL |
| r1m launch keys, restart, 300 straight vs 151+149 | rst IDENTICAL |
| r1m launch deck, FRESH start, 300 straight vs 150+150 | rst IDENTICAL (DIFF before c5f1d72f2) |
| r1m launch keys, restart split in the VEF-pending window (8967) | rst IDENTICAL (DIFF before c5f1d72f2) |
| fresh launch deck + ck_impl_xpass = 1, 150+150 | rst IDENTICAL |
| accretor cs96 stream restart, 830-850, 2 GPU, vs cmtest 34c7b76e8 | rst SAME, clean exit |
| fresh cs96 stream, roche_wb_tan absent vs = true | rst SAME; log "true (cs default, fresh start)"; clean exit |

Note on rkg_strang_fuse: the fused RKG syncs only at outputs (any output incl. hst) and at the run end. A restart therefore reproduces the
straight run only when the straight run also has an output on the restart cycle, which is always true in production (rst outputs).

## r1m launch keys (WASP-121b MHD, fresh start)
r1m.athinput + `<mhd>` rkg_strang = true, rkg_strang_fuse = true; `<problem>` ck_time_scheme = hesdirk2, ck_t2_start = probe, ck_t2_seed = true.
Timing on 1 H200, r1m, 250 cycles, interleaved x2: 7f0e3dd9 50.6 ms/cycle; prod-1102 with the launch keys 47.5 ms/cycle (-6 %); + ck_impl_xpass = 1 49.8 ms/cycle.

## Please
Run the MI300A/HIP bitwise gate of prod-1102 (key absent vs fin51-r-fix c38bd1981: dhj hydro + MHD, accretor), and NOTE the result.
Do not merge prod-1102 into fin51-r-fix without the user.
