# NOTE 2026-10-21 raven -> caltech: V12_READY — accretor FRESH on 077dda23 + v12lzc_be_gm; RSG E stays on 4de8168d

fin51-r-fix = 077dda23 (pushed) = 4de8168d + pgen-restore on the BE redo + finalize frees + accretor booking A + opt-in k1_reset.
A100 gates BITWISE vs 4de8168d on the default path: RSG E clip N32, HG global clip N32, hg145, accretor static.

## Accretor (V12_READY; races Raven run12 31291344, released now on the same bin + deck)
- Build ry_per_accretor H200 nofma at 077dda23 (same recipe as your 4de8168d accretor build). Post the md5.
- Deck: docs/handover/acc3d-prod-1018/acc3d_prod_v12lzc_be_gm.athinput.in (md5 5de0635e; substitute the 3 @FILES@ lines as for
  v11amwtr; Raven source deck md5 c6af572f). = v12lzc + time_scheme = be + implicit_gas_src = matched (USER-approved after the
  v12lzc NO-GO; GO screen: v11 speed, 0 redos, budgets closed, fast-cell count = v11). Also has coriolis_lz = true, implicit_cs2_clip = rad.
- FRESH from t=0 in a new run dir (static phase to 0.05, then stream; tlim 1e9), your usual in-job C/C2/R check first.
  Race: the first actual t=0 start wins; the other keeps running until told (USER).

## RSG E 5126654 — CORRECTION (USER, supersedes the first version of this NOTE): STAYS on 4de8168d
Do NOT switch RSG E to 077dda23 and do not build it for RSG. Keep BIN = your 4de8168d build (9a56b4a1) and the clip_rst5d deck
(a7e734c4) unchanged. All global runs (RSG E, HG global) on every site stay on 4de8168d; only the accretor moves to 077dda23.
