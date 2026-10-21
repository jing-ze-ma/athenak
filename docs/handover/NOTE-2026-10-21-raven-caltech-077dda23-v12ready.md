# NOTE 2026-10-21 raven -> caltech: V12_READY — accretor FRESH on 077dda23 + v12lzc_be_gm; RSG E 5126654 -> 077dda23 at start

fin51-r-fix = 077dda23 (pushed) = 4de8168d + pgen-restore on the BE redo + finalize frees + accretor booking A + opt-in k1_reset.
A100 gates BITWISE vs 4de8168d on the default path: RSG E clip N32, HG global clip N32, hg145, accretor static.

## Accretor (V12_READY; races Raven run12 31291344, released now on the same bin + deck)
- Build ry_per_accretor H200 nofma at 077dda23 (same recipe as your 4de8168d accretor build). Post the md5.
- Deck: docs/handover/acc3d-prod-1018/acc3d_prod_v12lzc_be_gm.athinput.in (md5 5de0635e; substitute the 3 @FILES@ lines as for
  v11amwtr; Raven source deck md5 c6af572f). = v12lzc + time_scheme = be + implicit_gas_src = matched (USER-approved after the
  v12lzc NO-GO; GO screen: v11 speed, 0 redos, budgets closed, fast-cell count = v11). Also has coriolis_lz = true, implicit_cs2_clip = rad.
- FRESH from t=0 in a new run dir (static phase to 0.05, then stream; tlim 1e9), your usual in-job C/C2/R check first.
  Race: the first actual t=0 start wins; the other keeps running until told (USER).

## RSG E 5126654 (not started)
- Build he_star_m1 H200 nofma at 077dda23 and switch run.cfg BIN/BIN_MD5 before it starts (default path bitwise vs 4de8168d).
  Deck stays your clip_rst5d (a7e734c4). If the build is not ready when 5126654 starts, 4de8168d is fine; switch at the next link.
