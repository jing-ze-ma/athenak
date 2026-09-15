---
name: session-state-2026-08-25
description: START HERE - job 190939 LEFT RUNNING toward tlim; the viper handoff doc's blow-up does NOT reproduce on orion CPU, and the exact input is pushed as cd2e8815
metadata:
  type: project
---

Entry point. Branch `polar-average-perf` at `cd2e8815`.
Predecessor: [[session-state-2026-08-24]].

## FIRST THING: check job 190939

`/orion/ptmp/jinma/Athenak/blowup/`, `p.shared`, 16 ranks x 7 threads, 12 h.
The **second link of a restart chain** -- 190900 walled out clean at
t = 6.629e5 = 2.174 rot after 11:59:22, and 190939 picked up from `rst/dhj.00066.rst`
automatically (`--dependency=afterany`). At the last look it was at 2.318 rot, no NaN.

```
squeue -j 190939; cd /orion/ptmp/jinma/Athenak/blowup
grep -icE 'nan|fatal' log.out.190939
grep -oE 'time=[0-9.e+-]+ dt=[0-9.e+-]+' log.out.190939 | tail -1
```
Rotation = 3.05e5 s. `tlim = 1.2e6` = 3.934 rot, reachable inside this slot. Markers:
general+grey dies at **2.509 rot** per the doc; reaching tlim is the doc's own bar for
"clean". **The user asked to be told when it dies.**

## The result of 2026-08-24, and it is a disagreement

**Read `docs/HANDOFF_dhj_ck_eos.md` in the repo first** -- the viper session's notes,
nine refuted hypotheses. Then: its fast reproducer (general EOS + correlated-k, limiter
OFF), which it says **dies at 0.674 rot**, ran clean past 2.17 rot on orion with dt not
merely steady but RISING. Details, and why the input is faithful, in
[[dhj-general-eos-ck-blowup]].

The exact file is committed and pushed as `inputs/mhd/deep_hot_jupiter_rt_blowup.athinput`
(`cd2e8815`, md5 `f52295ff8c321338369a91b76a5a1452`) so the viper session can run the
identical thing. **Only `problem/ck_table` and `problem/ck_data_dir` are orion-absolute
and need repointing.** The user was going to ask the viper session to re-run it.

The only known remaining difference is **CPU (Freya/orion, gcc 13 + OpenMPI 4.1) vs the
APU/GPU viper used**. Viper's HEAD was `c1d513c2`, the same commit, so it is not a
version difference.

## Two bugs found, NEITHER FIXED -- both need the user's go-ahead

1. [[rt-srclim-warn-rank-local]] -- the RT clip count is rank-local and the warning is
   SILENT when rank 0 is not the rank clipping. Distorts every clip count in the campaign
   record.
2. [[hst-cartesian-volume-on-spherical]] -- `.hst` weights cells by `dx1*dx2*dx3`, so its
   `mass`/`tot-E` are not physical integrals on this spherical mesh (8x off, measured).
   The handoff doc's "mass to 0.15 %" and "loses 67 % of its mass" may both be affected.

Fixing (1) is a ~5-line MPI_Allreduce plus a test; (2) is a `history.cpp` change or a
`user_hist_func`. Offered, not authorised.

## Conservation, measured

Mass +0.21 %, internal energy +0.56 % over t = 4.32e5 -- fine, not blow-up-like. Recipe
for measuring it (the hst cadence AND weighting both make hst useless here):
[[dhj-conservation-check]].

## Uncommitted in the repo (unchanged from 2026-08-24)

`data/exo_fms_ck/PROVENANCE.md`, `inputs/mhd/deep_hot_jupiter_rt_eos.athinput` (adds the
`<output3>` rst block, rewrites the stale `<meshblock>` comment),
`inputs/mhd/deep_hot_jupiter_rt_ideal_xe.athinput` (adds the whole `rt_ck`/`ck_*` block).
Rationale in [[session-state-2026-08-24]].

## Standing instruction from the handoff doc

The Gamma_1 lead touches `src/mhd/mhd.hpp:588-596`, and **the user said "let's not go
that way" on 2026-08-23**. Do not open it unless they raise it.
