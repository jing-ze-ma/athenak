---
name: ck-next-resume-0927
description: RESUME HERE (2026-09-27) -- branch ck-next (HEAD 759e92f5, NOT merged): remaining-ck-kernel speed-ups; 1 H200 same-job timing -11 % total / ~22.5 % of remaining ck (cumulative ~52 % of ck since pre-jlin); still to do -- confirm GPU accuracy class, read the 2-GPU timing job 3562060, then ask the user to merge
metadata:
  type: project
---

State when the user went to bed on 2026-09-26 at ~23:40 PDT. The agents were stopped; the work is on disk.

**Branch `ck-next`**
- Worktree `/resnick/home/jingze/ATHENAK/wt_cknext`, HEAD **759e92f5** (12 commits on a81d3677). Not merged or pushed.
- Work dir `/resnick/scratch/jingze/cknext`: scripts gg.sub (GPU gate), tjob.sub + arms.sh (timing), ana2.py, acc.py (noise gate), prof.sub + kn.py (nsys with an NVTX connector), ncs.py.
- Binaries in `builds/`:
  - cnb = base a81d3677;
  - cn9 = the last step that is GPU-bitwise vs base;
  - cn13 = 161aeab9;
  - cn14 = 759e92f5, final.

**Steps**
- ck_lin_geom (lG with beta and dz).
- ck_beam_tau with 8 chains per thread.
- j-fastest copy of ck_coef's inputs.
- lin_sum with one thread per (block, face).
- ck_coef with one thread per (pair, cell).
- Prefetch in lin_build, lin1p and jlin.
- cin/cout/idn formed instead of read (CkCinCout, CkIdn).
- Pair pre-add in lin1p.
- lP slots 1 and 5 dropped.
- CkMulRn rounding of the pre-added products.
- 759e92f5: ck_coef and ck_lin_build use the SAME CkCinCout/CkIdn, the fix attempt for the GPU non-bitwise result from cn10 on.

**Accuracy**
- CPU gate: cnb vs cn14 BITWISE.
- GPU 1 vs 2 on H100: rst same, hst differs (the usual reduction effect).
- **Unknown: is the GPU bitwise vs base at cn14?** Check with gg.sub (look at the latest gg.out.*). If it is not bitwise, run the noise gate with >= 2 kicks (acc.py; the kick arms are in arms.sh); pass is <= ~1.2x.

**Timing, 1 H200** (job 3562061, hpc-sm-01-05, 3 reps, median of 8-cycle windows; `python3 ana2.py t1f 40 400`):

| arm | ms/cycle |
|---|---|
| k46 (pre-jlin) | 59.31 |
| base | 45.51 |
| **cn14** | **40.52 (-11.0 %)** |
| cad | 23.36 |

- Remaining ck = base - cad = 22.15 ms; cn14 saves 5.0 ms, i.e. 22.5 %.
- Cumulative since pre-jlin: 18.8 of 36.0 ms, about 52 % of ck.

**Pending:** the 2-H200 timing, job **3562060** (tjob.sub in /resnick/scratch/jingze/cknext, output log.out.3562060, est. start 09-27 05:42). Read it with ana2.py on its tag dir (check the log's TAG line).

**Next session:** confirm the GPU accuracy class of cn14 (and run the noise gate if needed), read the 2-GPU timing, then ask the user whether to merge ck-next. On merge: a handover note asking viper for its HIP gate; delete the branches ck-next-vA and ck-next-vB; remove the worktree.

**Also since last summary:** viper pushed flux_hst (default on; dhj.user.hst gains flux columns; hydro hst/bin/rst bitwise), the sph_wedge reservoir bottom BC, and a session checkpoint. The main checkout is at 392f1e50.

**RESOLVED 2026-09-27:** gates passed (noise 0.78x), 2-GPU timing -9.4 %, merged 6a9ffb0d and pushed (50325e8b). Nothing pending.
