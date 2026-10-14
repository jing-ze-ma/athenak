# NOTE Caltech -> Raven: 4-GPU trials of the Linhao BSG and He giant FIN4 (USER 10-14), restarted from the 2-GPU prod rst

USER 10-14: try fewer, longer links with more GPUs, restarting from the newest rst of the running 2-GPU chains, and check whether a
GPU-layout change at restart gives any issue. The 2-GPU production chains are unchanged and keep running (bsgLIc 4487415 + 4487610-13,
hgF4 4489484 + 4489485-87); they remain the reference runs.

## Wide jobs (single long link each, 1 node x 4 H200, 4 ranks)
| run | job | MeshBlocks/GPU | wall | dir |
|---|---|---|---|---|
| Linhao BSG | 4490079 (-J bsgLIw4) | 4 | 30 h | bsgli_1014/wide4/run |
| He giant FIN4 | 4489939 (-J hgF4w) | 1 | 2 d | hegiant_1014/wide4/run |

- Same binaries / decks as prod (8e018272 md5 4a92728d; d505d21f md5 1440c053). No deck edits.
- At job start, run.cfg hard-links the newest VALID prod rst (rst_pick.py / rst_info.py; waits if < 120 s old) into the wide rst/;
  the normal hook restarts from it with its last_time keys. Link failure -> STOP (never a fresh start). Prod dirs are only read.
- A 2-node x 4 H200 BSG trial (4489938) was submitted and withdrawn before it started (user: not needed; 1 node fits in one link).
- Both wide jobs wait (afterany) for the layout test below; held again if it shows a 4-rank crash / NaN.
- Restarts across 2 -> 4 ranks are NOT bitwise by construction (reduction order); the wide runs are a separate branch of each run,
  for speed (s/cycle 2 vs 4 GPUs) and robustness, not a replacement of the bitwise 2-GPU chain.

## Layout test (laytest 4490288, 1 node x 4 H200, ~40 min, est. start 10-14 09:02)
From the same prod rst: 2 ranks vs 4 ranks, BSG 30 cycles (bsg3d.00002.rst, t 1.728e5, cycle 2051), He giant 50 cycles (newest valid
fin4 rst at job start). Reports rc, FATAL/NaN/DIVERGED/NC/NOT ADMISSIBLE, Picard, s/cycle, and 2-vs-4 diffs (bitwise or max rel
per variable). Results will follow in this NOTE.

## Other
- He giant FIN4 link 1 4489484 running ~1.1 s/cycle (FIN3 0.93); Linhao BSG link 1 at t 2.6e5 s, ~1.28 s/cycle, all clean.
- fin4-mem CPU bitwise TASK: running (debug QoS); its own NOTE follows.

## Layout test result (job 4490288, 6 min, 2 GPUs vs 4 GPUs from the same rst)

No crash, no NaN, rc 0 for all 4 runs. As expected, the 2-GPU and 4-GPU runs are NOT bitwise equal (different rank layout gives different sum orders).

- BSG, 30 cycles from bsg3d.00002.rst: dens rel 4.6e-10; eint 1.1e-7; velx 1.2e-7; vely/velz maxabs/max 2e-8; m1_e 1e-7; m1 f2/f3 maxabs/max ~1e-9. hydro.hst last row 7e-12; user.hst 2.5e-3 at col 17.
- He giant, 50 cycles from fin4 rst: dens rel 1.1e-9; eint 9.9e-8; m1_e 1.2e-7; vely/velz maxabs/max ~1e-6. **velx maxabs/max 4.8e-2** (one or a few cells; not localized yet, may be a near-floor or fast cell where round-off grows quickly). hydro.hst last row 2.4e-10; user.hst 2.1e-4 at col 19.

Verdict: a layout change only causes round-off-level divergence (apart from the one He giant velx cell to check), so the wide jobs continue. hgF4w 4489939 has started and seeded from fin4/rst/hegiant.00004.rst; bsgLIw4 4490079 is pending (Priority). The laytest rst hard links are removed.

Addendum (details from the test log):
- The He giant layout test started from hegiant.00000.rst (t=0, cycle 0), not from a fin4 rst.
- The He giant velx difference sits in the outer floor cells only (i >= 888 of 897, density 1e-20). That is floor noise, not a concern.
- FATAL, NaN, DIVERGED, NON-CONVERGED and NOT ADMISSIBLE counts are 0 in all four runs.
- Picard passes match between layouts: 6 for BSG, 11 for the He giant.
- s/cycle on 2 vs 4 GPUs (short runs, includes init): BSG 2.14 vs 1.53, He giant 1.39 vs 0.96.

## 10-14 morning: the 4-GPU runs are now the production runs; the 2-GPU chains are stopped (user)

Both 4-GPU runs overtook their 2-GPU parents overnight, so the user had both 2-GPU chains stopped. That covers the running link 1 plus the queued links: BSG 4487415 and 4487610-13, He giant 4489484 and 4489485-87. Their rst dirs are kept as a fallback (BSG prod/rst up to 00018, He giant fin4/rst up to 00012).

| run | job | t at check / tlim | s/cycle | dt | est. left |
|---|---|---|---|---|---|
| BSG 2 GPU (stopped) | 4487415 | 1.59e6 / 4.96e6 | 1.55 | 83.3 | - |
| BSG 4 GPU | 4490079 (30 h limit) | 1.70e6 / 4.96e6 | 1.09 | 83.3 | ~12 h |
| He giant 2 GPU (stopped) | 4489484 | 5.5e5 / 2.592e6 | 0.94 | 18.2 | - |
| He giant 4 GPU | 4489939 (2 d limit) | 8.4e5 / 2.592e6 | 0.48 | 19.0 | ~12 h |

On 4 GPUs: BSG is 1.42x faster per cycle and the He giant 1.95x. There is no FATAL or NaN in any log. Both 4-GPU runs should reach tlim inside a single job, and a report will follow when each ends.

Run dirs:
- BSG: /resnick/groups/carnegie_poc/jingze/bsgli_1014/wide4/run
- He giant: /resnick/groups/carnegie_poc/jingze/hegiant_1014/wide4/run
