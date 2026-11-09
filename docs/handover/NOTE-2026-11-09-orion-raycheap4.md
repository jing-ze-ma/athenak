# NOTE 2026-11-09 ORION -> viper: raycheap4 -- forcing-free time-order test P0-P3 (pulse2 transient decks)

Answer to `TASK-2026-11-09-orion-raycheap4.md` (corrected version, 42ad5569d). No code changes.

## Setup
- MOD = raycheap2-1109 `ff2ed77f1`, fresh build (gcc/13, openmpi/4.1, Release, MPI, `-ffp-contract=off`).
- Decks `tests_raycheap/pulse2/pulse2_<arm>_T<t_end>_c<cfl>.athinput` unchanged, 48 runs (4 arms x T 2.5/5/15 x cfl
  0.0375/0.01875/0.009375/0.0046875), 24 ranks each. All in one 6-node exclusive job (227961), 24 concurrent
  `srun -N 1 -n 24` steps, longest first. The whole job took 15 min; the critical path was P3_T15 finest (15 min, see below).
- All 48: rc=0, FATAL=0. conv.py on the last dump (all runs end at t_end exactly).

## Summary (conv.py, d = |E_k - E_k+1|_1/|E_finest|_1, 4 levels -> 3 d, 2 p)

| arm | T=2.5 d | p | T=5 d | p | T=15 d | p |
|---|---|---|---|---|---|---|
| P0 | [6.042e-03, 1.811e-03, 5.956e-04] | [1.74, 1.60] | [1.154e-03, 3.910e-04, 1.635e-04] | [1.56, 1.26] | [4.316e-07, 1.403e-12, 0.000e+00] | [18.23, inf] |
| P1 | [2.308e-02, 1.512e-02, 8.513e-03] | [0.61, 0.83] | [9.797e-03, 5.523e-03, 2.855e-03] | [0.83, 0.95] | [2.059e-03, 2.210e-04, 6.204e-06] | [3.22, 5.15] |
| P2 | [5.735e-03, 1.547e-03, 3.832e-04] | [1.89, 2.01] | [8.255e-04, 1.976e-04, 4.859e-05] | [2.06, 2.02] | [2.507e-06, 1.060e-07, 2.606e-08] | [4.56, 2.02] |
| P3 | [1.668e-02, 1.184e-02, 8.003e-03] | [0.49, 0.56] | [7.052e-02, 2.337e-02, 3.184e-03] | [1.59, 2.88] | [6.346e-02, 2.981e-02, 2.273e-02] | [1.09, 0.39] |

Reading (verdict at T = 2.5 and 5, as the TASK says):
- **P2 (vet_col frozen: time2_vet_col = lag, vet_col_every = 1e6): p = 1.89/2.01 (T2.5), 2.06/2.02 (T5). Clean 2nd order.**
  The transport stepping (hesdirk2 + the ramped inner flux) is 2nd order when the VET is not updated.
- **P0 (base: vet_col updated every step, time2_vet_col = extrap2): p = 1.74/1.60 (T2.5), 1.56/1.26 (T5), dropping
  with refinement.** Same counters as P2 (no BE steps besides the first, no fallbacks, no posfix). So the order loss
  comes in with the VET update/extrapolation itself, not with the stepping or the forcing.
- P1 (be): p = 0.61/0.83 (T2.5), 0.83/0.95 (T5) -> 1, the reference behaves.
- P3 (closure = eddington + implicit_hr_model, = raycheap3 C2e): p = 0.49/0.56 (T2.5), 1.59/2.88 (T5). Not a clean
  measurement: from cfl 0.009375 down (T2.5) and at every level (T5) hesdirk2 stage 2 is NOT ADMISSIBLE
  (min(E,T) = -1, i.e. negative E) a few times and the step is redone with backward Euler (stage fallbacks 1-4,
  backward-Euler steps 3-9, posfix 1-4 steps), so the scheme is partly BE and the count differs between levels.
  raycheap3 C2e (with forcing) gave p = 2.24/1.80; here the empty-field start with the steep front produces negatives.
- T = 15 (info): P0 d = 4.3e-7, 1.4e-12, 0 and P2 2.5e-6 ... 2.6e-8 -- by t = 15 the hesdirk2 runs are essentially at
  the relaxed steady state, the differences are at the solver tolerance (no time-order information). P1 (BE) still
  relaxing: d 2.1e-3 -> 6.2e-6. P3 T15 is spoiled by DT_RETRY dt cuts after t ~ 14.3 (reason NON-CONVERGED, resid
  1.8e-12 ... 2.8e-10 vs implicit_tol 1e-13) at all levels except c0.009375: 6, 4, 0, 36 retries; the finest ran 680
  cycles instead of 360 and took 902 s.

## Per-run counters (`time_scheme=hesdirk2:` line, rc, FATAL, NOT ADMISSIBLE count, wall)
`k1_redist ... time2_vet_col ... fallbacks` middle part of the counter line shortened to `...` (P0 time2_vet_col fallbacks = 0 in all 12 runs); P1's line is
`time_scheme=hesdirk2(be_only):`.
```
P0_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 19 s | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 29 s | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 54 s | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 102 s | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 38 s | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 55 s | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 106 s | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 198 s | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T15.0_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 77 s | hesdirk2: stage steps=44 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T15.0_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 145 s | hesdirk2: stage steps=89 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T15.0_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 276 s | hesdirk2: stage steps=179 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P0_T15.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 523 s | hesdirk2: stage steps=358 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 14 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=8 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 23 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=15 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 42 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=30 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 74 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=60 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 24 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=15 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 44 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=30 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 80 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=60 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 145 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=120 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T15.0_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 71 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=45 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T15.0_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 122 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=90 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T15.0_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 239 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=180 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P1_T15.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 407 s | hesdirk2(be_only): stage steps=0 backward-Euler steps=359 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 14 s | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 22 s | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 37 s | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 69 s | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 23 s | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 44 s | hesdirk2: stage steps=29 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 69 s | hesdirk2: stage steps=59 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 131 s | hesdirk2: stage steps=119 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T15.0_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 59 s | hesdirk2: stage steps=44 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T15.0_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 97 s | hesdirk2: stage steps=89 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T15.0_c0.009375 rc=0 FATAL=0 NOTADM=0 wall 189 s | hesdirk2: stage steps=179 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P2_T15.0_c0.0046875 rc=0 FATAL=0 NOTADM=0 wall 352 s | hesdirk2: stage steps=358 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P3_T2.5_c0.0375 rc=0 FATAL=0 NOTADM=0 wall 22 s | hesdirk2: stage steps=7 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=0 re-solves=0 cells=0 accepted-with-fallback stages=0
P3_T2.5_c0.01875 rc=0 FATAL=0 NOTADM=0 wall 31 s | hesdirk2: stage steps=14 backward-Euler steps=1 stage fallbacks=0 vet clips=0 ... time2_posfix: steps=1 re-solves=1 cells=28 accepted-with-fallback stages=0
P3_T2.5_c0.009375 rc=0 FATAL=0 NOTADM=1 wall 49 s | hesdirk2: stage steps=27 backward-Euler steps=3 stage fallbacks=1 vet clips=0 ... time2_posfix: steps=1 re-solves=1 cells=48 accepted-with-fallback stages=0
P3_T2.5_c0.0046875 rc=0 FATAL=0 NOTADM=2 wall 73 s | hesdirk2: stage steps=55 backward-Euler steps=5 stage fallbacks=2 vet clips=0 ... time2_posfix: steps=2 re-solves=2 cells=56 accepted-with-fallback stages=0
P3_T5.0_c0.0375 rc=0 FATAL=0 NOTADM=1 wall 49 s | hesdirk2: stage steps=12 backward-Euler steps=3 stage fallbacks=1 vet clips=0 ... time2_posfix: steps=2 re-solves=2 cells=68 accepted-with-fallback stages=0
P3_T5.0_c0.01875 rc=0 FATAL=0 NOTADM=1 wall 62 s | hesdirk2: stage steps=27 backward-Euler steps=3 stage fallbacks=1 vet clips=0 ... time2_posfix: steps=3 re-solves=3 cells=104 accepted-with-fallback stages=0
P3_T5.0_c0.009375 rc=0 FATAL=0 NOTADM=1 wall 93 s | hesdirk2: stage steps=57 backward-Euler steps=3 stage fallbacks=1 vet clips=0 ... time2_posfix: steps=2 re-solves=2 cells=76 accepted-with-fallback stages=0
P3_T5.0_c0.0046875 rc=0 FATAL=0 NOTADM=4 wall 150 s | hesdirk2: stage steps=111 backward-Euler steps=9 stage fallbacks=4 vet clips=0 ... time2_posfix: steps=4 re-solves=4 cells=112 accepted-with-fallback stages=0
P3_T15.0_c0.0375 rc=0 FATAL=0 NOTADM=8 wall 272 s | hesdirk2: stage steps=31 backward-Euler steps=37 stage fallbacks=6 vet clips=0 ... time2_posfix: steps=7 re-solves=13 cells=38692 accepted-with-fallback stages=0
P3_T15.0_c0.01875 rc=0 FATAL=0 NOTADM=6 wall 236 s | hesdirk2: stage steps=67 backward-Euler steps=29 stage fallbacks=6 vet clips=0 ... time2_posfix: steps=8 re-solves=17 cells=404604 accepted-with-fallback stages=0
P3_T15.0_c0.009375 rc=0 FATAL=0 NOTADM=8 wall 238 s | hesdirk2: stage steps=134 backward-Euler steps=46 stage fallbacks=8 vet clips=0 ... time2_posfix: steps=9 re-solves=24 cells=552300 accepted-with-fallback stages=0
P3_T15.0_c0.0046875 rc=0 FATAL=0 NOTADM=34 wall 902 s | hesdirk2: stage steps=263 backward-Euler steps=418 stage fallbacks=33 vet clips=0 ... time2_posfix: steps=33 re-solves=114 cells=3838948 accepted-with-fallback stages=0
```
DT_RETRY only in P3_T15 (c0.0375: 6, c0.01875: 4, c0.0046875: 36; all NON-CONVERGED, all after t ~ 14.3).

## Files
/orion/ptmp/jinma/raycheap4_1109: rc4.sbatch, rc4.out (job log), ana.sh, rc4_ana.out (full analysis output),
runs/*/run.log + wall.txt. bin dumps, sources and the build were deleted after conv.py (ptmp quota).

