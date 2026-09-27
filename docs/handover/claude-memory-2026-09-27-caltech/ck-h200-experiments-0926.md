---
name: ck-h200-experiments-0926
description: H200 ck speed-up experiments (2026-09-26 night): no lever reaches 25 % of ck; ck_beam_par -6.2 % total (11 % of ck, within noise) is the only win; host gaps only 1.5 %; ck_dif_dtau +106 %; RT_NB 2 slower; jlin rewrite is the main proposal
metadata:
  type: project
---

Agent report 2026-09-26, about 01:00 PDT. Branch `ck-h200` (worktree `wt_ckh200`; commits 2ed898a1 RT_NB2, 11d9a78c launch bounds). Not merged: experiments only. Work dir `/resnick/scratch/jingze/ckh200/`.

**Baseline.** 2 H200, nx1 256, production keys: 36.40 ms/cycle (1 GPU: 62.98). A cad-only arm (hydro + cad_lin) gives 15.18. So ck costs about 21.2 ms/cycle (58 %), and the 25 % target is 5.3 ms/cycle.

**Exp 1: host gaps are NOT a lever.**
- Steady-state idle is 1.5 % of the cycle (0.3 % in ck phases).
- The "40 % idle" in the first profile was startup, the first full call, and the output/finalize after the last cycle.
- The cudaEventSynchronize calls are Kokkos constant-memory launches, with 1-6 us gaps.
- Upper bound for a device-side convergence check or CUDA graphs: about 0.1 ms/cycle.
- ck cost split (1 H200):

| kernel | share of ck |
|---|---|
| store TsrtCkChain | 35 % |
| rt_chain_ck_jlin | 23 % |
| lin1p | 16 % |
| rt_apply | 7 % |
| lin_build | 7 % |
| lin_sum | 4 % |
| CkImplStep | 3 % |

**Exp 2: runtime keys** (2 H200; accuracy = deep eint rms 1-100 bar relative to a 1e-14 kick twin).

| arm | ms/cycle | change | accuracy |
|---|---|---|---|
| ck_beam_par | 34.15 | -6.2 % | 1.03x noise (nx1 76, 1 rotation); 1.15x (nx1 256, 0.25 rotation) |
| ck_dif_dtau=10 | 74.90 | +106 % | store kernel 3.3x slower; unexplained, needs profiling |
| both | 49.04 | +35 % | |

- ck_beam_par is not bitwise on GPU.
- The noise keys are `problem/seed=N problem/seed_amp=1e-14`. They must be in the input file.

**Exp 3: RT_NB and occupancy.**
- RT_NB=1 is invalid with nquad 2: the beam pairing needs both angles in one block.
- RT_NB=2 is +2.4 % slower.
- The store kernel has 164 registers, 18.75 % occupancy and 5.33 waves: it is latency-bound, not single-wave.
- Launch bounds: -0.3 % total (store -5 %, rt_apply +9 % from spills).

**Proposals** (together about 25 % only with the jlin rewrite):
1. Adopt ck_beam_par (-11 % of ck). User decision.
2. Rewrite or fuse rt_chain_ck_jlin (22 % of ck; 70 GB DRAM per call, L1 hit 8 %): maybe ~10 % of ck.
3. Launch bounds on the store kernel only (1-2 %).
4. rt_apply (250 registers, grid capped at 256 by the Kokkos reduce): about 3 %.

Excluded as not accuracy-neutral: FP32, fewer passes, every 8.
