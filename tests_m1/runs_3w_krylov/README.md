# Implicit-M1 Krylov on several GPUs: pipelined BiCGStab and a direct multi-rank halo

Date 2026-09-23, viper. Branch `m1-krylov` (worktree `/viper/ptmp2/jinma/wt_krylov`), base
89dc28d7 (m1-space2, fast path default on). Commit 988f11ae. Build and run tree:
`/viper/ptmp2/jinma/krylov_0923/` (`ref/` = git archive 89dc28d7, `new/` = this branch;
`bin/`, `cpu/`, `runs/`). Scripts are in `scripts/` here.

Motivation: `tests_m1/runs_3q_scscale/README.md` sect. 6 (sc-scale, bbf94802). The radiation
part of the cycle (27-38 ms) does not shrink from 1 to 8 GPUs; README_FAST sect. 4-5 blames
~588 launches per cycle and 2 blocking reductions per BiCGStab iteration, and
`implicit_halo_direct` only works when every neighbour is on the rank.

## 1. What changed

All new code is in the new file `src/rad_m1/rad_m1_krylov.cpp`. The only edits to
`rad_m1_implicit.cpp` are the two keys, one dispatch line in `ImplicitBiCGStabFused` and one
branch in `ImplicitHaloExchange`, so the branch should merge cleanly with m1-vimp (operator
work).

Two `<rad_m1>` keys, both default **false**. With both off nothing new runs and nothing new
is allocated.

**`implicit_krylov_pipe`** (needs `implicit_krylov_fuse = 3`): the pipelined
right-preconditioned BiCGStab of Cools & Vanroose (Parallel Computing 65, 2017, Alg. 4,
"p-BiCGStab").
- Auxiliary vectors w = A M^-1 r, s = A M^-1 p, z = A M^-1 s, t = A M^-1 w, v = A M^-1 z
  are carried by recurrences. The preconditioned copies are in iw's Krylov slots; phat,
  s, shat and rhat are in a new `kpw` array (4 vectors).
- Per iteration:
  - K1 is one fused vector update plus the reduction R1, (q,y) and (y,y).
  - P1 is M^-1 z, the halo and the operator (a plain par_for, without the reduction).
  - K2 is one fused update of x, r, rhat, w plus R2: (rt,r), (rt,w), (rt,s), (rt,z) and
    max|r|.
  - P2 is M^-1 w, the halo and the operator.
- **Neither reduction drains the device queue.**
  - Each is reduced into pinned host memory (`SharedHostPinnedSpace`).
  - The host waits on a HIP event recorded right after that kernel. By then the
    preconditioner (and, on one rank, the halo and operator) are already queued.
  - On several ranks the partial sums go into one `MPI_Iallreduce` of 6 doubles. It is
    posted before the halo exchange and completed after the operator launch.
  - Still 2 reductions per iteration, but both hidden. A single-reduction variant needs
    the (y,y) = (w,w) - 2a(w,z) + a^2(z,z) expansion, which cancels catastrophically near
    convergence, so it was not used.
- **Convergence** is tested on the recursive max|r| from R2, one P2 late. It is always
  confirmed on the true residual. A failed confirmation, or any breakdown, restarts from
  the true residual of x. Restart and fallback rules are those of the fused loop (3
  restarts, then line Jacobi).

**`implicit_halo_mpi`** (needs `implicit_halo_direct = true`): the implicit exchanges when
some neighbours are on other ranks.
- On-rank neighbours: the existing copy kernel. `hd_src` already holds -1 for off-rank
  neighbours.
- Off-rank neighbours: one precomputed cell list, one pack kernel, **one message pair per
  neighbour rank** (a private communicator) and one unpack kernel. The on-rank copy kernel
  is queued behind the pack and runs while the messages fly.
- Both sides enumerate regions sorted by (receiver gid, receiver direction), and each region
  in (k,j,i) order, so a message needs no header. This is the aggregation plan of
  vet_mb_agg (sc-scale).
- Pure copies of the same numbers, so the result is **bitwise** the general-machinery path.
- Same-level meshes without seams or poles only. Anything else prints a line and keeps the
  ordinary exchange.

## 2. CPU gates (`scripts/cpu_gate.sh`, `cpu_gate.log`)

Setup:
- gcc/openmpi build of ref (89dc28d7, md5 11fcefbf) and new (988f11ae source, md5 e0a8f010).
- Inputs: the runs_3p_fastdefault ones, `slab2d_def` (2-D seeded slab, 200 s) and
  `box3d_def` (84x32x32 in 4 blocks, 60 cycles). The two new keys are spelled out as false.
- 2 ranks: meshblock/nx2 = 16; 4 ranks: nx2 = 8 (2-D) or 1 block per rank (3-D).
- NON-CONVERGED = 0 and line-Jacobi fallbacks = 0 in all 36 runs.

**Switches off, bitwise vs 89dc28d7** (hydro and user hst identical):
- 2-D Eddington, 1, 2 and 4 ranks (`e_ref*`/`e_off*`);
- 2-D vet_sc full, 1 and 2 ranks (`v_ref*`/`v_off*`);
- 3-D Eddington, 1, 2 and 4 ranks (`e3_ref*`/`e3_off*`);
- 3-D vet_sc, 4 ranks (`v3_ref4`/`v3_off4`).

**`implicit_halo_mpi`, bitwise vs the current multi-rank path** (hst identical):
- 2-D Eddington, 2 and 4 ranks (`e_hm2`, `e_hm4`);
- 2-D vet_sc, 2 ranks (`v_hm2`);
- 3-D Eddington, 2 and 4 ranks (`e3_hm2`, `e3_hm4`);
- 3-D vet_sc, 4 ranks (`v3_hm4`).

**`implicit_krylov_pipe`, 2-D slab, 200 s.** Relative differences against `*_off` on 1 rank
(`scripts/stats.py`). `ctl` = off with rad_flux_inner changed by +1 ulp. `ph` = pipe +
halo_mpi.

| arm | F1top/Fin mean | KE1 end | KE2 end | dt mean | totE end |
|---|---|---|---|---|---|
| Eddington ctl | 5.2e-14 | 2.7e-9 | 2.3e-10 | 1.2e-12 | 2.7e-12 |
| Eddington pipe, 1 rank | 1.8e-13 | 2.2e-9 | 8.4e-9 | 2.3e-12 | 1.1e-12 |
| Eddington ph, 2 ranks | 2.3e-13 | 1.3e-9 | 5.8e-9 | 2.4e-13 | 3.5e-12 |
| Eddington ph, 4 ranks | 3.3e-13 | 6.5e-10 | 4.9e-9 | 3.0e-12 | 1.2e-12 |
| Eddington off, 4 ranks | 4.1e-13 | 5.5e-9 | 1.3e-8 | 2.9e-12 | 2.5e-12 |
| vet_sc pipe, 1 rank | 2.7e-13 | 9.0e-10 | 6.5e-8 | 3.9e-11 | 1.2e-12 |
| vet_sc ph, 2 ranks | 5.4e-13 | 2.1e-9 | 4.4e-9 | 5.8e-11 | 1.2e-13 |
| vet_sc off, 2 ranks | 2.3e-12 | 6.1e-9 | 1.2e-8 | 6.2e-11 | 9.9e-13 |

**3-D box, 60 cycles** (`scripts/cmp.py`, last row vs `e3_off`):
- pipe: totE 6.7e-13, KE1 1.1e-9, F1top 1.1e-12.
- ph on 4 ranks: 9.4e-13, 2.3e-10, 1.9e-12.
- ctl: 1.8e-12, 2.2e-10, 2.9e-12.
- vet_sc ph on 4 ranks vs off on 4 ranks: 1.9e-13, 1.3e-9, 8.1e-12.
- The net transverse momenta (hst cols 5-6) differ by up to 4.8e-3 of their own maximum.
  They are round-off noise: max|M2| = 4.5e9 against sqrt(2 KE2 M) = 4.1e21, a ratio of
  1e-12.

**Round-off and stability of the pipelined loop.** Inner iterations total, breakdowns
(restarts):

| case | fused (off) | pipe |
|---|---|---|
| 2-D Eddington, 1 rank | 25 611, 3 | 25 632, 1 |
| 2-D Eddington, 2 ranks (ph) | 25 956, 3 | 25 943, 3 |
| 2-D Eddington, 4 ranks (ph) | 27 132, 1 | 27 151, 3 |
| 2-D vet_sc, 1 rank | 25 621, 1 | 25 608, 4 |
| 2-D vet_sc, 2 ranks (ph) | 24 977, 6 | 24 963, 3 |
| 3-D Eddington, 1 rank | 2 018, 0 | 2 017, 0 |
| 3-D Eddington, 4 ranks (ph) | 2 018, 0 | 2 005, 0 |
| 3-D vet_sc, 4 ranks (ph) | 2 029, 1 | 2 038, 0 |

- The iteration counts agree to 0.1 %.
- Restarts: 0-4 per 2800-3400 solves in either loop. No fallbacks, NON-CONVERGED 0.
- The recursive-residual drift does not show at the working tolerance (Eisenstat-Walker,
  eta <= 1e-2). A restart is at most one extra start-up per solve.

**Restart bitwise** (pipe + halo_mpi, 2 ranks, 2-D slab):
- Method: run to t = 40 with a restart file at t = 20, then restart (`rs_full`/`rs_rst`).
- Identical: the final rst file, both final bin files, and all 20 hst rows after the restart
  (`scripts/hstrst.py`).

**Verdict: PASS.**

## 3. GPU cost

PENDING. Jobs submitted 2026-09-23 ~12:00, all with one binary `bin/athena_new_gpu` (md5
0c740d15) plus `head_*` arms on `bin/athena_ref_gpu` (89dc28d7, md5 68d89fc9):
- 11948842 `runs/f12`: 1 and 2 GPUs, apudev.
- 11948843 `runs/fw2`: weak, 2 GPUs, apudev.
- 11948845 `runs/f4`: strong and weak, 2 nodes, apu.
- 11948846 `runs/f8`: strong and weak, 4 nodes, apu.
- 11948847 `runs/prof`: rocprofv3 kernel + HIP runtime trace; launches and syncs per cycle
  from nlim 20 - nlim 10, 1 and 2 GPUs.

Setup:
- The box and layout of runs_3q_scscale (84x104x104, 8 blocks of 84x52x26; weak 0.9 M
  cells/GPU on 2/4/8 GPUs).
- vet_sc full, 64 rays, the FAST switches spelled out.
- 100 cycles; ms/cycle over cycles 10-100.
- Repeats a (in order) and b (reversed).
- `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`.
- Arms: off, hm (halo_mpi), pipe, ph (both), hyd (same box without `<rad_m1>`).

Analysis: `bash scripts/finish.sh` gives the timing tables, the profile counts and the GPU
exactness checks.
