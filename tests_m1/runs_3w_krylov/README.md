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

**Jobs.** All with `HSA_XNACK=1` and `HSA_NO_SCRATCH_RECLAIM=1`. Every run exited with
rc = 0. Summary: `bash scripts/finish.sh`.
- 11948842 (`runs/f12`, 1 and 2 GPUs, apudev);
- 11948843 (`runs/fw2`, weak 2 GPUs, apudev);
- 11948845 (`runs/f4`, 2 nodes, apu);
- 11948846 (`runs/f8`, 4 nodes, apu);
- 11948847 (`runs/prof`, rocprofv3).

**Binaries.** One binary, `bin/athena_new_gpu` (md5 0c740d15), for every arm. The `head_*`
arms use `bin/athena_ref_gpu` (89dc28d7, md5 68d89fc9).

**Setup:**
- The box and layout of runs_3q_scscale: 84x104x104 in 8 blocks of 84x52x26.
- Weak scaling at 0.9 M cells per GPU: 84x208x104 on 2 GPUs, 84x208x208 on 4, 84x416x208
  on 8.
- vet_sc full, 64 rays, the FAST switches spelled out.
- 100 cycles; ms/cycle is measured over cycles 10-100.
- Repeat a runs in order, repeat b in reverse order.

**Arms:**
- off = both keys false (the current path);
- hm = `implicit_halo_mpi`;
- pipe = `implicit_krylov_pipe`;
- ph = both;
- hyd = the same box without `<rad_m1>`.

**Checks:**
- **GPU exactness.** hst files are identical for head vs off (1 and 2 GPUs) and for hm vs
  off at 2, 4 and 8 GPUs, strong and weak.
- **pipe / ph vs off.** Max relative hst difference: user 2.0e-9 to 5.4e-9, hydro 3.3e-10
  to 1.1e-9. For scale, off on 1 GPU vs off on 2 GPUs differs by user 4.8e-9, hydro
  8.5e-10.
- **Solver health.** NON-CONVERGED 0 and fallbacks 0 in every arm. Breakdowns (restarts)
  0-2 per run.
- **Iterations.**
  - Picard passes per step are identical across arms.
  - Inner iterations per solve agree to 0.5 %: 15.26-15.30 (strong), 15.40-17.01 (weak).

**Strong scaling** (908 k cells). Values are ms/cycle, repeat a then b. rad = cycle minus
hydro, from the a/b means.

| GPUs | off | hm | pipe | **ph** | hydro | rad off | rad hm | **rad ph** |
|---|---|---|---|---|---|---|---|---|
| 1 | 58.9, 55.0 | - | 55.2, 55.5 | - | 23.1, 22.4 | 34.2 | - | **32.6** (pipe) |
| 2 | 55.6, 54.1 | 44.4, 45.2 | 55.2, 53.3 | **44.2, 42.8** | 16.5, 16.5 | 38.4 | 28.3 | **27.0** |
| 4 | 45.9, 44.9 | 39.0, 37.7 | 44.4, 43.6 | **36.3, 35.8** | 14.2, 14.1 | 31.3 | 24.2 | **21.9** |
| 8 | 39.6, 39.0 | 34.6, 37.7 | 38.7, 37.6 | **32.6, 35.9** | 12.4, 12.5 | 26.9 | 23.7 | **21.8** |

The head arms ran at 57.1, 55.7 on 1 GPU and 57.3, 54.5 on 2 GPUs, the same as off.

**Weak scaling** (0.9 M cells per GPU):

| GPUs | off | hm | pipe | **ph** | hydro | rad off | rad hm | **rad ph** | ph / hydro |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 58.9, 55.0 | - | 55.2, 55.5 | - | 23.1, 22.4 | 34.2 | - | 32.6 | 2.43 (pipe) |
| 2 | 84.5, 81.1 | 68.4, 64.7 | 84.8, 82.4 | **68.4, 66.4** | 26.3, 23.4 | 58.0 | 41.7 | **42.6** | 2.71 |
| 4 | 82.8, 85.1 | 69.3, 70.1 | 82.5, 82.8 | **68.8, 67.4** | 24.5, 23.5 | 60.0 | 45.7 | **44.1** | 2.84 |
| 8 | 90.8, 98.6 | 77.9, 77.1 | 92.2, 91.9 | **74.6, 75.6** | 25.0, 25.2 | 69.6 | 52.4 | **50.0** | 2.99 |

**Launches and host syncs per cycle** (`runs/prof`). rocprofv3 kernel + HIP runtime trace,
rank 0, counted as (nlim 20 - nlim 10)/10.

| arm | launches | hipStreamSync | hipDeviceSync | hipEventSync | queue drains (Stream + Device) |
|---|---|---|---|---|---|
| 1 GPU off | 795 | 268.7 | 49.0 | 78.2 | 318 |
| 1 GPU pipe | 845 | 70.0 | 49.0 | 175.6 | **119** |
| 2 GPUs off | 1185 | 265.6 | 195.8 | 78.2 | 461 |
| 2 GPUs ph | 1151 | 70.6 | 83.0 | 295.0 | **154** |

- hipEventSync is the wait on one kernel. The queue keeps running behind it, so it is not a
  drain.
- The 78 event syncs in "off" come from elsewhere in the code (the SC sweep), not from this
  work.

**What this says:**
- **`implicit_halo_mpi` is the multi-GPU lever.**
  - Strong scaling, radiation part: -10.1 ms (-26 %) on 2 GPUs, -7.1 on 4, -3.2 on 8.
  - Weak scaling, radiation part: -16.3 ms on 2 GPUs, -14.3 on 4, -17.2 on 8.
- **The pipe alone gains 1-2.6 ms**, 2-7 % of the radiation part: 1.6 ms on 1 GPU and
  0.6-2.6 ms on 2-8 GPUs. It removes 200 of the ~320 queue drains per cycle on 1 GPU and
  about 300 of 460 on 2 GPUs.
  - On top of hm it adds 1.3-2.3 ms strong and -0.9 to +2.4 ms weak, i.e. inside the
    repeat spread.
  - With event waits the device already stayed fed, so the drain count was not the
    bottleneck.
  - It costs about 50 more launches per cycle (start-up per solve).
- **Best (ph) vs hydro:**
  - Strong: the radiation part is 32.6 -> 27.0 -> 21.9 -> 21.8 ms on 1/2/4/8 GPUs.
  - It now does shrink with the GPU count (off: 34 -> 38 -> 31 -> 27). It still stops
    between 4 and 8 GPUs, at 1.75x hydro there.
  - Weak: the whole cycle is 2.7-3.0x hydro. It was 3.3-3.9x (off) and 3.3-3.6x in
    runs_3q_scscale.
  - The radiation part at 0.9 M cells/GPU is still 42-50 ms on 2-8 GPUs, against 32.6 ms
    on 1 GPU. The multi-rank overhead is 10-17 ms per cycle, down from 24-35.
- **Radiation is not yet below hydro anywhere.** The best ratio of the radiation part to
  hydro is 1.4 (1 GPU) to 2.0 (weak 8 GPUs).

## 4. Ranked TO-DO toward radiation < hydro, and group-batched multi-group

1. **Overlap the halo with the interior operator (multi-rank, 10-17 ms left).**
   - Every implicit exchange still blocks the host twice: the event after the pack, then
     `MPI_Waitall` before the unpack. The GPU idles through the message latency, about
     2 x 15 halos per solve.
   - Split `ImplicitStencilOp` into interior cells (launched before the wait) and a
     one-cell boundary shell (after the unpack).
   - Or go communication-avoiding: a 2-layer halo of zhat/what so that ONE exchange serves
     K1 + P1. That needs the precondition-and-operator pair on the overlap cells, which the
     block-local rbgs allows.
2. **Iteration count (15.3 inner per solve, 2.7 solves per step; 17.0 at weak 8 GPUs).**
   - The preconditioner is block-local, so it weakens as blocks shrink or multiply (weak
     8 GPUs: +11 % iterations).
   - A semicoarsening (x2-x3) multigrid V-cycle, or at least an rbgs sweep that crosses
     block faces through the same one-message halo, attacks both the count and the halo
     count.
   - Every iteration removed saves 2 halos + 2 preconditioners + 2 operators.
3. **Launches: 800-1150 per cycle.**
   - Remaining synchronous points outside the Krylov loop on 1 GPU: 70 hipStreamSync and
     49 hipDeviceSync per cycle. They sit in the Picard level, the SC sweep and hydro
     (fences, deep_copy); find them from the trace (`runs/prof/*/prof`).
   - With the pipe, the Krylov iteration has only event waits, so a HIP graph of
     K1-P1-K2-P2 (about 10 launches) per iteration is now feasible on 1 rank.
   - Fuse the two rbgs_fwd half-sweeps into one team kernel: both colours within a block
     are local.
4. **SC sweep, 6.6-7.5 ms at 2+ GPUs** (runs_3q_scscale sect. 5), unchanged here.
5. **Production:** turn `implicit_halo_mpi` on for every multi-rank implicit-M1 run. It is
   bitwise, and a candidate default. `implicit_krylov_pipe` is optional: exact to solver
   tolerance but worth only 1-2 ms.
6. **Group-batched multi-group** (the plan of runs_3q_scscale sect. 6, item 5), now with
   concrete carriers:
   - `implicit_halo_mpi` messages take the group index in the component slot (nq = number
     of groups). Still one message per neighbour rank.
   - The pipelined R1/R2 already carry a fixed 6-double MPI type. Widen it to 6 x groups for
     per-group dot products: still ONE `MPI_Iallreduce` per half-iteration, hidden behind
     the operator.
   - The kpw / iw Krylov slots gain a group index next to m. K1, K2, the preconditioner and
     the operator then launch once for all groups.
   - Since the multi-GPU cost is now latency (halos, launches) and not bandwidth, the first
     few groups should ride at well under linear cost. Measure with 2 and 4 groups before
     committing to a design.
