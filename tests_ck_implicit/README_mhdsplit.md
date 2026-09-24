# mhd-split: the operator-split hook (user_split_func) on the MHD task list

## The bug
- `ProblemGenerator::user_split_func` (src/pgen/pgen.hpp) was called only by `Hydro::RTStrangSplit`.
  Only `AssembleHydroTasks` adds that task to a task list.
- The MHD task list never called the hook. On `<mhd>`, the following did nothing, and no error was raised:
  - deep_hot_jupiter_rt `problem/ck_impl_once` (and so the c2 lever set and `ck_impl_every`);
  - `DhjRtFreezeSplit`, which is hydro-only and refused anyway;
  - red_giant `rt_strang` / `rt_once_per_cycle`.
- Found by the weak-scaling agent, job 11955965 (/viper/ptmp2/jinma/cs_weak_0924/README.md).

## The fix
- `MHD::RTStrangSplit` (src/mhd/mhd_tasks.cpp) mirrors `Hydro::RTStrangSplit`:
  - Strang mode: dt/2 in "before_timeintegrator" (stage 0) and dt/2 in "after_timeintegrator";
  - `user_split_once`: one call with the full cycle dt `pmesh->dt`, after the last RK stage.
  - Each call is followed by `MHD::RTOpSplitBvals`, a synchronous tail on u0:
    InitRecv, RestrictCC, PackAndSend, RecvAndUnpack, ApplyPhysicalBCs, Prolongate, ConToPrim, ClearSend/Recv.
  - Only u0 is exchanged. The split source does not touch b0, so the b0 ghosts that
    ApplyPhysicalBCs/Prolongate refill come out the same.
- `SaveMHDState` now depends on the pre-integrator split task, so it saves the state the stages start from.
  With no hook the task returns at once.
- Difference from hydro: hydro runs an extra, unconditional `ConToPrim` task after each split task (hydro_tasks.cpp).
  The MHD version runs ConToPrim once, inside RTOpSplitBvals, and only when a hook is enrolled.
  An MHD run without a hook therefore executes exactly the old task sequence.
- Hydro + MHD (ion-neutral), radiation and NR assemble neither list, so nothing calls the hook twice.
- New startup fatal: `MeshBlockPack::split_hook_tasks` is set by AssembleHydroTasks and AssembleMHDTasks.
  Both ProblemGenerator constructors fatal if user_split_func is enrolled and no list calls it.
  - Checked: red_giant rt_strang with `<hydro>` + `<mhd>` + `<ion-neutral>` stops with the new message
    (/viper/ptmp2/jinma/mhdsplit_0924/fatal).

## Gates
Everything is in /viper/ptmp2/jinma/mhdsplit_0924.
- Binaries: base = e89954e2 (rt-integration), new = this branch. Scripts: build.sh, cpugate.sh, cpugate2.sh, rstref*.sh, ab.sub.
- MHD state: bench/cs_mhd_prod4/rst/dhj.00127.rst, read in place.
  - t = 1.93151e7 s = rotation 63.3 (P = 2 pi / omega = 3.050e5 s), ncycle 1192186.
  - Input: the prod4 input plus the ck_impl keys at their defaults (mhd_ab.athinput).

### CPU gates (login node)
| gate | result |
|---|---|
| hydro dhj c2 (fast/wp_fast.athinput, sp, 6 cycles), base vs new | hst + rst BITWISE |
| hydro red_giant rt_strang (two_stream_sph_thick, n0 case, 6 cycles), base vs new | hst + rst BITWISE |
| MHD prod4 rst, production config (explicit ck, no hook), 3 cycles, base vs new | hst + rst BITWISE |
| MHD prod4 rst, T4 per stage (no hook), 3 cycles, base vs new | hst + rst BITWISE |
| MHD c2 (new) | one `### ck_implicit` report line per cycle (passes 3-4, 0 NOT-CONVERGED) |
| MHD c2 + ck_impl_every=4 (new), 8 cycles | full calls at ncycle 1192186 (first), 1192188 and 1192192 (ncycle % 4 = 0), linearised steps between |
| MHD restart, 6 straight vs 3 + restart + 3 | see below |

Restart: maximum relative difference of T in the final rst, against the straight run.
- Production (explicit ck): BITWISE.
- T4 per stage: 6.2e-6 (d 5.5e-7).
- T4 + once (the hook alone): 1.1e-5 (d 2.9e-7). The hook restarts like T4.
- c2: 0.12, rms 4.0e-4 (d 2.7e-3), in 14 cells at 1e-7..1e-5 bar.
  - The c2 levers (pred, xstep, jreuse) keep call-to-call state that the restart does not hold.
  - README_cadence says the same for hydro: c2 is not a bitwise restart.

### GPU A/B (apudev, 2 GPUs; jobs 11957490 and 11957491)
- Setup: MHD prod4 rst 127 to t = 1.93456e7 (0.1 rotation, about 1500 cycles); `fast/ab_ana.py <root> t4`.
- An earlier pair (11956962/63, 11957353) is kept in ab_old/ but is invalid.
  - ab.sub took T0 from `<output4>/last_time`, which is the next dump time, not the restart time.
  - t4 and t4x then hit the wall clock at a different time than the other arms.

| arm | rms day | max day | rms night | passes/call | non-conv | ms/cycle |
|---|---|---|---|---|---|---|
| t4 (T4 every stage, no hook) | 0 | 0 | 0 | 4.00 | 0 | 131.1 |
| t4x (tolerance 1e-10 only) | 2.2e-6 | 1.4e-4 | 4.8e-3 | 5.80 | 0 | 148.6 |
| o (T4 + once = lever 1) | **2.6e-4** | 1.6e-2 | 1.2e-2 | 4.00 | 0 | 78.2 |
| c2 (mhd-split) | **2.6e-4** | 1.6e-2 | 1.2e-2 | 3.00 | 0 | 52.2 |
| e4 (c2 + ck_impl_every=4) | 2.6e-4 | 1.6e-2 | 1.2e-2 | 3.01 | 0 | 34.1 |
| c2b (base binary: hook never called, RT off) | 1.9e-1 | 9.4e-1 | 8.9e-1 | - | - | 20.8 |

- On MHD, c2's day-side error is lever 1's: 2.6e-4. Hydro measured 3.0e-4 (README_fast).
- The base binary's c2 is the RT switched off: 0.19.
- On MHD, c2 is 2.51x cheaper than t4 and e4 3.84x.
