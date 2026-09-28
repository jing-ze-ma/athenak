---
name: dhj-run-to-run-nondeterminism
description: "REOPENED 2026-08-24: c39b794c fixed the outer-x1 bcc0 race, but GPU runs are STILL not reproducible run to run -- intermittent, correlated-k only so far, root cause unfound"
metadata: 
  node_type: memory
  type: project
  originSessionId: be8c40b1-14bc-4240-8f72-1fdd200e6113
  modified: 2026-08-19T15:59:31.919Z
---

## FOUND AND FIXED 2026-08-24: a missing `member.team_barrier()` in the x1 flux kernel

**The bug.** `src/hydro/hydro_fluxes.cpp` (and the identical site in
`src/mhd/mhd_fluxes.cpp`), general-EOS branch of the x1 direction:

```
switch (recon_method_) { ... PiecewiseLinearX1(member, m, k, j, il-1, iu, wder_, dl, dr); }
// reconstruction can undershoot; keep the interface pressure positive
par_for_inner(member, il, iu, [&](const int i) {
  dl(IDPR,i) = fmax(dl(IDPR,i), eos_.pfloor);   // reads a slot thread i-1 wrote
  dr(IDPR,i) = fmax(dr(IDPR,i), eos_.pfloor);
});
...
member.team_barrier();      // the only barrier -- and it is AFTER the floor
```

`PiecewiseLinearX1` writes `dl(n,i+1)` from the thread that owns `i`, so the floor loop
reads a scratch slot written by ANOTHER thread. Successive `par_for_inner` calls are not
implicitly synchronised, so whether the floor sees the reconstructed value or a stale one
is scheduling-dependent. **x2 and x3 do NOT have this**: those reconstructions write index
`i` from thread `i`. That is exactly why only the x1 flux ever differed.

**Scope.** The floor loop only exists under `if (nder > 0)`, i.e. **a general EOS only** --
`nder = eos_.IsGeneral() ? NDERIVED : 0`. Ideal-gas runs were never affected. This is why
the nondeterminism showed up in the tabulated-EOS work and nowhere else, and it is core
AthenaK code, not the pgen.

**The fix**: `member.team_barrier();` immediately before that `par_for_inner`, in both
`hydro_fluxes.cpp` and `mhd_fluxes.cpp`. Cost measured at **1.0 %**
(2.363e7 -> 2.339e7 zone-cycles/s, 400 cycles, same config).

**How it was found** (the chain that worked, after configuration-bisection stalled):
per-task checksums. A hook in `TaskList::DoAvailable` calls a function that prints the
wrapping integer sum of the raw bit patterns of `u0`, `w0`, `wder`, `wtemp` and the three
flux arrays after every completed task -- an EXACT, order-independent checksum, unlike a
floating-point sum which can cancel a difference away. Ten replicate runs, grouped by the
hash of their checksum stream: 8 agreed, 2 differed. The first differing line named it
outright -- `task=1` (`Hydro::Fluxes`) with `u`, `w`, `d`, `t`, `f2`, `f3` all bitwise
identical and **only `f1` different**. Identical inputs, different x1 flux.

**Verification** (clean Release build, 3 replicates each, 600 cycles, bitwise on every dump
AND the history):

| config | before | after |
|---|---|---|
| hydro + polar, no source terms | diverged at cycle 20 | **identical** |
| MHD + correlated-k + polar (production) | diverged at cycle 40 | **identical** |
| hydro, polar off, sources on | diverged at cycle 120 | **identical** |

Ten replicates under the checksum instrument also collapse to ONE stream, and that stream
is the majority (correctly floored) one from before the fix -- the barrier selects the right
branch, it does not invent a third answer.

**CORRECTION to what I wrote earlier today:** I claimed "at least TWO independent races",
because polar-off and `do_pole`-disabled builds looked clean. They were not clean -- the
same single bug was present and simply took longer to grow into the observable. One fix
makes every configuration reproducible. The earlier `do_pole` A/B was a false lead:
disabling it changed how fast a perturbation amplified, not whether one was created.

**Consequence for the science.** Every general-EOS GPU run in this campaign was
non-reproducible, and the interface pressure floor was applied or skipped at random in a
few cells per column. That is the source of the scatter in the dhj blow-up times
(0.370 / 0.485 / 0.584 / 0.685 / 1.464 rotations on identical input) -- see
[[dhj-ck-eos-blowup]]. Those comparisons should be re-run on a fixed binary before any of
them is trusted.

**State: COMMITTED as 6e600f12** on `polar-average-perf` (NOT pushed), 3 files:
`src/hydro/hydro_fluxes.cpp`, `src/mhd/mhd_fluxes.cpp`, and the new regression test
`tst/test_suite/nr/test_nr_geneos_repro_gpu.py`. All instrumentation was reverted; no
diagnostic code remains in `src/`. Verified binary:
`/viper/ptmp/jinma/claude_eos_gpu/athena_fixrel`.

**The regression test.** Runs the same binary three times on the shipped general-EOS dhj
input and requires bitwise-identical histories. `_gpu`, because a serial `par_for_inner`
runs in order and cannot expose the hazard. It builds its own binary (the setup is a user
pgen) and inherits the device from `tst/build/CMakeCache.txt`, so it targets whatever
`run_test_suite.py` was invoked for; it SKIPS when that cache shows no accelerator. No
correlated-k tables needed -- the shipped input defaults to grey.

**Validated both ways under pytest on an MI300A**: FAILED without the fix ("run 1 differs
from run 0"), PASSED with it. ~3 min including the build.

**Do not shrink the test grid.** Detection is very sensitive to configuration, measured:
production grid 64x64x128 (mb 16 or 8) = 3 of 3 replicate pairs differ; 16x16 and 32x32 =
0 of 3; a general-EOS 3D linear wave peaked at 2 of 5 and went to 0 of 5 with more blocks.
The race needs the pressure floor to bind on SOME interfaces but not all -- with pfloor
above every reconstructed value both orderings agree -- and enough resident teams for the
two inner loops to overlap. A smaller grid would pass while testing nothing.

To reproduce the bug for future work: `/viper/ptmp/jinma/claude_eos_gpu/athena_head`
(397b49b1, pre-fix) vs `athena_fixrel`, scripts `dhjtune*.sh`, `chk*.sh` there.

Also useful: viper's `python3` has no pytest; `module load python-waterboa/2025.06`
provides pytest 8.3.4, numpy and flake8.


## 2026-08-24 (II): DIRECT INSTRUMENTS -- the pack and unpack are BOTH exonerated

Three instrumented builds, each a temporary edit to `src/bvals/bvals_cc.cpp` that was
**REVERTED immediately after building** (`src/` is clean; the binaries live in
`/viper/ptmp/jinma/claude_eos_gpu`):

1. `athena_nopole` (`-DDHJ_POLE_PACK=0`) -- forces `do_pole = false`.
2. `athena_selfchk` (`-DDHJ_PACK_SELFCHECK=1`) -- **double-pack comparison**: wraps the
   "SendBuff" kernel in a 2-pass loop, saves every `recvbuf[n].vars` to a host mirror after
   pass 0, and compares after pass 1. Packing is a pure function of `u0`, so any difference
   is a race inside the kernel, with NO dependence on the dynamics.
3. `athena_ucnt` (`-DDHJ_UNPACK_COUNT=1`) -- **write counter**: an int array shaped like
   `u0`, zeroed before "RecvBuff", `Kokkos::atomic_add` at each unpack write, then a host
   scan for any cell written more than once (which would be two buffers targeting one cell).

| instrument | polar ON | polar OFF |
|---|---|---|
| double-pack comparison, 30 cycles | **0 differing buffer entries** | 0 |
| unpack double-write count, 4 cycles | **none** | none |

So **the pack kernel is internally deterministic and no ghost cell is written twice by the
unpack.** (The `athena_ucnt` runs end in a benign `Kokkos allocation "wcnt" is being
deallocated after Kokkos::finalize` abort -- the static View outliving finalize. The cycles
themselves completed; the result stands.) `hydro_bcs.cpp` has `default: break;` for x2, so
physical BCs never touch polar ghosts either -- no second writer there. AthenaK uses a
single execution space (no `partition_space`, no streams), so kernels cannot overlap.

**The observable, localised.** Fastest reproducer (hydro, polar, no source terms), first
divergence at cycle 20: **only `vely` and `velz` differ** -- dens, velx and eint are bitwise
identical -- across ~620 cells in the pole-adjacent blocks (lx2 = 0 and 1). Those are
exactly the two variables the polar exchange sign-flips (`signvar = -1` on IVY/IVZ), and
they are ~1e-10 here, i.e. near-cancellation quantities that amplify any perturbation.
Max |v_theta| is 7-8e-10 in the diverging AND in both reproducible runs, so the amplifier is
equally present when the answer is reproducible -- that argues against, but does not fully
kill, the reading that polar merely SUPPLIES the amplifier for a race living elsewhere.

**Where that leaves it.** The `DHJ_POLE_PACK=0` A/B says the `do_pole` branch is necessary
for the divergence; the double-pack check says that branch is not itself racy. Both can be
true only if the race lives in a kernel DOWNSTREAM that consumes polar ghost data -- and is
therefore worth hunting without reference to the pole at all.

**NEXT INSTRUMENT (not yet built):** checksum `u0` after every task in the stage, run twice,
and report the first task whose checksum differs. That names the guilty kernel directly
instead of bisecting configurations. Everything needed is in `driver.cpp`'s task loop.

## 2026-08-24 BISECTION -- localised to the `do_pole` branch of the boundary pack

Method: every row is TWO byte-identical invocations of the same binary, one GPU rank,
compared BITWISE on every bin dump and on the history. `athena_head` = 397b49b1 Release.
Scripts `det*.sh` in `/viper/ptmp/jinma/claude_eos_gpu`.

| config | cycles | result |
|---|---|---|
| MHD + c-k, polar ON (baseline) | 600 | diverges, first at cycle 40 |
| MHD + grey, polar ON | 600 | diverges |
| MHD, polar ON, `user_srcs=false` (RT off) | 600 | diverges |
| MHD, polar ON, user BCs -> outflow | 600 | diverges |
| MHD, polar ON, `bbot=1e-30` | 600 | diverges |
| HYDRO + c-k, polar ON | 600 | diverges |
| **HYDRO, polar ON, `user_srcs=false`** | 600 | **diverges at cycle 20 -- FASTEST REPRODUCER** |
| **MHD, polar OFF** (`use_polar_boundary=false`, x2 reflect, x2min/max given) | **3000** | **BITWISE IDENTICAL** |
| **HYDRO, polar OFF, `user_srcs=false`** | **3000** | **BITWISE IDENTICAL** |
| **HYDRO, polar ON, no srcs, built `-DDHJ_POLE_PACK=0`** | 600 | **BITWISE IDENTICAL** |
| MHD, polar ON, srcs ON, built `-DDHJ_POLE_PACK=0` | 600 | diverges -> a SECOND source |
| HYDRO, polar OFF, srcs ON (grey or c-k, split or not) | 600 | diverges -> the second source again |

`DHJ_POLE_PACK` was a temporary compile switch forcing `do_pole = false` at
`bvals_cc.cpp:104` and `bvals_fc.cpp:106`. **The source edit has been REVERTED**; the
binary survives as `/viper/ptmp/jinma/claude_eos_gpu/athena_nopole`.

**The amplification confound is EXCLUDED.** A "clean" run could in principle just be one
whose flow never amplifies a perturbation. Measured between dump 0 and dump 30 of the
diverging run and both clean ones: RMS d(rho)/rho = 5.149e+04 in ALL THREE, identical to
four digits, with comparable max |v_theta|. The reproducible runs are exactly as active.

**It is NOT an out-of-bounds access.** A full `-DCMAKE_BUILD_TYPE=Debug` build (Kokkos
bounds checking) runs the fastest reproducer 40 cycles with no violation.
Binary: `/viper/ptmp/jinma/claude_eos_gpu/athena_dbg`.

**So: at least TWO independent races.**
1. **In the `do_pole` branch of `MeshBoundaryValuesCC::PackAndSendCC`** (bvals_cc.cpp
   ~104-197, and presumably the bvals_fc twin). Disabling that branch alone makes an
   otherwise-identical, equally active hydro+polar run bitwise reproducible.
2. **A second one, not yet isolated**, that shows up as MHD+polar with the pack branch
   disabled, and as hydro+`user_srcs` with polar off. Candidates for the MHD half:
   `PolarAzimuthalAverageEr` (mhd_corner_e.cpp:460+, which READS as carefully written --
   per-block partials then an ascending-order combine) and `polar_local_sum_b`
   (bfield_bcs.cpp:339, NOT yet read).

**What I could NOT find by reading.** The do_pole branch differs from the normal path only
in the mirrored source index `jj = jl+ju-j`, the `signvar = -1` on IVY/IVZ, and the
destination slot -- and the destination remap is NOT the culprit, because the
`DHJ_POLE_PACK=0` build keeps the polar `dest` overrides from meshblock.cpp and is clean.
The buffer index arithmetic is algebraically identical to the normal path. I enumerated
the polar `dest` mapping for this mesh (1 x 4 x 8 blocks) at faces and x2x3 edges and
found no two senders writing the same (receiver, slot).

**NEXT STEP, not yet done:** instrument rather than read. Either (a) call the pack kernel
twice and compare the buffers -- packing is a pure function of `u0`, so any difference is
the race, with no dependence on dynamics; or (b) count writes per buffer slot with an
atomic counter and look for a slot written twice.


> **REOPENED 2026-08-24.** c39b794c fixed ONE race (the outer-x1 `bcc0` cross-thread read)
> but the code is still nondeterministic: two byte-identical invocations of the same binary
> on one GPU rank diverge. It is intermittent, appears only in correlated-k runs so far
> (2 of 3 pairs; grey 0 of 2 over 32 cycles), and first shows as B_phi differing in the
> outermost radial cells on meshblock corners at cycle 22. Evidence, jobs and the module
> A/B table are in [[dhj-ck-eos-blowup]], section "IT IS NOT THE TIMESTEP".


**The deep_hot_jupiter_rt GPU run is NOT bitwise reproducible.** Same binary, same input, ONE
MPI rank, one MI300A: two consecutive runs differ. Established 2026-08-19 with jobs 10957531
(200 cycles) and 10957548 / 10957572 (1 cycle). The t=0 dump is identical; divergence appears
within the FIRST cycle, so it is injected per-cycle, not accumulated chaos.

**Magnitudes after one cycle** (absolute, which is the only honest view here):
`dens` and `eint` differ by exactly 1-2 ULP; `bcc2` 1.5 ULP; `velz` 13; `vely` 23; **`velx` 53
ULP median**, max abs diff 6.6e-2 = 1.1e-3 of that field's rms. `bcc3` looks catastrophic in
relative terms (median 3%, max 7.6e20 ULP) but is an artefact -- at cycle 1 the azimuthal field
is still ~0 (rms 2.8e-7 against bcc1's 1.16), so it is near-zero over near-zero. Its absolute
difference is 2.4e-8.

**FIXED 2026-08-19 in commit `c39b794c`** (branch `polar-average-perf`, unpushed): the write was
split into its own kernel `usrboundaryx1_bcc_outer`, and the three out-of-range reads clamped.
Verified bitwise reproducible afterwards at 1 and 200 cycles in the production configuration.
**Runs made before `c39b794c` have a wrong outer-x1 magnetic boundary whenever
`bc_outer_maxwell = true`** -- against the fixed code, differences above 1% of field rms reach
radial index 45 of 64 and grow outward (96% of cells in the outermost shell for vely), while
density is untouched (no cell above 1% of rms). This includes the long runs 10955347/10955348.

**ROOT CAUSE (found by bisection):** a data race in
`src/pgen/deep_hot_jupiter_rt.cpp`, in the `bc_outer_maxwell` block of the outer-x1 user BC.
The kernel at line 1043 runs one thread per `(m,k,j)`; each thread WRITES
`bcc0(m,IB*,k,j,ie+i+1)` (lines 1102-1108) and then READS `bcc0(m,IBX,k,j+1,ie+i+1)` (line
1120) and `bcc0(m,IBX,k+1,j,ie+i+1)` (line 1122) -- cells written by OTHER threads in the same
launch. Whether a thread sees the new or the stale value depends on wavefront scheduling.

Bisection that proved it (1 cycle, same binary twice, `bench/polar_ab/`): `user_srcs=false`
still differs; resistivity removed still differs; `ix1_bc=ox1_bc=outflow` DETERMINISTIC;
`bc_outer_maxwell=false` DETERMINISTIC; inner-user-BC-only DETERMINISTIC; outer-user-BC-only
differs.

**Why it was never seen before:** on CPU the j loop runs in ascending order, so the read always
gets the stale value -- wrong, but deterministically wrong. Only a GPU exposes it. Two further
defects in the same lines: line 1116 reads `bcc0(...,ie+i+2)`, which the same thread writes on
the NEXT i iteration (read-before-write, so last cycle's value); and `j+1`/`k+1` at
`j=n2-1`/`k=n3-1` read one past the end of those dimensions.

**Fix shape:** compute the cell-centred ghost `bcc0` for all ghost cells in a separate kernel,
then do the Maxwell correction in a second kernel, so every read is of a completed write.

**Localisation:** radial indices 33-63 of 64 -- the OUTER HALF -- across 29 of 32 MeshBlocks and
every theta row. So it is NOT the polar boundary. It is also not `dt`: both runs print identical
times, and a perturbed dt would move every cell, not 0.02-0.6% of them. Leading suspect,
UNTESTED: the outer atmosphere is where dfloor/efloor fire at 8-11% of C2P calls
(see [[xe-resistivity-long-runs]]). No atomics exist in the MHD evolution path and there are no
WithoutInitializing allocations, so the usual causes are already excluded.

**It is pre-existing** -- reproduced with the pre-006bae07 binary snapshot in
`bench/xe_long/b3_e13/athena`. Not caused by the polar-average rewrite.

**Why it matters for verification:** any "is this change bitwise identical" test on this problem
has a noise floor. The polar rewrite's apparent effect matched the run-to-run noise variable for
variable (max abs diff / rms of 2.06e-10 for dens, 1.44e-05 for bcc3 in BOTH comparisons), so a
`cmp` of two output files proves nothing here. **Always compare a change against a same-binary
rerun, not against a single baseline run.**

**Repro:** `/viper/u2/jinma/ATHENAK/bench/polar_ab/` holds `onecycle.athinput`, `repro.sh`
(same binary twice, 200 cycles) and `oldone.sh` (pre-change binary, 1 cycle).
