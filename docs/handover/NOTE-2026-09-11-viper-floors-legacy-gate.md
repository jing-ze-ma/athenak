# NOTE 2026-09-11 (viper -> orion): the orion merge changes GPU answers at 1 ULP; gated as `floors_legacy`

Written on viper after pulling a10e367d (orion's polar-average-perf, 10 commits on top of c26b01ac).
HEAD = fork = **caad9247**. Read this before assuming any cross-merge A/B is "unchanged".

## 1. HEAD did not build on HIP (fixed, cdd7d2a5)

`src/hydro/hydro_newdt.cpp` (dt-collapse diagnostic, 1c6792b4) called
`dt_diag.template sync<HostMemSpace>()`. Kokkos rejects that when the DualView device is HIPSpace;
it only compiles where DevMemSpace == HostSpace (orion CPU). Replaced with `modify_device()` /
`sync_host()`, diagnostic path only. Please build GPU (or at least CUDA/HIP) before pushing device code.

## 2. The dhj RT run was NOT reproduced after the merge (A/B on apudev, 400 cycles)

Old = c26b01ac, new = a10e367d, same input, 2 ranks / 2 GPUs. Old-vs-old repeat is bitwise identical,
so GPU determinism holds and the difference is code.

| | cs (cs_ens/si/s01, general-EOS hydro) | sp (sp_mhd_prod3, MHD) |
| --- | --- | --- |
| first difference | cycle 12, one cell, 1 float32 ULP in vely | by cycle 104 (first dump) |
| rel. diff ~cycle 100 | 1.5e-12 dens, eint | 6.6e-12 dens, eint |
| rel. diff cycle 400 | 7e-7 dens, 6e-3 velocity | 1e-10 dens, 1e-3 bcc3 |

Growth after the seed is the known chaotic 1-ULP amplification; mass and total energy stay within
1e-10. Setup: /viper/u2/jinma/ATHENAK/bench/ab_orion/{cs,sp}/{old,new}.

## 3. Bisect: 8da093f5 (EOS floors commit) is the only source

20-cycle cs gate with per-cycle dumps, compared cell-by-cell with bin_convert against c26b01ac dumps
(headers differ because new defaults are recorded, so `cmp` is useless):

| build | content | gate |
| --- | --- | --- |
| a10e367d | all 9 + merge | changed, 1 cell at cycle 12 |
| c26b01ac + commits 1-3 | wb clamp, HLLC containment, wb_rmax | identical |
| a10e367d minus 7-9 | 1-6 + merge resolution | changed, same sequence |
| same minus 8da093f5 | 1-4, 6 + merge resolution | identical |

So the RT switches (4e2e4e0a), conduction (dccd2850), red-giant and diagnostics commits and the
hand-resolved merge are innocent for the dhj run. All new behaviour in 8da093f5 is switch-gated and
OFF in the dhj inputs; the change is code generation (hipcc contraction/reassociation) in rewritten
kernels. The hunk that actually mattered was the **general-EOS ConsToPrim** (general_hyd.cpp +
general_c2p_hyd.hpp), not the cs_raisev rewrite the diff reading had pointed at.

## 4. The gate: `EOS_Data::floors_legacy` (caad9247)

Set once in ReadEOS_Params: true iff `dfloor_keep_velocity`, `vceil`, `eos_floor_consistent`,
`efloor_from_ekin` are all off (the defaults). Under it the pre-8da093f5 kernels run verbatim:

- coordinates.cpp GnomonicEquiangleRaiseVel: legacy par_for launch, new parallel_reduce otherwise
- eos/general_c2p_hyd.hpp, eos/ideal_c2p_hyd.hpp: SingleC2P_*Legacy = the old inversions
- eos/general_hyd.cpp, eos/ideal_hyd.cpp ConsToPrim: the old 3-reducer launch; the switch-aware
  launch is `ConsToPrimFloors` in NEW files eos/general_hyd_floors.cpp, eos/ideal_hyd_floors.cpp
  (added to src/CMakeLists.txt)
- bvals/prolong_prims.cpp: calls the Legacy inversions under the flag

No default changed, no feature removed. Turning any switch on takes the new path (smoke-tested,
all four on, 20 cycles clean).

**Trap for anyone doing this again on HIP:** a legacy kernel in the SAME translation unit as the
new kernel still differed (16523 cells at cycle 12). hipcc inlines both bodies and changes the
contraction. Deleting the new kernel from the file restored identity, which is why the switch-aware
launches live in separate .cpp files. "Same source text" is not "same code" on GPU.

Gates after caad9247, vs the c26b01ac binary, 20 cycles, every cycle: cs ndiff 0, sp ndiff 0.

## 5. Consequences

- The default GPU answer is now bitwise the pre-merge one; sp_mhd_prod3 (built from 9a9396f7) can be
  restarted on caad9247 without a round-off seed.
- Any orion A/B taken between a10e367d and caad9247 on GPU carries the 1-ULP seed; CPU builds may or
  may not (not tested here).
- If orion turns the new floor switches on for production, the legacy path is bypassed and a fresh
  A/B against the switches-off run is needed anyway.
- Style: the repo script check_athena_cpp_style.sh curls the latest cpplint and now reports ~2200
  pre-existing hits at HEAD; cpplint with the repo CPPLINT.cfg gives 0 on the touched files.

Untouched: run/, tables. Leftover worktrees on viper: bench/wt_{c26b01ac,k3,k6,nEOS,nRT}.
Memory snapshot refreshed in docs/handover/claude-memory-2026-09-11/.
