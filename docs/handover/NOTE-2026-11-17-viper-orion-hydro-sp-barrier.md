# NOTE 2026-11-17 viper -> Orion (USER): CORRECTION, the hydro race is a missing team_barrier in the sp x2/x3 sweeps, not Kokkos scratch

This corrects item 1 of cdeca080a. The level-1 scratch explanation was wrong; the race agent found the real cause.

## Bug
- **Where:** `src/hydro/hydro_fluxes.cpp`, the x2 and x3 flux sweeps on spherical-polar meshes.
- **What:**
  - GridPiecewiseLinearX2/X3 reconstruct wl/wr into team scratch over [is-1, ie+1].
  - Then AddWbPrimFaceX2/X3 (wellbalance_static_reconst) and the sp x3 velocity rotation (polar_x3_shift = rotate) read-modify-write the
    same scratch entries over [il, iu], with **no member.team_barrier() in between**.
  - When [il, iu] differs from the reconstruction range (il = is unless fofc or sp_favg widens it to is-1), a lane reads entries another
    lane is still writing.
- **Effect:** a GPU-only, run-to-run nondeterminism, usually small. One MI300A test run produced 40 non-finite x3 fluxes, which then gave a
  BiCGStab breakdown in clean. CPU teams have one thread, so the CPU never shows it.
- **Not affected:** MHD, whose AddWb and rotation use the same [is-1, ie+1] range as its reconstruction.
- **How it was found:** a checksum bisect on 3 repeats of g6_sp1x1wb_base (sp, WB reconst). All inputs were bitwise, yet `hflux_x2` / x3
  differed at cycle 15, in sparse (k, j) rows.
- **What ruled out the alternatives:**
  - fences around par_for_outer: the runs still diverge;
  - the Kokkos HIP scratch fence patch: still diverge;
  - AMD_SERIALIZE_KERNEL=3: 2 of 3 runs still diverge, so the race is in-kernel.
  - The hydro scratch here is level 0 anyway.

## Fix
- One `member.team_barrier()` after the sp reconstruction in the hydro x2 and x3 sweeps, on fork branch clean-cc-race2-viper (6 repeats and
  gates running).
- viper will also put it on a branch off the production base, so production binaries can be rebuilt if the user decides.
- A full audit of every team-scratch kernel for the same class (missing barriers, single/broadcast, level-1 scratch, uninitialised scratch
  reads) is running on viper. We will send the findings.

## Exposure in productions
- **Exposed:** sp hydro with wellbalance_static_reconst, or sp hydro without fofc / sp_favg, on GPU.
- **Not exposed:** M35MMS (WB reconst off; fofc widens il to is-1, so the same lane writes and then reads).
- **Still being checked:** R136 and cubed-sphere runs.

## Request
- Check your staggered and production branches for the same pattern: any par_for_inner that writes team scratch, followed by a second
  par_for_inner over a different index range reading or modifying those entries, without team_barrier.
- The fix changes nothing on CPU.
