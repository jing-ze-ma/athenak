---
name: index-bugs-performance
description: Memory links for code bugs, GPU/CPU performance, kernels and build traps.
metadata:
  type: reference
---

## Bugs and performance

- **[FOFC on GPU VERIFIED](fofc-gpu-verified.md) — 48/48 HIP Release; Debug abort was a pre-existing dt-diagnostic over-read, fixed b18f3205**
- **[ORION MERGE round-off gated caad9247](orion-merge-roundoff-8da093f5.md) — `floors_legacy`; general-EOS c2p was the hunk**
- **[All-Mach solvers lhllc/ausmpup: algebra matches the papers](lowmach-solver-audit.md) — no CFL~M restriction**
- **[REFLECTING WALLS leaked mass under a blast: fixed](reflect-wall-mass-leak.md) — mirror the wall state**
- **[GPU nondeterminism fixed 6e600f12](dhj-run-to-run-nondeterminism.md) — missing team_barrier**
- [Polar MPI host-mirror bug fixed](polar-mpi-host-mirror-bug.md) — 3882e37f
- [Restart bug with tfloor_kelvin fixed](athenak-restart-tfloor-bug.md) — ba2f0943
- [dhj ideal input out-of-bounds fixed](dhj-ideal-input-oob.md) — 0d1f6f6a; use a Debug build
- [MHD flux stale x1 limits](mhd-fluxes-stale-x1-limits.md) — 3429f59f; latent
- [RT chain-parallel split kernel](rt-chain-parallel-split.md) — 38311a8a, live
- **[MAGNETIC ENERGY is not in the dumps](mhd-energy-not-from-dumps.md) — bcc is the mean of the fields**
- [rt_chain_ck SWEEP SPLIT measured 1.13x, shelved](rt-ck-sweep-split-measured.md) — patch in bench/prof/rt_new
- [Overall GPU profile of a dhj run](dhj-overall-gpu-profile.md) — RT 35.7%, fluxes 27.8%
- [Meshblock decomposition on one GPU](meshblock-decomposition-gpu.md) — 32 blocks beats 2 by 1.4x
- [nx1 ceiling from LDS](nx1-ceiling-lds.md) — MHD capped at nx1=264; hydro escape `hydro/scratch_level=1` 955c39df
- [par_for_outer team size](par-for-outer-team-size.md) — ~3% of wall, not worth it
- RT kernel history, superseded by rt-chain-parallel-split: [interleaving](inflight-rt-kernel-optimization.md), [band scaling](inflight-rt-band-scaling.md), [occupancy](rt-kernel-occupancy-limit.md)
- [IN FLIGHT: nx2=64 + GPU timings](inflight-nx64-timings.md) — paused 2026-08-16
- [AMR null-tree bug](amr-broken-on-branch.md) — fixed; the cubed-sphere refactor left `Mesh::ptree` unconstructed, and the same half-migration may bite elsewhere
- [Pre-existing branch breakage](branch-preexisting-breakage.md) — non-MPI build fails in bfield_bcs.cpp; hydro/mhd_linwave test scripts pass a `vflow` param that no longer exists
- [Event log MPI deadlock](eventlog-mpi-deadlock.md) — FIXED 27da6380: file_type=log hangs ANY clean MPI run (early return skipped the schedule advance -> Allreduce every cycle); likely upstream too
- [FOFC 1D segfault](fofc-1d-segfault.md) — pre-existing crash on the branch: FOFC + 1D hydro + PLM, not caused by the EOS work
- [FOFC COMPATIBILITY COMPLETE 09-11 incl. MHD on cs + static/dynamic WB tests (new pgen wb_atm; dynamic WB: fofc raises the Mach ceiling 39 -> >80) + MPI/SMR tests; COMMITTED b3345f55 + merged 7cee2eb7 on 09-12, not pushed; 49 fofc tests pass, bit-identical off; left: red-giant fofc run, GPU check, wiki, commit split from red-giant hunks](fofc-compatibility-plan.md)
- [hst uses Cartesian cell volume](hst-cartesian-volume-on-spherical.md) — .hst mass/tot-E are NOT physical integrals on the spherical dhj mesh
- [pcoord arrays unallocated in Cartesian](pgen-cartesian-coordinate-arrays.md) — which geometry arrays are placeholders on a Cartesian mesh, and the solar/cooling_convection audit for unguarded uses
- [RT source is semi-implicit](rt-source-semi-implicit.md) — FIXED 048dff30: exponential relaxation at rate 4E/e; cured the red giant dt collapse; the wall control is what ruled out the boundary
- [RT clip warning is rank-local](rt-srclim-warn-rank-local.md) — BUG, unfixed: the clip count is rank 0's share, and the warning is silent when rank 0 never clips
- [RT source in a transparent cell: RESOLVED](rt-transparent-cell-cancellation.md) — direct source in all 3 kernels + the REAL cause: the opacity lookup CLAMPED rho to the table edge (1e-14) and returned molecular opacity; table now extended to 1e-24 at the floor; background holds 411 K
- [two_stream_rt has THREE sweep kernels (grey/ck/generic)](two-stream-rt-three-kernels-trap.md) — a sweep edit must hit all three; red_giant uses rt_chain_grey; verify with rt_apply_debug. Also: <hydro>/wb_rmax added+verified
- [VIPER floors_legacy GATE 09-11 (caad9247): orion merge changed GPU answers 1 ulp via the general-EOS c2p; legacy kernels in separate .cpp under floors_legacy (all floor switches off); MERGED locally 7cee2eb7 on 09-12](viper-floors-legacy-gate.md)
- [VIPER MERGE LANDED 09-11 08:05: a10e367d on top of 9 topical commits; gates bit-identical both sides; not pushed](viper-merge-trial-2026-09-11.md)
