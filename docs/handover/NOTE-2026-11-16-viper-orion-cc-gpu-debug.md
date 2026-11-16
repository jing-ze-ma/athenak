# NOTE 2026-11-16 viper -> Orion (USER): sp/cs GPU divergence found and fixed: uninitialised vet_milne

Re TASK cc_gpu_debug (8e010e76d) and open problem 5.1 of HANDOVER-2026-11-16. The user asked viper to fix it.

## Fix: fork `clean-cc-gpufix-viper` = 073378e14 (one commit on origin/clean-cc 0eeac5fdb, 2 files, +11 -10)
1. `rad_m1.hpp`: the members that only `VetInit` sets get initialisers: `vet_nmu/nphi/nray = 0`,
   `vet_milne = false`, `vet_bc_bath = false`, `vet_x1per = false`, `vet_x1npass = 0`, `vet_axis_flux = false`,
   `vet_full = false`, `vet_eig_min = 0.0`, `vet_dump_every = 0`. Zero is what the CPU and A100 runs read (fresh
   pages), so CPU outputs and ref/ should be unchanged; Cartesian (vet_sc) overwrites them in VetInit as before.
2. `rad_m1_clean_closure.cpp` CleanFaceVector: `p2 = (dir == 2) ? p1 : x3f(m,k+1)`. The phi face of a 2-D
   (axisymmetric) cell read x3f one element past the end (memcheck, below); p2 is not used for dir 2, so bitwise.
Your clean-cc-gpufix audit can take or replace it; viper does not push origin/clean-cc (phase 1).

## Cause
`bool vet_milne;` (rad_m1.hpp:455) is set only in `VetInit` (rad_m1_vet.cpp:1780), which runs for the Cartesian
vet_sc closure only. On sp/cs (vet_gd, clean_rays = fv1) it is never set. On MI300A it read **true**, so
`VetBathLo()` returned -1 (`if (vet_milne || !VetBcEmit()) return -1`), and CleanFvSweep's inner face took
`t = emi ? iemit : ...` = -1: no bath enters the FV rays. g3_spthin cycle 0: CleanRays output `cln_fvi` sums to
5.1e-16 on MI300A vs 7.66e4 on A100, with every input of the sweep bitwise equal (cln_fvo, fvm, fvd, fvq, iw, u0;
cln_fva equal to round-off). Debug print at the sweep: MI300A `emi 1 iemit -1`, A100 `emi 1 iemit 4`
(implicit_ebath_x1min = 4). This explains why Cartesian was exact on both GPUs (VetInit sets it there), why MI300A
was deterministic (bitwise across reruns, with AMD_SERIALIZE_KERNEL=3 and HSA_XNACK=0), and why the sanitizers
missed it (a host read). Found by a checksum bisect of every clean View at each step of CleanClosure/CleanSolve
(temporary patch, not on the branch: /viper/ptmp2/jinma/ccfix_1116/debug_bisect.patch, env CLNDBG).

## Verification
Fixed binaries (incremental on the 32af6349e test trees, same two-file diff) are building on viper (MI300A) and
Raven (A100); G3/G4 on both GPUs vs ref/ follows in the next NOTE. The A100 rhse failure (eint <= 0 retries from
cycle 1, `reason=?` from rad_m1_clean.cpp:779) may be the same class (another uninitialised member read
differently on A100); the rerun will tell. If not, it is the next thing viper bisects with the same patch.

## D1-D4 (TASK cc_gpu_debug), A100 on 32af6349e unless noted
- D1 initcheck (nlim 2/3): g3_spthin, g4_cs_thin1, g3_rhse_sp2 **0 errors** (a host member is invisible to it).
- D1 memcheck (`--report-api-errors no`): spthin 0 errors; g3_rhse_sp2 192 x "Invalid __global__ read of size 8"
  in CleanFvSweep lambda instance 1 = `cln_fv_geo` -> CleanFaceVector x3f(m,k+1) at the 2-D phi face, 1 element past
  the allocation (fix 2 above). Without the flag the count is CUDA-aware OpenMPI pointer queries (noise).
- D1 racecheck spthin nlim 1: 4 shared-memory hazards, all in `hydro::Hydro::CalculateFluxes<hllc>` hflux_x1 team
  scratch (production hydro, not clean); not the cause here. Line numbers follow from a -lineinfo build (Raven job
  32406326, 0eeac5fdb + Kokkos bounds check).
- D2 bounds check: MI300A build of 32af6349e with Kokkos_ENABLE_DEBUG_BOUNDS_CHECK finished; A100 with the job above.
  Not needed for this cause; will be run once on the fixed code.
- D3 first differing cell: g3_spthin after 1 cycle differs in **every** cell (E 78 %, F1 100 %), from the first ray
  build (above), not a localised cell.
- D4 HSA_XNACK=0: no fault, bitwise unchanged (and AMD_SERIALIZE_KERNEL=3 bitwise unchanged).
- Also ruled out on MI300A: the sweep order (keys quantised + index tie-break: bitwise unchanged), the line solver
  (hybrid / pcr / thomas agree to 1e-16), fences before the host reads of unified-memory mirrors (no change).

## G1 cb score flag (handover 5.1)
`implicit_beam` lines of g1_cb_x0.5 are identical on A100 and MI300A:
`implicit_beam0: E_b=1 requested ang=-45 deg -> ray 1 (mu1 mu2 mu3) = (0.707107 -0.707107 0), in-plane ang=-45 deg,
1 ray(s), x2 patch [-0.109375, -0.09375]` and `implicit_beam1: ... ang=45 deg -> ray 0 (0.707107 0.707107 0),
in-plane ang=45 deg, ... x2 patch [0.09375, 0.109375]`. Please compare with the CPU run.log.
