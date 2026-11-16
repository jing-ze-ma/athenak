# NOTE 2026-11-16 viper -> Orion (USER): GPU race in the clean x1-line preconditioner (Kokkos HIP level-1 scratch), fixed

USER asked to tell you. FYI, plus one request at the end.

## Symptom
The clean scheme was not reproducible run to run on the real R136 state, on MI300A with 2 ranks.
- Same binary, same keys: most repeats were bitwise equal, some were not. Differences were 1e-4..3e-3 relative, from cycle 1, and once grew into 2 dt_retries.
- Independent of print code and of the kernel timer.
- `AMD_SERIALIZE_KERNEL=3 HIP_LAUNCH_BLOCKING=1` does NOT remove it, so it is inside a kernel, not a host/device fence race.

## Where (device checksum bisect of every BiCGStab vector)
- r136.00093 (cycle 487121): solve 19, BiCGStab iteration 2.
- y = M^-1 p, r0.v and alpha are bitwise equal in both runs.
- The first difference is z = M^-1 s from the x1-line preconditioner on bitwise-equal input.

## Cause
- Both x1-line kernels keep each line's rows in Kokkos level-1 team scratch: CleanCellLineSolveRows<kStore> and CleanCellLineFactorSolve (the latter from the opt-in clean_pc_factor). At nl = 738 that is 614 KB and 47 KB per team.
- Kokkos 4.6 HIP hands a released scratch slot to the next team of the same launch with no memory fence: kokkos/core/src/HIP/Kokkos_HIP_TeamPolicyInternal.hpp:369-404, hip_get/release_scratch_index. The CAS acquire has no fence, and the release is a plain store after __syncthreads.
- So a team can read rows a previous team left. The exact hardware path (L2 not coherent across MI300A XCDs, or late stores) is inferred, not proven.
- Proof: the fix below makes all repeats bitwise.

## Fix
Fork clean-cc-race-viper dad6ddabd, on 2503efe6c = the R136 cost keys; folded into the push-16 candidate fork clean-cc-push16 a97e368cb.
- Each team owns a fixed block of a new work array cln_lw for the whole launch.
- League = teams that run at once; each team strides over the lines.
- The arithmetic is unchanged. Files: rad_m1_clean_cell.cpp, rad_m1.hpp.

## Evidence (MI300A, R136 r136.00093, 2 ranks)
| arm | result |
|---|---|
| fixed, default keys | 6/6 bitwise |
| fixed, opt-in keys (pc_factor + pc_reuse=2 + conv_est + nl_norm=cell, efloor 1e-4) | 6/6 bitwise |
| fixed, pc_factor alone | 4/4 bitwise |
| old binary, pc_factor | every repeat differs |
| G1/G3/G4 (18 decks, incl. 3 x1line variants) | bitwise to the pre-fix binary |

## Same hazard outside the clean code (not fixed, not observed)
- hydro_fluxes falls back to level-1 scratch when (nvars + nder) * ncells1 * 8 * 2 bytes > 64 KB. With a general EOS (nder = 2):
  - R136 (nx1 738): 83 KB -> level 1.
  - HG (nx1 704): 79.5 KB -> level 1.
  - RSG (nx1 548): 62 KB -> borderline.
- No hydro divergence was seen: hydro_w dumps were identical in all 16 fixed repeats.
- Any other kernel that uses level-1 TeamPolicy scratch on HIP has the same exposure.
- A Kokkos-side patch is proposed but NOT applied anywhere: viper /viper/ptmp2/jinma/ccrace_1116/kokkos_hip_scratch_fence.patch. It adds __threadfence after acquire and before release, and an atomic unlock.

## Request (CPU, low priority)
On CPU the x1-line kernels run with one-thread teams (push 15 rule), so the fix should be bitwise there. It is in the push-16 TASK (fcc473cec, head a97e368cb). If you run any HIP or CUDA build on your side, check for level-1 scratch users in your stag2 code.
