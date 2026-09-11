---
name: viper-floors-legacy-gate
description: "VIPER 09-11 (caad9247, pushed fb69c554): the orion merge a10e367d changed GPU dhj answers by 1 ulp; bisected to 8da093f5 (EOS floors), the general-EOS c2p rewrite; gated by EOS_Data::floors_legacy (true iff dfloor_keep_velocity, vceil, eos_floor_consistent, efloor_from_ekin all off) running the OLD kernels verbatim from separate .cpp files"
metadata:
  type: project
---
Read docs/handover/NOTE-2026-09-11-viper-floors-legacy-gate.md. Under floors_legacy the pre-8da093f5
kernels run: GnomonicEquiangleRaiseVel legacy par_for, SingleC2P_*Legacy inversions, the old 3-reducer
ConsToPrim in general_hyd.cpp/ideal_hyd.cpp; the switch-aware launches moved to NEW files
src/eos/{general,ideal}_hyd_floors.cpp (ConsToPrimFloors) because a legacy kernel in the SAME
translation unit as the new one still differed on HIP. prolong_prims.cpp calls the Legacy inversions
under the flag. No default changed. Also on origin: HIP build fix (cdd7d2a5), implicit radial solve
split into per-cell/per-face kernels (77de5618, 1 ulp on GPU), EintFromCons subtracts magnetic
energy + rad_implicit_x1 for <mhd> (bc972cc8), cs seam geometry table (27ca5b13), dhj rt_use_cons.
**Why:** any GPU A/B across a10e367d..caad9247 is round-off seeded; production with the new floor
switches ON bypasses the legacy path and needs its own A/B. The red giant (I8, RG_fofc) runs with
switches ON, so it is on the NEW path; dhj is on the legacy path.
**How to apply:** when merging origin into the FOFC tree, the c2p vceil_test flagging (c0163568)
and the MHD vceil (step 4) must go into the *_floors.cpp / new kernels, and the FOFC test pass must
work on the legacy path too (fofc on + all switches off = legacy). Related: [[viper-hip-code-conventions]], [[fofc-compatibility-plan]].
