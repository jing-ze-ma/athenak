---
name: fixes-in-hydro-and-mhd
description: USER 09-27 - every M1/coupling fix must be done for hydro AND MHD, Cartesian AND spherical-polar; note rad_m1 currently couples to hydro only (no pmhd anywhere in src/rad_m1)
metadata:
  type: feedback
---
User 09-27: "remember to do all those fixes not only in hydro but also in mhd, both cartesian and sp."
**Fact (checked 09-27):** rad_m1 reads only pmy_pack->phydro (rad_m1_coupling.cpp:88 etc.); `pmhd` appears nowhere in src/rad_m1/. So M1 + MHD is not supported yet; the fixes cannot be applied to MHD until an M1-MHD coupling exists.
**How to apply:** every M1 fix (force_reference_work split, converged opacities, seam, sponge ...) gets a hydro + MHD x Cartesian + sp test matrix. Where MHD is not supported, say so and plan the coupling rather than silently skipping it.
