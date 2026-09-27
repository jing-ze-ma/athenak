---
name: box-restart-bug-fixed
description: 2026-09-13 box_convection RESTART corruption FIXED 929d752d (rg-box-convection, wt_rgbox, unpushed) - (1) gravitational potential phicc0/phi0 was filled inside the IC kernel AFTER the restart early return -> restarted runs had phi=0 (etotgrav made e_int wrong by rho g z; WB stencil read zeros); (2) BoxConvBC read w0 before the first c2p (the red-giant 3fd3a836 defect) -> walls came back as floor ghosts. Post-fix dens/ener/mom2/mom3 BITWISE for 100 cycles; mom1 residual 1e-7..1e-5 relative (~1e-3 cm/s) in ~50 cells, non-growing, OPEN. Any multi-link box chain needs a binary >= 929d752d. Test dir bench/bstar_fecz/rsttest.
metadata:
  type: project
---
Symptom: 1-D relax1d restart -> dt 1.49->0.60, rt_de_max clipping in 17 cells at cycle 0, |v1| 5e6. red_giant.cpp and
deep_hot_jupiter_rt.cpp on polar-average-perf share neither mechanism (checked). See [[restart-bc-reads-w0-bug]],
[[fecz-column-relaxation]].
