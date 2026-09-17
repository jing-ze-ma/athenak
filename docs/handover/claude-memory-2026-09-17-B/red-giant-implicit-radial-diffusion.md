---
name: red-giant-implicit-radial-diffusion
description: "IMPLICIT RADIAL RADIATIVE DIFFUSION IMPLEMENTED 09-10 21:00 (hydro/rad_implicit_x1, default off; per-column tridiagonal after the RK update, faces is/ie+1 stay explicit, conservation 3e-16, +4.4% cost, flag-off bit-identical): I1 (V9f from rst 5.0e5) PASSES the 5.7125e5 death (dt flat 30.6, tot-E unchanged) but the SAME evacuated-cell runaway reappears through the explicit ANGULAR faces at 5.7604e5 (m=85 i=342-344, rho 1/62 of the shell, T 9.6e4 -> 3e10). Next: conservative per-face heat-capacity cap rad_cap_ang=0.5 on x2/x3 (agent running, I2_implicit_cap)"
metadata:
  type: project
---
Code (uncommitted, working tree): conduction.hpp 93-115 (flag, imp_wrk 7-slot scratch),
conduction.cpp ctor 152-182 (fatal unless <hydro>, radiative, meshblock nx1 == mesh nx1),
RadFaceKappa 297-318 shared by explicit and implicit, explicit x1 kernel skips interior faces
when on, ImplicitRadialUpdate 559-810, NewTimeStep drops dt1; hydro_tasks.cpp task impcnd
between srctrms and sendu_oa (beta_dt = beta[stage-1]*dt); mesh.cpp collapse report says
"x1 implicit". Solver: e* from u0 (cs metric KE, etotgrav term), T* = Temperature(d,e*,wtemp),
alpha = 1/(rho cv), C_f = w_tau kappa_rad L(s*) A_f/dl_f frozen, Thomas per (m,k,j) column,
only u0(IEN) written. Reviewed by me: coefficients, units, telescoping all correct.
Gates: G2 flag off bit-identical to build_guard (data payload; only the embedded input text
differs); G3 vs explicit at 2.5e5: |dT/T| 1e-7 deep, 3e-4 at the photosphere, L_out 5e-7,
conservation 3e-16, cost 0.589 -> 0.615 s/cycle. Binary build_impl (athena_pin).
I1_implicit_x1 job 195813: dt 30.6345 through 5.8134e5, tot-E 9.8844e49, old culprit column
(gid 90 i=315) healthy (rho/median 0.18-1.5, smooth T, no checkerboard); new death at
5.760448e5 via dt2=1.4e-3 (angular explicit), culprit m=85 k=2 j=4 i=342-344 at the
photosphere, rho 3e-11 vs 1.9e-9 median, T 9.6e4 vs 3.4e3 -> 3e10 K. Cancelled at 5.84e5.
Why the cap and not implicit angular: angular faces cross MeshBlocks/ranks. The cap
f_f = min(1, 0.5/max(x_i,x_j)), x_i = dt/dt_cond,ang(i), is conservative, local, and makes
the angular operator unconditionally stable; a SYMMETRIC predicted-temperature flux is
UNSTABLE for x>1 (checkerboard amplification 1+4x(x-1)/(1+x)) so do not try that.
The physical question left open: WHY cells get evacuated 60x (grid-locked cube-edge/vertex
downdrafts, see [[red-giant-vertex-chimney]]). See [[red-giant-common-dt-collapse-571e5]].

## UPDATE 09-10 23:00: angular cap DONE, I2 REACHED tlim 9.0e5
rad_cap_ang (conduction.hpp 116-146, conduction.cpp 178-204 ctor, face_kcode 425-445,
radcapc2/c3/x/diag 554-708, x2/x3 kernels 715-758, report 763-785, NewTimeStep drops dt2/dt3;
hydro_tasks.cpp 405-409 sets stage_beta_dt). K_f includes 1/sin(alpha) on the cs (stricter).
Gates: cap unset -> byte-identical to implicit-only; both off -> byte-identical to build_guard;
cap=0.5 on the smooth 2.5e5 star: 0 capped, max x_i 7.4e-6 (angular dt_cond ~4e6 s!), byte-
identical. Cost +20.8% over explicit total (cap adds 13.7%: three extra kernels per stage).
I2 (job 195846): 5.0e5 -> 9.0e5 in 4937 s wall, dt 30.634 flat, ONE hydro dip (6.1 s at
5.83e5), tot-E 9.88440e49 flat. Capping only between 5.687e5 and 5.976e5 (401 calls, max 18
cells/62 faces, peak x_i 1.35e6 at the photosphere r 3.39-3.47e12, gid 85 first), then silent
for 3e5 s: the hot evacuated population (T up to 28x the shell median at 5.8e5) dissipated by
7.0e5 (max 3.3x, 0 cells >10x). Both fixes are production-ready pending the commit.
Binary: build_impl/athena_pin2 = I2_implicit_cap/bin/athena.
