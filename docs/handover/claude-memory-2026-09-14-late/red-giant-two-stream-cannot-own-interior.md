---
name: red-giant-two-stream-cannot-own-interior
description: "ANSWERED 09-10 19:00: the semi-implicit two-stream CANNOT own the optically thick interior (deposition throttled to src*t_rad, ~1e-11 of physical below tau~1e3; core L enters only via the conduction wall term; non-conservative). The 1.058 R killer cell is thin per cell (tau_cell~1) inside a tau>100 column, so the blend hands it to the diffusion operator whose kappa_rad is invalid there; cures = implicit radial diffusion (tridiagonal per column, whole radius in ONE MeshBlock so no MPI) or blend on local tau_cell"
metadata:
  type: project
---
Code facts (two_stream_rt.hpp anchors as of 09-10): icut = first face from the bottom with w<1
(rt_pre_cut ~988-1004); taublend forces int_at_cut=false (859) so no interior flux enters the
sweep; core L is the conduction wall term at conduction.cpp:381-387 (unweighted by w); the
semi-implicit relaxation de = deq(1-exp(-x)), x = src*bdt/deq, deq = src*t_rad, so de -> src*t_rad
when dt >> t_rad (1846-1870); de is applied straight to u0 (1999), no flux adjustment;
rt_chain* is a frequency chain (RT_NB=4), not spatial; sweep radial tier NN 72/136/264/520
(1322-1331), nx1=480 fits in one block. Max column tau at the wall 2.4e12. Restart path is clean
(BuildRadWeights every call). Large-tau layer coefficients are cancellation-free (1268-1272) but
the direct source (1274-1278) is a difference of nearly equal Planck functions at dtau~1e9.
Deeper-handover option that IS legitimate: rad_tau_lo=100, rad_tau_hi=1000 (still diffusion deep).
Killer-cell numbers: kappa*rho 7e-10, dx 2e9 -> tau_cell ~1.4, t_rad ~2e4 s (semi-implicit
valid there), while dt_cond used the unlimited kappa_rad. See [[red-giant-common-dt-collapse-571e5]].
