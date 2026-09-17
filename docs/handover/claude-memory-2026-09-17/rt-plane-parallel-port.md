---
name: rt-plane-parallel-port
description: 2026-09-13 the red-giant grey two-stream RT (src/utils/two_stream_rt.hpp) got a PLANE-PARALLEL mode (namespace flags two_stream_rt::rt_plane_parallel, rt_top_vacuum; X1V/X1F/DX1 helpers rebuild geometry from RegionSize; rt_pre_geom early return) wired into box_convection.cpp behind problem/rt_two_stream (cool_layer defaults to !rt_two_stream; one meshblock in x1, n1<=520). Commit 5f2be7d2 on rg-box-convection in bench/wt_rgbox (NOT pushed). Spherical path BITWISE unchanged (dhj ck + cs red_giant whole-run diffs byte-identical; 4/4 rad tests). Grey Eddington atmosphere - T error 3.5 % (256 cells, 2-pt quadrature), 2.1 % (512) = the two-stream's own angular error. HIP vs CPU 4e-15 (round-off) on the full 4e4 s run with the FIXED binary (job 11685737, 15-min apudev; a 1-D GPU run is latency-bound ~107 cycles/s)
metadata:
  type: project
---
Validation case inputs/hydro/box_convection_grey_atm.athinput + kappa_const.txt; run dir bench/rg_box/rt_check (compare.py /
compare9.py). Latent: the explicit radiative-diffusion dt limit ignores the tau-blend weight -> set hydro/rad_implicit_x1=true or
the step drops 13x. Repo style check fails on 2204 PRE-EXISTING cpplint violations on this branch (none from this change), so it
is not a usable gate. Used by [[fecz-box-projects]].

OPEN BUG 09-13 (~11:00): the plane-parallel two-stream + tau blend (rad_tau_lo/hi 10/100) does NOT conserve energy. B-star RT smoke
(bench/bstar_fecz/smoke_rt, job 11681669, 0.4 turnover): diffusion carries 1.005 F in the thick interior, emergent flux 1.09 F at t=0
(top correct), yet the box retains 25.6 % of F A t, 95 % as aT^4, HALF deposited continuously (12.8 % of F A) in the tau 10-100
handover -> 1.5e5 K isothermal corona (top cell rho x875, tau 3.9), dt 1.49 -> 0.63 s (sound speed), FeCZ heated 12 % and made
SUBADIABATIC (-3.8e-2). Whole box rings coherently (|<v1>|/rms 0.99), no convection (v_h 0.02 v_MLT). The spherical red giant
may share the defect undetected (13 % of L over 3e6 s = 7e-9 of its E). Diagnosis in progress (RT-port agent). Do NOT run
rt2turn until fixed. Analysis: smoke_rt/analysis/rt_profiles.py (adds tau, Eddington law, emergent flux, region energy budget).
