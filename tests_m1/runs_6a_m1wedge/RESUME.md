# RESUME: m1-wedge (sprhd_0926), gravity-bearing sp rad-hydro wedge with implicit M1 + vet_col

Brief (coordinator 09-26, redirect): follow docs/handover/TASK-2026-09-26-m1-real-wedge.md (rt-integration
206318c9): a gravity-bearing He-envelope sp wedge around the FeCZ with box_convection physics/tables; route
(a) red_giant M1 mode or (b) the sph_atm test pgen gains gravity + ic_profile + BCs; FIRST gate = static
grey hydrostatic + radiative balance with a point mass (test A); (B)/(C) optional; TASK "Gates" incl. GPU
lever timings and existing problems bitwise; CUDA-safe code; branch m1-wedge, doc docs/dev/m1_wedge_0926.md;
no merge into rt-integration, no push of rt-integration.

## Choice: route (b), as a NEW file
- src/pgen/tests/rad_m1_wedge.cpp: `m1_test = sph_wedge` (pgen_name rad_m1_beam, the h_none binary).
  Dispatched from RadiationM1Tests2 (rad_m1_tests2.cpp), declared in pgen.hpp, listed in src/CMakeLists.txt.
  Every other m1_test value is untouched (bitwise by construction; gated below).
- Why not red_giant: 4900 lines built around two-stream/conduction/MLT with active defaults (sponge on,
  required opac_table/teff/ptop, its own IC march); an M1 mode must neutralise dozens of switches.
  The test pgen already has every M1 hook; what it lacked (point-mass WB gravity = RedGiantGravity's WB
  branch, a column reader, the two-table reader of box_convection, closed x1 walls) is ~600 lines.

## Worktree / branch
- worktree /viper/ptmp2/jinma/sprhd_0926/wt, branch m1-wedge (from fork/rt-integration c1dc8226 = BASE).
- builds: scripts/build.sh <cpu|gpu> <tag> <rev> [problem] -> git archive snapshot in src/<tag>,
  kokkos -> /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos (4.6.2; the athenak submodule 4.4.00 has no
  GFX942_APU), binaries bin/athena_<tag>_<problem>_<cpu|gpu>.
  base = c1dc8226, d1 = 5dd4f3d9, d2 = 969ded08 (tables + seed).

## Inputs
- inp/gateA.athinput: grey point-mass atmosphere (ideal gas, opacity const kappa_t 1e4 (half absorption),
  c 300, a 3e5, GM 1, r 1..1.25 (96 cells), F_in 0.009 -> Gamma 0.3, tau(r_in) 4.4e3, rho_top 3e-5),
  wedge theta pi/2 +- 0.1, phi 0..0.2 (16x16), WB polytropic + phi_eff + force_reference wb_arad.
- inp/he.athinput: the He-box FeCZ star (V3edd) on a sphere: M 3.15 Msun, R0 = sqrt(GM/g0) = 3.2405e10,
  r = R0 + z for the V3edd box z range (-1.8149e8 .. 4.129e7, 84 cells, tau(r_in) 99), 2 deg x 2 deg wedge
  (= the box width, 32x32), table EOS gas only, Rosseland + Planck tables (HE_BOX_DATA), vet_col,
  hesdirk2, mg_gc, WB polytropic + phi_eff + wb_arad.
- IC file he/ic_wedge_he.txt from wt/tests_m1/runs_6a_m1wedge/ic/build_ic_sph.py --nsub 16 --pad 20
  (spherical grey RE + HSE column, Eddington closure, the box EOS/tables imported from
  bench/m1_stage2/ic/build_ic.py). Gamma max 0.751 at z = -9.23e7 (box: 0.7507).

## Runs (see sections below for results)
### CPU (login node, gcc/14 openmpi/5.0, 4 ranks unless noted)
- cpu/A1 (scripts/cpuA.sh d1 A1 6.0): gate (A) 96x16x16, t = 6 (~3.5 sound crossings, ~8 deep
  diffusion times), arms r3 (cfl 0.3), h9 (0.9), v4 (vet_col_every 4 + lag), h9v4, nowb (no WB,
  force_reference none, phi_eff off).  scripts/hst.py <dir> -> L(r)/L_in, interior Mach_rms, dE.
  RESULT: all arms Mach_rms (tau >= 1) 1.3e-4 -> 5e-5 at t = 6 (decaying ringing), L_bot/L_in 1.00039,
  L_top/L_in 0.9990-1.0002 oscillating; cfl 0.9 / every-4 differ from r3 by < 3 % of these numbers;
  nowb Mach_rms 3.7e-5 (the WB pair is not needed for the grey profile); NON-CONVERGED 0 everywhere.
  scripts/binv.py: interior max Mach <= 5e-3 (at r ~ 1.19-1.20, tau ~ 1-3), top max 1.3e-2;
  dt limited at r_in (radial cell / c_s), not at the top.
- cpu/Ares_d2 (scripts/cpuAres.sh d2): nx1 48/96/192 on 4x4 lateral, t = 3: Mach_rms(t=3) 4.0e-4 /
  1.2e-4 / 5.7e-5; L_bot-1 1.3e-3 / 3.9e-4 / 1.05e-4, L_bot offset 1.3e-3 / 3.9e-4 -> ~second order.
- cpu/gate_d2 (scripts/cpugate.sh d2): restart bitwise He (20+20 vs 40) TRUE, grey TRUE; T-S4 (vetevery
  wedge input, 32x32, 12 cyc) base vs d2 all bin/rst files identical; He 2 vs 4 ranks max rel hst diff
  2.2e-9 (tiny KE columns), <= 3e-12 elsewhere.  he_long: 300 cycles (t ~ 48 s).

### GPU (apudev, 2 GPUs unless noted; logs gpu/log/<job>.<id>.out)
- g1_d2 (11984593): He restart bitwise TRUE; 1 vs 2 GPU hst <= 4e-11 (KE 4.6e-9); long relaxation
  34 500 cyc, t 5 550 s: L_top/L_in 1 +- 2e-6, radial ringing Mach_rms 9e-4 -> 1.4e-3 (slow growth),
  lateral decays, no convection from the grid-scale seed.
- g2_d2 (11984595): lever timings at 84x64x64 from the IC (table in the doc); ref 18.1 ms/cycle.
- g3_d2 (11984594): test (C) hot spot; python3 scripts/anac.py gpu/g3_d2.
- gb (11984680/2/3): T-S4 / He box / dhj base vs n1 (b899f03c): all BITWISE, restarts bitwise;
  python3 scripts/gbcompare.py wedge box dhj.
- g4_n2 (11984805, apu 1 h, SUBMITTED, not analysed): He wedge 84x64x64, seed 1e-2 k = 4 over tau
  5-50.  Analysis: python3 scripts/hst.py gpu/g4_n2/run 20 1.6667  (Mach_rms and Mach_r columns:
  lateral KE = KE_int - KEr_int shows convection); for f in gpu/g4_n2/run/bin/*hydro_w*.bin; do
  python3 scripts/binv.py $f 3.23971e10 1.6667; done.  Then re-time the levers from its final rst
  (edit g2.sub: -r <rst> instead of -i, same arms).

## Branch state
- m1-wedge commits: 5dd4f3d9 969ded08 b899f03c f323a837 3b56c19d, merge of fork/rt-integration
  4979b3a1 = 38b2e2f3 (clean; merged CPU build + 20-cycle gate A smoke hst identical to d1).
- NOT pushed, NOT merged into rt-integration. Doc: wt/docs/dev/m1_wedge_0926.md. Inputs and scripts
  mirrored in wt/tests_m1/runs_6a_m1wedge/.
- On another cluster: build with -D PROBLEM unset (pgen_name rad_m1_beam, m1_test sph_wedge);
  rebuild the He IC with tests_m1/runs_6a_m1wedge/ic/build_ic_sph.py (needs bench/m1_stage2/ic/
  build_ic.py and the He tables; HE_BOX_DATA on Caltech), fix the absolute paths in he.athinput.
