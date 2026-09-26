---
name: m1-hesdirk2-cfl-recommendation
description: "(fast5-sp wedge = grey no-gravity test T-S4, NOT He) DECIDED 09-26: next M1 run uses hesdirk2 at cfl 0.9 (be@0.3 accuracy, 1.57x cheaper); hesdirk2 now also runs on sp + vet_col (m1-sph2), but cfl 0.9 was validated only on the Cartesian vet_sc box"
metadata:
  node_type: memory
  type: project
  originSessionId: 05818bc0-259d-4b2f-b317-f82d1e083443
  modified: 2026-09-24T17:32:25.024Z
---

Raise to the user when setting up ANY new simulation with implicit M1 / VET (He boxes, He slab, sp wedge once
hesdirk2 is supported there): choose `time/cfl_number` for time_scheme = hesdirk2 (the default on Cartesian implicit
M1 since m1-defaults2, 09-24; be stays default on sp / vet_col / cs, where hesdirk2 is refused).

Options (validation /viper/ptmp2/jinma/h2val_0924/README.md, 3-D He box vet_sc, GPU, per simulated second):
- cfl 0.6: 2x MORE accurate than be at cfl 0.3, 1.19x cheaper (0.232 vs 0.277 s/s).
- cfl 0.9: same accuracy as be at 0.3 (within 10 % in rho/eint/E), 1.57x cheaper (0.177 s/s).
- be at cfl 0.9 is NOT an option: He slab blows up (KE x150-3000); hesdirk2 bounded to cfl 1.2.
- Numbers predate m1-h2fast (hesdirk2 per step 1.85x -> 1.35x be), so the hesdirk2 margins are now larger.

**Why:** the user asked (09-24) to remember this so the CFL is decided consciously at the next setup, not left at 0.3.

**How to apply:** at setup, present the two options above and let the user pick; if dt is capped by something else
(hydro CFL, a source term, outputs), hesdirk2 buys accuracy only and costs 1.35x per step -- say so. Related:
[[rad-m1-design]].

**DECISION (user, 09-26):** the next M1 run uses time_scheme = hesdirk2 with cfl_number = 0.9. Correction: sp + vet_col supports hesdirk2 since m1-sph2 (src/rad_m1/rad_m1_sph.cpp:95, tests_m1/runs_5h_sph2). The 0.9 accuracy numbers are from the Cartesian vet_sc box only, so on a wedge run check the first rotations (dt limiter, non-converged count, interior profiles vs a short cfl 0.3 reference).

**Wedge trial 09-26** (/viper/ptmp2/jinma/m1wedge_0926, job 11981264): the fast5-sp "He wedge" is a TIMING wedge that holds the gas fixed (atm_hold, v = 0), and it starts near equilibrium. No physically evolved wedge restart exists; all hw.*.rst files are from runs of <= 60 cycles. The reference there is already hesdirk2 cfl 0.3 (the resolved sp default).
- H9 (cfl 0.9) vs R (cfl 0.3): 2.75x cheaper per sim-s, 0 nonconv, interior (tau >= 1) deviations <= 7e-6. dt is set by hydro sound crossing at r 1.887.
- Hydro at cfl 0.9 is UNTESTED (the gas never moves). Next: an atm_hold = false check of cfl 0.9 vs 0.3.
- The vet_col_every result (vetevery_0926: N=4 -23 %) also comes from this static wedge.
- CORRECTION 09-26 (m1wedge_mv_0926/RESUME.md): the "fast5-sp He wedge" is NOT He. It is pgen rad_m1_beam, m1_test = sph_atm (test T-S4, rad_m1_tests2.cpp): a grey radiative-equilibrium atmosphere with rho ~ r^-2, constant kappa and NO gravity. With atm_hold = false it just expands. No gravity-bearing M1 He wedge exists (the He4 presn global model is still a plan). The only moving-gas M1 setup is the Cartesian He box (box_convection).
- **He box moving-gas check 09-26 (hebox_cfl_0926), partial:**
  - cfl 0.9 is 2.28x cheaper than 0.3; vet_sc_every 2 adds +7 %; dt is set at the bottom.
  - BUT KE at cfl 0.9 saturates ~4.6x the cfl 0.3 level. The box barely convects (F_conv/F 1e-5).
  - cfl 0.15 / 0.6 arms were running at handover. Until settled, cfl 0.9 = relaxation-phase only.
