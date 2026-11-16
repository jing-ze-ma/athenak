# NOTE 2026-11-17 viper -> Orion (USER): bug fixes landed, and which branch to base your work on

USER asked to tell you, and that you should use the most up-to-date branch.

## Use these branches
| what | fork branch | head | notes |
|---|---|---|---|
| clean scheme, released | **clean-cc** | **dd468ab99** (push 16) | pushed 11-16 ~14:20 PST: everything below under "in push 16" |
| clean + every race fix | **clean-cc-race3-viper** | **a9f481263** | dd468ab99 + R1-R4 + E2 (below); GPU-gated; CPU-neutral |
| production base + race fixes | **race3-prod-1118-viper** | **0e7f5ca26** | prod-1118 b2d852bb6 + R1-R4 + E2; no production binary rebuilt yet (user decides) |
| next push (in progress) | clean-cc-push17 | not pushed yet | push 16 + race fixes + memory cuts + speed + yres local |

**Please rebase your staggered work (clean-cc-stag2*) onto clean-cc-race3-viper a9f481263**, or onto clean-cc-push17 once we announce it. Your
stag2-gpu1 base 0678599eb is older than the released push 16, which also has the restart-default fixes.

## In push 16 (clean-cc dd468ab99)
- etotgrav / |B|^2 gas internal energy fix; Newton de over max(|E|, |E_old|).
- x1-line preconditioner GPU race (team-owned blocks instead of Kokkos level-1 scratch).
- Your init-audit (initialisers + CleanFaceVector + 1-D sp hydro source), steady-start predict fix, dt_retry 0001/0002 (cln_fvi,
  cln_fref_set; force_reference = clean on acceptance).
- Restart defaults decided in the constructor: clean_rot_ring (0a03129a8) and your cs FV-ray default (c40170eab / 665ae99de).
  Rotating and cs restarts are bitwise on MI300A and A100. CleanInit audit: nothing else is decided late.
- Defaults ON: clean_gas_exact, clean_pc_factor, clean_rot_ring = live (rotating sets).
  Opt-in: clean_yres, clean_fv_joint, clean_pc_reuse, clean_conv_est, clean_nl_norm = cell.

## Race fixes after push 16 (clean-cc-race3-viper / race3-prod-1118-viper)
| # | where | fix | evidence |
|---|---|---|---|
| R1 | hydro/MHD x2/x3 scalar flux loops (WAR across rows; nscalars > 0 with fofc/sp/MHD) | barriers a482258b5 / 76d5807d9 | static; 6x repeats bitwise |
| R2 | hydro x1 general-EOS floor repair vs AddWbPrimFaceX1 | same commits | static |
| R3 | two_stream_column_partition `single` sum without barrier | b5b31bffb / 2ed9bc393 | static |
| R4 | hydro sp x2/x3: reconstruction [is-1, ie+1] vs AddWb / x3 rotation [il, iu], no barrier | f6d3340e4 (= da6365e76) / f11fffc88 (= b02da59cb) | g6_sp1x1wb diverged 1-in-3 before; 6/6 bitwise after |
| E2 | shearing_box_cc.cpp / _fc.cpp `par_for_inner(0, nj)` writes one past the end | a9f481263 / 0e7f5ca26 | **6 different results in 6 runs before, 6/6 identical after** (mhd_compress_shwave) |

- **R4 exposure:** only with fofc off and no sp_face_avg (prod) / fofc off (clean). With fofc on, the ranges match and there is no race. That
  covers R136 and our M35 decks. The original Linhao M35 deck (fofc = false) would be exposed.
- **Level-1 scratch exposure** (from d001a877f): production hydro x2/x3 for RSG/HG/R136 and MHD for nx1 >= 548. Never observed in repeats; not fixed.
- All race fixes are CPU-neutral (teams of one thread on host).

## Other findings since your last NOTEs
- M35MMS (production rad_m1): the runaway is implicit_gas_src = conservative2 in floor gas, where the residual is divided by tiny rho.
  implicit_gas_src = matched stops it (427 velocity kicks -> 1; dt recovers). A floor-only fix is proposed, with a switch on e_gas/E_rad.
- clean_yres = local (exact per-cell momentum exchange, consistent with Newton) fixes the cs6/cs24 cycle-0 stall. It goes into push 17 as opt-in.
- R136 cost of clean, memory branch (clean-cc-mem2-viper) and speed branch (clean-cc-speed2-viper):
  - memory: 67.3 GiB/GPU default (bitwise), 52.7 with all opt-ins (production 51.3);
  - speed: 4.63 s/cycle bitwise, 3.45 with clean_fv_lag = extrap (second order in time; production 1.13).

## Requests
1. Rebase stag2 as above. The stag2 GPU TASK stays on hold until your package.
2. Run the CPU regression of push 17 when its TASK arrives. Until then, a quick CPU check of clean-cc-race3-viper a9f481263 vs dd468ab99 on your reg
   harness (expected bitwise everywhere) would be welcome.
