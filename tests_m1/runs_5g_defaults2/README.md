# runs_5g_defaults2: hesdirk2 (+ implicit_vimp) and vet_col_surface_q become defaults

- **Branch.** `m1-defaults2` from rt-integration 1560d07a; worktree `/viper/ptmp2/jinma/wt_m1def2`.
- **Build and run tree.** `/viper/ptmp2/jinma/m1def2_0924`: `bin/` (base = 1560d07a,
  new = fe147c9f, the code of this branch; later commits change comments, tests and docs
  only), `cpu/` CPU gates, `gpug/` GPU gate (job 11965080), `gpuf/` fallback study (job
  11965081), `gates/` gates.py runs, `inp/` the inputs with keys added (`scripts/addkeys.py`).
- **Binaries (md5).** CPU box_convection base 6f3ed02d, new f8302bd3; CPU no-pgen base
  4cc333e3, new b004f515; GPU box_convection base 7b7cd4a3, new ffb9a91d. gcc 14 +
  openmpi 5, Release, MPI; GPU rocm 6.3, gfx942.
- **GPU environment.** apudev, `HSA_XNACK=1`, `HSA_NO_SCRATCH_RECLAIM=1`.

## 1. The changes

**`<rad_m1>/time_scheme`** (resolved at the top of `ImplicitInit`, rad_m1_implicit.cpp).
- The default is `hesdirk2` wherever hesdirk2 is accepted:
  - transport = implicit (not implicit_x1);
  - a Cartesian mesh: no spherical-polar (S1/S2), no cubed sphere, no polar boundary;
  - a closure other than tau / vet_col;
  - `<time>/integrator = rk2`.
- Everywhere else the default stays `be`, silently. No input that runs today becomes a fatal.
- A restart whose file lacks the key keeps `be`. The key was read only when named, so old
  files carry it only if their input named it.
- The resolved value is always echoed into the parameters, so it is written to the restart.
  Explicit input overrides.
- Evidence for the change:
  - `/viper/ptmp2/jinma/h2val_0924/README.md`: hesdirk2 at cfl 0.9 is 1.57x cheaper than
    be at cfl 0.3 at equal accuracy;
  - `tests_m1/runs_5f_h2fast/README.md`: 1.35x be per step.

**`<rad_m1>/implicit_vimp`.**
- It defaults on only where the resolved time_scheme is hesdirk2. It was validated there;
  under be it over-damps the P = 100 waves.
- It also needs every one of its own preconditions (otherwise it stays off, silently):
  - a multi-D mesh with implicit_solver = bicgstab;
  - nghost >= 2;
  - hydro with coupling, gas_feedback and dbg_gas_force.
- A restart whose file lacks the key keeps it off. The value is echoed; explicit input
  overrides.
- An input that names `time_scheme = hesdirk2` without `implicit_vimp` now gets vimp on.
- The accel levers (`implicit_vimp_fold` etc.) resolve after vimp, as before.

**`<rad_m1>/vet_col_surface_q`** (rad_m1.cpp).
- The default is true (`tests_m1/runs_5e_vetcol2`: T-S4 L1 2.1e-3 -> 5.2e-5 at n = 256;
  Milne 1.8e-3 -> 1.1e-4).
- A restart whose file lacks the key keeps false. Files written since 5e carry the key
  anyway, because it was always echoed.

**Tests and docs.**
- `tests_m1/gates/box3d_be_x.athinput` now names `time_scheme = be`, so the gates still test
  the be path. `slab2d_plm_vimp_be_x` already named it.
- `tst/test_suite/rad_m1/he_slab_m1.athinput` already named `time_scheme = be`. It now has a
  comment saying so. The REF numbers are unchanged; the tests stay be regressions.
- `docs/dev/switch_inventory_2026-09-24.md`:
  - the time_scheme and implicit_vimp rows now read DEFAULT-ON;
  - a new vet_col_surface_q row, also DEFAULT-ON;
  - decision 2 is updated.

## 2. Gates (all PASS)

| gate | cases | result | file |
|---|---|---|---|
| be named = bitwise base, CPU | box 1 rank; box 2 ranks overlap + faces; slab 1; slab vet_sc 2; slab2d_nd + be 2 ranks; box3d_nd + be 2; beam inputs cart_sym, marshak_cart, rw_cart, pp_np + be | 10/10 BITWISE | `RESULTS_cpu_gates.txt` (A_*) |
| be named = bitwise base, GPU (job 11965080, 2 GPUs) | slab2d_old + be 1 rank; box3d_nd + be 2 ranks; the same + vet_sc | 3/3 BITWISE (10 files each) | `RESULTS_gpu_gate.txt` |
| key absent (new) = base with time_scheme = hesdirk2 + implicit_vimp = true named, CPU | slab2d_nd 2 ranks; box3d_nd 2 ranks; box3d_nd vet_sc 2 ranks; cart_sym; marshak_cart; rw_cart (vimp named too, since it has hydro + coupling) | 6/6 BITWISE | `RESULTS_cpu_gates.txt` (B_*) |
| the same, GPU | box3d_nd 2 ranks | BITWISE | `RESULTS_gpu_gate.txt` (boxh2) |
| sp input without the key stays be | pp_np (sp): new vs base, no key | BITWISE, and no fatal. Base with hesdirk2 named is refused by the sp check, as expected | B_ppnd, B_pp_np |
| sp S1/S2/S5 + Cartesian vet_col, `vet_col_surface_q = false` named | the 15 runs_5e gate_sp cases | 15/15 BITWISE | `RESULTS_cpu_gates.txt` (C_*) |
| new vet_col_surface_q default | sph_atm_vc without the key (new) vs `= true` named (base) | BITWISE | C_sqdef |
| old restart without the key | base slab2d_nd (no key, be) to t = 100; continued to 200 by base and by new | new = base continuation, BITWISE; new = the base straight run (rst + 100/100 hst rows) | D_old |
| hesdirk2-default restart | new slab2d_nd 2 ranks, straight to 200 vs restart at 100 | BITWISE (rst, 100/100 hst rows) | D_h2 |
| `tests_m1/gates/gates.py` run + eval (new box binary) | 14 pairs | GATES: PASS | `RESULTS_gates_py.txt` |
| `tst/test_suite/rad_m1` (ATHENAK_M1_DATA = faces_0924/m1data) | 3 tests | 3 passed (333 s) | `RESULTS_pytest.txt` |

## 3. The gas-Newton bracketed-fallback flag of runs_5f

Job 11965081: 1 GPU, the runs_5f 3-D He box (84x104x104, 4 blocks, plm + vimp), 120
cycles, new binary with the hesdirk2 default. The input names the lever keys at their
defaults so that they can be overridden (`inp/box3d_fb.athinput`). Each arm is compared by
`tests_m1/gates/cmp.py` (`scripts/fbcmp.py`) with T, a tight reference: implicit_tol 1e-11,
implicit_lin_tol 1e-11, time2_lin_tol_fac 1. NON-CONVERGED = 0 in every arm.

| arm | Eddington fallbacks | vs T: hst_cons / rst_cons | vet_sc fallbacks | vs T: hst_cons / rst_cons |
|---|---|---|---|---|
| D (defaults: lin_tol_fac 10, one_pass_safety 30, one_pass 8) | 44 | 1.2e-9 / 1.1e-9 | 15 | 8.4e-10 / 1.7e-9 |
| L1 (time2_lin_tol_fac = 1) | 49 | 4.7e-10 / 6.4e-10 | 29 | 9.2e-11 / 8.1e-10 |
| S0 (time2_one_pass_safety = 0, i.e. 3) | 26 | 8.6e-10 / 1.1e-9 | 13 | 7.8e-10 / 2.1e-9 |
| L1 + S0 | 26 | 4.8e-10 / 5.7e-10 | | |
| OP0 (implicit_one_pass = 0) | 121 | 7.8e-10 / 1.3e-9 | | |
| T (tight) | 143 | 0 | 40 | 0 |

(About 2.7e8 cell-passes per arm. Full lines, including hst_dyn and rst_dyn, are in
`RESULTS_fallback.txt`.)

**Findings.**
- **time2_lin_tol_fac is not the cause.** With factor 1 there are more fallbacks, not fewer:
  49 vs 44 Eddington, 29 vs 15 vet_sc.
- **The rise from 26 to 44 (Eddington) comes from time2_one_pass_safety = 30.** For vet_sc
  the same lever changes the count only from 13 to 15.
- **The count grows as the solve gets tighter.** The tight reference has 143 / 40, and
  one_pass off has 121. The safeguard fires on cells whose gas residual has reached
  round-off, and a tighter solve (more passes) puts more cells there.
- **So a fallback is not a sign of an inaccurate solve.** The default arm is 1.2e-9
  (hst_cons) from the tight reference. D vs L1 is the same size (1.2e-9 hst_cons, 1.1e-9
  rst_cons). Every arm is within 2e-9 on the conserved quantities.
- **Decision: no default lowered.** No lever makes the result measurably worse. The count
  (1.6e-7 per cell-pass) is a round-off artefact of the safeguard's residual test, not an
  accuracy problem.

## Files

- `scripts/`:
  - `build.sh`;
  - `gate_cpu.sh` (A/B/C/D), `gate_gpu.sh`, `fb_gpu.sh`;
  - `gcmp.py` (bitwise), `rstcmp.py` (restart), `fbcmp.py`, `addkeys.py`;
  - `pytest_m1.sh`.
- `RESULTS_*.txt`: the tables above.
