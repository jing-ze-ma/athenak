# VET radiation test suite

A self-contained CPU test suite for the VET radiation scheme: one command, a PASS / FAIL
line per test, about 1 min (`--quick`) or 5 min (full) on 48 cores (`--jobs 12`).

```bash
# build a CPU binary (into build_vettest/) and run the quick subset
tests/vet_clean/run_tests.sh --quick
# or test a binary you built yourself
tests/vet_clean/run_tests.sh --exe build_cpu/src/athena            # full suite
tests/vet_clean/run_tests.sh --exe build_cpu/src/athena --only restart_sp,dt_retry
```

Requirements: python3 with numpy (nothing else), and for `--build` CMake, MPI (`mpicxx`) and
a C++17 compiler. The binary needs only the built-in problem generators, so the default
build (no `-D PROBLEM`) is enough. Every test runs on one MPI rank with
OpenMP threads, except `mpi_ranks` (2 ranks), which compares the dumps and `hydro.hst`
(the sums of `user.hst` depend on the order of the ranks' partial sums).

Options (`run_tests.py --help`): `--jobs N` tests at the same time (default: cores / 4),
`--star-exe PATH` a binary built with `-D PROBLEM=star_vet` (for `star_smoke`, `mlt_time` and
`wb_star`, which are skipped without it), `--threads T` OpenMP threads per run (default 4; on a shared node cap `--jobs`, the default
takes every core), `--workdir DIR` (default
`./vet_clean_runs`), `--keep` (keep the run directories of passed tests; failed tests are
always kept), `--record FILE` (write every measured number to a JSON file), `--keys
'block/key=value ...'` (extra input keys for every run, e.g. to test a new key on the whole
suite; a test's own settings win).

Exit status: 0 if all tests pass, 1 otherwise. A test whose run fails (non-zero return code)
is a FAIL with the path of its `run.log`.

## What the tests check

| test | deck | what | criterion |
|---|---|---|---|
| `conservation_sp` | `shell_sp` | closed rotating shell, 3-D spherical polar, 40 cycles (16 quick): total energy (gas + radiation), mass, total L_z, gas + radiation L_z (`dLz_gr`, no ledger) | dE/E, dM/M, dLz, dLz_gr <= 1e-13 |
| `conservation_cs` | `shell_cs` | the same on a cubed sphere, 12 cycles | as above |
| `mpi_ranks` | `shell_sp`, `shell_cs` | 2 MPI ranks against 1 on the same MeshBlocks, 6 / 4 cycles: the `floorcap_spin` keys on both meshes (every `CLEAN` counter and report line on) and `floor_pay` with `floor_audit_worst`; dumps and hydro histories bitwise. Needs an MPI build and launcher (`--mpirun`), else skipped | maxdiff = 0, `CLEAN FCONE_LZ` lines present |
| `conservation_lz` | `shell_sp`, `shell_cs` | the rotating shells with every L_z-exact cap / cut on (`clean_stag_fcone_lz`, `clean_lz_spin`, `clean_stag_fcone_emax`, `clean_wb_mcap`, `floor_from_ekin` with `vceil = 5e-3`, `floorcap_lz` on every cell), 40 / 12 cycles; restart 8 + 8 and a forced dt_retry (factor 1) with the same keys | dE, dM, dLz <= 1e-13 on both, dLz_gr <= 1e-13 on sp (cs: reported, as conservation_cs); restart and retry bitwise; no `CLEAN WARNING` |
| `floorcap_spin` | `shell_sp`, `shell_cs` | the `conservation_lz` keys + `floorcap_lz_spin` with `floorcap_rho = gas_rho` (cells that thin below the IC density are floor gas, capped at 5e-3; the rest take the rigid rotation), 40 / 12 cycles; restart 8 + 8, a forced dt_retry (factor 1), 1 vs 1 x 2 MeshBlocks in tile mode; the same without `floorcap_lz` for information | dE, dM, dLz <= 1e-13 on both, dLz_gr <= 1e-13 on sp; capped cells > 0; restart, retry and layout bitwise; no `CLEAN WARNING` |
| `floor_dt` | `shell_sp` | `hydro/dt_floor_courant = 0.6` with floor gas (`floorcap_rho = 100`), 16 cycles; restart 8 + 8 and a forced dt_retry | dE <= 1e-13; restart and retry bitwise; dt_gain reported; no `CLEAN WARNING` |
| `merge_lat` | `merge_sp`, `merge_cs`, `shell_sp` | merged lateral groups (`hydro/merge_lat_rc`, f = 4, 2, `merge_lat_mom = radial`): cell noise (density x 1e6, supersonic) 200 / 100 cycles on both meshes, radial flow, the radiating shell; `rigid` on the same noise | dM, dLz, dE, every group's change of M, E, sum V m_r, L and the lateral velocity of a radial flow <= 1e-13; restart, forced dt_retry, 1 vs 2 x 2 MeshBlocks (sp) bitwise; `rigid` must drift in L_z (>= 1e-8) |
| `restart_fixes` | `shell_cs`, `fofc_cs`, `merge_sp`, `shell_sp` | a restart file without `clean_cs_lat_wb`, `fofc_flag_exchange`, `merge_lat_mom`, `merge_lat_blockend`, `time2_gtime_slope` / `time2_gtime_pos` (written with the old value, the key's line then removed from it) | named with the old value: bitwise the run that stored it; not named: the NOTE "the fix applies", no `CLEAN WARNING`, a different run; the stored old value prints its own NOTE |
| `restart_sp` | `shell_sp` | 16 cycles straight vs 8 + restart + 8 | every dump and history row identical (difference 0) |
| `restart_cs` | `shell_cs` | 8 straight vs 4 + restart + 4 | identical |
| `dt_retry` | `shell_sp` | a forced step retry at cycle 4 with factor 1 (the step is repeated at the same dt) and with factor 0.5; the same redone as a failed step (`dt_retry_force_redo`): factor 1, and factor 0.5 on 1 vs 1 x 2 MeshBlocks (tile) with a restart just after the redo | factor 1: identical to the run without retry, also as a redo (the redo keeps the tolerance of any step); factor 0.5: runs on, one retry each; the redo converges (`CLEAN RETRY`), no `CLEAN WARNING`, layout and restart identical |
| `layout_tile` | `shell_sp` (16 x 16 lateral) | `<mesh>/layout_invariant_mode = tile`, `layout_tile = 8`, 8 cycles: 1 MeshBlock vs 1 x 2 MeshBlocks (along phi) and vs 2 x 2 (theta x phi, vacuum ray walls) | 1 x 2: identical (difference 0); 2 x 2: largest difference / largest value <= 1e-8 (see below); the 1 x 2 pair again with `dfloor_flux_limit` and `dfloor` just below the gas density (the limiter acts in every cell): identical, and its dE, dM, dLz <= 1e-13 |
| `layout_cs` | `shell_cs` (16 x 16 per panel) | tile mode on the cubed sphere: an optically thin rotating shell with turning rays, 8 cycles, one MeshBlock per panel vs 2 x 2; again with `clean_ray_min = 1` (most cells hold no rays) | identical (8e-3 of the largest lateral flux while a MeshBlock's high face took its last cell's ray sum, `RayFace`); dE, dM, dLz <= 1e-13; no `CLEAN WARNING` |
| `fofc_layout_cs` | `fofc_cs` | hydro shell on the cubed sphere whose fluxes FOFC replaces (density floor, cfl 0.8), 8 cycles: one MeshBlock per panel vs 2 x 2, 2 along r and 2 x 2 x 2; restart and a forced dt_retry of the 2 x 2 run | identical; E, L_z drifts of the 2 x 2 run, M drift of a fast shell without the density floor <= 1e-13 |
| `order_time` | `deckD` | time step halved at a fixed grid (3 levels quick, 5 full); successive-difference L1 of E, F_r, T | order of the finest pair >= 1.7; errors within 5 % of the reference |
| `order_stiff` | `wave_cart` | radiation-acoustic wave in the trapped regime (x = 100, P_rad/P_gas 1), time step T/64 ... T/512, strict solves (full only) | order of the finest pair >= 1.7; errors within 5 % of the reference |
| `order_space` | `wave_cart` | linear radiation-hydro wave, one period, N = 16 ... 128 cells with dt following the CFL | order of the finest pair >= 1.7 (full only) |
| `nlte_eps1e-1`, `nlte_eps1e-4` | `nlte_eps*` | two-level-atom non-LTE atmosphere (Jiang 2021, ApJS 253, 49, sec. 4.5): J vs the exact transport solution (`ref/nlte_exact_eps*.txt`) | L1 within 5 % of the reference value |
| `thermal_eq` | `thermal_eq` | uniform gas + radiation relaxing to equilibrium (Jiang 2021, sec. 4.1) vs the exact ODE | max relative error of E(t), T(t) <= bound |
| `mhd_equil` | `mhd_equil` | the relaxing uniform box with `<mhd>`, `eos = general` and a uniform field, every other `<mhd>` key at its default, 10 cycles; keys of `<hydro>` only under `<mhd>` | starts; gas + radiation energy and field <= 1e-13; 3 FATALs; no `CLEAN WARNING` |
| `mhd_layout` | `mhd_noise` | `<mhd>` with FOFC, Ohmic and ambipolar diffusion and the GS07 EMF on the cubed sphere, a spherical-polar wedge and a Cartesian box, 8 cycles: 1 MeshBlock (per panel) vs 2 along x1 and 2 x 2; 2 MPI ranks and restart of the 2 x 2 run | layouts <= 1e-12 (cs 5e-5 without `cs_seam_emf_exchange`); ranks and restart identical; E, L_z (box: momentum) drifts and div B <= 1e-13; FOFC fired |
| `mhd_vet` | `mhd_shell_cs`, `mhd_shell_sp`, `mhd_box` | `<mhd>` + `<rad_vet>` with Ohmic and ambipolar diffusion and FOFC, 8 cycles: closed shells and a periodic box; a second layout, 2 MPI ranks, restart, a forced dt_retry with factor 1 (full only) | layouts <= 1e-12, ranks / restart / retry identical; E <= 1e-12, M, L_z, div B <= 1e-13; no `CLEAN WARNING` |
| `mhd_geneos` | `mhd_shell_cs`, `mhd_shell_sp`, `mhd_box` | `mhd_vet` with `eos = general` and the ionizing table (`general_eos = table`), absorption on (full only) | as `mhd_vet` |
| `mhd_radwave` | `mhd_wave` | radiation-modified fast, slow and Alfven waves (Eddington closure, field at 45 deg to k) at two resolutions against the eigenmode of the gas + induction + moment system (`vett.radmhd_root`): P_rad/P_gas 0.01 - 100, thin to trapped, along x1, inclined 2-D and 3-D diagonal, Ohm + ambipolar, `eos = general` (full only) | phase speed <= 2e-3, damping per radian <= 4e-3 at the finer grid; no `CLEAN WARNING` |
| `mhd_wb` | `mhd_rhse_sp`, `mhd_rhse_cs`, `mhd_rhse_ct` | `wb_sp` / `wb_cs` under `<mhd>` and on Cartesian: no field, force-free field, `eos = general`; Cartesian `wellbalance_dynamic` (`wb_dirs = all`, `wb_arad_md`); a strong monopole (its J x B truncation is not withheld); restart (full only) | v / c_s <= 5e-8, monopole 0.1 <= 2e-5; restart identical |
| `wb_sp`, `wb_cs` | `rhse_sp`, `rhse_cs` | grey scattering envelope in hydrostatic and radiative equilibrium (spherical polar 2-D, cubed sphere); again with a rotating set (`_rot`: sp GL 11 x 12 turning about z every 2, cs GL 10 x 12 about a 3-D axis every 0.5; EMA ring, cross-fade, `vet_rotate_avg_ema_wb = true`) | max over the run of max |v| / c_s <= 1e-8, both runs |
| `shell_force_cs` | `shell_cs` | the opaque shell (`kappa_s = 300`) at rest and spherically symmetric, 8 cycles: the radiation pushes the gas along r alone, the same in every cell of a shell (`clean_cs_lat_wb`; `wb_cs` withholds the static force and cannot see this) | rms lateral / rms radial velocity and rms of the radial velocity about its shell mean / rms radial velocity: lat <= 3e-4, rad <= 2e-4 (6e-5 and 6e-5; 2.4e-3 and 6.4e-4 without the key) |
| `eddington_cs` | `shell_cs` | `clean_closure = eddington` with `vet_gd` set, 4 cycles (no rays are swept; the solve read their face sums and crashed), its restart pair, and `clean_consistency = false` with `vet_gd` on and off | drift of E, M, L_z <= 1e-13; restart identical; on and off agree to 1e-10 (3e-14); no `CLEAN WARNING` |
| `ray_clamp` | `shell_cs` | 4 radial cells: the near-tangent rays reach beyond the lateral band (`vet_gd_allow_clamp = false`), fixed set and rotating set (`clean_rot_ring = live`) | both stop with the clamp FATAL |
| `ray_cut_slab` | `star_small` | the slabs the fused sweep's two sets share below the ray cut (live ring, 2 ranks, full only, needs `--star-exe`): a restart with them against the unshared `clean_rot_fuse = false`, a stored cut far above the wanted one (every build falls back to two sweeps, counted, said once) and a forced `dt_retry` at factor 1 | dumps bitwise; the fallback said and counted only in its arm |
| `radwave` | `wave_cart` (80 x 4, x on 64 cells) | linear radiation-acoustic wave (amplitude 1e-3, at 1e-6 the thin wave ends below the solver tolerance; the exact eigenmode, 2 periods, cfl 0.3) at P_rad/P_gas, x = (10, 0.01), (100, 0.1), (10, 1), (1, 100), x = tau_cell c_eq / c the grid-scale trapping parameter; `mid_phase` is 3.17e-3 at N = 64 (bound 3e-3; 2.26e-3 at 80), and the thin phase error of 2e-5 that N = 64 used to show was a cancellation of time and space errors (1.1e-3 at cfl 0.15); phase speed and damping of a k = 1 fit vs the analytic root (`vett.radwave_root`); the x = 0.01 case also without `hydro/rad_signal_speed_diff`, x = 0.1 and 1 also with `hydro/rad_signal_speed_phase` (full only, ~1 min) | phase speed within 3e-3, damping within 5e-3 per radian, no growth; dt gain of the diff key within 2 % of 2.77; the phase key's dt >= the diff key's (2.15 at x = 0.1, 1 at x = 1) |
| `wedge_sym` | `shell_sp` | the shell at rest, spherically symmetric, GL 7 x 18 rays, 8 cycles, on a 20 deg wedge and (with `vet_gd_wall_sym = true`) 3.7 and 7.3 deg wedges: fphi = \|Σ F_phi\| / Σ \|F_r\|; the 3.7 deg wedge without the key; the 20 deg shell made theta-periodic with `vet_gd_theta_seam_sym = true`: fth = \|Σ F_theta\| / Σ \|F_r\|, and without the key; the theta-periodic 30 deg wedge with GL 11 x 12 turning about z (EMA ring, cross-fade): fth_rot, and GL 10 x 12 staggered (no mirror) | fphi, fth, fth_rot <= 1e-8; the runs without their key, and 10 x 12, stop with the FATAL |
| `large_step` | `shell_sp` | the closed rotating scattering shell of `conservation_sp` at c dt/dx ~ 200 to t = 60, 32 x 8 x 8 cells: the momentum a face flux gives the gas over one long step | kinetic energy of v_r `KE_r` within 10 % of 8.0e-4 (the value of 4 x finer radial cells and of half the step), `rho_min` >= 75 (start 100; 88 to 96 there), energy change of the gas `dEgas` <= 7e-4 (3.2e-4 measured) |
| `floor_pay` | `shell_sp` | `hydro/eos = general`, `tfloor` just above the gas temperature (the floor acts in most cells every cycle), `floor_pay = rad_column`, `floor_pay_skip`, `dfloor_flux_limit`, `floor_audit = 1`, 40 cycles: total-energy drift from t = 0; then 16 cycles with `floor_from_ekin = true` and `vceil = 1.2e-2`: energy, mass and L_z drift | dE <= 1e-13; floor energy paid > 0, unpaid 0; ekin dE, dM, dLz <= 1e-13; no `CLEAN WARNING` |
| `star_smoke` | `star_small` | the star3d_sp deck on 256 x 8 x 8 cells with the production scheme keys and the BSG tables, 20 cycles, ~15 s (needs `--star-exe`, a `-D PROBLEM=star_vet` binary) | every history column of `STAR_COLS` within 2e-4 of the stored history `ref/star_smoke.*.hst` (relative to its maximum); the two files are in git (`.gitignore` excepts them): re-record both by copying `bsg3d.{user,hydro}.hst` of the run when the scheme is meant to change; no `CLEAN WARNING` |
| `star_cs_ring` | `inputs/stars/star3d_cs` | the production RSG deck at 274 x 16 x 16 cells per panel (live ring, fused sweeps, general EOS, fofc): no start-up relaxation: 2 cycles on 6 ranks, restart to 3 on 6 and 3 ranks; a stale read in the scratch that lies over the rays' intensities (`mem_arena`) changes with the blocks per rank (needs `--star-exe` and MPI; full only; its runs set `GPU_MAX_HW_QUEUES=2`, which HIP needs with more than 2 ranks per GPU) | the two ends bitwise, no CLEAN WARNING |
| `wb_star` | `star_small` | the stars' well-balance mode (dynamic, `wb_x1`, `force_reference = wb_arad`) on the unseeded BSG column, 20 cycles, largest v/c_s: the production start (`prof`), the same with `merge_lat_mom = rigid`, and the column of `star_ic_file` with `star_gm_column`, `star_ic_balance` and one `star_ic_equil` round (`eq`) (needs `--star-exe`) | lateral velocity at r < 3.2e12 cm <= 1e-7 (5e-9; rigid 2e-4, must be >= 1e-6); radial velocity at r < 2.4e12 cm <= 3e-4 in the balanced column (9e-5), and >= 1e-3 in the supplied profile (2.3e-3), which is 0.6 - 0.9 % of g out of balance (the star's M(r) is 1 % above the point mass `star_gm`, constant mu) and linear between r = 3.36e12 and 3.65e12 cm, where the profile has no row: its photosphere starts at v/c_s 0.7 (0.1 in `eq`); no `CLEAN WARNING` |
| `mlt_time` | `star_small` | `problem/mlt_src_time`: the BSG column with F_MLT/F = 0.1 at the inner wall, ramped off by a cosine; the energy the solves add (`star_drive_ledger` [19], `P_esrc`) against the exact integral of the ramp at 4, 8, 16 steps, for `stage` and `start` (needs `--star-exe`) | order of the finest pair: `stage` >= 1.7 (2.81), `start` <= 1.3 (0.98: what the test detects); no `CLEAN WARNING` |
| `ceiling_metric` | `shell_sp`, `shell_cs` | the rotating shells (`gas_rot = 0.3`, `gas_vr = 0.2`, r 1-3), `hydro/eos = general`, `floor_from_ekin = true`, `vceil = 0.5` (below \|v\| in 1408 / 2928 cells) against the same run without the ceiling, cell by cell at t = 0 (cs: the velocity rebuilt from the contravariant v^a on the unit tangents, \|V\| the metric norm): v_phi kept (L_z), \|V\| = max(vceil, \|v_phi\|) where capped, else unchanged, v_theta 0, the radial part scaled in [0, 1]. A ceiling on Σ v_a² of mixed components or a split along the covariant phi component fails it on cs (dvphi 0.31). The cell-level C++ test `unit/ceiling_cell.cpp` (build line in its header) checks the same, the tfloor payment and the heat = metric KE change, at corner and ghost-extended points of all six panels, on sp, Cartesian and cs | dvphi, dcap, vtheta <= 1e-13; scale_ok |
| `floor_momentum` | `wave_cart` (8^3), `shell_sp`, `shell_cs` | the momentum ledger of the inversions' floors (`CLEAN FLOORMOM`): periodic box of sound waves (amplitude 0.2, negligible radiation) with `floor_from_ekin`, with `vceil` and with both, and with the energy exchange (`exch_ratio`), per cycle \|P(n+1) - P(n) - ledger\| / Σ\|p_gas\|; the t = 0 inversion of the rotating expanding shells against the dumps' p_r; restart, forced dt_retry, 1 vs 2 MeshBlocks | box and shells <= 1e-12; momentum removed in the box >= 1e-6; the same lines; no `CLEAN WARNING` |
| `key_fatal` | `shell_sp`, `mhd_noise`, `wave_cart`, `mhd_rhse_ct` | keys the code would ignore, 2 cycles: `hydro/reconstruct = ppm4` and `mhd/reconstruct = wenoz` on spherical polar, `mhd/reconstruct = ppm4` with `wb_dirs = all` (Cartesian), `ix2_bc = polar` without `use_polar_boundary`, `use_polar_boundary` on a Cartesian mesh, `hydro/wb_rho` without `wellbalance_dynamic`, `mhd/fofc_rsolver = hlle` without `fofc`, `mhd/polar_emf_diss` without a polar boundary, `time/dt_retry_views = all` on the cubed sphere (`shell_cs`); `reconstruct = plm` named vs not named | each stops with its FATAL; plm named and not named identical (difference 0) |
| `restart_keys` | `shell_sp` | a restart file storing `rad_vet/clean_tol = 5e-11` (default 1e-10), restarted with a `-i` deck that does not name it, one that does, and without a deck | one `CLEAN WARNING restart:` line with key, stored value and default; none when named; dumps identical to the restart without a deck |
| `momentum_box` | `wave_cart` (8^3, periodic 3-D) | total linear momentum, gas (`hydro_u`) + radiation (face momenta of `vet_face` / c^2), every cycle for 8 cycles; radiation wave along the xyz diagonal, amplitude 1e-3, `radwave_eig = false`, `nghost = 3`, tile 8: `dP_cyc` = largest component of \|P - P_0\| / Σ\|p_gas\| per cycle, `dP_max` its maximum (quick and full, ~10 s) | **a report, not a criterion**: passes when the run completes, the numbers are finite and there is no `CLEAN WARNING` |
| `p_ledger` | `wave_cart` (8^3, periodic 3-D), `shell_sp`, `shell_cs` | `hydro/p_ledger` (`CLEAN PLEDGER`): per step \|END\| / Σ V\|p\| of the Cartesian (box, amplitude 1e-3 and 0.1) and radial (shells) momentum ledger and every class's \|x\| dt / Σ V\|p\|; the shell_sp dumps and old history columns with and without the ledger | END and the classes that are no physical source (ntel, flr, mrg, gap; rad on the box) <= 1e-13; 8 cycles; `Pr_*` columns present; bitwise 0; no `CLEAN WARNING` |

### IDORT literature tests

The literature tests of the IDORT paper (Ma et al. 2025, sec. 4.1) and Jiang (2021, ApJS 253,
49) with the production keys, each against an exact or semi-analytic solution (decks
`idort_*`; the two-level-atom atmospheres are `nlte_eps1e-1` and `nlte_eps1e-4` above, the
fast-relaxation thermal equilibrium is `thermal_eq`). L1 is sum |q - q_exact| / sum |q_exact|.

| test | deck | what | criterion |
|---|---|---|---|
| `idort_beams` | `idort_beams` | two beams, 4 cells wide, crossing in vacuum at +-45 deg (256 x 64, x2 periodic, ring of 4 rays): E vs the exact straight beams (cell averages); `lum_top` the energy of the row 3 cells below x1max vs exact; `flux_err` the x1 energy flux (`vet_phi1`) through every row vs the injected flux; `F_nocone` the largest face momentum \|F\|/(c E_max) without `clean_stag_fcone` | L1 and `lum_top` within 5 % of 0.1772 and 1.193 (regression references, see below); `flux_err` <= 1e-5, `F_nocone` <= 1 |
| `idort_beams_mb` | `idort_beams` | the beams on 4 MeshBlocks vs 1: along x2 in tile mode (`layout_tile = 16`), along x1 in block mode | tile x2: identical; block x1: largest difference <= 1e-8 |
| `idort_shadow`, `idort_shadow_gl`, `idort_shadow_gd` | `idort_shadow*` | two beams at +-15 deg past an opaque ellipse (128 x 64; ring of 12 rays, the general 8 x 32 set, or with `vet_sc_set = gd` the global geodesic set of 162 rays, the beams snapped to its rays): E vs the exact attenuated beams | L1 within 5 % of 0.0138, 0.0126 and 0.0123; gd: 162 rays |
| `idort_cylinder` | `idort_cylinder` | homogeneous emitting cylinder (rho kappa = 10, a T^4 = 1) in vacuum, 128 x 128: E inside and outside, the energy flux Phi_x (`vet_phi1`, `F_out`) and the face momentum F_x (`vet_f1`, `FH_out`) outside along y = 0 vs the exact transport solution (`ref/idort_cylinder_exact.npz`) | each L1 within 5 % of 5.3e-4, 0.011, 0.0116, 0.203 |
| `idort_radshock_m2`, `_m5` | `idort_radshock_m*` | Lowrie & Edwards (2008) radiative shocks, Mach 2 and 5, 256 cells, started from the semi-analytic profile (`ref/idort_radshock_m*.txt`, its absolute path set by `run_tests.py`) with the Eddington closure of the diffusion-theory reference (`clean_closure = eddington`): rho, v, T_gas, T_rad at t = 2e-10 s, no shift | each L1 within 5 % of the reference (2.1e-4 ... 1.7e-3) |
| `idort_thermal_eq2` | `idort_thermal_eq2` | hot gas cooling into radiation (Jiang 2021 sec. 4.1), dt = 5 to 1e5 relaxation times: T(t) vs the exact ODE | max error within 5 % of 2.0e-3; total-energy drift <= 1e-13 |

`idort_shadow` and `idort_shadow_gl` are in `--quick`; the eight take about 45 s with
`--jobs 8`. The beams and shadow exact fields are computed in `run_tests.py`; the cylinder's
(`ref/idort_cylinder_exact.npz`, x > 0 on the 128-cell grid) is the transport integral over 64
Gauss polar angles x 4096 azimuths and the cylinder's periodic images. Notes:
- **Beams.** The 4-cell beams stay straight, at their angle and 0.8 cells wider at half maximum
  than the exact ones, but carry 1.19 x the exact energy (`lum_top`): the cap of the streaming
  share at the beam edges. They cross and go on: the field differs from the sum of the two
  single-beam fields by more than 0.05 E_b only in 40 of the 256 rows, around the crossings
  (1.3 % of the energy; the largest E is 1.13 x the exact). The excess falls as the beam
  widens: 1.76, 1.38, 1.19, 1.09 x at 1, 2, 4, 8 cells (L1 0.71, 0.35, 0.18, 0.093).
  Both are regression references. The energy flux through every row equals the injected flux
  to 6e-8 (`flux_err`).
- **Cylinder flux.** E outside the cylinder is within 1.1 % of exact transport and the
  energy flux Phi (`F_out`, `vet_phi1`, the flux the method paper measures) within 1.2 %
  (r F at x = 2-3.5 is 0.989 of exact). The face momentum F (`vet_f1`, `FH_out`) is 20 %
  above exact, a regression reference. The likely
  cause is the harmonic face opacity at the cylinder's sharp edge (not measured).
- **Shocks.** The decks name the production `time2_tableau = sdirk2`, with its default
  `time2_gtime_slope`: the L1 values are within 3 % of the method paper's bench table, which
  ran `hesdirk2` (Mach 2 and 5); the shock stays at the semi-analytic position (no shift
  needed).

### Classic radiation tests

The classic grey radiation tests (group `radiation`; decks `decks/rad_*.athinput`, ported from
the former `inputs/tests/rad_vet_*` decks, removed because they used the explicit transport, to
the production scheme: `transport = implicit`, `closure = sc`, the key block of the
other Cartesian decks, Marshak radiation faces).

| test | deck | what | criterion |
|---|---|---|---|
| `radshock`, `radshock_m5` | `rad_radshock_m2`, `rad_radshock_m5` | Lowrie & Edwards (2008, Shock Waves 18, 129) steady radiative shock, Mach 2 (subcritical) and 5 (supercritical), 256 cells, 2e-10 s from the semi-analytic profile, with the production closure (`idort_radshock_m*` runs them with `clean_closure = eddington`; same reference): L1 of rho, v, T_gas, T_rad | every L1 <= 1-2 %; L1 of rho, T_rad within 5 % of the reference (0.44 %, 1.1 %; Mach 5: 0.39 %, 0.74 %); E > 0 |
| `diffusion` | `rad_diffusion` | Gaussian pulse in a static pure scatterer at tau_cell = 1, 10, 1e3, 1e5 and in 2-D (64^2, tau_cell 1e3): variance growth vs 2 D t, L1 vs the exact Gaussian; a top-hat at tau_cell = 1e4 (quick and full, ~20 s) | \|growth/(2 D t) - 1\| <= 2e-3 (1e-4 in 2-D) at tau_cell >= 10; L1 within 5 % of the reference; top-hat within [bg, bg + amplitude] to round-off; tau_cell = 1 reported only |
| `marshak` | `rad_marshak` | non-equilibrium Marshak wave into a cold slab, constant c_v (no analytic solution: not Su & Olson) | E, T > 0 and <= the bath (maximum principle); E non-increasing into the slab; front x_f (T = T_b/2) within 5 % of the reference 0.535 |
| `advect_uniform` | `rad_advect_uniform` | uniform thick medium at 0.01 c at the coupling's fixed point (quick and full, ~3 s) | max relative change of T <= 1e-5, of E <= 1e-4 |
| `advect_pulse` | `rad_advect_pulse` | diffusing radiation pulse advected at 3e7 cm/s (Krumholz et al. 2007) vs the same pulse at rest shifted by v t | L1 <= 3e-3 at every dump; the last within 5 % of 7.75e-4 |
| `shadow` | `rad_shadow` | opaque clump in a thin gas illuminated by a 1740 K bath (Hayes & Norman 2003), 3 light-crossing times, vs the run without the clump | radiation left in the geometric umbra behind the clump <= 2 % (measured 0.61 %); cells whose sight lines miss the clump unchanged to 2 % |

Not ported: `beam`, `pulse1d` (free streaming of a set-up packet in vacuum: the gas is not
initialised by the problem generator and the sc formal solution has no source there, so
the packet does not move); `jump`, `advect_shear` (no such `vet_test` in the problem
generator); `equil`, `equil_table` (= `thermal_eq`), `atmosphere` (= `nlte_*`);
`g5_pulse2d` is part of `diffusion`.

## Star settings

The decks run the settings of the production star decks (`inputs/stars/`) where a test can
take them: cfl 0.3, `dt_retry` with the stars' factors, three ghost zones, tile mode, the
radiation signal speed with `rad_signal_speed_diff` and `_phase`, FOFC,
`dfloor_flux_limit`, `dfloor_keep_velocity`, the budget norms (`clean_clos1_budget`,
`clean_bud_maxnorm`), `clean_maxit` 200, `clean_wb_exch_local`, the default solver tolerances,
and on spherical-polar decks the 60 deg wedge with Gauss-Legendre 7 x 18 and the
interpolated, symmetrised walls (cubed sphere: 10 x 12). Not taken:

- everything in code units: the floors (`floorcap_*`, `floor_from_ekin` (general EOS only),
  `floor_pay`, `e_floor`, `vceil`), tables and the general EOS: these are covered by the
  tests that name them (`floor_pay`, `floor_dt`, `floorcap_spin`, `ceiling_metric`) and by
  `star_smoke`;
- `layout_tile`: the 60 deg wedge with 7 x 18 rays (the 2 x 2 difference is 7e-8, bound 1e-8)
  and the budget norms (8e-8), and `dfloor_keep_velocity` (it creates energy);
- `radwave`: the budget norms (one phase speed 3.1e-3, bound 3e-3) and, per case,
  `rad_signal_speed_phase`;
- `order_time`: the ray set (with 7 x 18 its coarse steps are pre-asymptotic, first order 0.9,
  and the finest F error doubles) and the wedge; `wedge_sym` sets its own rays and walls;
- decks with a fixed `dt` (`deckD`, `rad_diffusion`, `rad_marshak`, `thermal_eq`): `dt_retry`
  (a retry changes the cap of a fixed step);
- the theta-periodic seam key (the decks' theta walls are mirrors), the dynamic
  well-balancing (it needs the star problem generator: `wb_star`) and `clean_init = steady`
  (the decks start from their own analytic states); Cartesian decks have no ray set of the
  star kind.

Notes on the results of this release:
- **L_z.** `dLz` is the drift of `Lz_tot` (gas + radiation + the L_z ledger of the radiation
  solve). It is round-off (2e-16) by construction: the ledger books what the step does not
  conserve, so `dLz` checks the bookkeeping, not the conservation. The conservation error
  is the ledger's growth, `dLz_gr` (gas + radiation alone): 1e-15 on `shell_sp` (6e-7, from
  the face-flux cone cut, with `clean_stag_fcone_lz = clean_lz_spin = false`). On the cubed
  sphere the ledger stays at round-off; there `dLz_gr` (6e-5) is not a conservation measure,
  because the history's `Lrz` is the cell vector u0 F while the scheme's radiation momentum
  is its face momenta (`LzTotal`; ~0.6 % of the radiation's L_z apart).
- **Layout.** Tile mode gives bit-identical results for every layout of whole-tile
  MeshBlocks with one MeshBlock along x1 and x2 (theta): 1 x 2 (along phi) is identical.
  2 x 2 differs after 8 cycles by 8e-14 of the largest value (CPU): splitting theta puts ray walls (mirror or vacuum) inside the mesh,
  which is outside the tile guarantee, so the test bounds that difference instead of
  requiring 0.
  The `dfl` pair covers the stage-2 outflow limiter: under `time2_tableau = sdirk2` the u0
  ghosts held the stage-2 flux state while the interior held the gas state, so the two
  MeshBlocks of a face limited with different donors (before the fix: 1 x 2 differed by
  1.6e-2 of the largest value, dE 1e-8, dLz 3e-9).
- **Floors.** The floor energy of the initial conversion is paid before the t = 0 output, so
  `floor_pay` measures dE from t = 0 (`dE0`, the first step alone, is printed).
- **Narrow wedges.** On a phi-periodic wedge every ray set needs the width to be a multiple
  of its azimuthal period (GL n_theta x n_phi: 360/n_phi deg; LQ10: 90 deg; the default
  icosahedral set: 180 deg), else the run stops with a FATAL; `vet_gd_wall_sym = true` (with
  `vet_gd_wall_interp = true`) lifts it on narrow wedges. All suite decks meet the rule. A
  theta-periodic mesh needs `vet_gd_theta_seam_sym = true` (FATAL otherwise; `wedge_sym`).
  With the budget keys (`clean_clos1_budget`, `clean_bud_maxnorm`) a stage solve stops at the
  time-error tolerance, so the mirror symmetry holds only to that tolerance; the turning-set arm
  of `wedge_sym` therefore names them false, so it tests the mirror pairing under a `--keys` run.
- **Linear momentum.** The scheme does not conserve linear momentum (gas + radiation) in
  general (`p_ledger` books it by operator class): the face-flux cone cut (`clean_stag_fcone`), the well-balancing residual (`clean_stag_fwb`) and the
  floor velocity rules (`floorcap`, `dfloor_keep_velocity`, the velocity ceiling) change it
  outside the flux form. Where none of them acts the momentum exchange is exact; at amplitude
  1e-3 in `momentum_box` none acts and dP is round-off (below 1e-16).
  In a non-linear box (amplitude 0.1) the cut changes P by up to ~4e-3 of Σ|p_gas| per
  cycle. `momentum_box` therefore reports dP (printed, and kept with `--record`) and does
  not bound it; on spherical polar and cubed-sphere meshes Cartesian P is not a discrete
  invariant at all (local bases), so the test is Cartesian only.
- **Order.** Deck D's successive orders are 1.57, 1.61, 1.84 (E) over the five time steps:
  below 2 on the whole ladder. The criterion is on the finest pair.
  The exchange is not stiff in deck D and `order_space`; `order_stiff` is the stiff case.

The reference values and tolerances are in `ref/references.json` (measured on CPU with GCC 13
and `-ffp-contract=off`). Bitwise tests (restart, dt_retry, layout) need
no reference: they compare two runs of the same binary. `momentum_box` has only
`"equal": {"warnings": 0}`: its dP is recorded, not bounded.

`references.json` format, per test: `"max": {name: bound}`, `"min": {name: bound}`,
`"equal": {name: value}`, `"near": {name: [reference, relative tolerance]}`; `"last": true`
applies list-valued bounds (orders) to the last element only.

## The decks

`decks/*.athinput` are complete input files. Each starts with a comment on the problem and
names the test's physics and what differs from the defaults. A test changes only a few keys on top (for example `time/nlim`, `time/dt_fixed`,
`mesh/nx1`); `run_tests.py` writes the deck it actually ran as `in.athinput` in each run
directory. To run a deck by hand:

```bash
mkdir run && cd run && OMPI_MCA_hwloc_base_binding_policy=none OMP_NUM_THREADS=4 \
  ../build_cpu/src/athena -i ../tests/vet_clean/decks/shell_sp.athinput
```

## GPU builds

The same suite runs a GPU binary built WITHOUT `-D PROBLEM` (the built-in problem
generators; a `star_vet` binary cannot run these decks): `--exe build_tests/src/athena
--threads 1 --jobs 1` (on MI300A also `export HSA_XNACK=1`). The bitwise tests compare two
runs of one binary and must pass on a GPU too. The reference tests (`order_*`, `nlte_*`,
`thermal_eq`, `wb_*`) must agree with the CPU references within their tolerances; a much
larger difference points at the build (see the CUDA 12.2 warning in
`docs/vet_clean/BUILD.md`). On MI300A (ROCm 7.2) `layout_tile` finds 1 x 2 MeshBlocks (a phi split) 1.9e-11 from one
MeshBlock, where the CPU is bitwise: on a GPU the 1 x 2 case is bounded, not bitwise.

## Notebooks

`notebooks/` has one Jupyter notebook per test (conservation, thermal equilibrium, order,
well-balancing, non-LTE, momentum box, restart, radiation-acoustic waves, narrow wedges,
floors; `radiation.ipynb` for all the classic radiation tests): each runs its test through
this suite (or loads an existing run) and plots the result against the exact solution or the criterion. Only numpy,
matplotlib and jupyter are needed; see `notebooks/README.md`.

```bash
cd tests/vet_clean/notebooks && ATHENA_EXE=$PWD/../../../build_vettest/src/athena \
  jupyter nbconvert --to notebook --execute --inplace thermal_eq.ipynb
```

## Files

- `run_tests.sh`: builds a CPU binary if needed, then runs `run_tests.py`.
- `run_tests.py`: the tests and the PASS / FAIL logic.
- `vett.py`: readers for AthenaK `.bin` dumps and `.hst` files, the deck writer and the run
  helper (Python standard library + numpy).
- `decks/`: the input decks. `ref/`: reference values and exact solutions.
