# KEY_TABLE: the <rad_m1> input keys of the lean staggered clean scheme

Branch clean-stag2-lean-v207 (v20.7 + lean). Every key the source reads, its default for a fresh run, the comment at
the read and which decks name it. Restart files keep the values they recorded. Where a file lacks a key, the old
default applies, so restarts are unchanged (clean_face_opac / clean_stag_opac_update / clean_gas_newton_switch / clean_stag_fs_line /
clean_stag_vimp_cs / clean_stag_opac_drag / floorcap_mesh_vol / clean_stag_lat_sums: echoed in the log).

## 1. Fresh-run defaults = the production scheme (lean defaults flip)

These keys are common to every production deck (icrel BSG / M35 / T9 / T19, stag2/icrel/decks_start, t9base/decks).
A deck that leaves them out runs the production scheme.

| key | old default | new default | note |
|---|---|---|---|
| transport | explicit | **implicit** | set in recipe::CleanDefaults before the modules read it, so hydro / mhd rad_signal_speed agree |
| closure, vet_gd | m1, false | **vet_col, true** | scheme = clean on spherical polar / cubed sphere, unless the deck names closure or clean_closure |
| clean_conv_max | 1e-5 | **1e-4** | USER 11-21 |
| clean_stag_flim | true (lean flip 6f7375674; before: false) | **false** (recorded) | FDEF 11-21: flim cut the free-streaming face F (bench fs spread 0.05 -> 0.34); fwb + fcone bound F instead |
| clean_stag_fwb | false (read only when named) | **true** with the F solve and clean_consistency (recorded); a restart file without the key: false (echoed) | FDEF 11-21 (REAL3 a206152a6): the F rows well-balanced about the rays' state; free streaming kept |
| clean_stag_fcone | false (read only when named) | **true** with fwb (recorded); a restart file without the key: false (echoed) | FDEF 11-21 (REAL3 7c234abba): every cell |F| <= c E by construction; fcone without fwb fails fs |
| clean_stag_fwb_stage | - (new) | **true** with fwb (recorded); a restart file without the key: false (echoed) | FWBFAST 11-21: fwb's residual s of the stage-old state, once per stage round; A100 t9h cfl 0.9: cures the stage-2 divergence (t9r 26418) of the per-solve form; flz q booked dLz 1.4e-4 -> 2.5e-5 |
| clean_stag_fwb_t2 | - (new) | false | FWB2O 11-21: fwb_stage's s of the stage-time state E_old + h K (K1 / K2, CleanPredictE), fixed within the stage: deck D face F order 1.22/1.10/1.04 -> 1.93/1.97/1.99 (fwb off 1.86/1.93/1.96) |
| clean_face_opac | arith | **harmonic** (recorded) | USER 11-21; every mesh (the x1 face chi, distance-weighted harmonic mean) |
| vet_rotate_avg_ema_change_form | step (fresh, since v20.4) | **linear** (recorded) | USER 11-21; step costs the EMA ring's time order (deck D E/F p 1.66/1.39 vs 1.95/1.85); a restart keeps its recorded form: HG/RSG key files name linear |
| clean_stag_opac_update | false | **true** (recorded), every mesh | USER 11-21 (uniform settings on every mesh); oh3 with clean_stag_opac_passes = 3 (unchanged); acts only where the opacities depend on T |
| clean_gas_newton_switch, clean_gas_newton_switch_min | 0, 20 | **4, 6** (recorded), every mesh | USER 11-21; the old production's implicit_gas_newton_switch; acts only on a stalled solve (the cubed sphere's refresh needs it, W4h) |
| problem/he_sponge_form (he_star_m1) | exp | **lawson** (recorded) | USER 11-21 (SPONGE): the top sponge in the integrating-factor form of the SSP-RK stages, 2nd order in time; a restart file without the key keeps exp (NOTE) |
| problem/mlt_src_time (he_star_m1) | start | **stage** (recorded); start on a hesdirk2 run with another tableau | USER 11-21 (SPONGE): the MLT source at the implicit stages' time, 2nd order; with an adaptive closure the he_star restart block carries tgt^{n-1} (hdr[1] = 1: an older binary cannot read it); a restart file without the key keeps start (NOTE) |
| clean_stag_fs_line | pcr (clean_stag_fs_pcr), thomas on the cubed sphere (clean_stag_fs_pcr_cs false) | **pcr on every mesh** (recorded; the cubed sphere with clean_stag_fschur; device only) | USER 11-21 uniform settings (UNIFY S3). cs: round-off change (amplified at the H front, below the dt error), A100 hgs / rsgs_tf 2.9 / 2.1 % faster. The old keys are read when the new one is unnamed |
| clean_stag_vimp_cs | false (clean_stag_vimp off on the cubed sphere) | **true** (recorded) | USER 11-21 (UNIFY S4): the cubed sphere as spherical polar. Cost: rsgs_tf passes 37 -> 40 per 10 cycles (UNIFY-C s4/RESULTS.md) |
| clean_cs_lz_bnd_free | false (the radial boundaries' stress torque, booked in the L_z ledger) | **true** with clean_cs_lz_tel (recorded) | USER 11-22 (UNIFY-C): the booked torque was the outer boundary's one-sided lateral P n of the rotated ray sets, ~3e7 x the physical spin-down; L_z exact without a ledger |
| clean_stag_opac_drag | false | **true** (recorded) | USER 11-21 (UNIFY S9): a refresh rebuilds the drag share b and Q (CleanStagDrag); acts where clean_stag_opac_update and clean_stag_vimp run. t9h oh3: +1.3 % s/cycle, change ~1e-9 L1, 100x below the dt error, order unchanged |
| floorcap_mesh_vol | false | **true** with scheme = clean (recorded) | USER 11-21 (UNIFY S10c): the floorcap counters by the mesh volume, abs(dp) and L_z on every mesh (diagnostic; state bitwise) |
| clean_lz_ledger, clean_lz_ledger_k1 | clean_cs_lz_ledger / clean_cs_lz_k1 (cs), clean_sp_lz_bnd / _k1 (3-D sp) | **one name on every mesh** | UNIFY S10a (M1LzKey): a file or deck with an old spelling is read and written under it (old restarts bitwise) |
| clean_stag_lat_sums | false (central 4-cell lateral stress on sp / Cartesian) | **true** (recorded) with clean_stag_xdonor and no clean_stag_sbimp | USER 11-21 uniform settings (S8): the cubed sphere's lateral sums Y on sp / Cartesian (sp phi faces keep the torque form, L_z). Physics unchanged: bench WB rhse_sp2 / beams cb, cbg, ws66, shd = , deck D orders identical, t9h 300 cyc hst 4e-6..9e-4; cost +0.9 % A100 t9h, +3-6 % CPU deck D; a restart file without the key keeps false (echoed) |

These were already the fresh defaults: scheme = clean, clean_layout = staggered, clean_rays = sc, clean_xt_exact =
cons, clean_stag_xdonor = lag, clean_stag_xsrc = true, clean_stag_mg_color = mirror, clean_conv_kappa = 1.0,
clean_maxit = 800, clean_stag_thinf = true, implicit_pos_floor / _solve / _cons = true, implicit_pos_floor_gas_keep =
1, implicit_pos_gas / _cons = true, time2_face_fct = true (scheme = clean), zalesak_active = true, cs_vertex_mono /
cs_cross_mono (cubed sphere), and <hydro>/floor_pay = rad_column (fresh, eos = general with <rad_m1>).

## 2. Off by decision

| key | default | status |
|---|---|---|
| clean_pass_damp | 0 | **deprecated: do not use (USER 2026-11-21).** Do not name it in any deck or key file. The code stays until the USER decides. |
| clean_stag_thinf_held | false | off (USER 11-21: back to the v20.7 value on every mesh; LEAN6's fresh default true reverted; production names false from its next restart) |
| clean_stag_thinf_held_gate | false | off: the production decks do not name it. With it on by default, fresh t9h (A100, fwb off, 100 cycles) made 2.1e5 E-floor raises, against 0 with it off (LEAN6 bisect). |
| clean_stag_colflux | false | off (POS5: no gain; costs CW and dt on RSG) |
| time2_k1_xlim | off | off (POS5: HG dt 56 vs 98) |

## 3. Production keys that stay explicit (setup, not scheme)

These keys depend on the star or setup, so they stay in the decks:
- floors and caps: e_floor, floorcap_*, <hydro> tfloor / dfloor;
- well-balancing: force_reference(_work), wb_*;
- boundaries and fluxes: implicit_bc_*, implicit_flux_x1*, marshak_q;
- tables, units and opacity: arad, temp_unit_kelvin, rho_unit_cgs, kappa_unit, opacity, he_*;
- start: clean_init (steady | given), clean_init_dtfac, clean_init_n;
- ray set and rotation: vet_gd_gl_nmu / _nphi / _stagger, vet_gd_ls, vet_rotate_avg(_ema), vet_gd_rotate_dt, clean_rot_*,
  clean_cs_lz_tel (HG / RSG rotate, the icrel stars do not);
- GPU transport: vet_gd_halo_ipc.

Other keys the icrel decks name are old-production keys. Under scheme = clean the CleanPreflight table lists them
as having no effect, or the code ignores them (e.g. implicit_maxit, time2_pred_cap, vet_gd_twin*, report_newton_fb).

## 4. Every key (generated)

Generated by stag2/lean/tools/keytable.py.
- Default: the GetOrAdd / stag / pbool default expression ("named only" = the member's initial value applies, often
  set in the lines above the read).
- Decks:
  - P = icrel production decks;
  - H / R / X = hg / rsg / r136 production keys;
  - B = CPU bench decks;
  - G = gate decks.
- Meaning: the comment block above the first read (blank where the comment sits elsewhere).

| key | default | decks | read at | meaning (comment at the read) |
|---|---|---|---|---|
| advect_split | have_hydro && (thick_flux == M1_THICK_APHLL) | BG | src/rad_m1/rad_m1.cpp:1320 | (1c) the advective enthalpy-flux split of the E equation (design sect. 3 "Moving fluid"). It is what carries (4/3) E v where alpha -> 0 switches the transport flux off, so it is on by default exactly where alpha exists and there is a medium; with thick_flux = scaled / none it is available for compar |
| ap_form | "alpha2" | BG | src/rad_m1/rad_m1.cpp:1093 | which algebraic form of the ap_hll E-flux (design sect. 10; milestone 1c). `unified` needs no reconstruct-dependent switch, but it is MEASURABLY worse on every 1b gate that discriminates (T3 at tau_cell = 10, the resolution order, the Nyquist decay rate, T6), so the 1b pair stays the default: see ra |
| arad | 0.0 | PBG | src/rad_m1/rad_m1.cpp:1446 |  |
| athena_face_diss | named only | - | src/rad_m1/rad_m1_init.cpp:2095 |  |
| athena_face_diss_force | false | - | src/rad_m1/rad_m1_init.cpp:2098 |  |
| athena_face_mu1 | 0.5063576644955016 | - | src/rad_m1/rad_m1_init.cpp:2100 |  |
| athena_face_taucell | 4.0 | - | src/rad_m1/rad_m1_init.cpp:2099 |  |
| athena_face_taumin | 0.0 | - | src/rad_m1/rad_m1_init.cpp:2101 |  |
| be_legacy | global_variable::restart_run | - | src/rad_m1/rad_m1_init.cpp:145 |  |
| c_light | named only | PBG | src/rad_m1/rad_m1.cpp:49 | (1) parameters |
| cfl_rad | 0.4 | G | src/rad_m1/rad_m1.cpp:1070 | ---- end of the 3a hook |
| chat_over_c | 1.0 | PBG | src/rad_m1/rad_m1.cpp:55 |  |
| clean_bc | "halfrange" | B | src/rad_m1/rad_m1.cpp:181 |  |
| clean_bicg_stall_fac | 0.99 | - | src/rad_m1/rad_m1.cpp:243 |  |
| clean_bicg_stall_its | cln_linpass ? 0 : 50 | - | src/rad_m1/rad_m1.cpp:241 |  |
| clean_closure | "predict" | B | src/rad_m1/rad_m1.cpp:712 |  |
| clean_consistency | !curv // vg | B | src/rad_m1/rad_m1.cpp:333 |  |
| clean_conv | "budget" | G | src/rad_m1/rad_m1.cpp:766 | clean_conv = budget (default) / strict: budget measures every row in the cell's own budget for the step, S = the sum of its gross terms (old value, fluxes through the faces, exchange, work, beta terms), so that a change or residual is the fraction of the energy the cell exchanges in the step. Its to |
| clean_conv_est | named only | - | src/rad_m1/rad_m1.cpp:752 | clean_conv_est (read only when named): the Newton/Picard passes also stop once q/(1-q) times the change of a pass (q < 1/2 its ratio to the change of the pass before) is below clean_tol_nl; clean_nl_print: rank 0 prints the changes of every pass and the CLEAN XDONOR / CLEAN XT diagnostics (computed  |
| clean_conv_kappa | 1.0 | PHRX | src/rad_m1/rad_m1.cpp:773 |  |
| clean_conv_max | nw ? 1.0e-4 : 1.0e-5 | PHRXG | src/rad_m1/rad_m1.cpp:119 |  |
| clean_cs_lat_wb | false | - | src/rad_m1/rad_m1.cpp:524 | clean_cs_lat_wb (default off): a lateral face projects its cells' lateral stress sums after the rotation of the cell's radial direction into the face plane, so that a spherically symmetric P gives no lateral face force (CleanStagCsLat) |
| clean_cs_lz_abs | named only | - | src/rad_m1/rad_m1.cpp:498 |  |
| clean_cs_lz_bnd_free | nw | - | src/rad_m1/rad_m1.cpp:515 | clean_cs_lz_bnd_free (with clean_cs_lz_tel, default on for a fresh run, recorded): no stress torque through the radial boundaries (CleanCsTorque) |
| clean_cs_lz_tel | nw | PHR | src/rad_m1/rad_m1.cpp:510 | clean_cs_lz_tel (default on for a fresh staggered cubed-sphere run, recorded): its L_z by telescoping, the gas phi momentum of CleanCsTorque and the radiation's L_z that of its face momenta, so that the ledger holds the radial boundaries' torque alone (force_reference = clean withholds its phi row,  |
| clean_extrap_share | named only | - | src/rad_m1/rad_m1.cpp:728 | clean_extrap_share (read only when named): the moments at t^{n-1} and those of the end of the step share one array (a rejected step then redoes without the line) |
| clean_face_opac | nw ? "harmonic" : "arith" | B | src/rad_m1/rad_m1.cpp:213 |  |
| clean_fill_passes | 8 | B | src/rad_m1/rad_m1.cpp:175 |  |
| clean_gas_newton_switch | nw ? 4 : 0 | - | src/rad_m1/rad_m1.cpp:635 | clean_gas_newton_switch (default W = 4 on a fresh run, USER 11-21, recorded; 0 = off; W in [3, 64]), clean_gas_newton_switch_min (default 6 on a fresh run, else 20), clean_stall_fac (0.5): the old production's implicit_gas_newton_switch, a stalled solve's remaining passes on the bracketed gas root w |
| clean_gas_newton_switch_min | nw ? 6 : 20 | - | src/rad_m1/rad_m1.cpp:640 |  |
| clean_gas_tfloor | true | - | src/rad_m1/rad_m1.cpp:620 | clean_gas_tfloor (default true): the gas row's T bounded below by the EOS tfloor, and the gas eint written back from that root (E takes the rest) |
| clean_gd_face_fuse | named only | - | src/rad_m1/rad_m1.cpp:337 | clean_gd_face_fuse (read only when named; default true): sp sc rays, both lateral faces and the r face of a cell in one pass over the directions (cln_gd_fs, bitwise the three kernels cln_gd_lat x2 + cln_gd_rf) |
| clean_halo_depth | named only | G | src/rad_m1/rad_m1.cpp:792 | clean_halo_depth (read only when named; 0 = ng): clean_halo_direct sends at most this many ghost layers (old production's implicit_halo_depth); env CLEAN_POISON_GHOST=1 (gate only): the layers not sent are NaN |
| clean_halo_direct | true | - | src/rad_m1/rad_m1.cpp:788 |  |
| clean_halo_vars | named only | - | src/rad_m1/rad_m1.cpp:785 | clean_halo_vars (read only when named, 1..32, default 16): at most this many variables per ghost exchange of the clean solve (its scratch array holds that many) |
| clean_init | named only | PHRBG | src/rad_m1/rad_m1.cpp:804 | clean_init = given / steady (read only when named: an input without it keeps a byte-identical parameter dump). A start out of radiative equilibrium takes steady: the rays, a steady solution, run ahead of the given field's front. |
| clean_init_cnorm | 1.0e-10 | - | src/rad_m1/rad_m1.cpp:95 |  |
| clean_init_dtfac | eqk("he_ic_equil_dtfac", 1.0e6) | PB | src/rad_m1/rad_m1.cpp:825 |  |
| clean_init_n | static_cast<int>(eqk("he_ic_equil_relax_n", heq ? 40.0 : 50.0)) | PB | src/rad_m1/rad_m1.cpp:823 |  |
| clean_init_tol | 1.0e-8 | BG | src/rad_m1/rad_m1.cpp:827 |  |
| clean_layout | nw ? "staggered" : "cell" | PHRXBG | src/rad_m1/rad_m1.cpp:185 | clean_layout = staggered, the only layout (the cell layout was removed; a restart file without the key ran it) |
| clean_lin_cnorm | 1.0e-10 | - | src/rad_m1/rad_m1.cpp:96 |  |
| clean_lin_pass | nw // pl | HRX | src/rad_m1/rad_m1.cpp:130 |  |
| clean_maxit | 800 | PHRXBG | src/rad_m1/rad_m1.cpp:118 |  |
| clean_maxnl | 200 | HRXBG | src/rad_m1/rad_m1.cpp:117 |  |
| clean_maxpic | 12 | B | src/rad_m1/rad_m1.cpp:340 |  |
| clean_nl_norm | "global" | - | src/rad_m1/rad_m1.cpp:739 |  |
| clean_nl_print | named only | G | src/rad_m1/rad_m1.cpp:755 |  |
| clean_pass_damp | 0 | - | src/rad_m1/rad_m1.cpp:448 | clean_pass_damp (default 0 = off): from this pass of a solve on, a pass's change of E that turns back, l < 0, takes d/(1 - l) (CleanNewtonUpdate): period-2/3 pass cycles of the lagged exchange; a solve that converges before is unchanged |
| clean_pc_sweep | nw ? "sym" : "rb" | HRX | src/rad_m1/rad_m1.cpp:671 | clean_pc_sweep = sym (default of a fresh run, recorded; a restart file without the key keeps rb) / rb / rbr / rbrsym: the mg's smoothing on each level, a forward red-black sweep, + a second red pass (rbr), + a black-red post-sweep after the coarse correction (sym) (CleanStagVCycle) |
| clean_predict | true | B | src/rad_m1/rad_m1.cpp:176 |  |
| clean_predict_pos | named only | - | src/rad_m1/rad_m1.cpp:709 | clean_predict_pos = cut / smooth (read only when named; default cut): the predicted E of a depleting cell (CleanPredictE) |
| clean_ray_cut | named only | - | src/rad_m1/rad_m1.cpp:2588 |  |
| clean_ray_cut_count | 0 | - | src/rad_m1/rad_m1.cpp:2609 |  |
| clean_ray_cut_margin | 10.0 | - | src/rad_m1/rad_m1.cpp:2604 |  |
| clean_ray_cut_now | -1 | - | src/rad_m1/rad_m1.cpp:2608 | (the cut and its count in a restart file, CleanRayCut) |
| clean_ray_cut_r | 0.01 | - | src/rad_m1/rad_m1.cpp:2603 |  |
| clean_ray_cut_taumin | 0.0 | - | src/rad_m1/rad_m1.cpp:2605 |  |
| clean_ray_cut_wait | 10 | - | src/rad_m1/rad_m1.cpp:2606 |  |
| clean_ray_flux_in | named only | - | src/rad_m1/rad_m1.cpp:363 | (clean_ray_flux_in, the v20 A/B, is gone: the rays' state takes u0 F, the face momenta's, the only form) |
| clean_ray_min | 1.0e-3 | B | src/rad_m1/rad_m1.cpp:132 |  |
| clean_rays | "sc" | PHRXBG | src/rad_m1/rad_m1.cpp:134 | clean_rays = sc (vet_gd's own sweep; the FV rays fv1 / fv / lld are removed) |
| clean_retry_budget | false | - | src/rad_m1/rad_m1.cpp:782 | clean_conv_max (read above): the tolerance's ceiling, default 1e-4. The MeshBlock layout moves the answer to about it (a solve may stop a pass earlier or later), which a step dt_retry redoes amplifies (mpirej_retry, 1 vs 2 ranks: 1e-4 9.1e-2, 1e-5 9.2e-6) clean_retry_budget (default false): a step d |
| clean_ring_fields | named only | - | src/rad_m1/rad_m1.cpp:172 | clean_ring_fields (read only when named; default all): closure = the ring's offsets and Delta carry only the cell moments (CleanRingFields), the face sums are those of the swept sets b, b+1 (as production's ring) |
| clean_ring_jmin | true | - | src/rad_m1/rad_m1.cpp:584 | clean_ring_jmin (default on): a clean_rot_ring J ratio J_s/J_b is measured only where both sets' J are at least clean_ray_min E (a ratio of two ray-shadowed J is no ray pattern); elsewhere the offset is invalid (left out of the mix) |
| clean_rot_fuse | nw | HRX | src/rad_m1/rad_m1.cpp:149 | clean_rot_ring = off / live / jump (read only when named; live wherever the vet_gd rays rotate with a cross-fade, CleanInit): a rotating union swept as the current set (+ the next one), the older ones carried as moment offsets (CleanRingMoments) clean_rot_fuse (default on): the live ring's two sets  |
| clean_rot_ring | named only | HRXG | src/rad_m1/rad_m1.cpp:150 |  |
| clean_rot_ring_fill0 | named only | - | src/rad_m1/rad_m1.cpp:162 | (clean_rot_ring_fill0 = false, the fill's sweep of the older sets, is removed: the zero-offset fill is the only form) |
| clean_skip_time2_hist | named only | - | src/rad_m1/rad_m1.cpp:736 |  |
| clean_source_convex | false | - | src/rad_m1/rad_m1.cpp:591 | clean_source_convex (default off): the rays' vet_scatter source of the clean scheme (the state's own, CleanRays) with the absorption share kappa_P/chi at most 1, as the thermal source without vet_scatter: S = J + min(kappa_P/chi, 1) (B - J) between J and B (kappa_P > chi, a cool envelope's Planck me |
| clean_sphere_radius | named only | B | src/rad_m1/rad_m1.cpp:725 |  |
| clean_stag_colflux | false | - | src/rad_m1/rad_m1.cpp:492 |  |
| clean_stag_csvec_radial | true | - | src/rad_m1/rad_m1.cpp:360 | clean_stag_csvec_radial (true, the only form): the cubed sphere's cell vector exact for a radial field (CleanStagCsVector); false (the old least-squares fit) is gone |
| clean_stag_fcone | nw && cln_fwb | - | src/rad_m1/rad_m1.cpp:430 | clean_stag_fcone (default on for a fresh run with fwb, recorded; a restart file without the key: off): every cell's flux vector within its realisable cone, /F/ <= c E, by construction: the faces scaled towards 0 (CleanStagFCone) |
| clean_stag_flim | false | PHRX | src/rad_m1/rad_m1.cpp:381 | clean_stag_flim (default off; fwb + fcone bound F instead, and keep free streaming, which flim cuts): a face F^H above its donor's c E dissipated in the stress (conservative), the solution F within 5/4 (clean_stag_fwb: 1) of [-c E_R, c E_L], the cut booked (CleanStagFSolve); every mesh: the face F i |
| clean_stag_freal | cln_gflux && cln_fsolve | - | src/rad_m1/rad_m1.cpp:357 | clean_stag_freal (default on with the F solve and gflux): realisable by construction, the comoving part of sb at most 2 /beta/ hc kE E_k and the faces of a restart from the old scheme's file at most c E of their donor |
| clean_stag_fs_line | pcr ? "pcr" : "thomas" | - | src/rad_m1/rad_m1.cpp:233 |  |
| clean_stag_fschur | named only | - | src/rad_m1/rad_m1.cpp:291 |  |
| clean_stag_fsweeps | named only | - | src/rad_m1/rad_m1.cpp:266 | clean_stag_fsweeps N (read only when named, default 1): the staggered F solve's line solves up to N times, the line ends taking the neighbours' F of the last sweep |
| clean_stag_ftie | named only | - | src/rad_m1/rad_m1.cpp:248 |  |
| clean_stag_fwb | nw && cln_fsolve && cln_consist | - | src/rad_m1/rad_m1.cpp:388 |  |
| clean_stag_fwb_line | false | - | src/rad_m1/rad_m1.cpp:413 | clean_stag_fwb_line, clean_stag_fwb_ring (default off): fwb's s of the x1 faces less its line mean, of the phi faces (3-D spherical polar, periodic x3) less its ring mean, in proportion to the streaming share: no net radial force on a line, no torque on a ring, momentum and L_z exact without a ledge |
| clean_stag_fwb_ring | false | - | src/rad_m1/rad_m1.cpp:414 |  |
| clean_stag_fwb_stage | nw | - | src/rad_m1/rad_m1.cpp:398 | clean_stag_fwb_stage (default on for a fresh run with fwb, recorded; a restart file without the key: off): fwb's residual s of the stage-old state, once per stage (CleanStagFSolve), not of each solve's E (F_ray already is E_old's) |
| clean_stag_fwb_t2 | false | B | src/rad_m1/rad_m1.cpp:402 | clean_stag_fwb_t2 (default off): fwb_stage's s of the stage-time state E_old + h K, K the last slope (CleanStagFSolve): second order in time (no-op without fwb; FATAL with fwb but without fwb_stage) |
| clean_stag_kpos | named only | - | src/rad_m1/rad_m1.cpp:475 |  |
| clean_stag_kreal | false | - | src/rad_m1/rad_m1.cpp:491 |  |
| clean_stag_lat_sums | nw && cln_xd && !cln_sbimp | - | src/rad_m1/rad_m1.cpp:698 |  |
| clean_stag_mg_color | nw ? "mirror" : "parity" | PHRXBG | src/rad_m1/rad_m1.cpp:658 | clean_stag_mg_color = mirror (default) / parity: the mg's line colours, mirror-symmetric or (j + k) parity (CleanStagLineSolve) |
| clean_stag_opac_drag | nw | - | src/rad_m1/rad_m1.cpp:212 |  |
| clean_stag_opac_passes | 3 | - | src/rad_m1/rad_m1.cpp:207 | (old production implicit_opac_lag = pred, implicit_opac_lag_passes: the opacities of passes 0 .. N-1, then frozen: a change of the discretisation, order 2 for N >= 2) |
| clean_stag_opac_update | nw | B | src/rad_m1/rad_m1.cpp:204 |  |
| clean_stag_pc | named only | B | src/rad_m1/rad_m1.cpp:649 |  |
| clean_stag_pcr | true | - | src/rad_m1/rad_m1.cpp:192 |  |
| clean_stag_rayratio | false | - | src/rad_m1/rad_m1.cpp:579 | clean_stag_rayratio (default off): the bracket reads the rays only as ratios (their face flux over the donor's ray J, the model flux of E), times the state's E |
| clean_stag_sbimp | false | G | src/rad_m1/rad_m1.cpp:609 |  |
| clean_stag_stencil | !cln_sbimp | - | src/rad_m1/rad_m1.cpp:684 | clean_stag_stencil: the operator applied from its rows (CleanStencilRow), built once per pass (default, as production's implicit_op_stencil; on the cubed sphere the rows without the lateral stress, which the stored coefficients of its sums add, CleanStagCsLateral); not with clean_stag_sbimp |
| clean_stag_thinf | nw | PHRX | src/rad_m1/rad_m1.cpp:375 | clean_stag_thinf (default on): the staggered bracket's ray correction in the face half ranges, realisable (CleanStagBracket), not as the gated lagged flux DG |
| clean_stag_thinf_held | false | P | src/rad_m1/rad_m1.cpp:440 | clean_stag_thinf_held (default off, the v20.7 value, USER 11-21): the thin flux factor of a donor blended by theta = min(1, max(J, P)/(clean_ray_min E_old)), P its largest ray outflow: theta g(rays) + (1 - theta) h (its own half range); a ray-shadowed donor is no dam |
| clean_stag_thinf_held_gate | false | - | src/rad_m1/rad_m1.cpp:452 | clean_stag_thinf_held_gate (default off): theta of clean_stag_thinf_held times the old production's quasi-steadiness gate r^2, r = min(E_old/H, H/E_old), H = max(J, P) |
| clean_stag_tie_gid | nw && !pmy_pack->pmesh->multilevel | - | src/rad_m1/rad_m1.cpp:257 |  |
| clean_stag_vimp | cln_split | BG | src/rad_m1/rad_m1.cpp:597 |  |
| clean_stag_vimp_cs | nw | - | src/rad_m1/rad_m1.cpp:602 |  |
| clean_stag_xdonor | nw ? "lag" : "off" | PHRXBG | src/rad_m1/rad_m1.cpp:458 | clean_stag_xdonor = lag (default) / p1 / off: the off-diagonal stress flux an M-matrix term (CleanStagXDonor); clean_stag_xsrc (default on): also the face sources past their explicit share; clean_stag_xsweep (default off): a monotone line-Jacobi sweep after each linear solve |
| clean_stag_xdonor_freeze | 0 | - | src/rad_m1/rad_m1.cpp:444 | clean_stag_xdonor_freeze (default 0 = off): from this pass of a solve on, the xdonor cross flux C is held (CleanStagXDonor): its donor side and cap no longer switch from pass to pass (HG C2D: a 1-cell donor/cap flip, period-2 passes at 2e-4) |
| clean_stag_xpat | false | - | src/rad_m1/rad_m1.cpp:547 |  |
| clean_stag_xsrc | nw | PHRXBG | src/rad_m1/rad_m1.cpp:467 |  |
| clean_stag_xsweep | false | - | src/rad_m1/rad_m1.cpp:468 |  |
| clean_stall_fac | 0.5 | - | src/rad_m1/rad_m1.cpp:641 |  |
| clean_tol | 1.0e-10 | BG | src/rad_m1/rad_m1.cpp:94 |  |
| clean_tol_nl | 1.0e-11 | BG | src/rad_m1/rad_m1.cpp:131 |  |
| clean_tol_rel | named only | G | src/rad_m1/rad_m1.cpp:613 | clean_tol_rel = kappa (default 0.1, 0 = off) and clean_tol_rel_max (1e-4) (read only when named): the old production's implicit_tol_rel for the strict test, max(tol_nl, min(kappa r_0, clean_tol_rel_max)), r_0 the first pass's change |
| clean_tol_rel_max | named only | - | src/rad_m1/rad_m1.cpp:615 |  |
| clean_troot_rtol | 1.0e-13 | - | src/rad_m1/rad_m1.cpp:747 | (the exchange root stops at its residual noise, CleanExchangeTemperature) |
| clean_w_c0 | 0.1 | B | src/rad_m1/rad_m1.cpp:92 |  |
| clean_w_c1 | 0.5 | B | src/rad_m1/rad_m1.cpp:93 |  |
| clean_wb_fct | false | - | src/rad_m1/rad_m1.cpp:624 | clean_wb_fct (default false): the write-back's Zalesak limit on the written faces' energy flux, E > 0 by construction (CleanWbFct); its budget is the gas row's of clean_gas_tfloor with clean_xt_exact's work |
| clean_xt_cons | false | - | src/rad_m1/rad_m1.cpp:553 | clean_xt_cons (default off): the first iterate's gas T with A's G = 1 - beta^2/2, and after the first right-hand side T from the pass's gas row (its G and beta terms: one Newton update at X = E_k), so the first pass's A is the gas row's; clean_xt_norm (default off, study): the local exchange with E  |
| clean_xt_exact | nw ? "cons" : "off" | PHRXG | src/rad_m1/rad_m1.cpp:560 |  |
| clean_xt_norm | false | - | src/rad_m1/rad_m1.cpp:576 |  |
| closure | "m1" | PBG | src/rad_m1/rad_m1.cpp:1215 |  |
| coupling | !opac_zero | PBG | src/rad_m1/rad_m1.cpp:1368 |  |
| cs_cross_mono | named only | PHRXG | src/rad_m1/rad_m1_init.cpp:1822 |  |
| cs_cross_mono_k | named only | - | src/rad_m1/rad_m1_init.cpp:1825 |  |
| cs_precond_line | false | - | src/rad_m1/rad_m1_sph.cpp:168 |  |
| cs_seam_avg | named only | - | src/rad_m1/rad_m1_sph.cpp:447 |  |
| cs_vec_axis | named only | - | src/rad_m1/rad_m1_sph.cpp:503 |  |
| cs_vec_check | named only | - | src/rad_m1/rad_m1_sph.cpp:491 |  |
| cs_vertex_mono | named only | PHRXG | src/rad_m1/rad_m1_init.cpp:1804 |  |
| cs_vertex_mono_n | named only | - | src/rad_m1/rad_m1_init.cpp:1828 |  |
| dbg_cons_ledger | named only | - | src/rad_m1/rad_m1_init.cpp:1221 |  |
| dbg_eint_cells | named only | - | src/rad_m1/rad_m1_init.cpp:787 | DEBUG: dbg_eint_cells = N > 0 prints, in every solve whose gas write-back gives eint <= 0 somewhere, up to N of those cells per rank with their state. Read only when named; prints only; no number of the solve changes. |
| dbg_energy_tally | false | G | src/rad_m1/rad_m1_init.cpp:329 |  |
| dbg_gas_force | true | G | src/rad_m1/rad_m1.cpp:1373 |  |
| dbg_gas_force_trans | true | G | src/rad_m1/rad_m1.cpp:1375 |  |
| dbg_gas_heat | true | G | src/rad_m1/rad_m1.cpp:1374 |  |
| dbg_hydro_skip_t2 | named only | - | src/rad_m1/rad_m1_mr.cpp:83 |  |
| dbg_nonconv | named only | - | src/rad_m1/rad_m1_time2.cpp:405 |  |
| dbg_opac_patch | named only | - | src/rad_m1/rad_m1.cpp:1658 |  |
| dbg_opac_x1hi | named only | - | src/rad_m1/rad_m1.cpp:1661 |  |
| dbg_opac_x1lo | named only | - | src/rad_m1/rad_m1.cpp:1660 |  |
| dbg_opac_x2hi | named only | - | src/rad_m1/rad_m1.cpp:1663 |  |
| dbg_opac_x2lo | named only | - | src/rad_m1/rad_m1.cpp:1662 |  |
| dbg_rst_fp | named only | - | src/rad_m1/rad_m1_time2.cpp:487 |  |
| dbg_rst_fp_cycle | named only | - | src/rad_m1/rad_m1_time2.cpp:493 | dbg_rst_fp_cycle = C (read only when named): the first clean solve of cycle C also prints fingerprints at its setup and around each linear solve; dbg_rst_fp_dump: every fingerprint also writes vgd_cs to vgdcs_<tag>_r<rank>.bin (5 int extents, the data) |
| dbg_rst_fp_dump | named only | - | src/rad_m1/rad_m1_time2.cpp:496 |  |
| dbg_t2_pred | named only | - | src/rad_m1/rad_m1_time2.cpp:349 |  |
| dbg_tensor | "none" | G | src/rad_m1/rad_m1.cpp:1382 | DIAGNOSTIC, VET scaffolding (tests_m1/runs_3i_tensor): the Eddington tensor of the multi-D implicit solve is NOT taken from the cell's own flux. frozen computed once from the initial state (selected closure), never updated tilt frozen, with the tensor axis rotated by tilt*sin(2 pi x2/L2) about x3 ta |
| dbg_tensor_tilt | 0.3 | G | src/rad_m1/rad_m1.cpp:1384 |  |
| doppler_face | named only | - | src/rad_m1/rad_m1_init.cpp:1546 | <rad_m1>/doppler_face = donor / central (rad_m1.hpp dop_face; read only when named, default donor = bitwise) |
| doppler_vel | named only | - | src/rad_m1/rad_m1_init.cpp:1515 | <rad_m1>/doppler_vel = stage / vstar (rad_m1.hpp dop_vstar) |
| e_floor | (FLT_MIN) | PBG | src/rad_m1/rad_m1.cpp:1071 |  |
| f_source | "cell" | BG | src/rad_m1/rad_m1.cpp:1356 | (1d) the interface (well-balanced) form of the implicit F source, design sect. 13. Default `cell` reproduces milestone 1c bit for bit. |
| floorcap_ke | named only | - | src/rad_m1/rad_m1.cpp:1467 |  |
| floorcap_mesh_vol | nw | - | src/rad_m1/rad_m1.cpp:1487 | floorcap_mesh_vol (scheme = clean; default on for a fresh run, USER 11-21, recorded; a restart file without the key keeps off, echoed): the counters weighted by the mesh's cell volume, /dp/ and L_z of the gas momentum in the frame of y on every mesh (CleanVol, CleanLzArm); the state is the same |
| floorcap_report | named only | PX | src/rad_m1/rad_m1.cpp:1480 |  |
| floorcap_rho | named only | PX | src/rad_m1/rad_m1.cpp:1456 | floor-gas cap (rad_m1.hpp; read only when named) |
| floorcap_t | named only | - | src/rad_m1/rad_m1.cpp:1461 |  |
| floorcap_tmin_kelvin | named only | - | src/rad_m1/rad_m1.cpp:1464 |  |
| floorcap_v | named only | PX | src/rad_m1/rad_m1.cpp:1458 |  |
| force_reference | "none" | PBG | src/rad_m1/rad_m1.cpp:1550 | (2a) <rad_m1>/force_reference. `wb_arad` means: an EXTERNAL well-balanced scheme already delivers rho*arad_ref(z) to the gas momentum (box_convection's effective potential), so the coupling must hand the gas only the RESIDUAL radiative force. The reference array itself comes from the problem generat |
| force_reference_work | "auto" | PG | src/rad_m1/rad_m1.cpp:1574 |  |
| force_reference_work_pay | named only | - | src/rad_m1/rad_m1.cpp:1593 |  |
| gas_feedback | true | PBG | src/rad_m1/rad_m1.cpp:1369 |  |
| implicit_accel | "none" | BG | src/rad_m1/rad_m1_init.cpp:309 | ---- MILESTONE 3e: ANDERSON acceleration of the outer (Picard) iteration. `none` is not merely the default: nothing below is allocated and ImplicitAccelSave/Apply are never called, so an input file that does not name it is bitwise unchanged. |
| implicit_allow_multid | false | BG | src/rad_m1/rad_m1_init.cpp:1401 |  |
| implicit_am_cs | named only | BG | src/rad_m1/rad_m1_init.cpp:2454 |  |
| implicit_am_cs_cons_diag | named only | - | src/rad_m1/rad_m1_init.cpp:1306 |  |
| implicit_am_cs_deposit | named only | - | src/rad_m1/rad_m1_init.cpp:1325 | implicit_am_cs_deposit (rad_m1.hpp; read only when named) |
| implicit_am_cs_diag | named only | - | src/rad_m1/rad_m1_init.cpp:2462 |  |
| implicit_am_cs_energy | named only | - | src/rad_m1/rad_m1_init.cpp:1280 | amcsE-1111 implicit_am_cs_energy / _frac / implicit_am_cs_cons_diag (rad_m1.hpp) |
| implicit_am_cs_energy_frac | named only | - | src/rad_m1/rad_m1_init.cpp:1289 |  |
| implicit_am_cs_fold | named only | - | src/rad_m1/rad_m1_init.cpp:1311 | implicit_am_cs_fold (rad_m1.hpp; read only when named) |
| implicit_am_cs_kecap | named only | G | src/rad_m1/rad_m1_init.cpp:1264 | implicit_am_cs_kecap (rad_m1.hpp; read only when named) |
| implicit_anderson_beta | 1.0 | BG | src/rad_m1/rad_m1_init.cpp:1083 |  |
| implicit_anderson_m | 5 | BG | src/rad_m1/rad_m1_init.cpp:1082 |  |
| implicit_anderson_start | 1 | BG | src/rad_m1/rad_m1_init.cpp:1084 |  |
| implicit_bc_advect | !global_variable::restart_run | BG | src/rad_m1/rad_m1_init.cpp:2130 | runs_4f_drift: without it the Marshak end face carries only the comoving flux c q (E - E_bath), so the lab-frame radiation energy advected with the gas, A E = v (1 + chi) E, cannot leave (or enter) through the end: the end cell's E is raised by ~(A/(c q)) E at an outflow end (lowered at an inflow en |
| implicit_bc_x1max | "auto" | PBG | src/pgen/ry_per_accretor.cpp:3919 |  |
| implicit_bc_x1min | "auto" | PBG | src/pgen/ry_per_accretor.cpp:3918 |  |
| implicit_bcg_fallback | "line" | BG | src/rad_m1/rad_m1_init.cpp:976 |  |
| implicit_bcg_max_restarts | 2 | BG | src/rad_m1/rad_m1_init.cpp:975 | implicit_bcg_max_restarts / implicit_bcg_fallback (defaults 2 / line = the original, bitwise). After more than max_restarts breakdowns (a vanishing rho / omega / rhat.v, or a recurrence residual that met the tolerance while the true residual did not) the fused BiCGStab gives up. `line` then DISCARDS |
| implicit_bcg_rho_direct | named only | - | src/rad_m1/rad_m1_init.cpp:614 |  |
| implicit_bcg_sync | 1 | BG | src/rad_m1/rad_m1_init.cpp:371 | SPEED-UP 3: the host synchronisations of the BiCGStab loop. 0 = the original loop; 1 = the same recurrence with fused reductions and vector updates (3 blocking reductions per iteration instead of 5); 2 = 1 with alpha kept on the device, so rhat.v does not block either (1 rank only; with more ranks 2 |
| implicit_beam_n | named only | B | src/rad_m1/rad_m1_vet.cpp:5030 |  |
| implicit_blend | named only | BG | src/rad_m1/rad_m1_init.cpp:2924 |  |
| implicit_blend_fmode | named only | BG | src/rad_m1/rad_m1_init.cpp:2925 |  |
| implicit_blend_mode | named only | BG | src/rad_m1/rad_m1_init.cpp:2926 |  |
| implicit_bmom_half | true | BG | src/rad_m1/rad_m1_init.cpp:1786 | LIMIT 3 of the 3a findings: the DEFAULT is now `true`. It is a bug fix, not a tuning knob (the boundary face used to hand its whole momentum to one interior cell, which gave that cell 1.5 face-shares of radiative force: He column bottom-cell /v1/ 10.3 -> 0.47 v_MLT). The key is kept so that `false`  |
| implicit_cfl | -1.0 | BG | src/rad_m1/rad_m1_init.cpp:159 |  |
| implicit_clip_kecap | named only | G | src/rad_m1/rad_m1_init.cpp:1348 | implicit_clip_kecap (rad_m1.hpp; read only when named) |
| implicit_closure_lag | (md_full && !global_variable::restart_run) ? "step" : "pass" | BG | src/rad_m1/rad_m1_init.cpp:251 |  |
| implicit_closure_relax | 1.0 | BG | src/rad_m1/rad_m1_init.cpp:241 |  |
| implicit_closure_relax_thin | false | BG | src/rad_m1/rad_m1_init.cpp:245 |  |
| implicit_closure_thin_relax | global_variable::restart_run ? 0.0 : 1.5 | BG | src/rad_m1/rad_m1_init.cpp:274 |  |
| implicit_comoving_lat | named only | - | src/rad_m1/rad_m1_init.cpp:871 |  |
| implicit_cons_final | named only | - | src/rad_m1/rad_m1_init.cpp:1210 | (a) implicit_cons_final, dbg_cons_ledger (rad_m1.hpp): read only when named; the solver checks are made at the first solve prod-1112-cand: implicit_cons_final DEFAULT true (name false for the behaviour, bitwise). Not recorded when absent: a deck that cannot take it (no implicit_pos_floor here; line_ |
| implicit_conv_est | fixcl | BG | src/rad_m1/rad_m1_init.cpp:942 |  |
| implicit_conv_est_guard | named only | - | src/rad_m1/rad_m1_init.cpp:181 | (rad_m1.hpp): the guarded contraction estimate and the residual stop, read only when named (a run without them has the parameter dump of) |
| implicit_conv_est_nofloor | named only | - | src/rad_m1/rad_m1_init.cpp:187 |  |
| implicit_cs2_clip | named only | G | src/rad_m1/rad_m1_init.cpp:1930 |  |
| implicit_cs2_dacap | named only | - | src/rad_m1/rad_m1_init.cpp:1940 | implicit_cs2_dacap (read only when named; rad_m1.hpp) |
| implicit_cs2_dacap_vmin | 1.0e-4 | - | src/rad_m1/rad_m1_init.cpp:1948 |  |
| implicit_cs2_dbg | 0 | G | src/rad_m1/rad_m1_init.cpp:1900 |  |
| implicit_cs2_hll | "drag" | BG | src/rad_m1/rad_m1_init.cpp:1959 |  |
| implicit_cs2_shift | false | BG | src/rad_m1/rad_m1_init.cpp:1922 | consfix2: the hydro's dA as F0's explicit tendency (ImplicitCs2Shift) is opt-in; off, it stays in the solve's dA (clean order 2 on rw_p10t100 / rsw / rlat) |
| implicit_cs2_vimp | true | G | src/rad_m1/rad_m1_init.cpp:1905 | implicit_cs2_vimp: the vimp prediction of the conservative write-back velocity (ImplicitCs2Vimp, (b) in rad_m1.hpp). DEFAULT true (FIN3); false (diagnostic, bitwise implicit_cs2_dbg = 2) drops it: rw_p10t1 gas order 2 but rw_p10t100 E 1.75 at l 16-32, so not adopted. |
| implicit_cs2_vimp_npass | global_variable::restart_run ? 0 : 3 | G | src/rad_m1/rad_m1_init.cpp:1915 | implicit_cs2_vimp_npass: the prediction lags the neighbours' dA and the cross terms by one Picard pass, and the solves stop at 2 passes, so the E rows' dv misses the write-back's by O(dt) relative (implicit_cs2_dbg bit 32: 2.7e-3 / 1.4e-3 at rw_p10t1 l 1 / 2): a PER-STEP O(dt) term (rw_p10t1 hesdirk |
| implicit_cs2_wbsolve | true | G | src/rad_m1/rad_m1_init.cpp:1929 | implicit_cs2_clip (read only when named): gas (default: the gas is charged with the realisability-limited cell flux, so the limiter's change of F -- up to /F/ - cE in free-streaming cells -- is a momentum kick to the gas: ~100 km/s in 1e-8 floor gas, v12lz accretor) / rad (the limit's change stays w |
| implicit_cs_closure_metric | named only | - | src/rad_m1/rad_m1_init.cpp:868 | (read only when named): implicit_cs_closure_metric (cs): /F/, f and the unit n of a chi(f) closure with the panel metric, and A = v + D.v, beta.D.beta with D the MIXED tensor on face-normal components (the n = F//F/ closures; vet_col's n = r is unchanged by it) implicit_comoving_lat (sp, cs): the la |
| implicit_cs_e0_metric | true | G | src/rad_m1/rad_m1_init.cpp:860 | implicit_cs_e0_metric (cubed sphere; read only when named): the work array's transverse v and F are FACE-NORMAL components (v.n_xi, v.n_eta) on the skew panel basis, n_xi.n_eta = -cos_cell. Their Euclidean products miss the metric G^-1 = [[1, c], [c, 1]]/(1 - c^2) in beta.F, beta^2 and beta.P.beta o |
| implicit_det_reduce | det_def | G | src/rad_m1/rad_m1_init.cpp:1041 |  |
| implicit_dump_op | named only | - | src/rad_m1/rad_m1_init.cpp:720 |  |
| implicit_ebath_x1max | 0.0 | BG | src/rad_m1/rad_m1_init.cpp:2094 |  |
| implicit_ebath_x1min | 0.0 | BG | src/rad_m1/rad_m1_init.cpp:2093 |  |
| implicit_enthalpy | enth_up ? "upwind" : "plm" | PBG | src/rad_m1/rad_m1_init.cpp:1478 |  |
| implicit_eos_cache | false | BG | src/rad_m1/rad_m1_init.cpp:328 |  |
| implicit_eos_cache_check | true | BG | src/rad_m1/rad_m1_init.cpp:333 |  |
| implicit_eos_cache_check_every | named only | - | src/rad_m1/rad_m1_init.cpp:338 |  |
| implicit_eos_cache_nt | 2 | BG | src/rad_m1/rad_m1_init.cpp:332 | read only when named (an input without it is unchanged, also in the parameters a restart file records) |
| implicit_f0_init | "auto" | BG | src/rad_m1/rad_m1_init.cpp:1975 |  |
| implicit_face_mmat | named only | G | src/rad_m1/rad_m1_init.cpp:2477 | implicit_face_mmat (rad_m1.hpp; read only when named: the old scheme's face M-matrix, which scheme = clean never reads) |
| implicit_face_mmat_real | named only | - | src/rad_m1/rad_m1_init.cpp:2603 |  |
| implicit_face_mmat_rowsel | named only | - | src/rad_m1/rad_m1_init.cpp:2588 |  |
| implicit_face_mmat_trans | named only | - | src/rad_m1/rad_m1_init.cpp:2617 | DEBUG/ablation implicit_face_mmat_trans = false (read only when named): x1 only |
| implicit_face_opac | "arith" | G | src/rad_m1/rad_m1_init.cpp:2412 | <rad_m1>/implicit_face_opac = arith (default: the cell mean of M1FaceAvgX1, bitwise) / kirchhoff: every interior x1 face takes the Kirchhoff mean 1/(rho kappa_T)_f = int_{T_L}^{T_R} T^3/(rho_f kappa_T(rho_f,T)) dT / int T^3 dT, rho_f = sqrt(rho_L rho_R) (implicit_face_opac_ng Gauss points, default 8 |
| implicit_face_opac_n | !global_variable::restart_run | BG | src/rad_m1/rad_m1_init.cpp:1017 | implicit_face_opac_n (F2): rad_m1.hpp. defaults-1009 (user 10-09): default true on a fresh run; a restart whose file lacks the key keeps false |
| implicit_face_opac_ng | 8 | - | src/rad_m1/rad_m1_init.cpp:2421 |  |
| implicit_face_weight | global_variable::restart_run ? "equal" : "distance" | G | src/rad_m1/rad_m1_init.cpp:1436 | implicit_face_weight (see rad_m1.hpp): distance (default on the spherical-polar wedge) / equal (default elsewhere; the only choice on a Cartesian grid, which is uniform). On the wedge the key is recorded (GetOrAdd), so a restart keeps what the run started with. A RESTART whose embedded input predate |
| implicit_fast_kernels | ldef | BG | src/rad_m1/rad_m1_init.cpp:1651 |  |
| implicit_fclip_diag | named only | - | src/rad_m1/rad_m1_init.cpp:1861 | implicit_fclip_diag (rad_m1.hpp; read only when named) |
| implicit_final_resid_log | named only | - | src/rad_m1/rad_m1_init.cpp:176 | implicit_final_resid_log (DIAGNOSTIC, read only when named so that a run without it has the parameter dump of the run before the key) |
| implicit_final_resid_stop | named only | - | src/rad_m1/rad_m1_init.cpp:190 |  |
| implicit_final_resid_stop_rel | named only | - | src/rad_m1/rad_m1_init.cpp:194 |  |
| implicit_flux_x1max | 0.0 | PG | src/rad_m1/rad_m1_init.cpp:2092 |  |
| implicit_flux_x1min | 0.0 | PG | src/pgen/he_star_m1.cpp:674 |  |
| implicit_flux_x1min_mod_amp | named only | - | src/rad_m1/rad_m1_init.cpp:2046 |  |
| implicit_flux_x1min_mod_period | named only | - | src/rad_m1/rad_m1_init.cpp:2044 | (read only when named; test only): implicit_flux_x1min_mod_amp = A, implicit_flux_x1min_mod_period = P > 0: the inner face flux F0 (1 + A sin(2 pi t/P)), F0 = implicit_flux_x1min, set before each hesdirk2 stage solve at its time (Time2FormStage: t_n + dt), so the forcing of a time-order test is time |
| implicit_flux_x1min_ramp_t0 | named only | - | src/rad_m1/rad_m1_init.cpp:2055 |  |
| implicit_flux_x1min_ramp_w | named only | - | src/rad_m1/rad_m1_init.cpp:2053 | (read only when named; test only): implicit_flux_x1min_ramp_t0 = t0, implicit_flux_x1min_ramp_w = w > 0: the inner flux times 0.5 (1 + tanh((t - t0)/w)), a smooth switch-on (no start-up discontinuity), also at the stage times; the value at the start time is set here |
| implicit_flux_x1min_table | named only | - | src/rad_m1/rad_m1_init.cpp:2065 | (read only when named; absent = bitwise): the inner face flux from a (t, F or L) table (rad_m1.hpp ftab_on), set here at the start time and at the stage times by the hooks above |
| implicit_flux_x1min_table_t0 | 0.0 | - | src/rad_m1/rad_m1_init.cpp:2072 |  |
| implicit_flux_x1min_table_units | "flux" | - | src/rad_m1/rad_m1_init.cpp:2073 |  |
| implicit_g0_exchange | true (fresh) / restart file | BG | src/rad_m1/rad_m1_init.cpp:1178 |  |
| implicit_g0_limit | 1.0 | BG | src/rad_m1/rad_m1_init.cpp:1171 |  |
| implicit_gas_newton | full && !global_variable::restart_run | BG | src/rad_m1/rad_m1_init.cpp:326 | ---- MILESTONE 3g: the GAS-RADIATION energy coupling. implicit_eos_cache defaults false; neither key allocates anything nor is referenced when false. implicit_gas_newton: DEFAULT true since defaults-0927 (2026-09-27; runs_3g_newton_T VERDICT PASS, and every GPU timing and gate since 3g ran it; switc |
| implicit_gas_newton_switch | 10 | BG | src/rad_m1/rad_m1_init.cpp:1730 | implicit_gas_newton_switch = W (he_estall_1003). DEFAULT 10 on fresh runs (R2 default, user 10-03; it was off by default in 69e8f149 only because the switch with the mode-0 root find diverged the BSG envelope, bsg_1001/AB_SOLVER_1003.md, rst 00008: W = 10 died at cycle 14575, resid 0.31; with implic |
| implicit_gas_newton_switch_min | 20 | BG | src/rad_m1/rad_m1_init.cpp:1742 |  |
| implicit_gas_newton_u4 | false | BG | src/rad_m1/rad_m1_init.cpp:1134 | implicit_gas_newton_u4 (default false, opt-in study key): the gas Newton update of step (f) in a cell whose linearised row is EMISSION-dominated (4 c dt rho kappa_P a T^3 > rho c_v) moves u = T^4 by the Newton step, T' = (T^4 + 4 T^3 dT)^(1/4), instead of T by dT: then a T'^4 is exactly the emission |
| implicit_gas_src | gsdef ? "conservative2" : "matched" | PBG | src/rad_m1/rad_m1_init.cpp:1883 |  |
| implicit_gc_bands2 | 1 | - | src/rad_m1/rad_m1_init.cpp:545 |  |
| implicit_gc_bands3 | 1 | - | src/rad_m1/rad_m1_init.cpp:546 |  |
| implicit_gf_modes | 2 | - | src/rad_m1/rad_m1_init.cpp:555 |  |
| implicit_ghost_e_x1min | named only | - | src/rad_m1/rad_m1_init.cpp:2115 |  |
| implicit_ghost_x1min | named only | - | src/rad_m1/rad_m1_init.cpp:2111 |  |
| implicit_halo_depth | named only | - | src/rad_m1/rad_m1_init.cpp:3162 |  |
| implicit_halo_direct | fdef | BG | src/rad_m1/rad_m1_init.cpp:427 |  |
| implicit_halo_ipc | false | BG | src/rad_m1/rad_m1_init.cpp:653 | implicit_halo_ipc (m1-perf-0928, rad_m1_krylov.cpp): default false |
| implicit_halo_mpi | hmdef | BG | src/rad_m1/rad_m1_init.cpp:649 |  |
| implicit_halo_overlap | impl_halo_mpi && (global_variable::nranks > 1) && !global_variable::restart_run | BG | src/rad_m1/rad_m1_init.cpp:662 | implicit_halo_overlap (rad_m1_krylov.cpp, tests_m1/runs_3y_halo_overlap: round-off vs off, restarts bitwise). DEFAULT true since m1-accmerge wherever it is valid: implicit_halo_mpi on and more than one rank; otherwise false, silently. A restart whose file lacks the key (written before m1-accmerge, w |
| implicit_halo_ovl_faces | false | BG | src/rad_m1/rad_m1_init.cpp:670 | implicit_halo_ovl_faces (tests_m1/runs_4l_sync): default OFF since audit-m1-0928 (was true wherever the overlap is on, from 09-24 H200 weak scaling). On MI300A with moving gas it gains -0.1..+1.0 % (2 GPUs, one node) and -1.4..+0.8 % (2 nodes x 2 GPUs) per simulated second (He box, sph_wedge): no ga |
| implicit_hold | true | G | src/rad_m1/rad_m1_init.cpp:1568 |  |
| implicit_hr_cmem_fainit | named only | - | src/rad_m1/rad_m1_init.cpp:2282 |  |
| implicit_hr_cmem_gas_fa | named only | - | src/rad_m1/rad_m1_init.cpp:2287 | implicit_hr_cmem_gas_fa |
| implicit_hr_cmem_implicit | named only | - | src/rad_m1/rad_m1_init.cpp:2295 | implicit_hr_cmem_implicit (rad_m1.hpp) |
| implicit_hr_cmem_rst | named only | - | src/rad_m1/rad_m1_init.cpp:2284 |  |
| implicit_hr_q | named only | - | src/rad_m1/rad_m1_init.cpp:2727 |  |
| implicit_krylov_dev | named only | - | src/rad_m1/rad_m1_init.cpp:492 |  |
| implicit_krylov_dev_halo | named only | - | src/rad_m1/rad_m1_init.cpp:700 |  |
| implicit_krylov_fuse | fdef ? (cn_def ? 1 : 3) : 0 | BG | src/rad_m1/rad_m1_init.cpp:429 |  |
| implicit_krylov_pipe | false | BG | src/rad_m1/rad_m1_init.cpp:644 | multi-rank Krylov (rad_m1_krylov.cpp, tests_m1/runs_3w_krylov). krylov_pipe defaults OFF. halo_mpi (bitwise the ordinary exchange) defaults ON since m1-defaults wherever its preconditions hold: implicit_halo_direct, a same-level mesh (no SMR/AMR), no cubed-sphere seams and no polar boundary; otherwi |
| implicit_lag_restage | named only | - | src/rad_m1/rad_m1_init.cpp:1763 |  |
| implicit_lin_cnorm | 0.0 | BG | src/rad_m1/rad_m1_init.cpp:420 | implicit_lin_cnorm > 0 (m1-cnorm-fastpath, bsg_1001/CNORM_FAST_1003.md): the fast set stays on, with implicit_krylov_fuse = 1 instead of 3. The per-cell norm is taken in the reductions of ImplicitBiCGStab (fuse 0/1/2), which ImplicitBiCGStabTwo (fuse 3) does not have; fuse 1 is the p-update in the p |
| implicit_lin_ew_gamma | 0.9 | BG | src/rad_m1/rad_m1_init.cpp:957 |  |
| implicit_lin_ew_max | (fixcl && impl_bcg_sync >= 1) ? 1.0e-2 : 0.0 | BG | src/rad_m1/rad_m1_init.cpp:955 |  |
| implicit_lin_maxit | 200 | PBG | src/rad_m1/rad_m1_init.cpp:205 |  |
| implicit_lin_skip_true | named only | - | src/rad_m1/rad_m1_init.cpp:1217 |  |
| implicit_lin_tol | 1.0e-10 | BG | src/rad_m1/rad_m1_init.cpp:204 |  |
| implicit_lin_tolcap | named only | - | src/rad_m1/rad_m1_init.cpp:1212 |  |
| implicit_lin_tolcap_last | named only | - | src/rad_m1/rad_m1_init.cpp:1219 |  |
| implicit_lin_tolcap_stage | named only | - | src/rad_m1/rad_m1_init.cpp:1215 |  |
| implicit_line_solver | "pcr" | BG | src/rad_m1/rad_m1_init.cpp:350 | the x1 LINE SOLVE (preconditioner and line-Jacobi pass): thomas = one thread per column, serial recurrence (the original); pcr = one team per column, parallel cyclic reduction in team scratch (GPU). Same system; the answers agree to round-off, not bitwise. The gathered stack sweep (part_nblk > 1) is |
| implicit_lres_test | !fixcl | BG | src/rad_m1/rad_m1_init.cpp:941 |  |
| implicit_lz_diag | named only | - | src/rad_m1/rad_m1_init.cpp:2506 |  |
| implicit_lz_form | named only | - | src/rad_m1/rad_m1_init.cpp:2481 | implicit_lz_form (rad_m1.hpp, AM_CONSERVING_CS_1112.md sect. 10) |
| implicit_lz_mode | named only | - | src/rad_m1/rad_m1_init.cpp:2498 |  |
| implicit_lzx | named only | G | src/rad_m1/rad_m1_init.cpp:2517 | implicit_lzx (rad_m1.hpp, LZPHYS_1113.md); it closes the solve's whole L_z budget, so the implicit_am_cs deposit (a partial closure, after the write-back) is off inert (message, no slot) on a Cartesian mesh and with the gas fixed: no gas + radiation L_z budget to close there. Read only when named (t |
| implicit_lzx_diag | named only | - | src/rad_m1/rad_m1_init.cpp:2561 |  |
| implicit_lzx_weight | named only | - | src/rad_m1/rad_m1_init.cpp:2549 |  |
| implicit_marshak_face | mfdef ? "linear" : "cell" | G | src/rad_m1/rad_m1_init.cpp:1410 | implicit_marshak_face (m1-sp-order2, tests_m1/runs_5o_sporder2). DEFAULT linear on the spherical-polar wedge since m1-sp-order2b (tests_m1/runs_5q_sporder2b) with the fixed-tensor closures (eddington, vet_sc, tau, vet_col); cell with the lagged closures (m1, minerbo, kershaw: the runs_5h pp_np atmos |
| implicit_maxit | full ? 200 : 30 | PBG | src/rad_m1/rad_m1_init.cpp:203 | a lagged transverse coupling needs more outer passes than a pure column solve, so the ceiling is raised (not the count: the loop still stops at convergence). |
| implicit_mg_full_rows | named only | - | src/rad_m1/rad_m1_init.cpp:584 |  |
| implicit_mg_halo | true | BG | src/rad_m1/rad_m1_init.cpp:535 |  |
| implicit_mg_levels | 3 | BG | src/rad_m1/rad_m1_init.cpp:534 | default 3 since defaults-0927 (the recommended value, runs_5m_precond; was 2) |
| implicit_mg_sym | named only | - | src/rad_m1/rad_m1_init.cpp:581 |  |
| implicit_mg_x3line | true | - | src/rad_m1/rad_m1_init.cpp:573 |  |
| implicit_mg_x3line_halo | named only | - | src/rad_m1/rad_m1_init.cpp:590 |  |
| implicit_nc_dump | named only | - | src/rad_m1/rad_m1_init.cpp:782 | DEBUG (nc_cure_1002): implicit_nc_dump = N > 0 prints, in the last N passes of a solve that runs to implicit_maxit, the local state of each rank's worst-residual cell (optical depths, coupling stiffness, energy ratio, beta, f, the T-Newton flags). Read only when named; prints only; no number of the  |
| implicit_nofb_work | true | BG | src/rad_m1/rad_m1_init.cpp:1795 |  |
| implicit_od_cache | fdef | BG | src/rad_m1/rad_m1_init.cpp:428 |  |
| implicit_offdiag | "auto" | BG | src/rad_m1/rad_m1_init.cpp:229 | ---- MILESTONE 3b phase D: what happens to the off-diagonal Eddington terms, and how fast the closure is allowed to move between Picard passes. The defaults reproduce phase C on every 1-D / implicit_x1 configuration (there are no off-diagonal terms and no transverse closure there), and `operator` is |
| implicit_one_pass | 0 | BG | src/rad_m1/rad_m1_init.cpp:1660 | implicit_one_pass: default 0 (off) everywhere since audit-m1-0928 (8 + auto gained -2.8..+0.0 % on the moving He box and sph_wedge on MI300A, i.e. nothing, and changes results 5-22x the round-off spread; before, 0 only with force_reference = wb_arad, m1-keydefault-0927). Opt-in by name (e.g. 8). |
| implicit_one_pass_auto | !global_variable::restart_run | - | src/rad_m1/rad_m1_init.cpp:1675 | implicit_one_pass_auto (m1-onepass-auto): see ImplicitSolve. A kind of solve that goes auto_window check periods (auto_window * implicit_one_pass eligible solves) without one solve accepted after one pass is switched to one_pass = 0; after auto_reprobe such solves it is switched on again for one che |
| implicit_one_pass_auto_reprobe | 64 | - | src/rad_m1/rad_m1_init.cpp:1678 |  |
| implicit_one_pass_auto_window | 2 | - | src/rad_m1/rad_m1_init.cpp:1677 |  |
| implicit_one_pass_safety | 3.0 | - | src/rad_m1/rad_m1_init.cpp:1665 |  |
| implicit_op_check | named only | - | src/rad_m1/rad_m1_init.cpp:712 |  |
| implicit_op_check_tol | named only | - | src/rad_m1/rad_m1_init.cpp:714 |  |
| implicit_op_split_red | named only | - | src/rad_m1/rad_m1_init.cpp:468 |  |
| implicit_op_stencil | fdef && !x1per | BG | src/rad_m1/rad_m1_init.cpp:440 |  |
| implicit_op_team_red | false | BG | src/rad_m1/rad_m1_init.cpp:475 | default OFF since audit-m1-0928 (was on since 2026-09-25, m1-fast3: -12 % on the static T-S4 wedge). On moving gas on MI300A it gains 2.2-2.6 % (1 GPU) and 0.2-2.2 % (2 GPUs) per simulated second (He box, sph_wedge; 5 interleaved reps), under the ~3 % bar for a result-changing default (round-off). O |
| implicit_opac_lag | "none" | BG | src/rad_m1/rad_m1_init.cpp:1140 | implicit_opac_lag = none (default) / pred: pred evaluates the opacities in passes 0 .. N-1 only (pass 0 at the predicted T of the step predictor, or T^n on a cold start; N = implicit_opac_lag_passes, default 1) and freezes them (and drops the opac_newton face terms) for the remaining passes -- a cha |
| implicit_opac_lag_passes | (recipe::GetMode(pin) == recipe::Mode::none) ? 1 : 3 | BG | src/rad_m1/rad_m1_init.cpp:1149 | (r8 key set): DEFAULT 3 (was 1; recipe none keeps 1) |
| implicit_opac_newton | !rcp_none ? false : (impl_opac_update && impl_gas_newton && !global_variable::restart_run) | BG | src/rad_m1/rad_m1_init.cpp:1115 |  |
| implicit_opac_newton_guard | 0.5 | BG | src/rad_m1/rad_m1_init.cpp:1127 | implicit_opac_newton_guard (m1-opn-guard, 09-29): the Newton term above has the sign of the face flux times d kappa/dT and no bound; in fast optically thin cells (He presn wedge, across the Fe bump) it took a row diagonal of +802 to -9200, the row lost the M-matrix property and the solved E went neg |
| implicit_opac_newton_guard_mode | global_variable::restart_run ? 2 : 6 | BG | src/rad_m1/rad_m1_init.cpp:1388 | implicit_opac_newton_guard_mode: 0 = the diagonal test alone; + 1 = also the right-hand side, + 2 = also the sign of the neighbour entry (M1OpnGuardFace; He presn wedge gate: 7 floor clips with the diagonal alone, 0 with both), + 4 (m1-positivity) = also the row diagonal dominance (m1_impl_asm). DEF |
| implicit_opac_newton_slope_max | named only | BG | src/rad_m1/rad_m1_init.cpp:1398 |  |
| implicit_opac_newton_slope_off | named only | - | src/rad_m1/rad_m1_init.cpp:1395 | implicit_opac_newton_slope_off / _slope_max (read only when named; 0 = off). At an opacity cliff (He recombination: /d ln kappa/d ln T/ ~ 5-14 in the TOPS tables) the d(rho kappa_T)/dT term of a strongly coupled cell makes the E iterate cycle (He giant N897 photosphere, HEGIANT_PICARD.md: 77 passes, |
| implicit_opac_update | force_ref == M1_FREF_WB_ARAD | BG | src/rad_m1/rad_m1_init.cpp:1098 | implicit_opac_update: default true with force_reference = wb_arad (m1-keydefault-0927, docs/dev/ke_dt_0926.md 8), false otherwise |
| implicit_partition | "none" | BG | src/rad_m1/rad_m1_init.cpp:1697 | LIMIT 4 of the 3a findings is NOT implemented in 3a2: a column still has to live inside one MeshBlock along x1 (the fatal below). The option is parsed so that the input files and the gate scripts can already name it, and `gather` fatals rather than silently doing something else. |
| implicit_pcr_check | false | BG | src/rad_m1/rad_m1_init.cpp:361 |  |
| implicit_pcr_team | 0 | BG | src/rad_m1/rad_m1_init.cpp:360 |  |
| implicit_picard_log | 0 | BG | src/rad_m1/rad_m1_init.cpp:724 | the Picard pass count (bench/m1_picard_0923): a per-pass log, off by default |
| implicit_pos_floor | true (fresh) / restart file | PBG | src/rad_m1/rad_m1_init.cpp:1180 |  |
| implicit_pos_floor_cons | named only | PHRXG | src/rad_m1/rad_m1_init.cpp:1838 | implicit_pos_floor_cons (read only when named; rad_m1.hpp; DEFAULT true (USER 10-29); false = the old floor that creates energy) |
| implicit_pos_floor_gas_keep | 1.0 | PHRXBG | src/rad_m1/rad_m1_init.cpp:1852 |  |
| implicit_pos_floor_solve | impl_pos_floor && (!rst // pfs_rd) | PBG | src/rad_m1/rad_m1_init.cpp:1186 | (user 10-09): default = implicit_pos_floor on a fresh run (it needs it); a restart whose file lacks the key keeps false (c): a restart whose file lacks the key is no longer silently false: a WARN line says so, and <rad_m1>/implicit_pos_floor_solve_rst_default = true (read only when named) gives such |
| implicit_pos_floor_solve_rst_default | named only | - | src/rad_m1/rad_m1_init.cpp:1191 | prod-1112-cand: DEFAULT true (name implicit_pos_floor_solve_rst_default = false for the prod-1111 behaviour: such a restart keeps false) |
| implicit_pos_gas | true (fresh) / restart file | PHRXBG | src/rad_m1/rad_m1_init.cpp:1179 |  |
| implicit_pos_gas_cons | true | PHRXBG | src/rad_m1/rad_m1_init.cpp:1849 |  |
| implicit_pos_gas_frac | named only | - | src/rad_m1/rad_m1_init.cpp:1243 |  |
| implicit_pos_transport | named only | - | src/rad_m1/rad_m1_init.cpp:1588 |  |
| implicit_pos_transport_mem | named only | - | src/rad_m1/rad_m1_init.cpp:1597 |  |
| implicit_precond | pcdef | BG | src/rad_m1/rad_m1_init.cpp:507 |  |
| implicit_precond_fallback | "none" | G | src/rad_m1/rad_m1_init.cpp:998 | implicit_precond_fallback (default none = bitwise the old path). A linear solve whose BiCGStab stops at implicit_lin_maxit without meeting its test (not a breakdown: that is implicit_bcg_fallback) is redone ONCE with this preconditioner, from the same initial guess (the Picard iterate) and with the  |
| implicit_precond_float | false | BG | src/rad_m1/rad_m1_init.cpp:441 |  |
| implicit_predictor | fixcl ? "step" : "none" | PBG | src/rad_m1/rad_m1.cpp:179 | implicit_predictor (step unless named none): the stage solve starts from the predicted E (CleanFirstIterate), as the old production's step predictor |
| implicit_predictor_ecap | named only | - | src/rad_m1/rad_m1_init.cpp:1810 |  |
| implicit_predictor_frozen | false | BG | src/rad_m1/rad_m1_init.cpp:1076 |  |
| implicit_predictor_order | (ldef && impl_pred) ? 2 : 1 | BG | src/rad_m1/rad_m1_init.cpp:1683 |  |
| implicit_predictor_tcap | named only | PHRXG | src/rad_m1/rad_m1_init.cpp:1819 |  |
| implicit_rad_inertia | true | BG | src/rad_m1/rad_m1_init.cpp:1790 | the moving-medium fixes (rad_m1.hpp); every one off reproduces a61c1e75 the same default on the cubed sphere (the transverse deposit and the fixed point are the cs blocks of the write-back; cs-s4-1017 had false there) |
| implicit_rad_inertia_solve | named only | - | src/rad_m1/rad_m1_init.cpp:1792 | implicit_rad_inertia_solve (read only when named; rad_m1.hpp) |
| implicit_realisable_coupling | named only | BG | src/rad_m1/rad_m1_init.cpp:850 | implicit_realisable_coupling (gd_physfix_1005, read only when named; default false = bitwise the old code): the gas sees the REALISABLE radiation flux. In optically thin cells with c dt/dx >> 1 the Picard iterate of the cell flux is unrealisable (BSG far-thin top: mean /F/ = 10-20 c E, ~200 c E in t |
| implicit_recon | "dc" | BG | src/rad_m1/rad_m1_init.cpp:1462 | (scheme = clean: not read; the clean solve has its own face, and CleanPreflight lists the keys an input names) |
| implicit_recon_dgpass | named only | - | src/rad_m1/rad_m1_init.cpp:1775 |  |
| implicit_recon_lag | splag ? "step" : "picard" | BG | src/rad_m1/rad_m1_init.cpp:1760 |  |
| implicit_recon_npass | -1 | BG | src/rad_m1/rad_m1_init.cpp:1774 | MILESTONE 3c: freeze the deferred correction, and with it the plm limiter's choice, after this many Picard passes. The 3a2 finding is that the limiter keeps switching on a handful of cells and the iteration is a small limit cycle that never meets implicit_tol, so plm_dc costs implicit_maxit passes p |
| implicit_recon_w | -1.0 | BG | src/rad_m1/rad_m1_init.cpp:1706 |  |
| implicit_res_dmin | named only | - | src/rad_m1/rad_m1_init.cpp:831 | ...and the options that cut it. Since bench/m1_defaults_0923 they DEFAULT ON for the closures that do not read the iterate (eddington, vet_sc, tau: 1.6x on the GPU, statistics moved at the solver-tolerance level) and stay OFF for m1 / minerbo / kershaw, whose closure moves with the iterate: there th |
| implicit_res_floor | 0.0 | BG | src/rad_m1/rad_m1_init.cpp:1707 |  |
| implicit_res_rmax | named only | - | src/rad_m1/rad_m1_init.cpp:834 |  |
| implicit_resid_fatal | 1.0e-2 | PBG | src/rad_m1/rad_m1_init.cpp:903 | implicit_resid_fatal: a Picard solve that ends NON-CONVERGED with resid > this value (or a non-finite resid or lin_resid) is a diverged solve: FATAL before its state spreads (BSG production 10-01: one 200-pass solve with resid 1.2e6 at t = 2.91e5 s, the whole domain NaN ~15 cycles later). 0 = off. D |
| implicit_resid_fatal_masked | 1.0 | - | src/rad_m1/rad_m1_init.cpp:936 |  |
| implicit_solver | nd ? "bicgstab" : "line_jacobi" | BG | src/rad_m1/rad_m1_init.cpp:208 | (USER 10-19): default bicgstab under the new default set (utils/recipe.hpp); a restart whose input lacks the key / recipe none: line_jacobi |
| implicit_src_stable | !global_variable::restart_run | BG | src/rad_m1/rad_m1_init.cpp:947 | implicit_src_stable (stall_1002): the cancellation-free form of the gas-eliminated source row. DEFAULT true (user 10-02) for a fresh run; false (the plain form, bitwise) on a restart whose file lacks the key, with one line on rank 0. The value is recorded, so later restarts keep it; a named key alwa |
| implicit_stall_fac | 0.5 | BG | src/rad_m1/rad_m1_init.cpp:1743 |  |
| implicit_stall_trace | named only | - | src/rad_m1/rad_m1_init.cpp:795 | ...and with time_scheme = be, where Time2Init returns before reading it, the DEBUG (he_estall_1003): implicit_stall_trace = N > 0: in a solve still unconverged at pass implicit_stall_trace_p0 (default 40), fix each rank's worst cell of that pass and print its NCD line (tag 3) before step (c) and (ta |
| implicit_stall_trace_max | 30 | - | src/rad_m1/rad_m1_init.cpp:798 |  |
| implicit_stall_trace_p0 | 40 | - | src/rad_m1/rad_m1_init.cpp:797 |  |
| implicit_thin_freeze | 1.0e-2 | BG | src/rad_m1/rad_m1_init.cpp:882 | implicit_thin_freeze (m1-picard-aa, 095ca6ee): cells with c dt rho kappa_P below it at pass 0 keep their start-of-solve opacities for every Picard pass and take the frozen-opacity gas-T find. DEFAULT 1e-2 on fresh runs (user 10-04): on the He presn 128^2 wedge (he_mltpp_1002/FREEZE_AB.md, rst 00055, |
| implicit_timers | named only | - | src/rad_m1/rad_m1_init.cpp:800 |  |
| implicit_timers_mg_x3 | named only | - | src/rad_m1/rad_m1_init.cpp:587 |  |
| implicit_tol | 1.0e-8 | BG | src/rad_m1/rad_m1_init.cpp:160 |  |
| implicit_tol_rel | 0.0 | BG | src/rad_m1/rad_m1_init.cpp:169 | implicit_tol_rel = kappa (default 0 = off; opt-in study key). A solve that starts from the step predictor measures its first change r_0 = max /x_1 - x_pred//x: the predictor (an extrapolation of the last increments, O(dt^3) off under predictor_order 2) misses the fixed point by about the local trunc |
| implicit_tol_rel_max | 1.0e-4 | BG | src/rad_m1/rad_m1_init.cpp:170 |  |
| implicit_trans_fmax | 1.0 | BG | src/rad_m1/rad_m1_init.cpp:302 |  |
| implicit_trans_limit | "none" | BG | src/rad_m1/rad_m1_init.cpp:290 | ---- the TRANSVERSE realizability limiter. `none` is not merely the default: it takes the OLD expression for theta everywhere (see ImplicitTransTheta), so an input file that does not name it is bitwise unchanged. |
| implicit_tsolve_opac | true | G | src/rad_m1/rad_m1_init.cpp:734 | implicit_tsolve_opac (nc_cure_1002): the Newton-fallback gas-temperature root find with kappa_P(T), kappa_E(T) evaluated at the trial T (M1ImplTemperatureOpac) instead of at the lagged iterate; it removes the /dT//T = 0.365 limit cycles of shock-heated near-void cells (He presn wedge, /viper/ptmp2/j |
| implicit_tsolve_opac_mode | 1 | G | src/rad_m1/rad_m1_init.cpp:764 | implicit_tsolve_opac_mode (tsolve_root_fix_1003): 0 = the root find above; bit 1 (value 1) takes the root of the kappa(T) equation NEAREST the lagged T (M1ImplTemperatureOpacNear) instead of the halve/double bracket's; bit 2 (value 2) is a slope guard: a cell whose d ln kappa_P / d ln T at the lagge |
| implicit_tsolve_opac_slope_min | named only | - | src/rad_m1/rad_m1_init.cpp:775 |  |
| implicit_tsolve_opac_start | 20 | G | src/rad_m1/rad_m1_init.cpp:745 |  |
| implicit_tstart_old | named only | - | src/rad_m1/rad_m1_init.cpp:1807 |  |
| implicit_vimp | vdef | BG | src/rad_m1/rad_m1_init.cpp:1513 |  |
| implicit_vimp_bc | !global_variable::restart_run | BG | src/rad_m1/rad_m1_init.cpp:2138 |  |
| implicit_vimp_ewf | named only | - | src/rad_m1/rad_m1_init.cpp:1605 |  |
| implicit_vimp_fold | false | BG | src/rad_m1/rad_m1_init.cpp:1655 | implicit_vimp_fold: default OFF since audit-m1-0928 (gain 0.9-2.3 % per simulated second on the moving He box and sph_wedge, MI300A 1 and 2 GPUs; under the ~3 % bar for a result-changing default). Opt-in by name. |
| implicit_vimp_jscale | named only | - | src/rad_m1/rad_m1_init.cpp:1611 |  |
| implicit_vimp_kill | named only | - | src/rad_m1/rad_m1_init.cpp:1579 |  |
| implicit_vimp_lag_local | named only | G | src/rad_m1/rad_m1_init.cpp:1582 |  |
| implicit_wb_vcap | named only | - | src/rad_m1/rad_m1_init.cpp:1366 | implicit_wb_vcap (rad_m1.hpp; read only when named) |
| implicit_wb_vcap_measure | named only | - | src/rad_m1/rad_m1_init.cpp:1375 |  |
| implicit_wb_vcap_print | named only | - | src/rad_m1/rad_m1_init.cpp:1378 |  |
| implicit_wb_vkick_be_retry | named only | - | src/rad_m1/rad_m1_init.cpp:1257 |  |
| implicit_wb_vkick_reject | named only | - | src/rad_m1/rad_m1_init.cpp:1250 | implicit_wb_vkick_reject (rad_m1.hpp; read only when named: the old scheme's, which scheme = clean never reads) |
| implicit_work_fix | named only | G | src/rad_m1/rad_m1_init.cpp:1798 | implicit_work_fix (read only when named; rad_m1.hpp; DEFAULT true since USER 10-25; false = the old write-back) |
| kappa_e | kappa_p | BG | src/rad_m1/rad_m1.cpp:1165 |  |
| kappa_f | kappa_p | BG | src/rad_m1/rad_m1.cpp:1166 |  |
| kappa_p | 0.0 | BG | src/rad_m1/rad_m1.cpp:1164 |  |
| kappa_s | 0.0 | BG | src/rad_m1/rad_m1.cpp:1167 |  |
| kappa_unit | 1.0 | PG | src/rad_m1/rad_m1.cpp:1135 |  |
| marshak_q | 0.5 | PBG | src/rad_m1/rad_m1_init.cpp:1402 |  |
| mem_arena | named only | - | src/rad_m1/rad_m1_memarena.cpp:67 |  |
| opac_a | 0.0 | BG | src/rad_m1/rad_m1.cpp:1170 |  |
| opac_abs_kappa_s | named only | - | src/rad_m1/rad_m1.cpp:1154 |  |
| opac_abs_rho_max | named only | - | src/rad_m1/rad_m1.cpp:1152 | <rad_m1>/opac_abs_rho_max (read only when named, default 0 = off, bitwise inert): an ABSORPTION MASK of the table opacity for tenuous gas. Cells with a CODE density below it get kappa_P = kappa_E = 0 (no emission, no absorption) and the transport opacity <rad_m1>/opac_abs_kappa_s [cm^2/g] (electron  |
| opac_b | 0.0 | BG | src/rad_m1/rad_m1.cpp:1171 |  |
| opac_freeze | false | BG | src/rad_m1/rad_m1.cpp:1371 | debug switches (see rad_m1.hpp); every default is the production value |
| opacity | "const" | PBG | src/rad_m1/rad_m1.cpp:1113 | (1b) opacities, per unit mass, in code units. kappa_e defaults to kappa_p and kappa_f to kappa_p as well, so a single kappa_p gives a grey absorbing medium. |
| polar_ring_avg | false | G | src/rad_m1/rad_m1_ringavg.cpp:91 |  |
| polar_ring_avg_dist | true | - | src/rad_m1/rad_m1_ringavg.cpp:106 | the coarse phi distance (study key) |
| polar_ring_avg_ke | global_variable::restart_run ? "gas" : "rad" | - | src/rad_m1/rad_m1_ringavg.cpp:111 |  |
| reconstruct | "plm" | BG | src/rad_m1/rad_m1.cpp:2473 | reconstruction of (E, f_i). dc was added in 1b so that the thick-limit blend can be compared with the first-order scheme it reduces to (Berthon & Turpault / Bloch et al.); the three high-order kernels of src/reconstruct/ were added in 1c-B. |
| restart_pack_per_mb | named only | G | src/rad_m1/rad_m1_time2.cpp:635 |  |
| rho_ref | 1.0 | BG | src/rad_m1/rad_m1.cpp:1168 |  |
| rho_unit_cgs | 1.0 | PG | src/rad_m1/rad_m1.cpp:1134 |  |
| rsla_force | false | BG | src/rad_m1/rad_m1.cpp:1643 |  |
| rsla_vmax | -1.0 | BG | src/rad_m1/rad_m1.cpp:1642 |  |
| rsla_warn | 0.1 | BG | src/rad_m1/rad_m1.cpp:1641 | (2a) the RSLA start-up check of design sect. 2 |
| rst_state_ext | named only | - | src/rad_m1/rad_m1_time2.cpp:484 |  |
| scaled_prefactor | 20.0 | BG | src/rad_m1/rad_m1.cpp:1104 |  |
| scheme | scd ? "clean" : "default" | PHRXBG | src/rad_m1/rad_m1.cpp:85 |  |
| source_form | "full" | BG | src/rad_m1/rad_m1.cpp:1339 | the O(v/c) CONTROL of T4/T4b. Never a production setting: it is here to exhibit the spurious heating ~ (4/3) c rho kappa beta^2 E that the full source form cancels. |
| split_vel | "recon" | BG | src/rad_m1/rad_m1.cpp:1325 | (1c-B) which velocity builds the enthalpy flux at the face. `cell` is the 1c-A behaviour and is first order there; `recon` reconstructs v with the same method as (E, f_i), which is what makes A second order in a sheared flow. |
| subcycle | true | G | src/rad_m1/rad_m1.cpp:1072 |  |
| t_ref | 1.0 | BG | src/rad_m1/rad_m1.cpp:1169 |  |
| temp_unit_kelvin | 1.0 | PG | src/rad_m1/rad_m1.cpp:1133 | the unit conversions of the TABLE lookup, which is in kelvin and g/cm^3 (design sect. 4 / the stage-2 plan sect. 4). T[K] = T_code*temp_unit_kelvin: the same direction Conduction uses with Units::temperature_cgs(), which for the He box's <units> 1/1/1 with mu = 1 is mu m_u/k_B = 1.202724e-8 K per co |
| thick_flux | "none" | BG | src/rad_m1/rad_m1.cpp:1076 | thick-limit flux (design sect. 3, and the 1b correction to it recorded in rad_m1_closure.hpp: ap_hll is the BLEND form, which is what survives PLM). |
| time2_accept_fallback | named only | G | src/rad_m1/rad_m1_time2.cpp:610 |  |
| time2_adm_poslim | named only | G | src/rad_m1/rad_m1_time2.cpp:618 | <rad_m1>/time2_adm_poslim (rad_m1.hpp; read only when named) |
| time2_auto_be | named only | PHRXG | src/rad_m1/rad_m1_time2.cpp:622 | <rad_m1>/time2_auto_be = N (rad_m1.hpp; read only when named) |
| time2_be_only | named only | - | src/rad_m1/rad_m1_init.cpp:154 |  |
| time2_dbg_fail | named only | - | src/rad_m1/rad_m1_time2.cpp:402 |  |
| time2_dtr_poslim | named only | G | src/rad_m1/rad_m1_init.cpp:1585 |  |
| time2_enth_be | !global_variable::restart_run | PBG | src/rad_m1/rad_m1_time2.cpp:179 | time2_enth_be: with time2_enth_vel = central the stage solves build a_f as the mean of the two cells, but the backward-Euler solves of the same run (the start step, the time2_start slope probe, the rotation probe, the step after a failed stage) take implicit_enthalpy's own a_f (plm: the 4-point plm  |
| time2_enth_vel | named only | - | src/rad_m1/rad_m1_time2.cpp:181 |  |
| time2_exch_split | named only | - | src/rad_m1/rad_m1_time2.cpp:430 | time2_exch_split = none / psi (rad_m1.hpp; read only when named) |
| time2_face_fct | t2_ffct | PHRXBG | src/rad_m1/rad_m1_time2.cpp:563 |  |
| time2_face_fct_emin | named only | - | src/rad_m1/rad_m1_time2.cpp:568 |  |
| time2_face_fct_form | named only | - | src/rad_m1/rad_m1_time2.cpp:578 | (the one-pass donor form was removed: it fails the accuracy gate) |
| time2_face_fct_rounds | named only | - | src/rad_m1/rad_m1_time2.cpp:574 |  |
| time2_face_fct_sweep | named only | - | src/rad_m1/rad_m1_time2.cpp:571 |  |
| time2_floor_eint | named only | - | src/rad_m1/rad_m1_time2.cpp:614 | <rad_m1>/time2_floor_eint (rad_m1.hpp; read only when named) |
| time2_fref_wacc | true | PG | src/rad_m1/rad_m1.cpp:1607 | time2_fref_wacc (default true, USER 10-22; false = the old payment): the exact payment above also for the backward-Euler steps of a hesdirk2 run, which otherwise pay v' dt rho a_ref (B deck, first step: +76 energy units created, rE jump) |
| time2_inc_energy | named only | - | src/rad_m1/rad_m1_time2.cpp:413 | dbg_t2_admiss = N (DEBUG, default 0 = off, read only when named so the restart echo and every result stay bitwise): for the first N non-admissible hesdirk2 stages print which quantity went <= 0 (solved E, solved T, written-back gas eint) and, for each, the cell of its minimum with the stage OLD vect |
| time2_inc_exch | named only | - | src/rad_m1/rad_m1_time2.cpp:447 | time2_inc_exch = none / relax (rad_m1.hpp; read only when named) |
| time2_inc_exch_y0 | named only | - | src/rad_m1/rad_m1_time2.cpp:456 |  |
| time2_k1_pos_after_stiff | named only | - | src/rad_m1/rad_m1_time2.cpp:588 |  |
| time2_k1_pos_theta | named only | - | src/rad_m1/rad_m1_time2.cpp:591 |  |
| time2_k1_redist | named only | PHRXBG | src/rad_m1/rad_m1_time2.cpp:538 |  |
| time2_k1_redist_r | named only | - | src/rad_m1/rad_m1_time2.cpp:585 |  |
| time2_k1_reset | named only | - | src/rad_m1/rad_m1_time2.cpp:499 |  |
| time2_k1_reset_grow | named only | - | src/rad_m1/rad_m1_time2.cpp:597 |  |
| time2_k1_stiff | named only | PHRXG | src/rad_m1/rad_m1_time2.cpp:520 |  |
| time2_k1_xlim | named only | - | src/rad_m1/rad_m1_time2.cpp:508 |  |
| time2_k1_xpat | named only | - | src/rad_m1/rad_m1_time2.cpp:514 | <rad_m1>/time2_k1_xpat (default false, read only when named): Time2K1XLim's theta in Patankar form, positive by construction with no threshold (rad_m1.hpp); = pos: only where the unlimited loss exceeds 0.99 of the loser's holding (Time2K1XLim) |
| time2_one_pass_safety | rs ? 0.0 : 30.0 | BG | src/rad_m1/rad_m1_time2.cpp:368 | time2_one_pass_safety (tests_m1/runs_5f_h2fast): the implicit_one_pass safety factor of the stage solves. One-pass acceptance leaves a Picard error ~ tol/safety where the two-pass test leaves ~ q tol; the radwave G1 order at implicit_tol = 1e-11 needs the stage solves closer to the latter (safety 3: |
| time2_pred_cap | global_variable::restart_run ? 0.0 : 0.5 | PBG | src/rad_m1/rad_m1_time2.cpp:336 |  |
| time2_reject | named only | BG | src/rad_m1/rad_m1_time2.cpp:469 |  |
| time2_start | named only | PBG | src/rad_m1/rad_m1_time2.cpp:379 | time2_start / time2_start_frac (rad_m1.hpp; read only when named): be (default, bitwise) / slope. The Lie-split backward-Euler first step leaves an O(dt) start-up error that the second-order steps after it carry along as a free mode (h2order_1011: P = 10 linear wave, tau 100, v 1.1 1.2 1.0 -> see NO |
| time2_start_frac | named only | - | src/rad_m1/rad_m1_time2.cpp:392 |  |
| time2_tableau | named only | - | src/driver/driver.cpp:157 |  |
| time2_vet_extrap | named only | - | src/rad_m1/rad_m1_time2.cpp:200 |  |
| time2_vet_x2_fields | 9 | G | src/rad_m1/rad_m1_time2.cpp:296 |  |
| time_scheme | h2def ? "hesdirk2" : "be" | PBG | src/rad_m1/rad_m1_init.cpp:134 |  |
| transport | "explicit" | PBG | src/hydro/hydro.cpp:451 |  |
| vet_axis | "flux" | - | src/rad_m1/rad_m1_vet.cpp:2009 |  |
| vet_bc_bath | !global_variable::restart_run | G | src/rad_m1/rad_m1.cpp:2058 |  |
| vet_bc_diff | bdd | PBG | src/rad_m1/rad_m1.cpp:1806 |  |
| vet_bc_inner | "auto" | PG | src/rad_m1/rad_m1.cpp:2950 |  |
| vet_bc_reflect | "auto" | PG | src/rad_m1/rad_m1.cpp:2965 |  |
| vet_col_axis | "radial" | BG | src/rad_m1/rad_m1.cpp:1848 |  |
| vet_col_dump | named only | - | src/rad_m1/rad_m1.cpp:1858 |  |
| vet_col_dump_every | 1 | - | src/rad_m1/rad_m1.cpp:1860 |  |
| vet_col_fk_min | global_variable::restart_run ? (1.0/3.0) : 0.0 | BG | src/rad_m1/rad_m1.cpp:1732 | vet_col_fk_min (m1-sp-order2b): the lower clamp of f_K, default 0 on every mesh (a restart whose file lacks the key keeps 1/3, echoed) |
| vet_col_lat | named only | B | src/rad_m1/rad_m1.cpp:1865 | vet_col_lat (m1-vetcol-lat, rad_m1_vetlat.cpp): every key read only when named, so that an input without them writes the same restart and parameter dump as before |
| vet_col_lat_cut_thin | global_variable::restart_run ? 0.0 : 1.0 | PBG | src/rad_m1/rad_m1.cpp:2403 | vet_col_lat_cut_thin (VetLatCsFinish). DEFAULT 1 on a fresh run; a restart whose file lacks the key keeps 0 (the depth cut alone, bitwise) |
| vet_col_lat_every | named only | - | src/rad_m1/rad_m1.cpp:2389 |  |
| vet_col_lat_init_iter | named only | - | src/rad_m1/rad_m1.cpp:2392 |  |
| vet_col_lat_nmu | named only | - | src/rad_m1/rad_m1.cpp:2383 |  |
| vet_col_lat_npsi | named only | - | src/rad_m1/rad_m1.cpp:2386 |  |
| vet_col_lat_offdiag | named only | PBG | src/rad_m1/rad_m1.cpp:2411 |  |
| vet_col_lat_posfix | named only | G | src/rad_m1/rad_m1.cpp:2429 | vet_col_lat_posfix = drop (default, bitwise) / lag (read only when named): what the positivity trigger does to the lateral D_r,lat / Q operator terms (rad_m1.hpp vlat_posfix_lag) |
| vet_col_lat_taucut | named only | PBG | src/rad_m1/rad_m1.cpp:2395 |  |
| vet_col_lat_taumax | named only | - | src/rad_m1/rad_m1.cpp:2408 |  |
| vet_col_lat_taumin | named only | - | src/rad_m1/rad_m1.cpp:2405 |  |
| vet_col_ncore | 8 | BG | src/rad_m1/rad_m1.cpp:1724 |  |
| vet_col_nmu | 4 | BG | src/rad_m1/rad_m1.cpp:1726 |  |
| vet_col_nsub | 1 | BG | src/rad_m1/rad_m1.cpp:1725 |  |
| vet_col_reflect_amax | named only | - | src/rad_m1/rad_m1.cpp:1755 |  |
| vet_col_reflect_top | !global_variable::restart_run | BG | src/rad_m1/rad_m1.cpp:1745 |  |
| vet_col_source | named only | G | src/rad_m1/rad_m1_vet.cpp:1934 |  |
| vet_col_surface_face | !global_variable::restart_run | BG | src/rad_m1/rad_m1.cpp:1752 |  |
| vet_col_surface_q | !global_variable::restart_run | BG | src/rad_m1/rad_m1.cpp:1761 | vet_col_surface_q: DEFAULT true since m1-defaults2 (tests_m1/runs_5e_vetcol2: T-S4 L1 2.1e-3 -> 5.2e-5 at n = 256, Milne 1.8e-3 -> 1.1e-4). A restart whose file lacks the key keeps false; the resolved value is echoed. |
| vet_col_surface_qmax | 1.0 | BG | src/rad_m1/rad_m1.cpp:1764 |  |
| vet_col_surface_qmin | 1.0e-3 | BG | src/rad_m1/rad_m1.cpp:1763 |  |
| vet_eig_min | named only | - | src/rad_m1/rad_m1_vet.cpp:2033 |  |
| vet_gd | named only | PBG | src/rad_m1/rad_m1.cpp:318 |  |
| vet_gd_allow_clamp | named only | PBG | src/rad_m1/rad_m1.cpp:2107 |  |
| vet_gd_async | named only | - | src/rad_m1/rad_m1.cpp:1938 |  |
| vet_gd_band_colat | named only | - | src/rad_m1/rad_m1.cpp:2112 | (read only when named): the full sphere's band depth uses sin of this colatitude (deg) instead of sin(0) (rad_m1_vetgd_sph2.inc); default 30 |
| vet_gd_band_exit | named only | - | src/rad_m1/rad_m1.cpp:2104 |  |
| vet_gd_band_global | true | - | src/rad_m1/rad_m1.cpp:2140 |  |
| vet_gd_band_ring | named only | - | src/rad_m1/rad_m1.cpp:2206 |  |
| vet_gd_band_ring_refill | named only | - | src/rad_m1/rad_m1.cpp:2214 | vet_gd_band_ring_refill (read only when named, default false: speed over memory, USER; the raw band of clean sp): the outward pass fetches a shell's inward band from its owners again (one more exchange per shell) instead of keeping it in the save store (bitwise) |
| vet_gd_body_global | named only | - | src/rad_m1/rad_m1.cpp:2354 |  |
| vet_gd_cs_list_pad | named only | - | src/rad_m1/rad_m1.cpp:2183 |  |
| vet_gd_cs_xgroups | named only | - | src/rad_m1/rad_m1.cpp:2186 |  |
| vet_gd_fuse_shells | named only | - | src/rad_m1/rad_m1.cpp:2346 |  |
| vet_gd_geo | 3 | HG | src/rad_m1/rad_m1.cpp:1891 |  |
| vet_gd_gl_nmu | 9 | PHBG | src/rad_m1/rad_m1.cpp:1889 | (USER 11-17): the DEFAULT set of a FRESH start that names no set key, recorded for restarts; a restart whose input names none keeps HEALPix nside 3. (USER: <= 120 rays, rotation as in production; RAYSET2): the geodesic set vet_gd_geo = 3 (92 directions), averaged over 4 rotated sets (geodef7 below): |
| vet_gd_gl_nphi | np | PBG | src/rad_m1/rad_m1.cpp:1916 |  |
| vet_gd_gl_stagger | true | PBG | src/rad_m1/rad_m1.cpp:1917 |  |
| vet_gd_halo_compact | 2 | PBG | src/rad_m1/rad_m1.cpp:2163 | prod1104b default (USER 11-04): vet_gd_halo_compact = 2 on a FRESH vet_gd start (recipe::NewDefaults) when not named (bitwise; 29-40 % faster on cs vet_gd runs; the vet_gd setup still switches it off where it does not apply); recorded for restarts; a restart whose input lacks the key keeps 0; an exp |
| vet_gd_halo_compact_nan | named only | - | src/rad_m1/rad_m1.cpp:2323 | vet_gd_halo_compact_nan (read only when named; a gate): the compact exchanges put NaN on the band positions they do not send, so that a read of one shows |
| vet_gd_halo_exact | named only | - | src/rad_m1/rad_m1.cpp:2306 | vet_gd_halo_exact (read only when named; default 2 where it applies: MPI halo with vet_gd_halo_compact > 0 and lists, shell lists, no band_exit / async): true / 1 = lists cut from the masks, 2 = lists from the reads alone (both bitwise in the 1010 gates), false / 0 = the vgdspeed-1009 mask lists |
| vet_gd_halo_ipc | named only | PG | src/rad_m1/rad_m1.cpp:2331 |  |
| vet_gd_halo_ipc_fuse | named only | - | src/rad_m1/rad_m1.cpp:2343 | vet_gd_halo_ipc_fuse (read only when named; default true; CUDA lean IPC): an exchange's pack, flags, wait and scatter in one kernel (bitwise the three) |
| vet_gd_halo_lean | named only | - | src/rad_m1/rad_m1.cpp:2318 |  |
| vet_gd_halo_list_mb | named only | - | src/rad_m1/rad_m1.cpp:2290 | vet_gd_halo_list_mb (read only when named; default 4096; 0 = off): per-rank MB budget of the compact-halo lists (VetGdHlBuild); bitwise either way |
| vet_gd_halo_overlap | named only | - | src/rad_m1/rad_m1.cpp:2315 | vet_gd_halo_ipc (default below); env VGD_IPC |
| vet_gd_halo_pipe | named only | - | src/rad_m1/rad_m1.cpp:2349 |  |
| vet_gd_inner_reflect | named only | - | src/rad_m1/rad_m1.cpp:2277 | vet_gd_inner_reflect = diff / zeroflux (read only when named) |
| vet_gd_iter | named only | - | src/rad_m1/rad_m1.cpp:2357 |  |
| vet_gd_ls | 0 | P | src/rad_m1/rad_m1.cpp:1890 |  |
| vet_gd_nside | named only | - | src/rad_m1/rad_m1.cpp:1892 |  |
| vet_gd_overlap | named only | - | src/rad_m1/rad_m1.cpp:2296 | (read only when named): vet_gd_overlap = G (default 0 = off; bitwise) and its memory cap vet_gd_overlap_mb (default 1024), vet_gd_fuse_shells = H (default 1; H > 1 is a lagged exchange, NOT bitwise); rad_m1_vetgd.cpp VetGdSweep |
| vet_gd_overlap_mb | named only | - | src/rad_m1/rad_m1.cpp:2299 |  |
| vet_gd_pole_ring | named only | - | src/rad_m1/rad_m1.cpp:2124 | (read only when named): the polar ring buffers of the full sphere (exact whole-ring reads in the polar zone; default off: ~3 s/cycle more on the 4-H200 accretor with 8 rows, and the T6a pole test is unchanged by them) and their zone |
| vet_gd_pole_ring_colat | named only | - | src/rad_m1/rad_m1.cpp:2156 |  |
| vet_gd_raw | named only | - | src/rad_m1/rad_m1.cpp:2202 | vet_gd_band_ring (read only when named; default true where it applies, bitwise; VetGdInit): the band intensities lent per shell from a few wide slabs vet_gd_raw (read only when named): the clean sc rays' band on spherical polar keyed on the global cell (2226d2203); false = the exact-list band of pro |
| vet_gd_rebuild_every | named only | - | src/rad_m1/rad_m1.cpp:2273 |  |
| vet_gd_replace | named only | PG | src/rad_m1/rad_m1.cpp:2363 |  |
| vet_gd_rotate_align | true | HRXBG | src/rad_m1/rad_m1_vetgd.cpp:1371 |  |
| vet_gd_rotate_axis | named only | HG | src/rad_m1/rad_m1_vetgd.cpp:1421 | read only when named: the default is not added to the parameter dump |
| vet_gd_rotate_dt | named only | PHBG | src/rad_m1/rad_m1.cpp:1931 |  |
| vet_gd_rotate_every | named only | BG | src/rad_m1/rad_m1.cpp:1934 |  |
| vet_gd_rotate_mirror | named only | - | src/rad_m1/rad_m1_vetgd.cpp:1445 | rotwb2-1121 vet_gd_rotate_mirror (read only when named; checked in VetGdInit) |
| vet_gd_rotate_seed | 12345 | - | src/rad_m1/rad_m1_vetgd.cpp:1439 |  |
| vet_gd_rotate_seq | named only | - | src/rad_m1/rad_m1_vetgd.cpp:1435 | vet_gd_rotate_seq = r3 (default) / random (read only when named) |
| vet_gd_rotate_slope | fresh ? "correct" : "keep" | BG | src/rad_m1/rad_m1_vetgd.cpp:1377 | _slope DEFAULT correct on a FRESH run (Time2RotCorrect; acts only under time_scheme = hesdirk2, inert otherwise); a restart whose file lacks the key keeps keep (bitwise older runs). Echoed, so a later restart keeps it. |
| vet_gd_rotate_slope_frac | 1.0 - 1.0/std::sqrt(2.0) | G | src/rad_m1/rad_m1_vetgd.cpp:1389 |  |
| vet_gd_rotate_slope_reuse | true | G | src/rad_m1/rad_m1_vetgd.cpp:1391 |  |
| vet_gd_rotate_stretch | fresh ? 0.02 : 0.0 | BG | src/rad_m1/rad_m1_vetgd.cpp:1372 |  |
| vet_gd_rotate_t0 | 0.0 | BG | src/rad_m1/rad_m1_vetgd.cpp:1370 |  |
| vet_gd_rotate_xfade | false | HRXBG | src/rad_m1/rad_m1.cpp:1704 |  |
| vet_gd_rotate_xfade_fuse | named only | - | src/rad_m1/rad_m1.cpp:2021 |  |
| vet_gd_rotate_xfade_lean | named only | - | src/rad_m1/rad_m1.cpp:2012 |  |
| vet_gd_rotate_xfade_m | named only | - | src/rad_m1/rad_m1.cpp:2009 | (the combination checks: VetGdXfade, after every key is read) |
| vet_gd_rotate_xfade_pen | named only | - | src/rad_m1/rad_m1.cpp:2015 |  |
| vet_gd_seam_mask | named only | - | src/rad_m1/rad_m1.cpp:2258 |  |
| vet_gd_tan_operator | named only | - | src/rad_m1/rad_m1.cpp:2270 |  |
| vet_gd_tangential | named only | PHRXBG | src/rad_m1/rad_m1.cpp:2372 |  |
| vet_gd_theta_wall | named only | B | src/pgen/ry_per_accretor.cpp:6102 |  |
| vet_gd_thin_decades | named only | - | src/rad_m1/rad_m1.cpp:2261 |  |
| vet_gd_thin_parts | named only | - | src/rad_m1/rad_m1.cpp:2267 |  |
| vet_gd_thin_smooth | named only | - | src/rad_m1/rad_m1.cpp:2264 |  |
| vet_gd_thin_taumin | named only | - | src/rad_m1/rad_m1.cpp:2360 |  |
| vet_gd_twin | named only | PBG | src/rad_m1/rad_m1.cpp:2190 |  |
| vet_gd_wall_interp | named only | PBG | src/rad_m1/rad_m1.cpp:2231 |  |
| vet_interp | global_variable::restart_run ? "log" : "auto" | PBG | src/rad_m1/rad_m1.cpp:2070 |  |
| vet_kappa_es | 0.34 | G | src/rad_m1/rad_m1.cpp:1280 | cm^2/g -> code opacity per unit mass |
| vet_kappa_floor | 1.0e-5 | G | src/rad_m1/rad_m1.cpp:1281 |  |
| vet_mb_agg | named only | - | src/rad_m1/rad_m1_vet.cpp:2487 | vet_mb_agg (exact, default on): one message per neighbour rank and exchange; vet_mb_mom_fuse (exact, default off): the moments summed inside the ray launch; vet_mb_agroup = G (round-off, default 1 = off): the hybrid decomposition, rays split over the G consecutive ranks of a group (one node), space  |
| vet_mb_agroup | named only | - | src/rad_m1/rad_m1_vet.cpp:2493 |  |
| vet_mb_angles | named only | - | src/rad_m1/rad_m1_vet.cpp:3057 | PROTOTYPE vet_mb_angles: the angle-decomposed sweep (VetSweepAng) |
| vet_mb_force | named only | - | src/rad_m1/rad_m1_vet.cpp:1847 |  |
| vet_mb_halo | named only | - | src/rad_m1/rad_m1_vet.cpp:2586 | scaling options (read only when given, so that an input without them writes the restart file it always wrote): vet_mb_halo = K exchanges a band K w wide every K layers and recomputes the overlap redundantly (same arithmetic, same numbers); vet_mb_overlap sweeps the interior of the next layer while t |
| vet_mb_kernel | named only | - | src/rad_m1/rad_m1_vet.cpp:2592 |  |
| vet_mb_lag | named only | - | src/rad_m1/rad_m1_vet.cpp:2652 |  |
| vet_mb_mom_batch | named only | - | src/rad_m1/rad_m1_vet.cpp:2669 |  |
| vet_mb_mom_fuse | named only | - | src/rad_m1/rad_m1_vet.cpp:2490 |  |
| vet_mb_overlap | named only | - | src/rad_m1/rad_m1_vet.cpp:2589 |  |
| vet_mb_tblock | named only | - | src/rad_m1/rad_m1_vet.cpp:2626 | vet_mb_tblock = B (exact, default 1 = off; tests_m1/runs_4i_sctb/README.md): up to B launches per kernel launch, never across a band exchange, so B <= vet_mb_halo matters; without vet_mb_halo in the input, the default halo becomes max(3, B). vet_mb_tblock_tj / _tk: tile of the plane per team (0 = wh |
| vet_mb_tblock_direct | named only | - | src/rad_m1/rad_m1_vet.cpp:2641 |  |
| vet_mb_tblock_ovl | named only | - | src/rad_m1/rad_m1_vet.cpp:2638 |  |
| vet_mb_tblock_team | named only | - | src/rad_m1/rad_m1_vet.cpp:2635 |  |
| vet_mb_tblock_tj | named only | - | src/rad_m1/rad_m1_vet.cpp:2629 |  |
| vet_mb_tblock_tk | named only | - | src/rad_m1/rad_m1_vet.cpp:2632 |  |
| vet_milne | false | - | src/rad_m1/rad_m1_vet.cpp:1874 |  |
| vet_nmu | 4 | B | src/rad_m1/rad_m1_vet.cpp:1872 |  |
| vet_nphi | 8 | B | src/rad_m1/rad_m1_vet.cpp:1873 |  |
| vet_path | "none" | G | src/rad_m1/rad_m1.cpp:2091 |  |
| vet_ring | named only | B | src/rad_m1/rad_m1_vet.cpp:2130 |  |
| vet_ring_groups | 1 | B | src/rad_m1/rad_m1_vet.cpp:2132 |  |
| vet_rotate_advance | true | HRXBG | src/rad_m1/rad_m1_vetgd.cpp:1309 | vet_rotate_advance (read only when named; default false = bitwise). (USER 11-05): DEFAULT true on a FRESH vet_gd run (recipe::NewDefaults: fresh start or <job>/recipe = prod1019) when the key is not named; recorded, so the run's restarts keep it; a restart whose input lacks the key keeps false (logg |
| vet_rotate_avg | 4 | PHBG | src/rad_m1/rad_m1.cpp:1930 |  |
| vet_rotate_avg_ema | 128.0 | PHRG | src/rad_m1/rad_m1.cpp:3157 |  |
| vet_rotate_avg_ema_change | named only | - | src/rad_m1/rad_m1.cpp:3169 |  |
| vet_rotate_avg_ema_change_form | "linear" | HRG | src/rad_m1/rad_m1.cpp:3178 | (linear, the fresh default again since USER 11-21: step costs the ring's time order, E/F p 1.66/1.39 vs 1.95/1.85 on deck D; a restart keeps its recorded form) |
| vet_rotate_avg_ema_change_j | named only | - | src/rad_m1/rad_m1.cpp:3172 |  |
| vet_rotate_avg_ema_dema | named only | - | src/rad_m1/rad_m1.cpp:3192 |  |
| vet_rotate_avg_ema_diag | named only | - | src/rad_m1/rad_m1.cpp:3196 |  |
| vet_rotate_avg_ema_lc | named only | - | src/rad_m1/rad_m1.cpp:3206 |  |
| vet_rotate_avg_ema_tau | named only | - | src/rad_m1/rad_m1.cpp:3203 |  |
| vet_rotate_avg_ema_wb | named only | - | src/rad_m1/rad_m1.cpp:3161 | vet_rotate_avg_ema_wb acts on the EMA ring only (no silent no-op) |
| vet_rotate_avg_weight | named only | - | src/rad_m1/rad_m1.cpp:3129 | vet_rotate_avg_weight = j (default: the sets weighted by their formal-solution J, i.e. the normalised moments of the union of the K sets) / equal |
| vet_sc_kernel | named only | - | src/rad_m1/rad_m1_vet.cpp:2193 |  |
| vet_sc_set | named only | - | src/rad_m1/rad_m1_vet.cpp:2144 | vet_sc_set = gd (read only when named): the vet_gd direction set (the product vet_gd_gl_nmu x vet_gd_gl_nphi [_stagger], or vet_gd_ls / vet_gd_geo) with its pole on x1. On a Cartesian block this x1-layer sweep is vet_gd's shell sweep (the layers are its shells) and vet_gd's z rotation turns the set  |
| vet_scatter | true | PBG | src/rad_m1/rad_m1.cpp:1221 |  |
| vet_scatter_form | "ma" | PBG | src/rad_m1/rad_m1.cpp:1227 |  |
| vet_scatter_j | "relaxed" | PBG | src/rad_m1/rad_m1.cpp:1239 | vet_scatter_j default = relaxed (J = E*, the end-of-step pair of B*). |
| vet_scatter_kappa_e | (fresh && gen_eos_tab) ? "eos" : "const" | PG | src/rad_m1/rad_m1.cpp:1267 | eos where a tabulated general EOS supplies x_e (fresh runs) |
| vet_seg | sgd | PBG | src/rad_m1/rad_m1.cpp:1828 |  |
| vet_seg_cross | false | - | src/rad_m1/rad_m1_vet.cpp:2008 | vet_seg_cross (default false): where the ray crosses a row or layer boundary of the upwind column before the x1 face and that cell is a material edge to the stencil, the last part of the upwind half-segment takes the chi and S of the cell it crosses (VetSegCross; with vet_interp != lin) |
| vet_source | nd ? "relaxed_esc" : "gas" | PBG | src/rad_m1/rad_m1_vet.cpp:1933 | (USER 10-20): DEFAULT relaxed_esc under the new default set (utils/recipe.hpp: a fresh run, or <job>/recipe = prod1019; it equals gas on the non-LTE tests and keeps relaxed where the exchange is stiff); a restart whose input lacks the key, or recipe = none, keeps the old vet_sc default gas vet_col_s |
| vet_source_noesrc | fresh | BG | src/rad_m1/rad_m1.cpp:1711 | vet_source_noesrc for vet_sc too (rad_m1_vet.cpp). (#11): DEFAULT true on a fresh run, as on vet_col / vet_gd (the VET source of a backward-Euler solve reads the physical E^n, not E^n + dt (chat/c) esrc); a restart whose file lacks the key keeps false (echoed). Bitwise wherever esrc = 0. |
| vet_source_noesrc_h2be | named only | P | src/rad_m1/rad_m1.cpp:2228 |  |
| vet_tensor | global_variable::restart_run ? "uniaxial" : "full" | B | src/rad_m1/rad_m1_vet.cpp:2025 |  |
| vet_x1_npass | 64 | - | src/rad_m1/rad_m1_vet.cpp:2051 | NOT a slip (rad_m1.hpp vet_mirror_npass): vet_x1_npass is the one pass-count key of both multi-pass x1 sweeps, the periodic one (read again below) and the cavity mirrored at both x1 ends (vet_bc_reflect = mirror, default 4, >= 2) |
| vet_x1_periodic | named only | - | src/rad_m1/rad_m1_vet.cpp:2063 |  |
| vet_x1_periodic_inflow | global_variable::restart_run ? "boundary" : "source" | - | src/rad_m1/rad_m1_vet.cpp:2111 | vet_x1_periodic_inflow: the first pass's entry intensity. source : BOTH directions enter with the face source S (the intensity of an infinite uniform medium): a uniform state gives I = S on every ray and H = 0 exactly, whatever npass (default on fresh runs) boundary : the ordinary boundary intensiti |
| vet_x1_ptol | fixed ? 0.0 : 1.0e-6 | - | src/rad_m1/rad_m1_vet.cpp:2092 |  |
| zalesak_active | true | PHRXBG | src/rad_m1/rad_m1_init.cpp:2583 |  |
| zalesak_cap_rounds | named only | G | src/rad_m1/rad_m1_init.cpp:2575 | zalesak_cap_rounds (inert without implicit_face_mmat; read only when named) |
| zalesak_colsweep | named only | - | src/rad_m1/rad_m1_init.cpp:2578 |  |
| zalesak_halo_skip | named only | - | src/rad_m1/rad_m1_init.cpp:2568 | (rad_m1.hpp): halo exchange of the Zalesak R only when it can change |
| zalesak_local_rounds | named only | - | src/rad_m1/rad_m1_init.cpp:2571 |  |
