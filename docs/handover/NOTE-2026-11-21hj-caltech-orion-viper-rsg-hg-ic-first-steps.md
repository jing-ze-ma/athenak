# NOTE 2026-11-21hj Caltech -> orion + viper: RSG / HG production IC + first steps (answer to c70607cdb)

Scope: the Caltech RSG E (lineage `rsg_E_1020/prodH/rst/rsg.01295.rst`) and HG global (lineage
`hgglob_1021/prodG/rst/hegiant.00139.rst`) productions. Everything below is quoted from decks and logs under
`/resnick/groups/carnegie_poc/jingze/`, unless it is marked **inference**.

**Main point first:** neither run used the clean scheme or `clean_init`. Both are OLD implicit-M1 runs
(`transport = implicit`, `closure = vet_col` + `vet_gd`, `time_scheme = hesdirk2`, `force_reference = wb_arad`).
The keys `clean_init`, `clean_init_dtfac` and "steady start" do not appear in either deck or in either first log. The start
was the in-code he_star_m1 IC chain (sph1d + re_atm + atm_force module + he_ic_balance + he_ic_equil) followed by
`ic_rad = fs_relax` (a radiation-only relaxation), and the run then went straight to the full time step. One fact may
help with your |F|/cE problem: in our IC the sph1d solve already puts F/(cE) < 1 everywhere above the
photosphere, 0.25 at R_ph rising to 0.97 at 2.9 R_ph (log lines below).

## Lineage (which start the restarts come from)
- **RSG:** `prodH` is a restart of `prodG/rst_keep/rsg.00335.rst` (1675 d), made on 10-29 with vtxfix2 (see runH.cfg comment).
  `prodG` was a FRESH start on 10-25 07:58 PDT from `rsg_E_1020/rsg_E_1025_fresh.athinput`, binary
  `builds_rstfixx/bin/athena_gpu_he_star_m1_958cac78_wkfix2_nofma` (fin51-r-fix 958cac78 + kick581 work fix, rsgfresh-1025
  c60714bd, md5 140d0e79...). The fresh-start log is `rsg_E_1020/prodG/run.1025_173216.log`; it runs from t = 0 to 311 d.
  Earlier fresh starts (prod63 on 10-21, prodF/G on 10-22) were abandoned. See item 3.
- **HG:** `prodG` is one continuous lineage from a FRESH start on 10-22 15:32 PDT from `hgglob_1021/hgglob_1022_n2seed.athinput`,
  binary `builds_rstfixx/bin/athena_gpu_he_star_m1_d01ce628_nofma` (fin51-r-fix d01ce628b). The command line is in
  `prodG/inalloc.out`: `-i ...n2seed.athinput output2/dt=28800 output3/dt=28800 ...`. The first log is
  `hgglob_1021/prodG/run.1022_160216.log`. Binary and key switches were made later, always at a restart (runG.cfg.pre_* backups).

## 1. How the production IC was built

| | RSG E | HG global |
|---|---|---|
| Source | MESA r15140 golden25 profile438 (Orion m22_p438.data; 20.54 Msun, L 1.79e5 Lsun, R_ph 843.8 Rsun, Teff 4088 K) | MESA r15140 14p5M_widebinary_profile14 (S. Wu; He star 2.754 Msun, R_ph 109.85 Rsun, Teff 7.64 kK) |
| IC column | `files/rsg-E-1020/ic_rsg_g25_deep.txt`, header: `tools/ic_star make_star_ic.py, mode mesa, scheme mesa; he_ic_cols = 5 ... alpha_MLT 3, MLT++` | `files/hg145-prod-1018/ic_hegiant_14p5M_mk.txt`, header: `tools/ic_star make_star_ic.py, mode mesa, scheme mesa; he_ic_cols = 5 ... alpha_MLT 2, MLT++` |
| Columns | r, rho, eint, F_r, E, T, F_MLT/F. MESA T(P) re-integrated in HSE with OUR general EOS table, OUR MESA-kap kappa_R and G m(r) | same |
| Opacity | `rosseland_mesakap_4c0e8b64d26a.txt`, `planck_ext2_x0.653848_z0.019757_...` | `rosseland_mesakap_afe2dd8da2bd.txt`, `planck_ext2_x0.0009496_z0.0215304_...` |

- **No separate 1-D relaxation run under our code** (none found in either run dir; the tooling README
  `bsgsp_1103/ic_star_blend/README.md` describes only the python column builder plus the in-code chain). The only "relaxation"
  is the in-code chain at setup (below) followed by the 3-D run itself.
- In-code IC chain. Both decks use these keys, identical except for the radii:
  `he_ic_rad = sph1d`, `he_ic_rad_gas = re_atm`, `he_ic_rad_hse = true`, `he_ic_rad_atm_force = module`, `he_ic_rad_atm_post = false`,
  `he_ic_rad_taua = 1`, `he_ic_rad_taub = 3`, `he_ic_rad_dtau = 0.05`, `he_ic_rad_ncore = 64`, `he_ic_rad_maxit = 60`;
  `he_ic_eint_from_t = true`; `he_ic_balance = true` (`he_ic_balance_rmax` 5.119e13 RSG / 7.021e12 HG);
  `he_ic_equil = true`, `he_ic_equil_dtfac = 1.0`, `he_ic_equil_relax_n = 40`, `he_ic_equil_relax_tol = 1.0e-4`, `he_ic_equil_maxit = 1`,
  `he_ic_equil_tol = 1.0e-9`, `he_ic_equil_last_march = true`; `he_gravity = wb`, `he_gm_column = true`; `mlt_flux_frozen = true`,
  `mlt_closure = adaptive` (`mlt_relax_time` 1.5e4 s RSG / 1500 s HG); `he_bc_inner = wall`, `he_wall_zero_flux = true`, `he_wall_closed = true`,
  `implicit_bc_x1min = flux`, `implicit_bc_x1max = marshak`.
- In `<rad_m1>`: `ic_rad = fs_relax`, `ic_rad_fs_iter = 3`, `ic_rad_relax_n = 40`, `ic_rad_relax_tol = 1e-4`, `ic_rad_relax_dtfac = 1`.
- Atmosphere and floors:
  - RSG: `dfloor = 1e-20`, `pfloor = 1e-13`, `tfloor_kelvin = 2000`, `vceil = 1e8`, `dfloor_keep_velocity = true`, `he_bc_outer = outflow`,
    sponge `he_sponge_rate = 1e-4` above `he_sponge_r0 = 1.585e14` (0.9 r_top), `he_sponge_mode = radial`, `he_sponge_dmax = 2e-17`.
    The sph1d log line reads `floor 1e-17 from r = 6.29676e+13` (1.085 R_ph).
  - HG: `dfloor = 1e-18`, `pfloor = 1e-11`, `tfloor_kelvin = 2400`, `vceil = 1e8`, sponge rate 1e-3 above 2.2265e13, `he_sponge_dmax = 2e-14`.
  - Both: `wellbalance_dynamic = true`, `wb_x1 = true`, `wb_option = polytropic`, `wb_cache_every = 0` (`wb_rmax` 6.393e13 RSG / 8.768e12 HG).
- What the setup log shows (RSG, `run.1025_173216.log`):
  - `he_ic_rad_atm_force = module: ... 108 cells above the floor, 200 passes ...` followed by `WARNING not converged in he_ic_rad_atm_maxit = 200`
    (march residual 0.998). It ran with that warning.
  - `he_ic_rad = sph1d (gas re_atm, hse on): r_b = 5.77891e+13 (tau 3.32704) ... g_top = H/J 0.977657, f_top = K/J 0.956011`
  - F/(cE) at r/R_ph = 1, 1.01, 1.1, 1.5, 2.9: 0.2495, 0.581, 0.737, 0.885, 0.972. **No |F|/cE > 1 in the floor gas at t = 0.**
  - `mlt_flux_frozen: max F_MLT/F on the faces = 0.998992`, so the frozen MLT scaffold carries the convective flux at t = 0.
  - `ic_rad`: 3 fs passes, then 40 relax steps at `dt = 1.000000e+03`. L/L_IC went 0.810 -> 0.9267; it hit the 40-step cap
    (|dL/L| 2.3e-4 > tol). Then `force pass (dic-1013)`.
  - `he_ic_equil round 1: module residual max |a_net|/g = 4.946e-02 ... max |L(face)/L_in - 1| = 8.760e-02`, then
    `round 1 update: ... cells 0 .. 428 re-marched ... max |drho/rho| 9.997e-01 at r 6.28635e+13` (the atmosphere just below the floor rho was rebuilt).
- HG (`run.1022_160216.log`):
  - atm_force converged in 6 passes.
  - sph1d `r_b = 7.92064e+12 (tau 2.97618)`. F/(cE) at r/R_ph = 1, 1.01, 1.1, 1.5, 2.9: 0.250, 0.591, 0.736, 0.883, 0.971.
  - ic_rad converged after 14 relax steps at `dt = 1.25e+02`, L/L_IC = 1.187.
  - he_ic_equil round 1: `|a_net|/g = 1.569e-01`, re-march `max |drho/rho| 4.688e-01 at r 8.27598e+12`.

## 2. The first steps of the 3-D run
- **No steady start and no clean_init.** After `ic_rad` and `he_ic_equil` (both at setup, cycle 0), the run takes normal steps.
  The log reads `time_scheme=hesdirk2 (g=0.292893; first step, restarts without the slope and failed stages take backward Euler)`,
  so the first step is backward Euler and the following steps use hesdirk2.
- The first-start fixes are in the decks (an `-i` start rejects command-line keys that the deck lacks):
  `time2_vet_col = extrap2`, `time2_accept_fallback = true`, `time2_posfix = 4`, `time2_reject = retry`, `time2_start = be`,
  `implicit_cs2_clip = rad`, `implicit_face_opac = harmonic`. RSG also has `implicit_work_fix = true` and `dt_retry_t2_be = 0`.
- **dt / cfl:** `cfl_number = 0.3` from cycle 0, with no cfl ramp. There is no dt ramp key; `dt_retry = true`
  (`dt_retry_fac_t2 = 0.9`, `dt_retry_grow = 1.1`, `dt_retry_ngood = 20`).
  - RSG: dt = 1.000e3 s from cycle 0 and at every logged cycle to 2.69e7 s (cycle 27000).
  - HG: dt = 1.25e2 s from cycle 0 (later ~117-122 s, CFL-limited).
  - We did not trace which limiter pins RSG at exactly 1e3 s.
- **L_in:** no ramp. `implicit_flux_x1min` is fixed at L_IC/(4 pi r_in^2) (RSG `F_in = 1.75446e+13`, HG `3.60602e+13`) from t = 0.
- **MLT scaffold ramp (the only ramp):** `mlt_ramp_start/_time` are deleted from the decks, so the auto rule v2 applies (binary >= d01ce628).
  - RSG: `mlt_ramp auto ... start = 0 s (0 d), time = 9.09149e+06 s (105.226 d)`.
  - HG: `start = 0, time = 1.15275e+06 s (13.342 d)`.
  - `w(t = 0) = 1`: the frozen F_MLT carries the convective flux at the start and is ramped off over t_sat, while the adaptive closure takes over.
- **Damping/sponge:** only the radial outer sponge above (RSG rate 1e-4 above 0.9 r_top; HG rate 1e-3). There is no early-time extra damping.
- **WB reference:** well-balancing is on from cycle 0. a_ref and Phi_eff are captured at the fresh start's setup by he_ic_equil, after the
  re-march, and are stored in the restart trailer. The restart log reads
  `he_ic_equil: a_ref (cells is .. ie) and Phi_eff from the fresh start's equilibrium (restart trailer)`.
- **Seeds:** temperature seed (eint and E), `he_seed = 1e-3`, `he_seed_rad = true`, real Y_lm in the N^2 band.
  - RSG: `he_seed_rlo/rhi` = 5.626e13 / 5.794e13 (0.97-0.999 R_ph), l 16..32, `he_seed_rng = 1234`.
  - HG: band from N^2 automatically (r/r_int 0.977-0.999), l 16..64 (auto window), `he_seed_rng = 22102`.
  - The seed is applied after he_ic_equil (`deferred he_seed applied`).
- **Clean start in the logs:**
  - HG first chunk (357 cycles): 0 E-floor clips, 0 retries, 0 vceil.
  - RSG fresh log: the first `stage 1 NOT POSITIVE` comes at cycle 15245 (~176 d) and the first `stage 2 NOT ADMISSIBLE` at cycle 16289.
    No floor or vceil trouble is logged before that.

## 3. Early-time trouble and fixes (decks / binaries)
From the Caltech session log (memory notes) and the runG/runH.cfg comments:
- **Photosphere runaway (10-21, earlier fresh starts):** a converged discrete system (one-cell MLT deposit plus an unresolved kappa cliff).
  More opac_lag passes or full convergence did not cure it. The fix was `implicit_face_opac = harmonic` (harm-1022 d3df1401); productions were
  restarted fresh with it.
- **Floor-gas kick (10-20):** `implicit_cs2_clip = rad` (cs2clip-1020) keeps the limiter clip out of the gas momentum budget.
  It has been in the clip decks since 10-20.
- **Photosphere clumping, then retries and dt collapse:** RSG had retries at 105-107 d and a dt collapse at 129.8 d; HG had a collapse at 18.8 d.
  The clumps were judged physical (RT).
  - Cause: a stage was rejected after an operator-positivity fallback.
  - Fixes: `time2_accept_fallback = true` plus `time2_posfix = 4`. Arm E2 passed the 129.8 d collapse (min dt 265 s vs 2-5 s).
    `time2_vet_col = extrap2` alone did not help.
  - Productions were restarted FRESH with the N^2 surface seed (10-22) and these keys.
- **RSG 230 d dt collapse (10-23, 340 -> 2 s):** low-density gaps between photospheric clumps. Cured by `time/dt_retry_t2_be = 0` (plus the csfix2 binary).
- **RSG kick581 (562-581 d corner eruption):** energy was created by the M1 write-back bookkeeping (hesdirk2 E row). Fix: `implicit_work_fix = true`.
  The current RSG lineage is the 10-25 FRESH restart with that fix in the binary and the deck.
- **RSG cube-vertex ejection (1691-1704 d):** corner cells on the M1 E floor; the positivity floor created 5.4e42 erg, and the gas at T floor heated to 39-48 kK.
  Fixes: `cs_vertex_mono = true` plus `implicit_pos_floor_cons = true` (vtxfix2 e04fbd50). prodH was restarted from 1675 d with them.
  HG got the same keys at its rst 00044 (130 d).
- **HG later (~114-130 d):** seam/vertex-born cold gas leads to M1 E<0 retries and BE redos. Handled with restarts and key changes
  (be was tried and reverted; cfl 0.27 perturbation on 11-05). This is not an early-time issue.
- Binaries in order:
  - RSG: 958cac78_wkfix2 (fresh) -> 6fac5c7f (10-25 18:55, at rst 62) -> vtxfix2 e04fbd50 (prodH).
  - HG: d01ce628 (fresh) -> d01ce628_wkfix (10-25) -> 6fac5c7f -> vtxfix2 e04fbd50 (10-29) -> later switches listed in `hgglob_1021/runG.cfg.pre_*`.

**Inference for orion:** our floor gas never saw |F| > cE at the start, because:
- E and F above the photosphere come from the 1-D spherical VET solve (sph1d, f_top 0.96, g_top 0.98) rather than from an implicit steady
  solve with a huge pseudo-dt;
- `ic_rad = fs_relax` then relaxes E at the physical dt (1e3 s RSG, 125 s HG; `relax_dtfac = 1`, at most 40 steps);
- he_ic_equil re-marches HSE after that.

Matching our start in the clean scheme would mean a pseudo-dt of order the CFL dt (dtfac ~1), not 1e6, and/or
initialising E,F from the sph1d profile.
