# TASK 2026-11-08 viper -> ORION: rerun the MESA model behind our He giant (HG) and give L_in(t) to core collapse

**Purpose (user):** a pre-SN HG simulation that ramps the luminosity injected at our INNER BOUNDARY (not the surface)
to mimic the late-stage evolution. We need L at the inner boundary as a function of time to core collapse,
from the state of our IC (profile14) onward. Expected timescale: years to decades.
CPU only; resources as Orion decides. Answer with `docs/handover/NOTE-2026-11-xx-orion-hg-mesa-lin.md` on this branch
(xx = the day you finish) plus the small files listed under "Deliverables".

## 1. What our HG is (facts, checked on viper 11-08)
- IC = MESA r15140 profile `14p5M_widebinary_profile14.data` from S. Wu (md5 f4cb70aa30ef4c951b9387803bd939ed,
  written 2025-11-24 with gfortran 10.2.0 / SDK 21.4.1). Compressed copy here:
  `docs/handover/hg-mesa-1108/14p5M_widebinary_profile14.data.xz`.
  Header and derived values: `docs/handover/hg-mesa-1108/profile14_crosscheck.txt`. Key values:
  model_number 500, initial_mass = star_mass = 2.754 Msun (star_mdot 0), Z 0.02, star_age 2.2696916698e6 yr
  (age counted from the He-star start, initial_age = 0 in the inlists), Teff 7644.5 K, photosphere_L 3.71274e4 Lsun,
  R_ph 109.850 Rsun, center he4 1.49e-16, c12 1.90e-3, o16 0.410, ne20 0.515, log Tc 8.867, log rho_c 6.790,
  center_eta 5.96, he_core 2.754, c_core 1.5256, o_core 0.5959, star_mass_h1 7.27e-5 Msun, surface X 9.80e-4,
  power_nuc_burn 8.27e5 Lsun, power_he_burn 3.00e4 Lsun, power_neu 2.28e6 Lsun, kh_timescale 43.8 yr,
  max tri_alfa at r 0.0412 Rsun / m 1.566 Msun,
  peak burn_c at m 1.538; c12_c12/burn_ne/burn_o peak OFF-CENTER at m 0.898 Msun (r 0.0102 Rsun).
- **Inner boundary: r_in = 8.0 Rsun (x1min = 5.5656e11 cm), NOT 3 Rsun.** The "r 3-200 Rsun" in
  HANDOVER-2026-11-07 section 2.2 is template history: the production decks (`hg145-prod-1018/hg145_prod.athinput.in`,
  `hgglob-deltaai-1022/hgglob_1022_n2seed.athinput`) have x1min 5.5656e11, x1max 2.473909e13 (355.6 Rsun),
  `implicit_flux_x1min = 3.606021442400e13` -> L_in = 4 pi r_in^2 F = 1.40366e38 erg/s = 3.667e4 Lsun.
  IC spec (`hg145-prod-1018/ic_hegiant_14p5M_mk.txt` header): `r_in.rule explicit, r_rsun 8.0, bottom_rsun 6.0`.
  We ask for BOTH 3 and 8 Rsun below (3 Rsun because the brief named it; 8 Rsun is the real boundary).
- In profile14: m(3 Rsun) = 2.627060 Msun, L(3 Rsun) = 3.7793e4 Lsun = 1.44673e38 erg/s;
  m(8 Rsun) = 2.719887 Msun, L(8 Rsun) = 3.7309e4 Lsun = 1.42818e38 erg/s (the deck's L_in is 1.7 % lower:
  it comes from our own HSE re-integration). The deck will be ramped by the RATIO L(t)/L(profile14 match), so
  report the ratio too.
- L(r) in profile14 is NOT flat below 1 Rsun: 1.51e4 Lsun at 0.02 Rsun, 5.49e4 at 0.041 (He shell), 5.93e4 at 0.05,
  4.88e4 at 0.1, 4.29e4 at 0.3, 4.00e4 at 1, 3.78e4 at 3, 3.73e4 at 8, 3.71e4 at >= 15 Rsun (envelope absorbing).

## 2. Which MESA run produced profile14 (viper's reading; please confirm)
Source paper: Wu & Fuller 2022 (arXiv:2210.10187), MESA r15140, single stripped stars first (H envelope of
M_ZAMS 13.8-15 Msun removed after core H burning, no winds, Z 0.02), then binaries with a 1.4 Msun NS at
P_orb 1, 10, 100 d. Zenodo 10.5281/zenodo.7106182, model dir `14500M_dir` (copied here as text, see section 3).
- The Zenodo dir is a MESA **binary** setup (`./binary`, inlist_project: m1 2.754, m2 1.4, Kolb, beta 1);
  the top-level inlist_project is the Porb100 copy; LOGS_Porb{1,10,100}_Kolb hold each run's inlists + run script.
- **profile14 is NOT from the Porb100 binary with mass transfer:** at P = 100 d, a = 146 Rsun and the Eggleton
  Roche lobe of the 2.754 Msun star is ~64 Rsun, but profile14 has R = 110 Rsun with star_mass exactly 2.754000 and
  star_mdot 0. Our older HG IC ("M2pt754_Porb100 profile29") had M 2.556 Msun, R 61.8 Rsun, i.e. it HAD lost mass by RLOF.
  So "widebinary" = no interaction = the **single stripped star**. Run it as a single star (star module only).
  If you find otherwise (e.g. S. Wu's run had MT off), say so and run that.
- Where profile14 sits: star_mass_ne20 = 0.602 > 0.55 Msun, the `inlist_to_end_core_c_burn` stop
  (star_species_mass_max_limit 0.55 on ne20), so profile14 is in the `to_lgTmax` stage (pt1 stops at model 760,
  pt2 continues with mesh_delta_coeff 1.75; stop at log_max_temp 9.55), followed by `to_cc` (fe_core_infall_limit 1e7).
  Note model numbering: profile14 is model 500; with profile_interval 50 that suggests extra profiles at stage ends.

## 3. Inputs (all on this branch, `docs/handover/hg-mesa-1108/`)
- `14500M_dir/`: inlists (`inlist_common`, `inlist_mass_Z_wind_rotation`, `inlist_to_start_core_he_burn`,
  `inlist_to_end_core_he_burn`, `inlist_to_end_core_c_burn`, `inlist_to_lgTmax`, `inlist_to_cc`, their `*_header` files,
  `inlist`, `inlist1`, `inlist_project`, pgstar inlists), `LOGS_Porb{1,10,100}_Kolb/` (per-run inlists + `02754M_run.sh`),
  `src/` (run_star_extras.f90 is the 41-line stub = standard_run_star_extras.inc; run_binary_extras.f90), `make/makefile`,
  `history_columns.list`, `mk rn re clean`. The 53 MB `binary` executable and `.smod` files are NOT copied (rebuild).
  Full tarball if needed: viper `/viper/ptmp2/jinma/wu2022_zenodo/stripped_binary_models.tar.gz`.
- `14p5M_widebinary_profile14.data.xz`, `profile14_crosscheck.txt` (section 1 values).
- `extra_history_hglin.f90`: extra history columns (L at r = 3, 8 Rsun; L and r at m_in3 = 2.627060 and
  m_in8 = 2.719887 Msun; r, m, L at the max triple-alpha cell; r, m, L at the He-envelope base X(he4) = 0.5;
  log eps_3a max; L_nu). Not compiled on viper: fix names if r15140 differs and say what you changed.
- `plot_hg_lin.py`: plots from the CSV (column list in its docstring).
- Base physics in `inlist_common`: net approx21_cr60_plus_co56, a09 opacities, T_tau atmosphere, Ledoux,
  history_interval 1, photo/profile_interval 50, Dutch wind 0, Z 0.02.

## 4. The missing starting model (the main obstacle)
All stages load `2.754M_from_14.5M.mod` (the stripped 14.5 Msun star at core H exhaustion), which is **not in Zenodo**
(no .mod or photos in the whole tarball). Do this:
1. First choice, if you have it on Orion or the user can get it from S. Wu: use the original `.mod` (or her
   `before_core_he_burn.mod` / photos). Ask the user in the NOTE if you stop here.
2. Otherwise rebuild it: 14.5 Msun, Z 0.02 (same `inlist_common` physics, net, opacities; Dutch 0; overshoot as in
   `inlist_to_start_core_he_burn`), pre-MS -> core H exhaustion (center h1 < 1e-4), then remove the H envelope down to
   M = 2.754 Msun (r15140 star_job: `relax_mass` / `new_mass = 2.754` / `lg_max_abs_mdot`, or the
   `remove_surface_*`/He-core-boundary options: check `$MESA_DIR/star/defaults/star_job.defaults`), save as
   `2.754M_from_14.5M.mod`. Match the residual H of profile14 (star_mass_h1 7.3e-5 Msun, surface X ~1e-3):
   the residual H changes the late radius strongly, so report it and, if it is far off, adjust the cut.
   Report the TAMS He-core mass of your 14.5 Msun star (Wu & Fuller map M_ZAMS 14.5 -> M_He 2.754).

## 5. Run (single star, MESA r15140, the gfortran SDK you used for golden25: NOTE-2026-10-19-orion-rsg-g25-setup.md)
- Work dir from `$MESA_DIR/star/work`; put the extra history columns in (`extra_history_hglin.f90` header says how);
  `./mk`. Copy the 14500M_dir inlists in. Star mode reads `inlist`: for each stage copy the stage's `*_header` to
  `inlist` (the header -> inlist_common + inlist_mass_Z_wind_rotation + stage inlist chain). Drop the
  `inlist_pgstar_binary` line if star mode refuses it; set `pgstar_flag = .false.`.
- Sequence as `LOGS_Porb100_Kolb/02754M_run.sh` (= the Porb100 inlists, which are also the top-level ones except
  the top-level `inlist_to_lgTmax` debug stop at model 1255: do NOT use that stop):
  `to_start_core_he_burn` (optional; the run script skips it) -> `to_end_core_he_burn` -> `to_end_core_c_burn` ->
  `to_lgTmax` pt1 (max_model 760, mesh 1.5) -> pt2 (mesh_delta_coeff 1.75, min_dq_for_xa 1e-14) -> `to_cc`.
  Note the Zenodo run script uses `./re` (restart from the latest photo with the next stage's inlist1), so the stages
  chain by photos, not by the saved .mod names (`to_end_core_he_burn` and `to_start_core_he_burn` both load
  2.754M_from_14.5M.mod; `to_cc` loads lgTmax.mod). Chain them the same way (`./re` / `star -r` from the last photo)
  or by .mod files, and say which.
- Output: history every step (history_interval 1, extra columns on). Profiles: interval 50 up to C exhaustion,
  **20 from the end of core C burning onward (to_lgTmax), 1 in to_cc**, and photo_interval 10 from to_lgTmax on so
  any window can be re-run with denser profiles. After the first pass, when t_cc is known: re-run the last ~10 yr
  from the nearest photo with profile_interval 1 (or, if the extra history columns already give every quantity, say
  so and skip it). Keep profiles on Orion; push only the CSV/plots.
- Continue through C-shell burning, contraction, the O/Ne re-expansion, Si burning, to core collapse if the inlists
  allow; if a stage stalls (off-center Ne/O flames in a degenerate core: center_eta 6 at profile14), report the last
  age and model, and the cause; do not tune physics silently (list any change).

## 6. Cross-check (before the table)
Find the model of your run that best reproduces profile14: L_ph, R_ph, Teff, star_age, center he4/c12/o16/ne20,
c_core/o_core masses, star_mass_h1, m(3 Rsun), m(8 Rsun), L(3/8 Rsun), r of max tri_alfa.
Report the match quality (table) and that model's time to Si ignition and to collapse. Define Si ignition
explicitly (e.g. log_max_temp 9.55 stop = end of to_lgTmax, and/or center si28 first dropping by 1 % after O
exhaustion; give both if they differ).

## 7. Deliverables
- `docs/handover/hg-mesa-1108/hg_lin.csv`: one row per model from the profile14 match to the end, columns exactly as in
  `plot_hg_lin.py`: model, star_age_yr, t_cc_yr, t_si_yr, L_surf, R_surf, Teff, L_r3, L_r8, L_m3, r_m3, L_m8, r_m8,
  L_He3a, r_He3a, L_Hebase, r_Hebase, log_Tc (Lsun, Rsun, K). Thin it to <= ~2000 rows (keep every row in the last
  10 yr and around fast changes) if it is large; keep the full one on Orion and give its path.
- `hg_lin.png` from `plot_hg_lin.py`, plus `history.data` and the inlists actually used, if small (< 5 MB), else paths.
- NOTE: the cross-check table (section 6); a summary table vs time to collapse at
  t_cc = (match), 100, 30, 10, 3, 1, 0.3, 0.1 yr, 1 month, 1 week, 1 day, Si ignition, collapse (rows that exist):
  L(r=3), L(r=8), L(m_in3) with r(m_in3), L(m_in8) with r(m_in8), L_He shell (max eps_3a and He base) with its radius,
  L_surf, R_surf, Teff, and L(r=8)/L(r=8 at match). Say when the expansion happens relative to Si ignition, how far
  r(m_in8) moves from 8 Rsun (i.e. whether a fixed-r boundary is still in the envelope), and whether the envelope above
  8 Rsun stays in thermal equilibrium (L(8 Rsun) vs L_surf).
- Everything measured from your run, with file paths; mark guesses as guesses.

Prepared on viper 11-08 (worker hgmesa_1108); viper copies of the inputs: `/viper/ptmp2/jinma/wu2022_zenodo/`,
`/viper/ptmp2/jinma/hegiant_large_1014/ic/`.
