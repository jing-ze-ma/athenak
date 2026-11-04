# Reproduction run (Ma+2026 BSG): initial condition, grid and column check — 10-06

Status: the IC file is written and the CPU column check is done. **HOLD** until the RT setting is decided:
no 3-D inputs, smokes or chain were made. Work dir: `/viper/ptmp2/jinma/bsg_1001/repro/` (`ic/`, `ic3d/`,
`grid/`, `col/`). The SK -67 133 work stays parked as it is on disk (see the note in `SK67_PREP.md`).
Source of every number: the scripts and logs named in each section.

## 1. The IC: `repro/ic3d/ic_ma2026.txt` (he_ic_cols = 5, 7 columns, fmlt = 0)

Builder: `repro/ic/build_ic_ma2026.py` (a copy is in `ic3d/`). Run: `python3 build_ic_ma2026.py ic_ma2026_full.txt`
(16 s, login node). Format: r rho eint F_r E=aT^4 | T fmlt.

- **Interior** is L. Ma's `ma2026_files/BSG/MESA_profile_combined_cgs.txt`, used **as is**: no re-balance, and T and
  rho exactly as in the file. The one duplicate row at 19.994 Rsun was removed. Rows from 18.0 Rsun (this covers the
  3 ghosts below 20 Rsun) up to 46.767 Rsun.
  - eint = 1.5 rho k T/(0.62 m_u); E = a T^4; L = 4 pi r^2 F_r = 5.581461e38 (the median over the file).
- **Photosphere.** The file has none: its last step jumps from 48.32 to 52.49 Rsun, where T is still 63 kK. The paper
  integrates the radiative-only profile outward "to the photosphere" and, in footnote 5, adjusts "the inner boundary by
  a small amount" for convergence. I did the same:
  - Outward radiative-only integration: dT/dr = -3 kappa_R rho L/(16 pi a c r^2 T^3) and dP_tot/dr = -rho GM/r^2.
    Point mass of 20 Msun; ideal gas mu 0.62 plus radiation; kappa_R from the padded Ma table with the code's lookup
    (bilinear in logs). The photosphere is where L = 4 pi r^2 sigma T^4.
  - Started from the file's own (T, rho) at 48.318 Rsun, no photosphere is ever reached. T stalls at 62.7 kK, and the
    gas pressure runs out at 60.15 Rsun. The file's 52.49 Rsun row (63.4 kK, 2.6e-11) lies on this same failed
    branch, so Linhao's integration most likely stopped the same way.
  - A photosphere exists only above a critical density adjustment (a separatrix). The adjustment needed at the start
    row is +33 % at 48.32 Rsun, -1.9 % at 46.02 and -7 % at 44.0 Rsun. It changes sign between 46.55 and 46.77 Rsun
    (`ic/probe2-4.py`).
  - **Used:** start at the file row r = 46.76706 Rsun with rho x (1 + 3.4637e-3); the critical value is 3.46359e-3.
    The adjustment is shot on the supercritical side so that the attached atmosphere gives tau_R(r_photo) = 2/3
    (tau 0.675 in the written file, whose floor rows are coarse). The file rows above 46.767 Rsun are replaced by the
    integration (1e8 cm steps). There is no density inversion between 46 and 53 Rsun.
  - **r_photo = 53.030 Rsun, T_photo = 15488.5 K, rho_photo = 3.38e-11, P_tot = 215.4, P_gas = 70.2**, and
    kappa L/(4 pi c GM) = 0.43 at r_photo. Every near-critical start between 45.6 and 47.7 Rsun puts the photosphere
    at 52.4-53.7 Rsun. The paper quotes a 1-D photosphere of ~50-51 Rsun: that is MESA's, while ours is the
    radiative-only profile. This is question 1 below.
- **Atmosphere** (the paper, literally): T = T_photo and rho = rho_photo exp(-(r - r_photo) g(r) rho_photo/P_photo).
  - g = GM/r^2 is the local value. P_photo is the **total** pressure. The paper writes g/P_photo; I restored the factor
    rho_photo so the exponent is dimensionless.
  - Scale height at r_photo: P/(rho g) = 0.470 Rsun. This form is not hydrostatic.
  - Floor: 1e-20, our dfloor; the paper does not state its floor. It is reached at 72.05 Rsun and held constant to
    83 Rsun.
- **F_r** = L/(4 pi r^2) everywhere, including the atmosphere. The paper sets F_r = 0 above the photosphere for its
  intensity IC (see section 3).
- **How he_star_m1 loads it.** he_ic_balance defaults to false and is set false here, so nothing re-integrates or
  re-balances the IC.
  - The pgen resamples the column in log onto he_nfine = 262144 fine points and interpolates that to the cells. E is
    taken from the file (cols 5).
  - F_r becomes both the initial M1 flux F1 and the wb_arad reference force kappa_R F_r/c. That is a scheme choice of
    ours (the x1 well-balanced pair), not a re-balance.
  - The column must lie inside the tables (log rho -20..-3.98, log T 4.19..6.11: it does) and cover the ghosts.

## 2. Grid: `repro/grid/design_N267.txt` (`design_plateau_ma.py 267 34.0 0.45 80.0 20.0 0.153`)

- **r: 20-80 Rsun, nx1 = 267.**
  - dr is 0.153 Rsun on the plateau 34.0-55.86 Rsun and at most 0.1531 over the whole 35-55 Rsun. Elsewhere it is
    0.3056-0.3060.
  - The tanh transitions are 0.45 Rsun wide; the largest ratio between neighbouring cells is 1.236. 52 cells lie below
    35 Rsun, 131 in 35-55 and 84 above 55 Rsun.
  - Keys: x1min 1.391400e12, x1max 5.565600e12, f_stretch_r_p_amp -6.8085e-01, xa 1.723088670756e-01,
    xb 7.035565721502e-01, w 5.507821105970e-03.
- **Angular grid: 256 x 256 everywhere**, the same pi/3 x pi/3 wedge (theta pi/3-2pi/3, phi 0-pi/3).
  - That equals their refined level. Their 128 x 128 outside 35-55 Rsun cannot be had without SMR.
- **Blocks:** 16 blocks of 267 x 64 x 64. Total 17.50 M cells, against 10.49 M for Ma (1.67x) and 13.50 M for arm 2
  (206 x 256^2, 1.30x).
- **Cells per H_p,tot = P_tot/(rho g) on this IC:**
  - 22.1 at the Fe kappa peak (41.3 Rsun, kappa 1.59, H_p 3.38 Rsun);
  - at least 6.6 in 20-35 Rsun;
  - 2.9 at r_photo, where H_p is 0.44 Rsun (arm 2: 2.35);
  - at least 2.7 in 35-55 Rsun, near 52.8 Rsun.

  The paper's ">10 cells per scale height in 35-55 Rsun" therefore cannot hold for H_p,tot near the photosphere, even
  on their own dr = 0.153. This is question 9 below.
- **Cost.** The column dt is 123.1 s, against 124.3 s for the arm-2 column. So the cost per physical time is about
  1.30 x 124.3/123.1 = **1.31x arm 2**. That is a CPU-column estimate, not a GPU measurement. Memory per GPU rises by
  1.30x per block, which needs a GPU smoke (on hold).

## 3. The paper's setup (sect. 2.2) item by item, mapped to our keys

| item | Ma+2026 | ours | status |
|---|---|---|---|
| domain, wedge | r 20-80 Rsun, theta pi/3-2pi/3, phi 0-pi/3 | same | same |
| grid | uniform 196 x 128^2 base, one SMR level in 35-55 Rsun | stretched r (above), 256^2 everywhere | closest available |
| EOS, gravity | ideal mu 0.62, gamma 5/3; point mass 20 Msun (GM) | eos ideal, units/mu 0.62, he_gm 2.6542488e27 | same |
| opacity | TOPS kR, kP; kappa_a = kR - kappa_e; kP exchange; clamp at the edges | the padded Ma tables (he_opac_table / he_planck_table), no he_opac_logd_min; kappa_F = kR, kP exchange | same |
| IC interior | radiative-only MESA re-integration, as is | the file as is, he_ic_balance = false | same |
| IC atmosphere | isothermal T_photo, rho_photo exp(-(r - r_photo) g/P_photo) to the floor | literal form, P_tot, local g, floor 1e-20 (sect. 1) | same, with the assumptions listed as questions |
| IC radiation | two-zone I, E = aT^4, F_r = L/4 pi r^2 below r_photo, 0 above | E = aT^4; F1 = F_r = L/4 pi r^2 everywhere | **not possible without code**: F_r = 0 above r_photo needs a separate initial-M1-flux column or key, because the file's F_r is log-interpolated and is also the wb_arad reference force |
| velocity, seed | zero, no seed | zero, he_seed = 0 | same |
| inner BC hydro | ghosts = IC rho, P (unscaled); first ghost v_r solved so the HLLC contact speed is 0 (zero mass flux); v = 0 elsewhere | he_bc_inner = wall (default) + he_wall_closed = true: ghosts = the IC scaled to the edge cell, mirrored v_r, exactly zero mass flux through the wall face (wall_closed_ix1) | closest available. The exact form needs code: unscaled IC ghosts (as in he_bc_inner = inflow) combined with the wall_closed_ix1 zero-flux face |
| inner BC radiation | two-zone I in the ghosts with F = L/4 pi r^2 | implicit_bc_x1min = flux, implicit_flux_x1min = L/(4 pi r_in^2) = 2.2942e13 | closest (the same flux) |
| outer BC hydro | copy rho, P, v_t; v_r copied if > 0, else 0 | he_bc_outer = noinflow (default): ghosts are the IC scaled to the edge cell (= a copy on the constant floor at 80 Rsun), v_r copied out, 0 in | same in effect |
| outer BC radiation | outgoing I copied, incoming I = 0 | implicit_bc_x1max = marshak, flux 0 (vacuum) | closest (the moment analogue) |
| sponge | none stated | none in the check (he_sponge_* removed) | same; see the drift below for whether the top needs one |
| floors | not stated | dfloor 1e-20, tfloor 5e3 K, pfloor 1e-12 | question |
| CFL | 0.3 hydro, 0.025 radiation (explicit transport sub-steps) | cfl_number 0.3; implicit M1 (no radiation CFL) | same / not applicable |
| Riemann solver | HLLC | rsolver hllc | same |
| reconstruction, integrator | not stated | plm, rk2 | question |
| RT | implicit discrete ordinates, 120 angles | M1, closure vet_col in the check (the RT setting is undecided) | not possible (a different method) |

Our scheme keys with no counterpart in the paper: wb_option polytropic + wellbalance_dynamic, force_reference wb_arad,
fofc, etotgrav, rad_signal_speed and the R2+B solver keys of arm2r2b_4n.

## 4. CPU column check: `repro/col/chk` (input `col/chk.athinput` from `col/mk_chk.py`)

- **Setup.** The 4 x 4 thin column with the production solver keys, Ma IC as is, he_ic_balance false,
  mlt_flux_frozen false, no sponge, he_bc_outer noinflow, nlim 100.
  - Binary: `sk67/bin/athena_he_cpu_11243930` (the existing one; it reads the padded tables).
  - Analysis: `col/ic_check.py`, `col/ana.py`, `col/logsum.py`.
- **Result.** rc 0; 0 FATAL, 0 nan; 100 cycles to t = 1.23e4 s; **dt 122.9-123.1 s**. There were 440
  NEWTON-FALLBACK lines, all at r > 53 Rsun, most of them in 54-61 Rsun. FOFC acted at 72-80 Rsun (the floor).
  Mass drift was 2.7e-15 and L_top/L_in was 1.014.
- **Cycle 0 against the IC file**, max |dump/IC - 1|:
  - rho: 2.1-2.9e-7 in every region (20-25, 25-35, 35-53, 53-80 Rsun);
  - eint: 2.1-2.7e-7;
  - E: up to 2.7e-7;
  - v = 0.

  This is the same 2-3e-7 level as the resampling seen in the SK67 check.
- **The transient** (max over the region, between dumps about 2000 s apart; deep < 35 Rsun, env 35-53, atm above
  r_photo with rho > 1e-17, floor):

| t [s] | deep \|v\| | env \|v\| | atm \|v\| | floor \|v\| | atm dlnrho, dlnT per ~2000 s |
|---|---|---|---|---|---|
| 2.1e3 | 0.35 km/s | 0.06 | 2.0 | 1.8 | 1.3e-2, 0.37 (the atmosphere's T adjusts at once) |
| 6.0e3 | 0.83 | 0.43 | 6.0 | 5.2 | 2.8e-2, 3.0e-3 |
| 1.2e4 | 1.39 | 1.54 | 12.2 | 10.5 | 5.6e-2, 4.5e-3 |

- **Reading the transient.** Everything falls inward and the speeds grow linearly. The deep maximum is at 22.9 Rsun,
  with an acceleration of about 11 cm/s^2, i.e. about 1 % of g. That is the ~1 % HSE mismatch of the file with our
  EOS (MA2026_FILES.md section 3). The env maximum is at 53.0 Rsun.
- The literal atmosphere falls at about 100 cm/s^2 (about 0.6 g at 56 Rsun). It is neither hydrostatic nor in
  radiative equilibrium. dt did not change in 100 cycles.
- This is the start transient of the paper's setup. It is documented here, not removed.

## 5. Questions for Linhao (short list for the user to send)

1. The 1-D photosphere. Your profile ends with one step from 48.32 to 52.49 Rsun (T 63 kK, still inside the star).
   - Which r_photo, T_photo, rho_photo and P_photo did the 3-D IC use?
   - Did the radiative-only integration reach tau = 2/3, and which "inner boundary" adjustment (footnote 5) did you
     make?
   - Ours, from a +0.35 % density adjustment at 46.77 Rsun, gives 53.03 Rsun, 15.5 kK and 3.4e-11 g/cm^3. The paper
     says about 50-51 Rsun.
2. The atmosphere formula: is it rho_photo exp(-(r - r_photo) g rho_photo/P_photo)? Is P_photo the total or the gas
   pressure? Is g the local GM/r^2 or the value at r_photo?
3. The density floor value, and the pressure and temperature floors if any.
4. Reconstruction (PLM/PPM?) and the time integrator (vl2/rk2?); the number of ghost layers.
5. Inner BC:
   - are rho and P in all the ghost layers the unscaled IC values?
   - is the zero-contact-speed v_r applied only in the first ghost layer?
   - what ghost T is used for the intensities?
6. Intensity IC: confirm F_r = 0 above r_photo (with E = a T_photo^4 there).
7. Any sponge, damping or velocity ceiling (also during the first days)? Any reduced speed of light?
8. Was the IC adjusted for HSE with mu = 0.62 and the point mass? With our EOS the file is about 1 % out of balance.
9. The SMR region: exactly 35-55 Rsun in r for both refined directions? Which scale height gives ">10 cells per scale
   height"? With P_tot/(rho g) we find 2.7 cells near 53 Rsun at dr 0.153.

## 6. Next (on hold: the RT setting)

- 3-D input from `col/chk.athinput`: the 256 x 256 wedge, 16 blocks of 267 x 64 x 64, no seed, no sponge, and the
  closure per the RT decision. Then a Raven gpudev or apudev smoke (memory: 1.30x arm 2 per block) and the chain
  (job prefix bsg_repro, run dir /raven/ptmp/jinma/bsg_raven_1002/repro_*).
- Optional, if wanted (code needed): F_r = 0 above r_photo for the M1 IC, and unscaled inner ghosts with the
  wall_closed face.
- Linhao's answers may move r_photo, the atmosphere and the floor. The IC rebuild then takes about 20 s.
