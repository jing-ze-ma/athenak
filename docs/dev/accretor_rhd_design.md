# Accretor RHD scoping: general EOS + table opacity + implicit M1/VET for ry_per_accretor

Status: design note only (2026-10-08, Caltech). No code, builds or jobs. Base: branch
accretor-adiab-1008 (= origin/accretor-1007 + 65ec7626 `env_t_relax_stream`).
Line numbers refer to this tree: `src/pgen/ry_per_accretor.cpp` (= RP), `src/pgen/he_star_m1.cpp` (= HS),
`src/rad_m1/*` (= M1).

## 0. Summary and recommendation

- **Build it inside `ry_per_accretor`**, as a third thermodynamic mode of the envelope stage
  (`problem/inner = envelope` + a new `problem/thermo = rhd`). Copy the needed pieces of `he_star_m1`
  into the accretor; do not turn `he_star_m1` into a binary-star pgen (section 5).
- **Use `closure = vet_col` only. `vet_gd` cannot run on a 4-cell theta slab** (section 4): its lateral
  halo is capped at min(nx2, nx3) = 4 cells, but it needs 10-50. It also treats theta walls as periodic.
- **Physical units**: add a `<units>` block (Rsun, km/s, and a fixed physical density unit). Mdot then
  becomes a per-run input (`rho_stream`); the star's IC no longer depends on it (section 2).
- **The hot 300 km/s ambient cannot be left to the radiation.** It is optically thin (tau_es ~ 5e-4)
  and has ~350x less thermal energy than the local radiation field. Its temperature is set entirely by
  whatever Planck opacity the floor cells get, and an edge-held table value cools it in ~1e-6 s.
  Floor/ambient cells need kappa_abs = 0 (scattering only), plus the existing T relaxation. This is new
  code in the module (section 2.3).
- **Cost**: ~0.5-1.0 s/cycle on 1 H200 (central 0.7), against 0.03 s/cycle hydro-only. That is
  ~1 H200-day per half orbit and ~3 H200-weeks for 10 orbits (2.5M cycles). The Picard count at
  dt ~ 1.25 s is the main unknown. One calibration smoke measures it (section 6).
- **Order**: EOS interface (bitwise gate) -> general EOS without radiation -> M1 on the static envelope
  (flux = L gate) -> stream on -> cost decision (section 7).

## 1. Ideal-gas / isothermal assumptions in ry_per_accretor and their EOS replacement

Stage 1 (`inner = surface`, absorbing wall, `<hydro>/eos = isothermal`) stays as it is. RHD is
refused there (a new fatal check). Everything below is the envelope stage (`inner = envelope`).

| where (RP line) | assumption | through the EOS interface |
|---|---|---|
| 290-366 `RyPerBC` (stage 1) | `iso_cs`, exp(-dPhi/cs^2) ghosts | unchanged; fatal if `<rad_m1>` exists and inner = surface |
| 505, 521-523 inner ghost | T = (gamma-1) e/rho, floored at 1e-3 c_ph^2 | `eos.Temperature(da, ea)`; floor = `eos.tfloor` |
| 524-537 inner ghost fill | isothermal HSE walk `exp(-(Phi_g - Phi_a)/T)`; `IEN = rho T/(gamma-1)` | `WBAdvance(eos, mode, ...)` along Phi_wb (pattern: HS 1640-1677, the hse top). Or: the IC column scaled to the edge cell (HS 1682-1704). `IEN = e_g + KE + rho Phi`. Rad ghosts via `radm1::M1FillGhost` (HS 1676/1703) |
| 547-568 outer ghost (no window) | same isothermal walk in the TRUE Roche Phi | same `WBAdvance`, or keep the isothermal walk with T = eos T of the edge cell (the ambient is isothermal by design) |
| 553-557 stream window | `tg = env_cs_stream^2`, `e = rho tg/(gamma-1)` | `e = eos.EnergyFromTemperature(dg, T_don)` (new key in kelvin). Rad ghost E = a T_don^4, F_r = (4/3) E v_r (thick advective limit) |
| 597, 616 WB gravity | `p = (gamma-1) e` | `eos.Pressure(d, e)`, exactly as HS 1460-1462 |
| 604-639 thermal relaxation | target `P/rho = cref2` (EnvC2 polytrope / c_ph^2 / c_amb^2), `e_target = rho ct/(gamma-1)` | **off** in the envelope and the stream under RHD (`env_t_relax_env = env_t_relax_stream = 0` already means off, line 637 + 65ec7626). **Kept** for the ambient class only, target T_amb in kelvin -> `eos.EnergyFromTemperature(d, T_amb)` |
| 640-653 floor sponge | none (kinetic energy only) | unchanged |
| (missing) | no radiation-force work | add the `force_reference_work = split` term of HS 1464-1487 (`rho a_ref v_r` into IEN) |
| 779-806 EnvSetup keys | c_ph, c_amb, c_stream in km/s (mu 0.62 implied) | kelvin keys (`env_t_ph`, `env_t_amb`, `env_t_stream`), converted through the EOS. Keep c_ph^2 = P/rho(photosphere) as the *geometric* scale height for the r_top / equipotential-top search (820-866, 857); that is geometry, not thermodynamics |
| 926-937 ambient rho_amb(r, phi) | isothermal HSE at c_amb in the Roche Phi | unchanged (fully ionized at 6.8 MK, ideal to <1e-3); P/rho from the EOS at T_amb |
| 940-951 `cref2_` | n = 3 polytrope `EnvC2(psi)` | replaced by the IC column T(psi) (section 1.1); kept only to classify "inside" (630) |
| 1021-1101 IC (`ryper_ice`) | analytic n = 3 polytrope (`EnvRho`/`EnvC2`, 183-191); per-face isothermal walk (1067-1074) = the rest state of `wb_option = isothermal`; `IEN = d T/(gamma-1)` (1093) | the column from tables (1.1); discrete balance with the code's own walk under `wb_option = polytropic` (the `he_ic_balance` solve, HS 717-800, per phi column); `IEN` from `EnergyFromTemperature` |
| 1126-1136 checks | `eos.is_ideal` required | allow `eos = general`; require `eos_radiation = false` (already enforced by M1, `rad_m1.cpp` 506-545) |
| `<hydro>/wb_option = isothermal` (input) | isothermal background | `polytropic` (the radiative envelope is not isothermal; this is the He giant / AG Car choice, and `he_ic_balance` requires it, HS 727-730). `env_wb_depth` (980-993) is unchanged |

Coriolis, the Roche gravity arrays (`gr_`, `gp_`, `gpl_`), the stream ballistics, `WrapPhi`, the flux
integrals (657-713) and the sponge do not touch thermodynamics and stay as they are.

### 1.1 The envelope IC under RHD

`he_star_m1` reads a 1-D column file (r, rho, eint, F_r [, E, ...]; HS 337-500, built by
`make_ic_*.py --hse own` with the code's EOS) on a fine radial grid. It then sets
Phi_eff = Phi - int a_ref dr (HS 631-664) and balances the cells discretely (HS 717-800). The accretor
cannot use r as the column variable, because the Roche equipotentials bulge (photosphere 9.00 Rsun at
phi 90, 9.485 at phi 0). So:

1. **Column in psi = Phi_s - Phi_wb** (the variable the accretor already uses, RP 949).
   Build it offline (python, `scripts/make_ic_accretor_rhd.py`, new), with the run's EOS and opacity
   tables. Use an Eddington-grey T(tau) atmosphere on the photospheric equipotential, then integrate
   inward: radiative diffusion at constant L, HSE with P_gas + P_rad, and r_eq(psi) from the phi = 90
   cut. Output columns: psi, rho, eint (or T, see `he_ic_eint_from_t` on branch he-ic-eint-from-t),
   F_eq, E.
   Alternative: the Wade+2026 MESA gainer profile at the onset of Case A, remapped to psi (the user
   would have to supply it).
2. **Pgen**: T(psi), rho(psi) -> cells through `pwc`. Then add the reference acceleration
   a_ref = kappa_t F/c (face form, HS 701-715) to `phicc_wb`/`phi_wb_x1f` as -int a_ref dr along each
   phi column. This needs a 2-D (r, phi) Phi_eff, which `EnableWBEffectivePotential` already allows
   (RP 911-922 fills it per (m,k,j,i)).
3. **Discrete balance**: the `he_ic_balance` march with T fixed, done per phi column (the accretor's
   `ryper_ice` already loops over (m, k) columns, RP 1054). The ambient takeover rule (RP 1083) is
   unchanged.

**Physical numbers (n = 3 analytic model, for orientation).** At r_in = 6.3 Rsun:
psi ~ 1.45e5 (km/s)^2, T ~ 2.8 MK, rho ~ 2.6e-3 g/cc, P_rad/P_gas ~ 0.15 (beta ~ 0.87).

**Physics caveats.**
- A radiative envelope in Roche geometry cannot be in exact radiative equilibrium (von Zeipel), so
  expect a weak transient.
- A 16 Msun MS star at log L 4.73 may carry weak HeII (~50 kK) and Fe-bump (~200 kK) convective zones.
  Check the column for Schwarzschild instability, and decide whether 2-D r-phi convection there is
  acceptable (no MLT scaffold is planned).

## 2. Physical units

### 2.1 Fixing the density unit

Today the run is scale-free in rho (RP 18-22). The only physical input is the ratio
env_rho_ph = rho_ph/rho_peak(r_out) = 0.055 (Mdot 1e-4 Msun/yr). Opacity (kappa rho) and
aT^4 vs rho kT break the scaling:
- Every Mdot becomes its own run.
- The stellar IC depends on rho_ph in g/cc.

Recommended `<units>`:

| key | value | note |
|---|---|---|
| length_cgs | 6.957e10 | Rsun (unchanged) |
| time_cgs | 6.957e5 | Rsun/(km/s) (unchanged) |
| mass_cgs | rho_u x Rsun^3 = rho_u x 3.367e32 | rho_u a FIXED physical unit |
| mu | 1 | general EOS keeps mu_ref = 1 (eos.cpp 220-225) |

Two choices:
- **(a) rho_u = 5.11e-8 g/cc** (the env13 stream peak). Every code number of env13 stays the same
  (rho_stream 1, env_rho_ph 0.055, dfloor 1e-8). A hydro-only env13 restart is unit-compatible for
  rho and v (IEN must be recomputed).
- **(b) rho_u = rho_ph = 2.82e-9 g/cc.** The star is fixed, and `rho_stream` = Mdot/(1e-4) x 18.1
  carries Mdot.

**Recommendation: (b) for an Mdot scan, (a) for the first runs** (re-use of env13 restarts and
numbers). Derived code constants (case a):
- pressure unit 511 dyn/cm^2
- T_unit (mu_ref 1) = m_u (1 km/s)^2/k = 120.3 K
- arad = a T_unit^4/P_unit = 3.1e-9
- c_light = 2.998e5

The pgen already checks units against the tables (HS 528-541); copy that check.

### 2.2 Floors and ambient in physical units (case a)

| quantity | code | cgs | comment |
|---|---|---|---|
| dfloor = rho_amb | 1e-8 | 5.1e-16 g/cc | below every opacity table (log R ~ -17 at 6.8 MK; tables stop at log R -8) |
| ambient at photospheric equipotential | 1e-7 | 5.1e-15 g/cc | T_amb = 300^2 km^2/s^2 x mu m_H/k = 6.8e6 K |
| photosphere | 0.055 | 2.82e-9 g/cc | 28.3 kK; log R ~ -3.9; kappa ~ 0.5-1; dtau per plateau cell (dr 2.5e-3 Rsun) ~ 0.3 |
| stream peak at r_out | 1 | 5.11e-8 g/cc | 26.2 kK; tau across the 0.735 Rsun width ~ 2-3e3 |
| r_in | ~5e4 | ~2.6e-3 g/cc | ~2.8 MK |

Tables needed:
- **EOS** (X 0.70, Z ~ 0.014-0.02): log rho -17 .. 0, log T 3.6 .. 7.3. The ambient (6.8 MK) is
  above the He giant's eos_logt_max 6.8.
- **Opacity**: Rosseland + Planck on the same grid, X 0.7. `data/stellar_opac` has
  rosseland_gs98_x0.7_z0.014; the BSG used TOPS x0.7 z0.008; a Planck table must be built
  (`data/stellar_opac/planck_tools`).
- Set `he_opac_logd_min` (~ -17) with `he_opac_extend_hold` (HS 222-246, 512-519).
- `tfloor_kelvin` (eos.cpp 177-231) ~ 5e3 K.

### 2.3 The hot ambient under radiation (main physics risk)

Optical depths and energies of the ambient (case a):
- tau_es over 4.5 Rsun ~ 5e-4: **optically thin**. True free-free kappa_P ~ 1e-16 cm^2/g, so it is
  **thermally decoupled**.
- Gas thermal energy ~ 7 erg/cc. The radiation energy density near the star is ~ 0.5 a T_eff^4
  ~ 2.4e3 erg/cc, i.e. ~350x larger.

Whatever kappa_P the floor cells get therefore sets T_amb. With an edge-held table value
(kappa_P ~ 1e-2 assumed), the emission time at 6.8 MK is ~ 3e-7 s. The ambient would collapse to
T_rad ~ 2-3e4 K in one step and fall onto the star (it is supported only by its 300 km/s pressure).

Required:
- A **new `<rad_m1>` key** (e.g. `opac_abs_rho_min`, or a pgen-filled per-cell multiplier). It zeroes
  kappa_P/kappa_E (absorption) below a density, or for cells the pgen flags as ambient, and keeps
  kappa_es = 0.34. Default off, bitwise. Today there is no such hook: only `dbg_opac_patch`, a debug
  box, and `opac_freeze`.
- The existing T relaxation to T_amb stays for the ambient class (RP 630-639). It is an energy
  source/sink that the radiation does not see, so it must be booked in the energy budget.

Radiation force on the ambient: Gamma_e = kappa_es L/(4 pi G M c) ~ 0.09, harmless.

The ambient is a stand-in for the real B-star wind (Mdot ~1e-7 Msun/yr, line-driven). Grey M1
cannot drive it. Open decision: keep the artificial hot ambient, or test a cooler, denser one.

## 3. Radiation: how he_star_m1 hooks rad_m1, and what the accretor needs

The module is built by the mesh from `<rad_m1>` (`pmbp->pradm1`); its tasks are in
`rad_m1_tasks.cpp`; restarts write/read its four moments (`outputs/restart.cpp` 82-207,
`pgen/pgen.cpp` 237-245). The pgen only configures it:

| duty | he_star_m1 | accretor (new) |
|---|---|---|
| sanity | HS 296-331: `<hydro>`, sp/cs, etotgrav + WB x1 dynamic, `opacity = table`, `force_reference = wb_arad`; one block along x1 (365-367) | same; one x1 block is already required (RP 1043-1052) |
| tables | `HsReadOpacityTable` (192-271), `SetOpacityTables` + unit check (503-541) | copy verbatim |
| range check | column inside the EOS and opacity tables (580-618) | over the IC column + the stream state + T_amb (the ambient is OUTSIDE the opacity grid by design: check only the cells with kappa_abs on) |
| bottom flux | `implicit_bc_x1min = flux`, `implicit_flux_x1min` = column F(r_in) checked (620-630; M1 `rad_m1_implicit.cpp` 1450-1455) | same; L = 5.37e4 Lsun (log 4.73) over 4 pi r_in^2. One scalar for the whole face: OK, the r_in = 6.3 equipotential is nearly spherical (tide ~ r^3) |
| top | `implicit_bc_x1max = marshak` (+ `vet_col_surface_q`) | marshak (vacuum) on the whole face, including the stream window (below) |
| force reference | Phi_eff + `SetForceReference(aref)`, `fref_wsplit_ok` (683-715) | the 2-D Phi_eff of 1.1 |
| reference work | `HeStarGravity` split term (1464-1487) | add to `RyPerSrcEnv` |
| x1 ghosts of E, F | `M1FillGhost` (1676, 1703); the implicit solve imposes its own face BCs | the same in `RyPerBCEnv`; stream window: E = a T_s^4, F = (4/3) E v |
| wall flux fix | `he_wall_zero_flux` (1489-1513) | the accretor uses `wall_closed_ix1` (RP 812), which has the same effect for free slip; keep it |
| restart | everything rebuilt from the input (header 44-47) | already so (EnvSetup rebuilds everything, then returns on restart at RP 1019) |
| history | L_bot..L_top, E_rad, Picard (1718-1812) | **the accretor hst is already at NHISTORY_VARIABLES = 22** (RP 748-753; outputs.hpp 19): drop the 9 rate columns or raise the global limit (touches every build; check bitwise) |

**Stream inflow radiation.** The M1 face BC is one type per x1 face; there is no per-cell mask. Leaving
the window as Marshak (vacuum) is acceptable:
- The stream's radiative diffusion time across its half width (~2e3 s = 3e-3 code) is ~6x shorter
  than its infall time r_out -> surface (~1.3e4 s = 0.019 code).
- So the stream's thermal state is set radiatively on the way in, not by the BC.

Donor irradiation is the real omission at the top. The donor (log L 4.84) gives the gainer's facing
side ~20 % of its own flux. Adding it needs a per-column incoming flux at r_out (new module code). It
is optional, a later stage.

**Rotating frame.** No radiation Coriolis is needed:
- The M1 source terms use the rotating-frame gas velocity.
- The dropped frame terms are O(Omega r/c) <= 6e-4 (Omega r <= 185 km/s at r_out). They only act on
  free-streaming radiation, whose momentum the gas does not feel.
- Coriolis stays explicit on the gas (RP 619-620) and does no work.

**etotgrav.** The gas-radiation exchange changes the gas internal energy; with etotgrav the conserved
energy also carries rho Phi (RP 893-908). This is the same combination he_star_m1 runs (HS 317-322).
The 09-28 double-count guard (no explicit centrifugal work) is unchanged.

**Hydro signal speed.** `rad_signal_speed` defaults on with implicit transport (`hydro.cpp` 243).
- At the photosphere and deep inside it changes dt little (beta ~ 0.87).
- In the impact shock the post-shock gas is radiation-pressure dominated: ram pressure ~6e7 dyn/cm^2
  gives T ~ 4e5 K and P_rad/P_gas ~ 5. Expect dt -20..-40 % there (stream v + ~200 km/s on
  dr 2.5e-3 Rsun).

## 4. VET on this geometry

| check | vet_col | vet_col_lat | vet_gd |
|---|---|---|---|
| one MeshBlock along x1 | needed (vetcol 418-419); accretor has it (640 = 640) | needed | needed (vetgd 650) |
| 3-D | multi-D, implicit transport (416) | sp only (vetlat 161-167) | `three_d` (vetgd 269); nx2 = 4 is 3-D: OK |
| outer x1 | Marshak or reflect | - | **Marshak only** (vetgd 270-275) |
| theta extent / walls | none: column-local, no communication (vetcol 105-112) | physical x2 faces handled (vetlat 960-966) | theta AND phi treated as periodic wraps (vetgd 33-37, 634-660); a reflecting theta wall is not represented |
| lateral reach | none | the hydro ghost band; oblique reads beyond it clamped (counted) | halo capped at min(nx2, nx3) (vetgd 316-334) |

**vet_gd is out.** On the accretor grid r dtheta = r dphi = 0.003 r, and the near-tangent chord
sqrt(r_{l+1}^2 - r_l^2) gives a needed halo of:
- ~10 cells in the photospheric plateau (dr 2.5e-3 at r 9);
- ~50 cells at r_in (dr 0.077 at r 6.3).

With nx2 = 4 the cap is 4, so nearly every oblique read would be clamped. For comparison, the He giant
input records that 32x32 lateral blocks already fail vet_gd (hegiant_scout128_N445_fresh.athinput,
meshblock comment). vet_gd would need nx2 >~ 64, i.e. a real 3-D wedge (16x the cells).
vet_col_lat with ng = 2-3 ghosts is in the same position: mostly clamped reads.

**Use closure = vet_col.** It has no theta assumption, and the memory rule says vet_col is enough for
the He star and dhj. With nx2 = 4 and reflecting theta, the slab is a theta-invariant band. M1's theta
flux is zero at the walls, which is consistent with that.

Where vet_col's radial-column assumption breaks:
1. **The stream above the photosphere.** It is an optically thick ribbon 24 deg off radial, flanked by
   thin ambient. Its phi-flank losses are carried by the M1 F_phi with the tangential factor
   (1 - f_rr)/2 from a radial column that sees the ribbon as a radial slab. This is wrong at the flanks
   (thin limit); the thick interior is diffusion and fine.
2. **The impact shock / hot spot.** Strongly non-radial flux, at the column scale of the stream width.
3. **Columns shadowed by the stream.** vet_col sees the occulter only along its own radial line.

Energy is conserved in all three; f errors only redirect flux.

**A 2-D modeling limitation, independent of closure.** The theta-invariant band gives the stream no
vertical (z) radiative losses: a real stream tube loses through ~2x more surface. Stream cooling is
underestimated by roughly that factor.

## 5. Reuse vs new code: recommendation

**Option A: extend ry_per_accretor (recommended).**
- The accretor-specific machinery is ~900 lines and all of it is needed: the Roche potential with
  phi dependence (159-180), PhiWB with spin and the equipotential cap, the ballistic stream and window
  BC, Coriolis, the RK-weighted flux integrals and history (657-765), the ambient classification,
  `env_wb_depth`, and r_meas.
- he_star_m1's machinery is r-only: the fine radial column `HsLinInterp(cphi, ...)`, BCs from the
  scaled 1-D column, face-luminosity history at fixed r. ~1200 of its ~2000 lines (MLT scaffold,
  adaptive closure, seeds, wind IC, cs support) are not needed.
- What the accretor needs from he_star_m1 is ~300 lines, mostly self-contained helpers:
  - the table reader;
  - the range and unit checks;
  - the column -> E = aT^4, kappa_t, a_ref code;
  - the reference-work term;
  - `M1FillGhost` calls;
  - the L/E_rad/Picard history pieces;
  - the `he_ic_balance` march.
- Copy them first. Factor into a shared header (e.g. `src/pgen/m1_star_util.hpp`) only after both pgens
  pass their bitwise gates, so he_star_m1 production builds are not touched.

**Option B: he_star_m1 variant + Roche + stream (rejected).** It would have to import exactly the
accretor list above, and replace its r-indexed column, BCs and history with (r, phi) versions. That is
more code, and it puts the He giant / AG Car production pgen at bitwise risk.

**New module code (both options)**, default off and bitwise:
- the absorption-opacity mask for floor/ambient cells (2.3);
- optional: a per-column incoming flux at the outer x1 face (donor irradiation);
- possibly: NHISTORY_VARIABLES (3).

## 6. Cost estimate (1 H200)

| run | cells | s/cycle | s per cell-cycle |
|---|---|---|---|
| accretor env13, hydro only | 5.24e6 | 0.030 (33 cyc/s) | 5.8e-9 |
| He giant N445, M1 + vet_gd + general EOS (Picard ~18 at N897 smoke) | 7.29e6 | 0.97 | 1.33e-7 |
| AG Car B, M1 + vet_gd + general EOS (Picard 7-24) | 7.86e6 | ~1.6 | ~2.0e-7 |

Radiation is ~95 % of the M1 runs' cost. Scaling per cell gives **0.7 s/cycle** for the RHD accretor
(range 0.35-1.05), if dt stays ~1.8e-6 code (1.25 s):

| span | cycles | wall time (0.7 s/cycle) | range | H200-h | at $1.9/h |
|---|---|---|---|---|---|
| half orbit | 1.27e5 | 25 h | 12-37 h | 25 | ~$47 |
| 10 orbits | 2.54e6 | 21 d (~21 x 24 h links) | 10-31 d | ~490 | ~$930 |

Hydro-only reference: 64 min per half orbit.

What moves it:
- **Picard / BiCGStab counts (dominant).** dt = 1.25 s is ~15x smaller than the He giant's ~18 s.
  The stiff zones are the photosphere (thermal time ~1 s, diffusion number ~400 per step) and the
  impact shock; the deep interior is trivially diagonal (D dt/dx^2 ~ 1e-6) and the ambient explicit-
  like (c dt/dx ~ 20). Small steps on a slowly changing state may give 3-6 Picard passes instead of
  15-20, i.e. 0.3-0.4 s/cycle. A radiation-dominated impact shock may push it up instead.
- **dt.** `rad_signal_speed` in the impact zone may cut it 20-40 % (3).
- **vet_col only** (no vet_gd): cheaper build than the He giant's.
- **Thin x2 (4 cells + 2x2 ghosts).** Halo and ghost overhead per active cell is ~2x that of a 64-wide
  block: worse cache use in the implicit stencil sweeps.
- **Levers if too slow** (no accuracy sacrifice per the user rule): nx3 2048 -> 1024 if radiation
  smooths the azimuthal structure (to be shown); 2 GPUs (8 phi blocks split well). Multi-rate
  radiation (`rad_m1_mr.cpp`) exists but is opt-in only, never recommended.

**Calibration smoke** (after stage S3 below; ~4 H200-h total, Caltech testing):
1. **Envelope only, M1, fresh**, 300 cycles, `time/nlim`. Measure s/cycle, Picard mean/max, BiCGStab
   iterations, NON-CONVERGED count, and L(r_meas) drift.
2. **Stream impacting.** Either fresh with the stream for ~15k cycles (the stream needs ~0.02 code time
   = ~11k cycles to reach the photosphere; ~3 h at 0.7 s/cycle), or a hydro-only env13 restart at 0.1
   orbit converted by a remap script (rho, v kept in unit choice (a); eint from T via the table;
   E = aT^4, F = 0; pattern: hegiant `he_remap_giant.py`), run 2000 cycles. The converted start
   readjusts in a few hundred cycles; that is fine for cost, not for science.

Read the steady Picard count and dt from run 2. Then decide: 1 or 2 GPUs, nx3, and how many orbits.

## 7. Staged build order with gates

All new keys default off. Every stage keeps env13 bitwise (hst identical, bin identical apart from the
parameter dump) and is compared GPU vs CPU.

| stage | content | gates |
|---|---|---|
| S0 | IC column script (psi-column from tables; or the MESA remap); table files (EOS X 0.7 to log T 7.3, Rosseland + Planck X 0.7, extended to log rho -17) | column inside the tables; HSE residual < 1e-6; L(r) const to 1e-4; Schwarzschild check reported |
| S1 | EOS interface in the accretor: every `gm1` site of table 1 through `eos.*` (still `eos = ideal`) | env13 + env13_nostream 200 cycles: **bitwise vs 65ec7626** (if the EOS calls round differently, <= 1e-14 rel and the user decides). Old keys bitwise |
| S2 | `eos = general` (table), `<units>`, kelvin keys, `wb_option = polytropic`, discrete-balance IC with the analytic polytrope T(psi), relaxation as today, no radiation | static envelope (no stream) 0.05 orbit: max v < 10 km/s, MR swing < 1e-5 (T2 of NOTE-2026-10-07); stream on 0.1 orbit without dt collapse; Mdot_in = analytic within 1 % |
| S3 | `<rad_m1>`: tables, column IC (S0), Phi_eff + a_ref, split work, flux bottom, marshak top, vet_col, relaxation off inside, kappa_abs mask in the ambient, history columns; **no stream** | (i) startup range checks; (ii) static envelope holds 0.05 orbit as S2; (iii) L(r_meas) and L(r_out) = L_gainer (5.37e4 Lsun) within 5 % after the transient; flux-weighted T_eff 28.3 kK within 3 %; (iv) gas + rad + grav energy changes = boundary fluxes + relaxation booking, to 1e-6 rel per orbit; (v) ambient keeps T_amb and rho_amb (no collapse); (vi) restart bitwise (rst at N vs straight run); (vii) CPU vs GPU per the H200 validation rule; calibration smoke 1 |
| S4 | stream on: window eint from T_don via the EOS, rad ghosts, stream relaxation off | 0.1 orbit without dt collapse or NON-CONVERGED; Mdot_in matches; L(r_out) - L_* vs the released accretion energy G M Mdot (1/R_acc - 1/r_out), booked against the enthalpy advected through r_meas; j_acc vs hydro-only env13 at the same Mdot (the science comparison); calibration smoke 2 |
| S5 | cost decision with the user; optional: donor irradiation (outer incoming flux), Mdot scan (unit choice b) | - |

Code touched per stage:
- S1-S4: ry_per_accretor.cpp only.
- S3: one key in `rad_m1_opacity.cpp`/`rad_m1.cpp` (kappa_abs mask); possibly outputs.hpp
  (NHISTORY_VARIABLES).
- S5: `rad_m1_implicit.cpp` (outer incoming flux).

Rough size: ~600-800 new lines in the pgen, ~50 in the module.

## 8. Risks and open decisions

**Top risks**
1. **The hot ambient vs radiation (2.3).** Without an absorption mask it collapses within one step. With
   the mask it is an artificial, radiatively inert, relaxed region next to a radiatively cooling stream.
   The stream/ambient interface was already the failure point of q10.
2. **Cost and convergence in the impact zone.** The radiation-pressure-dominated shock (P_rad/P_gas ~ 5)
   and the ~1 s photospheric thermal time may need many Picard passes or hit `implicit_resid_fatal`.
   The cost range is 0.35-1.05 s/cycle, i.e. 10 to 31 H200-days for 10 orbits.
3. **Closure and geometry.** vet_gd is impossible on nx2 = 4. vet_col is wrong at the stream flanks and
   the hot spot. The theta-invariant slab removes the stream's vertical radiative losses. All three bias
   the stream's cooling, i.e. the quantity RHD is meant to add over the isothermal runs.

Secondary:
- the IC transient (von Zeipel; possible weak sub-surface convection);
- EOS/opacity table coverage of 5e-16 g/cc and 6.8 MK;
- the history-column limit.

**Open decisions for the user**
1. **IC source:** in-house radiative psi-column from our tables (S0 script) or the Wade+2026 MESA
   gainer profile (needs the file).
2. **General EOS vs ideal gas mu 0.62 + table opacity** (the BSG precedent, simpler). General EOS adds
   HeII ionization just below the photosphere; small effect, small extra cost.
3. **Density unit / Mdot:** unit choice (a) (env13 numbers and restarts reusable) or (b) (fixed star,
   Mdot scan by rho_stream); which Mdot values.
4. **Ambient:** keep the 300 km/s hot ambient with an opacity mask + relaxation, or test a cooler,
   denser ambient.
5. **Outer radiation BC:** vacuum Marshak (default), or add donor irradiation (new module code).
6. **Budget:** full 640 x 4 x 2048 for how many orbits (~1 H200-day per half orbit); reduced nx3; 1 vs 2
   GPUs; where it runs (Caltech is testing only).
7. **Accept vet_col-only**, or a 3-D wedge (nx2 >~ 64, ~16x cells) to make vet_gd possible.

## 9. S3a status (2026-10-08, Caltech, branch accretor-rhd-1008)

Code (ry_per_accretor + one module key): `<rad_m1>` on the envelope (inner = envelope, thermo =
general, env_ic = column). Tables problem/env_opac_table + env_planck_table (module lookup also
recomputes kappa_t along the column). Radiation-force reference: Phi_eff = Phi_wb + G(psi),
G = -int kappa_t F/c dr_eq along the column (F ~ g on equipotentials: G is a function of psi),
a_ref = -dG/dr per cell (face difference), split work in RyPerSrcEnv; the IC march, the inner
ghost walk, gpl_ and wb_phimax all use Phi_eff. Flux bottom (checked against F_col(r_in) to
1e-3), Marshak top, M1FillGhost (copy / vacuum). `<rad_m1>/opac_abs_rho_max` (+ opac_abs_kappa_s):
no absorption below that density, electron scattering kept (applied inside M1TableOpacities, so
at every lookup). History: the 9 rate columns become Lmeas Ltop Erad Etot KE Ebnd Erel Pic Picmax
in rad mode only (NHISTORY_VARIABLES untouched). Input template docs/dev/accretor_rhd/s3a.athinput.in.

Gates (dev tables TOPS X 0.7 Z 0.008, the Z 0.02 TOPS query is pending): env13, env13_cool,
S2 bitwise; restart bitwise (40 vs 20 + 20); GPU = CPU to round-off (10 cycles, nx3 256).
Static run 0.02 orbit (5855 cycles): Lmeas/L 0.963-0.966 (rising slowly), T_eff 28.05 kK, mass
through r_meas 7e-13 Menv, photosphere moved <= 1 cell, v <= 0.5 km/s in the envelope, budget
drift 1e-6 of int L dt after the transient, 0 NON-CONVERGED; 0.49 s/cycle (vet_col 0.32 s = 64 %),
Picard mean 5.6 (2-12), max 14; dt 1.57e-6 (1.09 s).

The 3.5 % luminosity deficit is geometric, not numerical: with T = T(psi) and F ~ g the
equatorial band of the Roche-distorted star radiates <(g/g90)(r/r90)^2/cos alpha> = 0.9605 of
the phi = 90 column's L at the photosphere (0.9875 at r_in), docs/dev/accretor_rhd/
band_luminosity.py. The band cannot shed flux to the poles (2-D slab), so it heats on the
envelope's thermal time until it radiates the imposed L_in (T_eff(phi 90) then +1 %).

### 9.x User decisions 2026-10-08 (after S3a)
- **L_in = band equivalent** for the 2-D equatorial slab: the inner flux is lowered from the phi=90 column's L = 4.68e4 Lsun
  by the von Zeipel band factor (0.9605 with the dev tables; recompute with the final tables via band_luminosity.py), so the
  band radiates what it receives and T_eff stays 28.3 kK. The current flux-BC consistency check against F_col(r_in) must
  allow this factor (a named key, e.g. problem/env_lin_band_factor, default 1 = bitwise). **For a 3-D run (full theta),
  go back to the full column L: the band factor applies only to the theta-thin slab.**
- **Opacity tables:** X 0.70, Z 0.02 Rosseland + Planck requested from viper
  (TASK-2026-10-08-viper-accretor-opacity-tables); rebuild the S0 column and rerun S3a with them.
- **Stream window radiation:** no injection; vacuum/Marshak window (P_rad/P_gas ~0.7 % at the stream peak, 5.1e-8 g/cc,
  26.2 kK; only the low-mass wings reach order unity).

### 9.y S3a final (2026-10-09): Z 0.02 tables, L_in = band equivalent
- **Tables:** viper's X 0.70 Z 0.02 GS98 ext2 pair (Rosseland ad516c7e, Planck 99d668f6; md5 checked against
  docs/handover/accretor-tables-1008/MD5SUMS). Run copies and the column in rhd/s3/final/.
- **S0 column** (make_ic_accretor_column.py, same EOS dump, defaults): L 4.6822e4 Lsun (phi 90 column, unchanged: L =
  4 pi R_acc^2 sigma T_eff^4), rho_ph 6.02e-10 g/cc (dev Z 0.008: 6.86e-10), Gamma_edd(ph) 0.350, Fe bump at r 8.5:
  kappa 2.44, P_rad/P_gas 1.41, Gamma 0.557 (dev: 1.71, 0.79, 0.39); HSE residual 2.2e-6. Schwarzschild-unstable
  (radiative column, convection not included): r 8.961-8.993 (T 31-46 kK, HeII/H, max nabla - nabla_ad 0.047) and
  r 8.584-8.790 (T 113-203 kK, Fe bump, max 0.069; dev: 8.637-8.786, 0.040). The physical column is kept (user).
- **Band factor:** problem/env_lin_band_factor = 0.960504 (band_luminosity.py at psi = 0). It is purely geometric (T =
  T(psi), F ~ g), so the new tables do not change it. implicit_flux_x1min = 0.960504 x 1453030.1277 = 1395641.2498
  code (L 44973 Lsun at r_in, phi 90). The column's flux and the force reference G(psi) are left at the column's
  values: at phi 90 the photospheric flux in the band-equivalent state is F_col (the band deficit comes from the
  other longitudes); deeper, where the band ratio is 0.9875, F(phi 90) is ~2.7 % below F_col, a 0.4 % error in
  the reference force at Gamma 0.14, absorbed by the hydrostatic adjustment. **3-D: factor 1 and the full flux.**
- **Gates (9a706823):** env13 and env13_cool vs 76afb8f4, S2 vs 6e314189: hst, bin and cycle/dt lines identical;
  the S3a dev input without the new key vs ee94815b (= 1377cdb4 code): hst, cycle/dt and bin data identical.
- **Static run (stream off, 0.02 orbit, 5840 cycles, 1 H200, rhd/s3/final/static):** Lmeas/L_in 1.0027-1.0065 after
  t = 2e-4 (last 1.0056), Ltop/L_in 1.0055; T_eff(phi 90 equivalent) 28.34 kK (+0.14 %); photosphere (T_rad =
  T_eff) moved +2e-3 Rsun at phi 90 (< 1 cell), +2e-3 at phi 0; mass through r_meas 1.2e-12 Menv; envelope rms
  |v| <= 0.4 km/s (KE share 3.6e-6 above r_meas); budget drift after the transient -1.3e-6 of int L_in dt;
  0 NON-CONVERGED, 0 floor clips; Picard mean 5.7 (last quarter 3.2), max 15; dt 1.53e-6 (1.06 s); 0.46 s/cycle
  (vet_col 0.32 s = 69 %). The two unstable zones had not developed visible convection in 0.02 orbit (rms v
  0.19 km/s in 8.6-8.8).
- **Diagnostic caveat (pre-existing, same in the dev run):** above the photosphere the CELL-centred F1 (bin output)
  gives 4 pi r^2 F1 = 0.70-0.80 of the face flux that the history (Lmeas, Ltop) and the energy budget use; the face
  flux is the conserved one. Below the photosphere they agree to 1 %.

## 10. S3b: the physical Plaskett stream (2026-10-08/09, Caltech, accretor-rhd-1008)

Estimates: docs/dev/accretor_rhd/stream_physics.py (EOS dump + the Z 0.02 Rosseland/Planck tables);
output in the run directory rhd/s4/stream_physics.out.

### 10.1 Parameters: literature vs env13

| quantity | literature (reference) | env13 / run | note |
|---|---|---|---|
| initial masses | 18 + 16.2 Msun, q_i 0.9 (Wade+2026 Sect. 5); Fig. 3 caption 18.2 + 16.8 | M_d 18, M_a 16 | Wade's own text and caption disagree; Sect. 6.4: gainer "M ~ 16 Msun, R ~ 9 Rsun" when it starts accreting |
| initial period | 3.81 d (Sect. 5); "a 3.7-d orbit" (Discussion) | 3.69 d (a 32.5575 Rsun) | 3.81 d gives a 33.27 Rsun (+2.2 %): every geometric quantity (d_L1, r_out, grid, column, band factor) would move; open decision |
| Mdot at onset | ~1e-4 Msun/yr (Sect. 6.4) | 1e-4 | same |
| gainer Teff, L; donor Teff, L, R | not in the text; read from the HR figure (env13: 28.3 kK, log L 4.73; donor 26.2 kK, log L 4.84, R = R_L 12.8) | same | figure read, not tabulated; the Zenodo MESA files would settle it |
| mass transfer | conservative (assumed by Wade+), non-rotating models | spin 1 | the run's spin 1 is our choice |
| L1 Mdot factor | Ryu+2025 Eq. 4: isothermal 0.721 - 0.149 tanh^2(0.522 log q) = 0.721; adiabatic 0.649 + ... = 0.649. Ryu+2025 tie the ADIABATIC case to optically thick overflow (photosphere outside the Roche lobe) and the isothermal one to optically thin overflow | 0.721 (isothermal) | Plaskett at 1e-4 Msun/yr overflows from tau ~1.5e3 (10.2): Ryu's adiabatic case applies (Q 0.649, FWHM 1.2x narrower, peak 1.1x); env13 used the isothermal case |
| L1 structure | Ryu+2025: isothermal: M 1.1-1.3 at L1, vertical HSE, in-plane expansion at 0.2-0.3 of the overflow speed; adiabatic: M 0.9-1, FWHM 1.2x narrower, peak 1.1x the analytic profile | Gaussian widths sigma_y 0.932 c/(Omega sqrt B), sigma_z c/(Omega sqrt C) at T = 26.2 kK | Ryu+2025 simulate the L1 region only (no radiation); nothing in either paper fixes the stream T |
| stream T | none | isothermal at the donor Teff 26.2 kK | see 10.2: the L1 gas comes from tau ~1.5e3 in the donor atmosphere |
| r_out state | none (derived) | rho 5.11e-8 g/cc, width 0.735 (arc), v (-122.9, 55.5) km/s, T 26.2 kK | ballistic orbit and pressure-supported widths (env13 RESULTS.md sect. 2-3) |

### 10.2 What the real stream looks like (estimates)
- **env13's stream is extremely optically thick** with the real tables: at r_out kappa_R 83 cm^2/g (H-ionisation opacity
  peak at 26 kK, 5e-8 g/cc), tau across one sigma 2e5 (in plane) and 1.3e5 (vertical); diffusion time 2.6e3 flights
  (flight L1 -> r_out 3.5e4 s, r_out -> impact ~1.2e4 s). Inside the nsig 3 window the column outside any point
  has tau >= 13; tau = 1 is reached only at 3.5 sigma (rho 1e-10 g/cc), outside the window. Thin gas cools in
  ~1e-5 s (kappa_P 2.8e3), i.e. it is always at radiative equilibrium.
- **An isothermal stream at the donor Teff is not consistent with Mdot 1e-4.** Optically thick overflow (Kolb & Ritter
  1990 form with Ryu's factor, Mdot = Q 2 pi/(Omega^2 sqrt(BC)) int c_T dP_gas over the donor's grey atmosphere, EOS +
  Rosseland table, mean donor g 3.0e3) needs the overflow to reach tau ~1.5e3: T_L1 ~ 151 kK, rho_L1 5.2e-8 g/cc,
  P_rad/P_gas 1.2, c_T 45 km/s, sigma_y,z ~1.15 Rsun (2.4x env13's). With g reduced near L1 (x0.3): 250 kK. The
  isothermal-at-Teff estimate puts 6.2e-7 g/cc at L1 (12x denser, 2.4x narrower).
- **Along the flight the core is neither isothermal nor adiabatic:** t_diff/t_flight is 3.6 at L1 and 0.4 at r_out for the
  adiabatic state. Bounds at r_out (pressure-supported widths, Mdot conserved): adiabatic 77 kK, rho 8.5e-9 g/cc
  (0.17 code), arc width 1.81, P_rad/P_gas 1.0; one-zone adiabatic + diffusion model 65 kK, rho 1.6e-8 g/cc (0.31 code),
  arc width 1.33, sigma_z 0.81, P_rad/P_gas 0.33, tau across 4e3. The surface is at the thin radiative-equilibrium T,
  19.4 kK with both stars, 16.9 kK with the gainer only (vacuum r_out).
- **Consequence for the no-injection decision (9.x):** its premise (P_rad/P_gas 0.7 % at the peak) holds only for the
  26.2 kK isothermal stream. The physical stream has P_rad/P_gas 0.3-1 at r_out; injected as gas only, its gas
  would hand ~40 % of its internal energy to the radiation field in the first cells. This is a decision for the user.
- **Thin wings and M1:** the window carries only the tau >= 13 part. Where gas is thin (stream surface, window edges,
  hot ambient interface) the implicit coupling (BE, Newton gas update) puts the gas at the local radiation temperature
  within one step (t_cool ~1e-5 s << dt ~1 s); vet_col sets the closure along radial columns, so the stream's
  sideways (phi) radiation is closed with the radial Eddington factor (wrong at the flanks, design sect. 8 risk 3);
  the slab has no vertical losses (tau_z 1e5 makes that a small error for the core, but not for the surface layer).
  The hot ambient is masked (electron scattering only), so it does not exchange energy with the stream.
- **Donor irradiation** (point-source dilution of the donor disc, 26.2 kK, R 12.8): F_irr/F_gainer = 0.23 at the
  substellar point (phi 0), 0.21 at the impact longitude (14 deg), 0.09 at 45 deg, 0 beyond ~75 deg. It raises the
  facing photosphere T by up to (1.23)^(1/4) = +5 % and the stream's surface equilibrium T from 16.9 to 19.4 kK.
  Implementing it needs an incoming-flux term in the module's x1max Marshak BC (rad_m1_implicit.cpp, phi-dependent),
  i.e. module code; not done. Recommendation: not needed for the first impact/accretion runs (the impact ram pressure
  is ~2e3 P_ph and the core is opaque), needed before quoting the facing-side photosphere or the stream's cooling.

### 10.3 S4 calibration smoke (2026-10-08/09, 98b98ab4, 1 H200 unless stated; rhd/s4/)
Restart from the S3a-final static state (t = 0.00916) with problem/stream = true; stream radiation not injected,
env_t_relax_stream = 0 (enforced). Two stream states at r_out (same ballistic v, Mdot 1e-4 at the midplane):
**env13** (26.2 kK, rho 5.11e-8 g/cc, arc width 0.735; P_rad/P_gas 0.007, tau per arc sigma 2.2e5) and **hot**
(the one-zone estimate of 10.2: 64.6 kK, rho 1.56e-8 g/cc = 0.306 code, width 1.328; P_rad/P_gas 0.33, tau 4.7e3;
gas only, so it violates the no-injection premise; run as a sensitivity/cost case).
- **Stream approach** (both): front r 13.06 / 12.58 / 12.01 / 11.35 / 10.59 at t = 0.0117 / 0.0142 / 0.0167 /
  0.0192 / 0.0217, max |v| 180 -> 365 km/s; Picard mean 12-13 (static 3-6), max 15-17; dt 1.50-1.53e-6 (unchanged:
  r-limited); 0.63-0.68 s/cycle (static 0.46).
- **env13 stream: Picard DIVERGED at t = 0.02411 (cycle 15636)**, resid 3e8 after 200 passes, worst cell
  (k 51, i 450): r 9.585, phi 9.05 deg, i.e. the MASKED HOT AMBIENT (rho 3-5e-8 code = 2e-15 g/cc, below the
  absorption mask 1e-6 code; T_gas 6.6-6.9 MK, T_rad 22 kK, v_r -60..-80 km/s) being compressed ahead of the
  stream's leading edge just outside r_meas (dt had dropped to 1.1e-6). Reproduced bitwise from the t = 0.01916
  restart. With rad_m1/implicit_res_dmin = 1e-6 (masked cells out of the stopping test) the rest converges but
  the masked cells themselves diverge (masked_resid 1.4 > implicit_resid_fatal_masked 1, same cycle). This is
  design risk 1 (2.3, the hot ambient under radiation) materialising; not fixed (module work).
- **hot stream: impact reached, no NON-CONVERGED** (to t = 0.02816, wall limit): front at r 9.69 (t 0.0242) and
  9.13-9.24 (t 0.0267-0.0282, below the local photosphere); impact shock at r 9.47, phi 14.5 deg: rho 0.11 code,
  T 103 kK, P_rad/P_gas 3.7 (10 at t 0.0267); MR (mass in through r_meas) 0.014 code by t 0.028. In the impact
  phase dt min 9.3e-7, mean 1.29e-6 (t 0.022-0.025), 1.38e-6 (0.025-0.028), 1.44e-6 after; Picard mean 13.5-14,
  max 15; 0.62-0.67 s/cycle.
- **2 H200** (2 ranks, restart t 0.01916, hot, 1450 cycles): 0.365 s/cycle vs 0.672 on 1 H200 over the same cycles
  (1.84x; same dt sequence); static fresh start 0.256 vs 0.42 (1.64x). **nx3 1024** (static, 1 H200): 0.234 s/cycle
  (0.557x of 2048), dt unchanged (r-limited).
- **Cost extrapolation** (impact-phase dt 1.35e-6 code, 0.65 s/cycle at nx3 2048 on 1 H200; P_orb 0.458266):
  | config | s/cycle | half orbit (wall) | 10 orbits (wall) | 10 orbits GPU-days |
  |---|---|---|---|---|
  | 2048, 1 H200 | 0.65 | 30.6 h | 25.5 d | 25.5 |
  | 2048, 2 H200 | 0.353 | 16.6 h | 13.9 d | 27.7 |
  | 1024, 1 H200 | 0.362 (static ratio) | 17.1 h | 14.2 d | 14.2 |
  | 1024, 2 H200 | ~0.23 (assumed 1.6x, unmeasured) | ~10.7 h | ~8.9 d | ~17.8 |
  The impact shock lowers dt further as it deepens (min 9.3e-7 seen); treat these as lower bounds.
