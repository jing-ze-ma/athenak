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
