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

## 11. S5: 3.81 d, physical ambient, hot optically thick stream (2026-10-09, Caltech, accretor-rhd-1008)

User decisions 10-09: (1) physical ambient (floor gas at radiative equilibrium, absorption on, no
mask), (2) the hot one-zone stream with its radiation injected (E = a T^4, F = (4/3) E v; reverses the
10-08 no-injection decision), (3) P = 3.81 d. Derivation (inputs, approximations, numbers):
docs/dev/accretor_rhd/plaskett_s5_1009.md; runs rhd/s5/.

### 11.1 Code (keys default off, bitwise)
- `problem/env_amb_mode = radeq` (default hot): IC atmosphere ends at rho_amb, ambient = rho_amb at T_col;
  requires env_t_relax = 0, no opac_abs_rho_max, hydro/dfloor = rho_amb. The pgen applies the floor itself
  after the sources (added gas at the cell's v and T): mass booked in hst column 8 `Mfl` (replaces Jstr),
  energy in Erel. The existing floor sponge stays.
- `problem/stream_rad = true` (default false): window M1 ghosts E = a T_s^4, F = (4/3) E v; in the implicit
  solve (face-flux x1 BC, ghosts not read) the window columns of the outer face are an optically thick
  INFLOW face (RadiationM1::SetX1maxBathColumns): comoving flux 0, lab inflow (1 + chi) v E_s; elsewhere
  Marshak/vacuum. The module stores the lab advective radiation flux through r_out (ibadv_x1max), booked in
  Ebnd. First version (29d5a6c4) used E_s as a Marshak bath, c q (E - E_s): an isotropic 65 kK wall
  shining ~450x the advected radiation into the thin ambient; it diverged in the first step (rhd/s5/
  smoke_bath_fail). Final: 3f37ebf8.

### 11.2 Status at the mem-1009 merge (interim, 3f37ebf8)
- Gates (29d5a6c4, new keys unset): env13, env13_cool vs 76afb8f4, S2 vs 6e314189, S3a-final input vs
  9a706823: hst, bin, cycle/dt identical. GPU = CPU to round-off (nx3 256, 10 cycles, radeq + stream_rad,
  stream on from t = 0; both 29d5a6c4 and 3f37ebf8).
- Static (radeq, stream off, 0.02 orbit, 1340 cycles, rhd/s5/static): Lmeas/L_in 1.0031-1.0060 (last 1.0049),
  T_eff(phi 90 eq) 28.33 kK, photosphere moved +1-2e-3 Rsun (< 1 cell), MR 4e-14 Menv, budget drift -1.5e-6
  of int L dt after the IC transient, mass budget closes (Mfl 1.2e-9), floor gas T_gas = T_rad to 1e-4,
  0 NON-CONVERGED, Picard 2.1; dt 7.04e-6 (4.6x S3a's 1.53e-6: the hot ambient set dt), 0.49-0.58 s/cycle.
- Stream smoke (rhd/s5/smoke, restart t 0.00946): approach fine (Picard 11-15, dt 7.03e-6, 0.95 s/cycle,
  budget 4e-7, compressed ambient T_gas/T_rad = 1); the stream core cools 63 kK (r_out) -> 52 (r 12) -> 35
  (r 11) -> 23 kK (front). Picard DIVERGED at t = 0.02428 (resid 0.58) at r 9.62, phi 6 deg, where the
  stream head meets the top of the atmosphere near the L1 cap: the same place and time as the S4 env13
  failure. Before it, the whole column from the photosphere top (r 9.465, f_rr -> 1) to the stream head
  collapses radiatively 20 -> 9.5 kK in 2e-4 code. Not cured by cfl 0.15 or Anderson acceleration.
  Also not cured by implicit_flux = blend / blend + implicit_flux_faces = all (mem-1009 code; both diverge in
  the first step after the switch); ap_hll/berthon/closure_lag pass are not available with transport =
  implicit on this mesh.
- **Gas-only control (stream_rad = false, same hot one-zone gas, radeq ambient; rhd/s5/smoke_gas,
  3f37ebf8): through the first impact, 0 NON-CONVERGED** to t = 0.0300 (6824 cycles): the window gas hands
  its radiation to the dark window at once (T 65 -> 15.6 kK at r_out, P_rad/P_gas 0.01), the cold stream
  hits at phi ~15 deg: shock T 84 kK (t 0.0284) -> 110 kK (t 0.0300) at r 9.43, P_rad/P_gas 9 -> 3.4, rho
  0.15 code; MR 0.0206 code by t 0.030 (Min 0.40); budget -2e-6 of int L dt, mass budget closes; Picard 12.4
  in the impact phase, dt 1.52-1.64e-6, 0.65 s/cycle; compressed ambient ahead of the stream T_gas/T_rad =
  1.000 (0.997-1.007). So the physical ambient converges where the S4 masked ambient diverged; what fails is
  the radiation-rich (injected) stream head meeting the atmosphere.
- **mem-1009 merge (ce3cc213):** env13, env13_cool, S2: hst, bin, cycle/dt identical to the refs; S3a-final
  input vs 9a706823: bin (hydro + m1) and cycle/dt identical, history sums differ in the last digit (rad-mode
  reductions, e.g. Lmeas row 0 5.165626921119900e4 vs ...899e4, hydro momenta 1e-16 relative); the injected
  smoke diverges at the same cycle with the same residual (bitwise).
- **Cost (impact phase, gas-only stream, dt 1.6e-6; half orbit 0.2366 code = 1.48e5 cycles):** 2048 x 1 H200
  0.655 s/cycle -> 26.9 h per half orbit, 22.4 d per 10 orbits; 2048 x 2 H200 0.356 s/cycle (1.84x) -> 14.6 h,
  12.2 d wall (24.4 GPU-d); nx3 1024 x 1 H200 (static ratio 0.51, dt radial-limited, unchanged) -> 13.7 h,
  11.4 d; 1024 x 2 H200 (1.6x assumed) ~8.6 h, ~7.2 d wall. Approach phase: dt 7.0e-6 (4.6x S3a) at
  0.95-1.05 s/cycle (Picard 11-15). The injected stream ran at Picard ~15 (vs 12 gas-only).

### 11.3 Option (a): transverse T profile (user 10-09), a844c97d
- `problem/stream_t_profile = diffusion` + `stream_t_surf` (default uniform, bitwise): T^4 = max(T_c^4
  cos(pi/2 erf(|x|/sqrt2)), T_eq^4), x = dphi/sigma (lowest cooling mode of a slab in optical depth at
  constant opacity; T_eq 19239 K, both stars); gas and injected radiation follow the local T. Window
  T 64.8 / 53.9 / 33.5 / 24.2 / 19.2 kK at 0/1/2/2.5/3 sigma; P_rad/P_gas 0.40-0.52 inside 2.5 sigma (1.0 at 3);
  injected radiation 2.8x lower (0.37 of the gas enthalpy flux). Input adds implicit_pos_floor_solve = true.
- Gates (a844c97d): env13, env13_cool, S2 identical to the refs; S3a-final and an S5 uniform stream_rad input
  (nx3 256, stream on, 20 cycles) identical to ce3cc213 (hst, bin, cycle/dt).
- Smoke (rhd/s5/smoke_prof, restart of the static state): approach as before (Picard 7 -> 15, dt 7.03e-6,
  budget 1.7e-7, mass closes); core 62.6 (r_out) / 50.5 (12) / 33.0 (11) / 22.8 kK (10). **Picard DIVERGED at
  t = 0.02432** (resid 3.4e14) at r 9.50, phi 8.5 deg: the same head/atmosphere contact. Diagnostic dumps
  (crash_prof): between t 0.02406 and 0.02431 the column from the photosphere top (r 9.45, where f_rr -> 1
  and the cell-centred L collapses) up through the arriving head (rho 3e-5..1e-3 code, v -420 km/s) cools
  21 -> 10 kK with T_gas = T_rad; the head above radiates inward (F < 0). Stopped per instruction; module
  work (option d) not started.

## 12. New binary: fin51-r-fix merged, dt_retry + ROL2U on the S5 stream smoke (2026-10-15, Caltech, accretor-rhd-1015)

User 10-15: "continue the accretor work with the new binary"; tests only (budget ~15 H200-h, used 0.97).
Runs: /resnick/groups/carnegie_poc/jingze/accretor_1015 (gates/, smoke/, JOBS.txt); builds nofma HOPPER90
(bsg_1009/build_inc_nofma.sh, BUILD_ROOT builds_accr): fdb18517 md5 4e89a437..., f960ed04 md5 c2ed9dda...,
reference a844c97d (old accretor code, same tree and flags) md5 8b004c6e....

### 12.1 Merge (eb880b72, f960ed04, fdb18517)
- origin/fin51-r-fix b16b66db merged into accretor-rhd-1008 7fedd194. Two conflicts, both unions: rad_m1.hpp
  members (per-column outer bath iebath_hi_col/ibadv_x1max + littests-1011 boundary beams ibm_*) and the
  ImplicitSolve face flux (the inner beam face first, then the Marshak branch with the per-column inflow face).
- dt_retry visitors regenerated (tools/dtretry_gen.py, clang 16.0.6): + hyd wb_phimax, m1 otab.amask_rho/kes,
  iebath_hi_col_on, iebath_hi_col, ibadv_x1max; nothing else changed against the fin51-r-fix files.
- ry_per_accretor per-cycle state: stage_ctr_ and the 13 flux/budget registers acc0_/acc1_ advance in every
  stage, so they are registered with dtretry::RegisterValue (a failed attempt was booked twice in Min..Ebnd/Erel).
  ibadv_x1max (read by the hydro-stage budget from the previous solve, like f0x1) is in the m1 minimal list
  (driver/dt_retry.cpp; empty elsewhere). Setup-only arrays (rbath_, tcol_, rhoamb_, cref2_, gpl_) and
  output-time state (hist_prev_, rad_pn0_/ps0_) need nothing.

### 12.2 Gates (1 H200, debug QoS, 20 cycles unless noted)
| gate | result |
|---|---|
| keys absent vs a844c97d: env13 | state bitwise (bin, cycle/dt); hst: the new code writes no extra final row at nlim (t = 0 row identical) |
| keys absent: S3a-final, S5 uniform, S5 profile | NOT bitwise: fresh-run defaults changed since 10-08 (implicit_gas_src conservative2, implicit_rad_inertia, implicit_nofb_work, implicit_f0_init auto, vet_bc_diff face, vet_seg jump, vet_source_noesrc, implicit_face_opac_n, implicit_pos_floor_solve). S5 profile fresh: dt 7.05e-6 -> 2.4-2.7e-6 from cycle 7 (stream on from t = 0) |
| same, those 9 keys pinned to the old values | S3a-final bitwise (bin, cycle/dt); S5 uniform / profile: cycle/dt identical, bin differs at round-off (max rel 1.2e-7 in E, eint after 20 cycles; Lmeas/Ltop first differ at cycle 9); source not isolated (code change in the merged rad_m1, no key) |
| restart (S5 profile deck, 20 straight vs 10 + rst + 10) | bin bitwise, hst state columns bitwise (budget registers restart at 0 as before) |
| forced retry at cycle 5, dt_retry_force_fac 1 | hst, bin, cycle/dt bitwise to no retry |
| same + dt_retry_views = all + poison of every changed non-minimal View except m1.iw, m1.fref_wacc (19 Views) | bitwise |
| forced retry fac 0.5, minimal vs + dt_retry_views_add = m1.iw,m1.fref_wacc | bitwise (they are scratch) |
Poison notes: fref_wacc is reset by acc = g0 acc + w with g0 = 0 at stage 1, so NaN survives as 0*NaN (scratch,
not a bug); iw likewise is fully rewritten (fac 0.5 test above), poison hits a NaN-propagating read.

### 12.3 S5 option (a) stream smoke on the new binary (restart of rhd/s5/smoke_prof rst 00005, t 0.02366)
Keys: time/dt_retry = true + the full ROL2U set (rad_m1/implicit_tol_rel 0.1, implicit_tol_rel_max 1e-4,
implicit_opac_lag pred, implicit_opac_lag_passes 2, implicit_gas_newton_u4 true). New fresh-run defaults stay
off on this restart (old values kept). Binary f960ed04 (= fdb18517 bitwise on the gates; only the budget
checkpoint differs). smoke/smoke_a, job 4591413.
- **Passes the first impact and runs to t = 0.0300 (cycle 7370), no FATAL, NON-CONVERGED 0.** 21 dt retries,
  all `eint<=0` (5-216 cells), all in t 0.02432-0.02443 (the 10-09 divergence time), one cut each except chains
  of 2-3 (dt down to 1.8e-7), cap released after ~10-20 cycles.
- Picard mean 4.15 per solve (10-09: 12-15), max 15; tol_eff mean 1e-4 (tol_rel); BiCGStab 13.5 inner its.
- 0.51-0.56 s/cycle (10-09 approach 0.95-1.05); vet_col tensor builds take 1332 of 2159 s (333 ms each, 62 %).
- Column at the 10-09 failure cell (k 50, i 440; r 9.55, phi 8.9): T 21.5 kK (t 0.02406) -> 17.8 kK (t 0.02426,
  head arriving, 10-09: 10 kK and divergence) -> 19.7 kK under the head (t 0.02446, rho 1-3e-4 code, v -425
  km/s), T_gas = T_rad to 4e-4.
- Impact: hottest dense gas r 9.43 phi 15: 66 kK (t 0.0261) -> 119 kK (t 0.0300), rho 0.116 code, P_rad/P_gas
  5.4 (gas-only control: 110 kK, 3.4); MR 0.0322 code by t 0.030; Lmeas/L_in 1.01 -> 1.11.
- Budgets: mass closes (residual -1.2e-9); energy residual 1.0e-5 of int L_in dt (max 2.2e-4 after the transient;
  10-09 approach 1.7e-7): the looser tol_rel solve.
- Variant without ROL2U (dt_retry only, smoke_b_norol): fails at the same cycle 3642 (t 0.024316): 4 cuts to
  1.0e-7, then dt 6.3e-8 < time/dt_min 1e-7 FATAL at cycle 3654. So ROL2U is what carries the impact; dt_retry
  absorbs the 21 eint<=0 events on top. The smaller-CFL variant was not needed.

## 13. Bigger domain: the stream as a volume source, donor Roche-lobe sink (2026-10-15, Caltech, accretor-bigdom-1015)

User request 10-15 (option a): put r_out well beyond L1 and inject the stream inside the domain, not as the r_out
boundary window. Branch accretor-bigdom-1015 off accretor-rhd-1015 (af45fa46); keys default off (old decks bitwise).

### 13.1 Grid: the old grid kept face for face, then a log tail (mesh/f_stretch_r_tail_*)
- New StretchRPoly option: `f_stretch_r_tail_nin` (old nx1 640), `f_stretch_r_tail_x1max` (old x1max 13.7926),
  `f_stretch_r_tail_ij` (junction face 520). Faces 0..520 are the old map evaluated at xi_old = i/640 (agree with
  the 640-cell grid to 5e-15); above r_j = 10.4157 the faces grow by a constant ratio up to x1max:
  r = r_j (x1max/r_j)^((xi - xi_j)/(1 - xi_j)). The old grid is fully coarse (dr 0.0278) from face ~515, so r <= 10.4
  (the whole envelope, the 8.756-9.598 plateau, the WB region r < 8.605, r_meas 9.573) is unchanged.
- x1max = 1.5 d_L1 = 24.33996, nx1 = 840: tail ratio 1.002656 vs the old grid's 1.002675 at the junction (dr steps
  down by 0.45 % once, then grows smoothly). dr = 0.0366 at r 13.8 (old 0.0282), 0.0645 at r_out (dr/r 0.00265, the
  phi cell is r dphi = 0.00307 r: cells near square).
- Cost: 840/640 = 1.31x cells per cycle; dt is set near the star (unchanged grid there), so ~1.31x per unit time
  (measured in 13.4).
- Restart from the old 640-cell rst is not possible (different mesh) -> fresh static phase.

### 13.2 Stream nozzle (problem/stream_inject = volume, stream_r, stream_ncell)
Choice: an internal Dirichlet region (a "nozzle") rather than an explicit mass-source term S = Mdot/V.
- The nozzle is stream_ncell = 2 cells above the x1 face nearest stream_r = 13.7926 (the old r_out; face 626 at
  13.7974, nozzle r 13.7974-13.8708), inside the same phi window |dphi| < 3 sigma. After all other sources of every
  stage these cells are reset to exactly the old ghost state: rho_s exp(-x^2/2), the ballistic (v_r, v_phi) at
  stream_r (unchanged: phi 3.410 deg, v_r -126.46, v_phi 60.68 km/s), T(phi) of the diffusion profile, and with
  stream_rad E = a T^4, F = (4/3) E v.
- Why: the face below the nozzle then sees the same states the old r_out face saw (PLM slopes vanish across the
  two equal cells, as across the two ghosts), so the inflow and its radiation are the old ones by construction,
  including the Riemann-flux response if the gas below pushes back. An explicit source would fix Mdot but spread
  the mass over the source cells and need an injection velocity mixing rule; the reset needs neither. The stream
  is supersonic (Mach ~4, v_r + 3 c_s < 0), so the nozzle does not leak outward.
- Booking: the mass, inertial Lz, gas energy and radiation energy the reset adds per stage are Min, Jin and
  (negative) the hydro / radiation out registers of Ebnd, so Mdom and Etot budgets close as before. Expected Min
  rate = sum rho_s e^{-x^2/2} |v_r| r^2 dtheta dphi = 19.80 code/time (old run start0: 19.80).
- Outer boundary: no window; outflow everywhere (v_r = max(v_r, 0)), problem/outer_rho = copy (now honoured under
  thermo = general: ghost rho and T copied; the hse walk in the true Roche potential would pile e^cap densities
  into ghosts near the donor). Radiation: Marshak with a dark per-column bath (= the old top outside the window).
  Mout = all of r_out.
- problem/env_r_top_max = 13.7926 caps the WB top radius: beyond L1 the Roche potential falls below the envelope's
  top cap again inside the donor lobe, and an uncapped Phi_wb would give psi > 0 there (envelope T(psi), envelope
  radiation IC near the donor). With the cap Phi_wb is flat above the old r_out: everything r <= 13.79 is identical.

### 13.3 Donor sink (problem/donor_mask = roche) and donor irradiation (problem/donor_irr)
- Mask = the exact equatorial Roche lobe of the donor: Phi(r, phi) < Phi(L1) and x = r cos phi > d_L1 (the
  equipotential, no sphere approximation needed: RochePot is already there and the band is the orbital plane; the
  condition x > d_L1 picks the donor's lobe, L2/L3 lie outside r < 24.34). It touches L1 at r = 16.23 and reaches
  |phi| = 7.9 / 13.7 / 20.3 deg at r = 18 / 20 / 24.34.
- Reset after every stage to the floor ambient: rho_amb, at rest in the rotating frame, T = the IC's ambient
  T_col (the capped column top). Removed mass and gas energy (incl. KE, rho Phi) are booked in hst columns
  Mmsk, Emsk (replacing Mwal, Jwal, which are identically 0 with the closed free-slip wall). Radiation untouched.
- Why a reset sink: the simplest robust choice; the gas has v = 0 there so no CFL limit; inflow from L1 is
  absorbed exactly like a donor surface (no momentum/energy feedback), and a reset at floor density never creates
  pressure gradients that push gas out.
- donor_irr = true (default false; needs donor_mask and <rad_m1>): the lobe also holds E = a T_d^4, F = 0
  (isotropic black body, comoving) every stage, T_d = problem/donor_teff [K], default 26200 K (the donor T_eff of
  the Wade+2026 MESA track in the S5 deck comment; no luminosity key). The light leaves the lobe surface through
  the M1 transport (a transparent region at E = aT^4 radiates ~ sigma T^4 per area into the half space); the
  radiation energy the reset adds is booked in Ebnd as radiation in.

### 13.4 Test (1 H200, deck accretor_1015/bigdom/bigdom.athinput)
NOT RUN on Caltech (user 10-15: handed to Raven, TASK-2026-10-15-raven-accretor-bigdom.md on bsg-files-1009).
Status: code c22b696e builds (Caltech GPU nofma, H200; binary builds_accr/bin/athena_gpu_ry_per_accretor_c22b696e_nofma)
but has NOT been executed: no regression gate (old decks bitwise with the new keys absent) and no run of the
nozzle, mask or donor_irr yet. Deck + the S5 column and opacity tables: docs/handover/accretor-bigdom-1015/
(replace @HDIR@). Plan: fresh static phase (stream false) to t 0.00946, then restart with problem/stream=true to
t 0.025 (first impact ~0.024-0.026), ROL2U keys and time/dt_retry as in orbit1. Reference (start0/orbit1, 640
cells, boundary window): Min rate 19.795 code/time (analytic 19.80), dt 7.03e-6 before impact, 0.51-0.56 s/cycle
during the approach; ballistic at R_acc phi 16.5 deg; impact hottest gas r 9.43 phi 15 deg near t 0.026.
Pass checks: no FATAL; dMin 19.8 +- 1 %; stream path as in start0 (bin snapshots at t 0.0189, 0.0237); Mmsk small
and the mask quiet; dt not set in the nozzle or the mask; cost per cycle ~1.31x.

### 13.5 User decisions 10-15 (84b55d97): L1 nozzle, donor luminosity, extra history columns
**Stream from L1 (problem/stream_inject = l1; supersedes the 13.2 nozzle at 13.79 for production, which stays
available as stream_inject = volume).** Lubow & Shu (1975) state at L1, the same reset-nozzle mechanism as 13.2:
- Nozzle = stream_ncell (2) cells whose OUTER face is the x1 face nearest d_L1 (face 687, r 16.2165; nozzle r
  16.130-16.2165), phi window |phi| < 3 sigma_phi about phi = 0.
- State: isothermal sonic flow along -x, v = -c_T e_x (v_r = -c_T cos phi, v_phi = +c_T sin phi, rotating frame),
  c_T = 45.45 km/s (the launch speed of t_don, problem/stream_l1_v); uniform T_L1 = 155.3 kK (problem/stream_l1_t,
  the stream_physics_1009 L1 state); E = a T_L1^4, F = (4/3) E v.
- Width: the in-plane hydrostatic Gaussian of the L1 saddle, sigma_y = c_T/(Omega sqrt(B)), B = A - 1, A = (1-mu)/d_a^3
  + mu/d_d^3 = 7.994 (a = 1 units), B = 6.994: sigma_y = 3.4228/2.6446 = 1.294 Rsun (sigma_phi 0.0798 rad = 4.57 deg
  at d_L1), narrower than the 1.536 Rsun arc width at 13.79 (problem/stream_l1_width overrides).
- Mass rate: the peak density follows from problem/stream_mdot = 19.795 (the orbit1 Min rate) as the exact discrete
  sum of rho c_T cos(phi) A over the window cells of the nozzle's inner face (the face flux of a supersonic/sonic
  state is its upwind flux). rho_L1 = 0.675 code (3.45e-8 g/cc) vs the 3-D one-zone 5.25e-8 g/cc: the band carries
  the old window's rate per unit height, not the 3-D Gaussian. P_rad/P_gas ~2 at L1 (radiation is in M1, the EOS is
  gas only); the gas then expands along the ballistic path (the r = 13.79 window state of 13.2 was this path's one-zone
  result).
- Booking: Min = mass added by the reset (closes the mass budget); the new column **Mnzf** = mass through the
  nozzle's inner face (RK flux integral) = the stream that leaves L1 toward the star; Min - Mnzf = what leaves the
  nozzle backward/sideways (expansion over the saddle; into the mask = Mmsk). The pass check is Mnzf rate = 19.795.
- Mask vs nozzle: nozzle cells have r < d_L1, so x = r cos phi < d_L1 and the mask (x > d_L1) never contains them
  (the kernel also gives the nozzle precedence). The lobe's cusp touches the nozzle's outer face at phi = 0 only; no
  gap was added (if the back-leak into the mask is large, a mask offset x > d_L1 + gap is the knob to add).

**Donor irradiation by luminosity (problem/donor_lum [Lsun], default 10^4.84 = 69183, the log L of the S5 deck header
from the Wade+2026 track).** T_d from sigma T_d^4 = L_d / (4 pi R_L^2), R_L the volume-equivalent Roche-lobe radius
(Eggleton 1983, q = M_d/M_a = 1.125: 12.943 Rsun) -> T_d = 26.0 kK (vs the 26.2 kK HR read). Area definition: the
band is an equatorial slice of the 3-D donor, so it must carry the PER-AREA flux of the whole donor; L_d is the
full-sphere luminosity and 4 pi R_L^2 the full donor surface (the same band <-> sphere convention as Lmeas/Ltop,
which scale the band by 4 pi / band solid angle). problem/donor_teff overrides when given (fallback). On in the
handover deck (intended on for production). Column **Eirr** = cumulative radiation energy the lobe reset adds.
Caveat: E = a T_d^4 isotropic in a transparent (floor-density) lobe emits ~ sigma T_d^4 per area through its edge
only approximately (M1 closure at the edge); check the band's donor luminosity from Eirr/Ltop in the test.

**History**: with a nozzle or the mask the user hst has 24 columns (23 Mnzf, 24 Eirr; NHISTORY_VARIABLES is 24);
older decks keep 22 and their dt_retry registration list.

### 13.6 Regression gate and smokes on Caltech (10-15, 1 H200, debug QoS; accretor_1015/gates/g6-g10)
**Regression gate g6 (84b55d97 vs the orbit1 binary fdb18517, new keys absent, 20 cycles): BITWISE** on all five:
env13, S3a-final, S5 uniform, S5 profile (fresh, nx3 256 as gate 1) and the orbit1 lineage restart (smoke_a rst
00009 + the orbit1 ROL2U/dt_retry keys): hst files identical, bin files identical as whole files, cycle/time/dt lines
identical. The later commits (5fc471d5, 82197717) change only code under accx_on_ (nozzle or mask keys), so this
carries over; not re-run.

**L1 nozzle, first attempts (g6-g9): dt collapse at cycle 2-3** in every variant (mask on/off, stream radiation on/off,
T_L1 155 or 64.8 kK): the sonic L1 state expands BACKWARD (v + 3 c_s > 0) into the near-vacuum beyond L1; cells
there reached 1e5-1e9 km/s with garbage energies (the 13.79 nozzle is supersonic, v + 3 c_s < 0, and runs). Two
fixes, both active only with the new keys:
- 5fc471d5: the radeq floor keeps momentum (v = m/rho_amb) and caps |v| at env_v_floor_max (old: v = m/du, which
  blows up when du << rho_amb). Not sufficient alone.
- 82197717: with stream_inject = l1 the sector BEHIND the nozzle (r above its outer face, |dphi| < 3 sigma + 1
  sigma) is part of the donor sink (that space is the donor's atmosphere, which is not modelled).

**L1 smoke g10 (82197717, the handover deck, donor_irr on, stream on from t = 0, no static phase, 24 min):**
1563 cycles to t 0.0078, no FATAL, 1 NON-CONVERGED, no dt retry. Nozzle face rate Mnzf 21.5 (t 0.0002) -> 20.84
(t 0.0074), still falling: +5 % over stream_mdot = 19.795 (the sonic face's HLLC flux exceeds rho c_T A; the
design rate counts only the upwind part). Min - Mnzf ~ 0 (the back sector sink takes 1.4e-4 in total, Mmsk).
dt 7.05e-6 -> 4.0e-6 by t 0.006 (location not measured: no bin output; static-first ordering not run).
0.84 s/cycle vs 0.51-0.56 for orbit1 during its approach (840 cells, donor_irr). Eirr grows (2.3e4 code by
0.0074); Ltop 6.0e5 -> 1.12e6 (the donor's light leaving through r_out).
**Volume nozzle at 13.79 (g8, 5fc471d5, from t = 0): 217 cycles stable**, Mnzf rate 19.4 over the first 3e-4 (fill
transient), dt 1e-6 to 3.8e-6 with eint<=0 retries (start0's fresh stream-on-at-t=0 also dropped to 2.4e-6).
Open: whether to calibrate stream_mdot (about 18.8 for a 19.795 face rate at t ~0.007) or accept the sonic-face
excess; where dt is set after the stream switch-on.

### 13.7 Pinned stream rate: supersonic nozzle at r 14.5 (user 10-15; supersedes the L1 nozzle for production)
User: the stream rate must be PINNED at stream_mdot = 19.795. The sonic L1 face (13.5) delivers +5-8 % (HLLC
star-region flux), so the nozzle moves inward to where the stream is clearly supersonic and the face flux is the
upwind rho v_r A. **stream_path_1015.py** (docs/dev/accretor_rhd/, output stream_path_1015.out; run with the spack
python 3.10 + py-scipy 1.11.1 + py-numpy 1.25.2): the one-zone core of stream_physics_1009 (adiabatic gas + trapped
radiation, diffusion loss) along the TRUE ballistic path from L1 (launch c_T 45.45 km/s) with the true transverse
curvature; at 13.79 it gives T 65.7 kK, arc width 1.568, rho 0.231 code (stream_physics_1009's linear-ramp model:
64.8 kK, 1.536, 0.241). Mach of |v| with c_tot (adiabatic, gas + radiation, EOS table) and c_gas (gas only):

| r | t from L1 | phi deg | v_r | v_phi | T K | Prad/Pgas | c_gas | c_tot | M_gas | M_tot | arc width |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 16.20 | 0.0005 | 0.00 | -45.5 | 0.6 | 107042 | 1.23 | 48.7 | 67.0 | 0.93 | 0.68 | 1.604 |
| 16.00 | 0.0049 | 0.05 | -46.7 | 5.9 | 104550 | 1.14 | 48.2 | 65.0 | 0.98 | 0.73 | 1.607 |
| 15.75 | 0.0100 | 0.22 | -50.9 | 12.5 | 98752 | 1.03 | 46.8 | 61.6 | 1.12 | 0.85 | 1.674 |
| 15.50 | 0.0147 | 0.49 | -57.2 | 18.9 | 92072 | 0.91 | 45.2 | 57.9 | 1.33 | 1.04 | 1.749 |
| 15.25 | 0.0188 | 0.82 | -65.2 | 25.1 | 85927 | 0.81 | 43.7 | 54.5 | 1.60 | 1.28 | 1.791 |
| 15.00 | 0.0223 | 1.21 | -74.4 | 31.2 | 80747 | 0.72 | 42.3 | 51.7 | 1.91 | 1.56 | 1.793 |
| 14.75 | 0.0255 | 1.62 | -84.3 | 37.3 | 76482 | 0.65 | 41.2 | 49.3 | 2.24 | 1.87 | 1.765 |
| **14.50** | **0.0283** | **2.06** | **-94.8** | **43.3** | **72966** | 0.58 | 40.2 | 47.3 | 2.59 | **2.20** | **1.720** |
| 14.25 | 0.0308 | 2.52 | -105.7 | 49.4 | 70024 | 0.53 | 39.4 | 45.7 | 2.96 | 2.55 | 1.667 |
| 14.00 | 0.0330 | 3.00 | -116.9 | 55.5 | 67534 | 0.49 | 38.7 | 44.3 | 3.35 | 2.92 | 1.612 |
| 13.79 | 0.0348 | 3.41 | -126.5 | 60.7 | 65743 | 0.46 | 38.1 | 43.3 | 3.68 | 3.24 | 1.568 |

Chosen: **inner face at stream_r = 14.5** (closest to L1 with M_tot >= 2; radial |v_r|/c_tot 2.00, |v_r|/c_gas 2.36,
so the code's gas-only HLLC flux is upwind). Deck: stream_inject = volume, stream_r 14.5, stream_mdot 19.795 (the
peak density 0.2726 code follows from the discrete face sum; rho_stream unused), env_t_stream 72966 K (diffusion
T profile, floor 19239 K, as before), stream_width 1.720 (arc). The pgen's ballistic state at 14.5 (phi 2.06 deg,
v_r -94.8, v_phi 43.3) is used. Nozzle cells r 14.5106-14.5878.
Stabilisers (keys, defaults): problem/stream_back_sink = true for any nozzle (the sector behind it, |dphi| < 4 sigma,
is sink): NEEDED also at 14.5: without it the gas behind the nozzle (v_r + 3 c_gas > 0 there) expanded into the
vacuum at 4e3-7e3 km/s and set dt (dt_report, g11: dt 2.9e-7, DT_RETRY eint<=0). problem/floor_keep_momentum =
true (radeq floor: v = m/rho_amb, |v| <= env_v_floor_max): kept on; not separately re-tested at 14.5 (timebox).
problem/dt_report = true prints the hydro dt-setting cell at every hst output.
g11 (14.5, back sink OFF): dMnzf 19.10 over the first 1e-4 (startup), then **19.796** (t 1e-4 to 2e-4).
**g12 (7fb237a8, the handover deck + dt_report, back sink ON, stream on from t = 0, 1 H200):** to t 0.0021 (interim,
job still running at the push): dMnzf **19.802 -> 19.798** per 1e-4 interval (pinned to 0.04 % from the first
interval), Min - Mnzf constant 0.0162 (the nozzle's initial fill), Mmsk 4.6e-6; no DT_RETRY, no NON-CONVERGED;
**dt 7.047e-6 flat**, set by the envelope (cell r 8.787 in the fine plateau, at rest), not by the nozzle or the mask.
(The L1 nozzle's dt 4e-6 is therefore L1-specific.) First impact: flight from 14.5 to the star ~0.045 - 0.0283 =
~0.017 after the stream reaches the nozzle state... i.e. stream on at 0.00946 -> impact near t ~0.027-0.030.

## 14. 3-D version: theta range +-20 deg (2026-10-16, Caltech, accretor-bigdom-1015)

User request 10-16: a 3-D version of the bigdom accretor (nozzle, donor sink/mask, donor irradiation, band-to-sphere
bookkeeping) and a half-resolution benchmark for Raven. Code 7ef6d8af.

### 14.1 Scope: what is 3-D now, and what is NOT
- **3-D mode** switches on automatically when the theta range reaches beyond 2 deg from the orbital plane (the 2-D
  band is +-0.35 deg; no key). It changes only the nozzle, the donor mask, the back sink and the Mnzf register;
  every 2-D path is the old code (gate 14.5: bitwise).
- **NOT 3-D yet (required before any 3-D production):** the forces, IC and WB are still the band's: gravity is the
  discrete gradient of Phi(r, phi) IN THE PLANE (arrays (m,k,i), no theta component, no theta dependence), Coriolis
  is the in-plane pair (a_r += 2 Omega v_phi, a_phi -= 2 Omega v_r; no sin/cos theta factors, no theta component),
  and the envelope IC, the WB background, the ambient (rhoamb_, tcol_, cref2_) are theta-invariant. The star is
  therefore a theta-extruded copy of the equatorial band: no vertical gravity confines the stream or the splash,
  and the tidal/centrifugal terms at latitude are the equatorial ones (centrifugal off by cos^2 20 = 0.88 at the
  walls). The +-20 deg deck is a COST AND PLUMBING benchmark (cells, dt, s/cycle, nozzle/mask/irradiation), not
  physics. The physics step: Phi(r, theta, phi) on (m,k,j,i) for g_r, g_theta, g_phi, the full rotating-frame
  Coriolis, and the IC/WB/ambient on the 3-D equipotentials (EnvICGeneral, PhiWB, tcol_ with j).

### 14.2 Nozzle Gaussian in z (problem/stream_sigma_z, stream_mdot = the full stream)
- State in the nozzle cells (r in the stream_ncell shell above stream_r, |dphi| < 3 sigma_phi, |z| < 3 sigma_z,
  z = r_f cos theta with r_f the nozzle's inner face, so all nozzle cells of a theta row share one state, as all
  cells of a phi column do in 2-D): rho = rho_s exp(-(x^2 + (z/sigma_z)^2)/2), x = dphi/sigma_phi; v = the ballistic (v_r, v_phi),
  v_theta = 0 (vertical hydrostatic stream); T = the diffusion profile in the elliptical radius sqrt(x^2 + (z/sz)^2);
  stream_rad E = a T^4, F = (4/3) E v as in 2-D.
- **sigma_z = 1.0054 Rsun**: the same one-zone core as the 2-D deck (stream_path_1015.py at r 14.5: c_iso =
  sqrt(P_tot/rho) = 39.19 km/s, Phi_zz = 8.622 Omega^2 there, sigma_z = c_iso/(Omega sqrt(Phi_zz)); the identity
  Mdot = rho_pk |v| 2 pi sigma_perp sigma_z with the script's rho 1.265e-8 g/cc, |v| 104.19, sigma_perp 1.564 gives
  the same 1.0054). For reference sigma_z at L1 is 1.21 (Lubow-Shu/Ryu), at 13.79 0.867 (S5 window model).
  3 sigma_z = 3.02 Rsun = +-12.0 deg at r 14.5: inside the +-20 deg walls.
- **stream_mdot (3-D) = 254.843 code = 1.0e-4 Msun/yr** (6.3028e21 g/s x 6.957e5 s / 1.7206e25 g). The peak density
  follows from the exact discrete sum over the 3-D window cells of the nozzle's inner face (r_f^2 (cos th_j - cos
  th_j+1) dphi), so the supersonic nozzle delivers exactly stream_mdot (as in 2-D). The pgen prints the central
  slice: the rate per unit height at z = 0 and its 2-D band equivalent (height r_f x 0.0122718).
- **Consistency with the 2-D deck (19.795):** at the same peak density the 2-D band rate is the central slice,
  Mdot_band = Mdot_3D r_f dtheta_band / (sqrt(2 pi) sigma_z). With sigma_z 1.0054 and 1e-4 Msun/yr this is
  **18.01**, 9 % below 19.795; the 3-D peak density is then 0.2480 code = the one-zone core at 14.5 (0.2476). The
  2-D value 19.795 is the 13.79 window rate (S5 one-zone with sigma_z 0.867 there) carried unchanged to the 14.5
  nozzle, i.e. the 2-D nozzle at 14.5 runs 10 % denser than the one-zone core there (0.2726 vs 0.2476). Two exact
  alternatives if the 3-D run must reproduce 19.795 in its central slice: stream_mdot = **280.1** (1.099e-4
  Msun/yr) with sigma_z 1.0054, or sigma_z = **0.9146** with 254.843. The deck uses the physical pair (1e-4, the
  one-zone sigma_z); USER DECISION.

### 14.3 Theta boundaries: reflecting at +-20 deg
Chosen: reflect (hydro and M1, as the band). Why: (i) the mass/energy budgets (Min, Mout, Macc, Mmsk, Ebnd) count
only the x1 faces, so reflecting theta walls keep them closed without new registers; (ii) with the theta-invariant
forces (14.1) there is no vertical gravity, so an outflow (zero-gradient) wall would let the unconfined
atmosphere/splash drain or inflow with no restoring force; (iii) the stream stays > 3 sigma_z from the walls
(nozzle +-12 deg; sigma_z shrinks toward the star as Phi_zz grows), so only the splash and the envelope circulation
at latitude see the walls. Revisit with the 3-D gravity: then an outflow wall above the stream (or a wider range)
is the alternative.

### 14.4 Mask, irradiation, back sink, bookkeeping in 3-D
- Donor mask = the exact 3-D Roche lobe: RochePot3(r, theta, phi) < Phi(L1) and x = r sin(theta) cos(phi) > d_L1
  (centrifugal about the z axis through the CM: -Omega^2 [(x - x_cm)^2 + y^2]/2). At latitude the lobe shrinks as
  the true equipotential does.
- donor_irr: unchanged (E = a T_d^4, F = 0 in the 3-D lobe; T_d from donor_lum over 4 pi R_L^2).
- Back sink (stream_back_sink): the sector behind the nozzle also limited to |z| < 4 sigma_z.
- Mnzf: the nozzle's inner face inside |dphi| < 3 sigma_phi AND |z| < 3 sigma_z (the face outside the window is
  ambient).
- Band-to-sphere: Lmeas/Ltop scale the domain sum by 4 pi / (domain solid angle) (= 4 pi / (2 sin 20 deg x 2 pi) for
  +-20 deg); unchanged code. env_lin_band_factor stays 0.963286: the star in this code is the equatorial band
  extruded in theta (14.1), so the bottom flux of every theta column is the band's.

### 14.5 2-D regression gate (Caltech, 1 H200, debug QoS; accretor_1015/g3d/gate)
7fb237a8 (the bigdom code head) vs 049c1f39, the handover deck as run (b145: nozzle 14.5, mask, donor_irr, stream on
from t = 0, dt_report), 20 cycles: **BITWISE** - user and hydro hst, all hydro_w and m1 bins, both rst files
identical as whole files; cycle/time/dt and dt_report lines identical (only the wall-clock field differs). Same
result for the first 3-D commit 7ef6d8af.

### 14.6 3-D half-resolution benchmark (deck docs/handover/accretor-bigdom-1015/b3d_half.athinput)
- Grid: nx1 420 (the bigdom 840 halved: old 640-cell map at half resolution, f_stretch_r_tail_nin 320 / ij 260, every
  other face of the 2-D grid; nozzle cells r 14.472-14.627), nx2 114 (theta 69.96..110.04 deg, dtheta = dphi =
  2 pi/1024), nx3 1024: **49,029,120 cells**; 4 MeshBlocks 420 x 114 x 256 on 4 H200 (one block in x1 for the nozzle
  and vet_col). (The 320 x 114 x 1024 = 37,355,520 of the old domain at half resolution.)
- Keys as the 2-D handover deck plus stream = true from t = 0, stream_mdot 254.843, stream_sigma_z 1.0054;
  env_lin_band_factor 0.963286 kept (14.4).
- **Sanity run (049c1f39, 4 H200, debug QoS, 200 cycles to t 0.00282; accretor_1015/g3d/b3d):** no FATAL, 0
  NON-CONVERGED, 0 dt retries, Picard mean 2.8 (max 5). dt **1.408e-5** flat (2 x the full-resolution 2-D 7.05e-6),
  set by an envelope cell at rest (r 8.786, the fine plateau), not by the nozzle or the mask. Nozzle face rate
  Mnzf **254.93 -> 254.89** per 2e-4 interval after the fill (+0.02-0.04 % of 254.843); peak rho 0.249341, central
  slice 18.011 in 2-D band units (14.2). z profile at the nozzle (bin, phi cell 0.034 rad): Gaussian in z_f with
  sigma_z 1.0054 to 0.2-2 % (first commit, cell-centre z; the r_f fix of 049c1f39 makes it exact), cut at
  |z| = 3 sigma_z, floor outside. Mask Mmsk 7.4e-5 by t 0.0028, smooth; Eirr grows linearly, Ltop 1.44e5 -> 1.90e5
  Lsun-equivalent (donor light leaving r_out; Lmeas 4.45e4 -> 4.54e4 as in 2-D). GPU memory 24.6 GB per H200.
- Cost: **0.80 s/cycle** steady (cycles 130-200), 1.0-1.9 s/cycle during the first ~120 cycles (switch-on
  transient). P_orb 0.473169 code / dt 1.408e-5 = 33,600 cycles/orbit -> **~30 H200-h per orbit** at 0.80 s
  (~45 at 1.2 s/cycle), pre-impact. Per GPU 15.3 M cell-cycles/s, consistent with the 2-D run (6.9 M cells, 0.51-0.56
  s/cycle, 13 M/s). Full resolution (8x cells, 2x cycles): ~480 H200-h per orbit on the big domain, memory ~0.8 TB
  (>= 8-16 H200). These are pre-impact, quiet-envelope numbers: the 2-D orbit1 impact needed dt retries and more
  Picard passes, and the true 3-D forces (14.1) are not in yet; the user's planning numbers (~300 half / ~2,500 full
  per orbit on the old domain, x1.3 for the big one) are 5-10x above this measurement.

## 15. v4: donor light scaled to the 3-D dilution, gravity-darkened donor, stream-lit back sink, clean sinks (2026-10-16, Caltech, accretor-bigdom-1015)

Raven (NOTE-2026-10-16-raven-accretor-bigdom.md, ADDENDUM 1) showed that in the 2-D band the donor light (13.3,
lobe held at a T_d^4) over-irradiates the gainer: E_rad near the star ~5x at all phi, |F| ~ cE mostly azimuthal,
~30x the 3-D dilution L_d/(4 pi D^2) ~1.5e5 code at the facing side; this drives the thin atmosphere (r 9.1-9.6)
to 200-940 km/s and dt down 2-3x (donor_irr = false: dt 7.03e-6 flat, 47 km/s). User decisions 10-16:
scale the donor light down to the physical strength (keep it on); gravity-darken the donor; the back-sink sector
holds the stream's own light; zero the sink momentum after the radiation step. **Status: implemented (5b0d4bc1),
builds (GPU nofma, md5 064f2963), NOT RUN.** The project was handed to Raven on 10-16 (user) before the Caltech
static check and the bitwise gate ran; Raven owns the tests (TASK-2026-10-16-raven-accretor-bigdom-HANDOVER.md).

### 15.1 Gravity darkening (problem/donor_gdark = beta, default 0.25; donor_gdark_tfloor, default 0.3)
- The lobe surface is the L1 equipotential of the rotating-frame Roche potential the code uses (RochePotXYZ =
  RochePot3 in Cartesian). g_eff = |grad Phi| (analytic, RocheGeffXYZ). Von Zeipel for a radiative envelope:
  F ~ g_eff, T_loc^4 = C g_eff^(4 beta), floored at (tfloor T_mean)^4.
- Each mask cell holds E = s a T_loc^4 (F = 0), with g_eff taken at the SURFACE point where the ray from the donor
  centre through the cell leaves the lobe (DonorSurfT: march in 0.005 a steps, capped at the L1 distance, then
  40 bisections), not the cell's own |grad Phi| (which goes to 0 at the donor centre and is not a surface value).
  A ray through the neck (no crossing) gives g = 0 -> the floor.
- Normalisation: int sigma T_loc^4 dA = L_d (problem/donor_lum) over the FULL 3-D lobe surface (rays on a 720 x 48
  quarter grid about the L1 axis; the lobe is symmetric in y and z), also in the 2-D band (the band holds the
  equatorial cut of the 3-D donor). T_mean^4 = L_d/(sigma A_lobe). beta = 0 keeps the old uniform T_d
  (4 pi R_L^2 Eggleton normalisation) exactly.
- Values (a 33.2596, q 1.125, L_d 69183 Lsun): A_lobe 2133 Rsun^2 (4 pi R_L^2 2105), in-plane perimeter 86.27 Rsun,
  T_mean 25934 K; in-plane T_loc (python twin of the host integral, 360 x 24 rays) by angle alpha from the L1
  axis: 7.8 kK (0.25 deg, the 0.3 T_mean floor), 16.9 kK (5 deg), 21.4 kK (15), 24.2 kK (30), 26.2 kK (60, 90),
  25.5 kK (120), 24.2 kK (180, the back; outside the domain, r_out 24.34 < a). The code prints its own table
  (alpha 0..180 and the pole) at start-up; check it against these.

### 15.2 donor_irr_scale (multiplies the held lobe E, uniform or gravity-darkened)
- problem/donor_irr_scale given: used as is. Absent: **1 in 3-D**; **auto in the 2-D band: s = s_geo x s_m1**.
- s_geo = (L_3D/(4 pi D^2)) / (L_2D/(2 pi D)) = int T_loc^4 dA / (2 D oint T_loc^4 dl): the band with reflecting
  theta walls is a donor infinitely extended in z, whose light dilutes as 1/D (cylinder) instead of 1/D^2 (sphere);
  for a sphere s_geo = R/D. D = a - problem/donor_irr_rcmp (default 9.6, the facing photosphere) = 23.66.
  s_geo = 0.523 (uniform), 0.636 (beta 0.25). Target 3-D flux at D: 1.54e5 (uniform) / 1.52e5 code (printed).
- s_m1 = problem/donor_irr_m1fac (default 1): the band M1/vet_col transport's excess over ideal cylindrical
  dilution. **Not calibrated.** Raven measured the facing flux ~30x the target with s = 1 (uniform); ideal
  transport would give only 1/s_geo ~2x, so the remaining ~15x is M1 numerics (held isotropic cell emits ~cE/2,
  beams do not spread, vet_col radial columns cannot represent the lateral light; fix 3). If that holds,
  s_m1 ~ 1/15 and s ~ 0.035-0.04 (Raven's stopgap: ~1/30, T_d,eff ~11 kK). Calibrate from a static arm with s = 1:
  s_m1 = target / (measured facing donor flux at r 9.6, minus the donor_irr = false arm) / s_geo.

### 15.3 Back-sink sector (nozzle -> L1): the stream's own light (user 10-16)
- Was: held a T_d^4 like the lobe (doubled Eirr in the stream phase). Now: within the nozzle width (|dphi| < 3 sigma,
  3-D also |z_f| < 3 sigma_z) it holds E = a T^4 with the NOZZLE's T profile (core env_t_stream, diffusion profile
  StreamT, floor stream_t_surf), F = 0 (the sink gas is at rest); outside it (3-4 sigma) its radiation is left to the
  transport (dark, like the ambient). Only with stream_rad = true. The T(r) of the one-zone model
  (stream_path_1015.py) is not in the code, so the nozzle profile is held unchanged back to L1. Not scaled (it is the
  stream's light, not diluted donor light). Booked with the nozzle radiation (register 3 -> Ebnd).

### 15.4 Sink momentum (problem/sink_zero_v, default true)
- New pgen hook user_m1post_func (pgen.hpp), called by RadiationM1::HydroConToPrim after the M1 coupling and before
  the fluid refresh. RyPerM1Post zeroes the momentum of every sink cell (lobe and back sector, not the nozzle),
  removes the kinetic energy from E_tot and books it in Emsk. Active only with donor_irr on or the back sink live
  (stream on). The hydro dt is computed before the M1 stage, so dt is unchanged. Removes Raven's r 14.6 row artefact
  (v_r up to 1e6 km/s at floor density) from outputs and Emsk (untested).

### 15.5 Gate status
- Not run. Expected bitwise vs 08625ed7 (binary 049c1f39) on the 2-D v3 deck with donor_irr = false and the new keys
  absent: every new path is behind donor_irr or the live back sink, and sink_zero_v is read without being added to
  the parameter dump (5b0d4bc1). The hook is called (empty) in every M1 run of this pgen.


## 16. Back-sink sector light switch: problem/stream_sector_light (2026-10-17, Raven, accretor-sector-light-1017)

Why: the v4 stream test on Raven (NOTE-2026-10-17-raven-accretor-v4.md) collapsed dt to 2-7e-7 through the impact.
The suspected driver is the light held in the back-sink sector (section 15). It is the stream's own light, UNSCALED:
core 72966 K, sector-mean T^4 = (56.2 kK)^4. That is ~132x the scaled donor surface value 0.22 (24.2 kK)^4.
Before this change the sector light could only be switched together with stream_rad.

Key `problem/stream_sector_light`. It is read only when named, so the PAR_DUMP of older decks is unchanged:
- `full` (default): the v4 rule, bitwise unchanged.
- `scaled`: E x donor_irr_scale (needs donor_irr). The band's 2-D -> 3-D factor of the held donor light
  (s_geo x s_m1) is applied to the stream's held light too: the same transport, and the same cylindrical-vs-
  spherical dilution. For the stream this is an UPPER bound. Its vertical extent is smaller than the donor's:
  int T^4 dz / (2 D T_c^4) = 2.12 sigma_z / (2 x 5.1) ~ 0.21 for sigma_z 1.0054, against s_geo 0.636. With
  s_m1 that gives ~0.07 instead of 0.22.
- `surface`: the sector follows the donor-lobe rule (donor_irr light: scaled, gravity-darkened, booked in Eirr).
  This is the v3 rule (v3: lobe light at T_d on the whole sink, sector included). It needs donor_irr.
- `dark`: no light is held in the sector, which is left to the transport like the ambient. The stream's light then
  enters only through the 2-cell nozzle.

**Gate (Raven A100, gpudev 31249502): BITWISE.** The binary (62119e93, md5 0ba9d066) was compared with v4 5b0d4bc1,
with the key absent, over 20 cycles:
- S, the v4 stream switch-on (cal deck from the v4 static rst, stream on; exercises the sector): hst, hydro_w and m1
  bins, and the cycle/dt lines are identical. With `stream_sector_light = full` named, hst and cycle/dt are
  identical; the bins differ only by PAR_DUMP.
- v4 gates A (base deck, donor_irr off) and G2 (s5 + stream + ROL2U): identical.

**Confirming tests.** All runs are on Raven A100, from v4/cal/static/rst/ryper.00001.rst (s = 0.22, t 0.00946); see
/raven/ptmp/jinma/accbig_1016/NOTE.md, RESULTS 10-17b. dt is min/median over each window.

| arm | dt 0.0095-0.015 | dt 0.02-0.0265 | contact 0.027 | dt after contact | thin-atm max\|v\| t 0.013 / 0.020 / ~0.031 | Ltop t 0.0107 / 0.0285 | E r 9.4-10 median t 0.0142 / 0.0284 |
|---|---|---|---|---|---|---|---|
| v3 | 5.6e-6 / 6.7e-6 | 1.5e-6 / 2.9e-6 | smooth | 1.2-1.6e-6 | 376 / 717 / 1117 | 3.7e5 / 7.8e5 | - / 61 |
| v4 full | 2.9e-6 / 4.4e-6 | 3.1e-7 / 7.6e-7 | 49 retries over the run | 4.5-5.5e-7 | 660 / 1429 / 2394 | 1.2e6 / 4.1e6 | 151 / 308 |
| scaled | 6.7e-6 / 7.0e-6 | 2.9e-6 / 5.0e-6 | 1 retry (t 0.0251) | 1.5-1.6e-6 | 193 / 468 / 1068 | 3.0e5 / 7.5e5 | 23 / 58 |
| surface | 7.03e-6 flat | 7.02e-6 flat | eint<=0 burst: FATAL at dt_min 1e-7; 61 retries, recovered at dt_min 1e-8 | 1.56-1.65e-6 | 80 / 139 / 469 | 6.9e4 / 1.0e5 | 5.0 / 7.4 |
| dark | 7.03e-6 flat | 7.02e-6 flat | same as surface: FATAL at 1e-7; 34 retries, recovered at 1e-8 | 1.57-1.63e-6 | 77 / 128 / 474 | 6.3e4 / 8.9e4 | 4.5 / 6.5 |

- The impact matches v3 in every arm. Dense infall reaches r 9.4-9.65 at t 0.0272-0.0284, and dense-gas heating
  starts at t ~0.031 at phi 13.6-13.8 deg (v3: 14.0).
- Floor-gas e/rho stays at 2.0-2.3e3 in every arm (v4 full: up to 1.8e6).
- **Conclusion: the unscaled sector light caused the v4 dt collapse.**

## 17. Ray-traced external donor irradiation: problem/donor_irr_mode = ray (2026-10-17, Raven, branch accretor-raydonor-1017)

USER 10-17 (night): fix option (A) for the vet_col lateral-beam artefact (accbig_1016 NOTE "Fix 3"; vet_col's per-column
tensor diag(f_K, (1-f_K)/2, (1-f_K)/2) and the diagonal-only M1SphCurv turn the donor's lateral light into an azimuthal
beam wrapping the star: night/facing dE ~1). The donor's light is taken OUT of the M1 field. It becomes an external
beam, traced along straight rays with exp(-tau) and the star's exact shadow, and deposited where it is removed. Only
what the gas re-emits (and the scattered part) is M1 radiation. Code: src/pgen/ry_per_accretor_raydonor.hpp (included
in ry_per_accretor.cpp); the hooks in ry_per_accretor.cpp touch only the donor-irradiation path.

### 17.1 Keys (read only when donor_irr_mode is given; default m1 = sect. 15 bitwise)
- problem/donor_irr_mode = m1 | ray (needs donor_irr = true, donor_mask = roche, <rad_m1>; ray: the 2-D band only).
- donor_ray_scale (1): factor on the surface intensity (test knob). donor_irr_scale is NOT applied in ray mode.
- donor_ray_npsi (24), donor_ray_nel (4, <= 8): azimuths and Gauss-Legendre elevations per column.
- donor_ray_every (1): the beam is recomputed every N cycles (first hydro stage).
- donor_ray_step (0.5): march step / min(dr, r dphi). donor_ray_grow (4). donor_ray_taumax (40). donor_ray_force (true).
- donor_ray_deposit = rad (default) | gas; donor_ray_opac = total (default) | rosseland.
- donor_ray_dump_dt (0 = off): rank 0 writes ryper.raydonor.NNNNN.bin (q, f_r, f_phi, F_r, F_phi, chi on the global
  r-phi grid; reader tool accray_1017/tool/rdump.py).

### 17.2 Method
- Source: the exact 3-D Roche L1 lobe, each surface point a grey isotropic emitter I = (c/4 pi) a T_loc^4 with the
  donor_gdark law of sect. 15 (same C, same floor, L_d over the full lobe). The 3-D solid angle is integrated
  directly, so there is no band dilution factor. A point or extended source in the M1 sense is NOT justified here:
  from the facing photosphere (D 23.7 to the donor centre) the lobe spans up to ~46 deg, and the gainer (R ~9) is not
  small against d_L1 = 16.2.
- Directions, per (r, phi) column P: npsi azimuth midpoints over the lobe's in-plane angular span seen from P (from a
  720-vertex polygon of the in-plane cut), x nel GL nodes in elevation on [0, e_max(psi)], doubled for -e. e_max is
  where the 3-D ray grazes the lobe, from the lobe height z_s(x, y) sampled at 32 points along the in-plane chord; the
  hit point of each ray (first s with z_s >= s tan e) gives g and hence I. This geometry is static and tabulated at
  setup (float W = I cos e w dpsi, cos e; 354 MB per A100 for the 840 x 2048 band on 4 ranks).
- Band = z-uniform (reflecting theta walls), so an elevated ray crosses the gas of its in-plane projection over a path
  longer by 1/cos e: tau_3D = tau_proj / cos e. At the photosphere this is the plane-parallel slab of the sphere (the
  band cylinder's tangent plane is the sphere's at the equator), so horizon and limb are those of the 3-D star.
  Beyond the in-plane lobe entry the ray is in the masked lobe (floor, tau ~ 0); beyond r_out it is vacuum; below r_in
  it is opaque.
- Attenuation (long characteristics): one march per (column, psi), starting at the cell's donor-side face, through a
  GLOBAL extinction map chi(k, i), the theta mean of <rad_m1>/opac from the last M1 solve. The first update of a run
  calls Opacity() on the current state. The map is gathered over ranks by MPI_Allreduce every update; every entry has
  one owner, so the gather is exact and deterministic.
  - donor_ray_opac = total (default): chi = max(rho kappa_T, rho kappa_E).
  - donor_ray_opac = rosseland: rho kappa_T.
  - Step: donor_ray_step x min(dr, r dphi). It doubles per step up to donor_ray_grow (4) cells while chi h < 1e-3.
  - The march stops at tau > 40.
- Deposition (photon conserving): per ray I dOmega exp(-tau_up/cos e)(1 - exp(-chi L/cos e))/(L/cos e), with L the
  in-plane chord of the cell.
  - donor_ray_deposit = rad (DEFAULT): the whole removed power goes into the M1 E at the last hydro stage, with F
    unchanged. The implicit gas-radiation coupling gives the absorbed part to the gas at the iterate's kappa(T).
  - donor_ray_deposit = gas: eps = kappa_E/chi of it heats u(IEN) explicitly, every stage, RK weighted. This is
    UNSTABLE: in the 10-17 smoke, the rho ~5e-6 thin atmosphere (kappa_P ~2000 cm^2/g) went from e/rho 2e3 to 8e5 in
    one step and cannot cool through its now tiny kappa_P(T); dt fell from 7e-6 to 6e-7.
  - With rad, the beam's local energy density (|F_inc|/c <~ 0.25 code) is < 1 % of the gainer's own E at the
    photosphere (~3-10), so the thin-gas T barely depends on where the deposit goes.
- Momentum: q/c x rho kappa_T/chi along the beam, i.e. the flux-mean (transport) opacity, as the M1 module couples the
  gainer's own flux. The work v.f goes to u(IEN). With the Planck mean (smoke) the facing thin atmosphere got ~28x
  gravity inward (v_r -290 km/s in 60 cycles).
- Booking: the energy given (deposit + work) is radiation IN, in Eirr (hst 26) and Ebnd.
- The lobe and the stream_sector_light = surface sector hold no light (as donor_irr = false): the lobe becomes a
  transparent sink for the gainer's light (no reflection effect).
- Update every donor_ray_every (1) cycles at the first hydro stage. The beam is held per volume between updates.
- 3-D wedge: FATAL. The z-uniform projection no longer holds, so each (cell, direction) needs its own 3-D march
  (12.5 M cells x ~200 directions, ~100x the band cost per update). Options: a coarse ray grid with interpolation,
  update every ~100 cycles, or short characteristics per direction.
- Diagnostics: setup self-check (unattenuated F_inc on the line of centres, table vs 16x finer quadrature); ray dumps
  (donor_ray_dump_dt).

### 17.3 Reference: the gravity-darkened donor's facing flux (accray_1017/tool/ref_flux.py)
This is an independent surface integral over the 3-D lobe, sum I cos dA u/d^2, with the same gdark law and the same
normalisation. Unattenuated F_r on the line of centres:

| r | beta 0.25 | beta 0 | sphere L/(4 pi D^2) |
|---|---|---|---|
| 9.6 | 7.09e4 | 1.46e5 | 1.52e5 |
| 10.0 | 7.18e4 | 1.51e5 | 1.58e5 |
| 10.4 | 7.25e4 | 1.57e5 | 1.63e5 |

**The gravity-darkened L1 face gives 0.47x the uniform-sphere value.** The old calibration target 1.52e5 was the
uniform sphere.
