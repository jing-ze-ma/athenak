# NOTE orion -> raven (2026-10-19): the RSG g25 setup of the Arepo2 + IDORT runs

Reply to `TASK-2026-10-19-orion-rsg-g25-setup.md` (80c93a5d2d). Small files are in `docs/handover/rsg-g25-orion/`
(listed at the end). Large files stay on Orion; their paths are given below. Units are cgs. "R" means R_ph late in the
run (5.83e13 cm) unless a different R is named. L_c = 6.854917e38 erg/s is the core heating rate (= MESA L).

Code and run:
- Code: Arepo2-IDORT-dev, IDORT as a module (`src/modules/idort`), repo `/orion/u/jinma/AREPO/Arepo2-IDORT-dev`.
  Pushed head 0ab6a8294a. Setup files are in `examples/idort/G25_setup/`; method docs in
  `documentation/modules_idort.md` and `examples/idort/README.md`.
- Production run: `long_v1`, at `/orion/ptmp/jinma/idort_dev/runs_g25/runs/long_v1`, build 36f0bdb1e0, run 0 -> 3.725e6 s.
- Reference: the old code Arepo-IDORT (IDORT-dev) production "cons12". It runs on Fritz (job 4325157); an Orion rerun is
  at `/orion/ptmp/jinma/idort_dev/cons12_orion`. Same IC and almost the same physics; differences are in point 8.
- Analysis: `/orion/ptmp/jinma/idort_dev/STATUS_phys_cmp.md` ("long_v1 analysis"), `/orion/ptmp/jinma/idort_dev/phys_cmp/v1/`.

Items marked **[uncertain]** could not be checked from files on Orion.

---

## 1. Stellar model

| item | value |
|---|---|
| source | MESA r15140 (gfortran), `golden.zip` run 25, `LOGS1/profile438.data` (= `m22_p438.data`, md5 29277949b03aa45ecd346436117fa7aa) |
| MESA run | binary-mode MESA: m1 = 22.0, m2 = 17.6 Msun, P_orb = 1e6 d, so effectively a single star. `inlist_cluster` only sets new_Z = 0.02, new_Y = 0.28, kap Zbase = 0.02. The base inlists are not in golden.zip **[uncertain: the rest of the inlist]** |
| phase | core He burning RSG, age 8.386e6 yr, X_c(He4) = 0.441, M_He-core = 9.40 Msun |
| why this profile | among the core-He-burning RSG profiles, it has the largest fractional radius of the convective-envelope base (memory note rsg-ic-profiles) |
| M (MESA) | 20.540 Msun |
| L (MESA) | 1.7907e5 Lsun = 6.855e38 erg/s |
| R_phot (MESA) | 843.80 Rsun = 5.869e13 cm |
| Teff (MESA) | 4088 K |
| core treatment | **point mass + replaced core, no inner boundary** (details below) |
| composition | **uniform**: X = 0.653848, Y = 0.326375, Z = 0.019757 |

Core treatment:
- Point mass M_pm = 2.11248e34 g (10.62 Msun), held fixed at the box centre.
  - Gravity on the gas is exact (N^2) with an Arepo spline softening, epsilon = 6.2878e11 cm.
    The kernel radius h = 2.8 epsilon = 1.7606e12 cm ≈ R_c, so the point mass is exactly Newtonian beyond R_c.
- Gas fills the core down to r = 0.
- Core radius R_c = 1.7611e12 cm = 25.3 Rsun = 0.030 R_MESA (the IC builder's `fraction_cut` 0.03).
- Inside R_c the MESA structure is replaced by the "c1s" core model (`coreB_side_g25_c1s.hdf5`; the attributes and profile
  are in `rsg-g25-orion/core_c1s_profile.txt`). Its description (HDF5 attribute `note`):
  - a model with a smooth target Brunt-Väisälä profile N^2(x) = A x^2/(x^2 + xk^2), with A = 5.5e-10, xk = 0.25,
    blended into the "c1" Gamma ramp over x = 0.5-1;
  - kappa_eff = Gamma 4 pi c G m_eff / L, with L = L_c q(r); Gamma_in = 0.15, Gamma(R_c) = 0.55;
  - M_pm is chosen so that M_gas(0) = 0 while m_tot(R_c) is unchanged;
  - gas mass inside R_c = 1.117e33 g (0.56 Msun);
  - the structure is identical to the "c1" model from R_c to R_b = 2 R_c = 3.52e12 cm.

  The script that built it (`core_optB/c1s/c1smodel.py`) is not on Orion **[uncertain: details beyond the attribute note]**.
- The core is radiative and convectively stable by construction:
  - kappa_eff replaces the opacity for r < R_c (point 4);
  - the heating L_c is deposited as gas energy with a cubic-spline kernel of radius r_h = 0.2 R_c = 3.522e11 cm
    (`SPHERICAL_SOURCE`, point 6).
- No other central source and no inner boundary.

Composition:
- The 6-species `NuclearComposition` in every cell: H 0.65385, He4 0.32638, C12 0.00215, N14 0.00380, O16 0.00773,
  Ne20 0.00610. These are the MESA envelope surface values, applied to the whole star including the core region.
- Nothing changes it (no sources).
- The EOS uses the cell's X (one table node at X = 0.6538476353670762). The opacity tables are for X = 0.70, Z = 0.02
  (single composition, see point 4).

Total mass in the box: gas 1.98170e34 g (9.963 Msun) + point mass = 4.0942e34 g (20.58 Msun).
This is 0.2 % above MESA's 20.540 Msun **[uncertain: the cause of the 0.2 %, not investigated]**.

## 2. Initial condition

Construction, in three steps:
1. **Base IC.** `rsg-g25-orion/create-ics-golden.py` (Orion copy of arepo-snap-util `create-ics.py` with a command line;
   venv `/orion/u/jinma/AREPO/venv-ics`). Run as
   `create-ics-golden.py --name RSGg_golden25_p438_hr --profile m22_p438.data --masslimitfactor 1e-4 --grid-ndir 128`.
   - It uses `stellar_ics` (Ohlmann et al. 2017 style): `GiantStarMesaController`, core cut at `fraction_cut` 0.03 R
     replaced by a softened point mass + n = 3 profile, OPAL EOS5 Z = 0.02.
   - Cells are placed on HEALPix shells (`mode='combined'`, `lower_mass_limit` 1e-4, `factor_gridtoradius` 6), with a
     background of rho 1e-17, T 200 K.
   - Three nested background boxes are added: 2x (64^3), 10x (32^3) and 50x (32^3) the star box of 3.503e14 cm, so
     BoxSize = 1.7515e16 cm.
   - Output: 6,092,886 gas cells + 1 point mass (2.1319e34 g). File:
     `/orion/ptmp/jinma/Arepo/ICs_golden/ic/RSGg_golden25_p438_hr/rhoin-RSGg_golden25_p438_hr.ic.hdf5`.
2. **Production IC.** `rhoin-RSGg_golden25_p438_hr_bgr2_coreB_opal_c1s.ic.hdf5` (md5 6137633e6101916bd0105aeb7a95fbb1, 755 MB),
   path `/orion/u/jinma/AREPO/RSGg_g25_IDORT/ics/`. It was made from the base IC elsewhere, not on Orion.
   I compared the two files directly:
   - same IDs and generator coordinates (bitwise);
   - density changed only for r < 3.52e12 (the c1s core; ratio 0.93-1.83);
   - u changed everywhere in the interior: ratio 0.32 in the core, 0.72 at 1e13-5e13 cm, 1.0 near the surface. The base u
     included radiation; the production u is the gas-only OPAL u ("opal"; the radiation is carried by IDORT);
   - CSM/background ("bgr2") changed from rho = 1e-17 to rho ≈ 1.1e-18 (median 5.87e13-3e14) falling as about r^-2,
     reaching the 1e-20 floor beyond ~1.5e15; u = 4.84e10 erg/g (T ≈ 630 K) in all of it;
   - point mass 2.1319e34 -> 2.11248e34 g.

   The names bgr2/coreB/opal/c1s are read from the data; the scripts are not on Orion **[uncertain]**.
3. **Module conversion.** `rsg-g25-orion/convert_ic.py` changes only two things:
   - BoundaryRole: the outermost 32^3-lattice layer becomes WALL_SOLID (5768 cells), the next layer WALL_FLUID (5048 cells);
   - Header BoxSize is set.

   Output: `/orion/ptmp/jinma/idort_dev/runs_g25/ics/g25_full.hdf5` (731 MB). `Masses` holds the density
   (`MESHRELAX_DENSITY_IN_INPUT`: m = rho V of the first mesh).

Relaxation, damping and seed noise:
- **No relaxation and no damping.** cons12 had RelaxTime 0 and RELAXOBJECT/DAMP_STAR_CORE not compiled; the module run
  has none either. The run starts directly from the mapped IC with v = 0 everywhere, so the first ~1e6 s are the
  relaxation transient (point 7).
- **No seed noise.** The perturbations come from the HEALPix-shell discretisation and the Voronoi mesh.
- Radiation start (`IdortInitMode 2`): a steady-state radiation solve at fixed gas, so snapshot 000 holds the converged
  field. Escaping L at t = 0 is 4.64e38 (0.677 L_c).
- f-hybrid start: the share f a T^4 is moved into the gas energy, |dT/T| ~ 1e-15.

IC zones (direct cell counts from the IC; rho = IC density):

| r [cm] | cells | rho range | note |
|---|---|---|---|
| < R_c = 1.761e12 | 54,361 | 1.5e-5 - 3.96e-4 | c1s core (rho_c 3.96e-4) |
| R_c - 1e13 | 145,548 | 4.3e-7 - 1.35e-5 | |
| 1e13 - 5e13 | 525,432 | 4.1e-9 - 3.9e-7 | convective envelope (MESA) |
| 5e13 - 5.87e13 | 3,020,664 | 1e-17 - 3.7e-9 | fine photosphere layer, dx down to 3e11 |
| 5.87e13 - 6.5e13 | 14,680 | 8e-18 - 1e-17 | |
| 6.5e13 - 3e14 | 2,091,392 | 3.8e-19 - 8e-18 | CSM, r^-2 |
| 3e14 - 9e14 | 177,936 | 4.3e-20 - 3.8e-19 | |
| > 9e14 | 62,873 | 1e-20 - 4e-20 | nested 32^3 lattices |

Profiles of the IC (after the steady radiation start) and of the late state, in 0.02-dex shells, are in
`profile_ic_t0.txt` and `profile_late_t3.667e6.txt`. Columns: r, n, rho, T_gas, T_rad, L_lab, L_comoving, v moments, kR,
kP, f, dx, mean cell mass, M_r. The MESA profile (selected columns) is in `mesa_g25_p438_selected.txt`.

## 3. Domain, resolution, gravity

Domain and mesh:
- Box: cube of 1.7515e16 cm, star at the centre. The faces are 8.76e15 cm = 150 R from the centre.
  The mesh is periodic, but the faces are separated by the vacuum layer (point 6), and gravity is not periodic.
- Moving Voronoi mesh (`VORONOI_CELL`), regularisation CellShapingSpeed 0.5, CellMaxAngleFactor 2.25.
- Cells: 6.09M (IC) -> 8.6M (1e5 s) -> 12.5M (1e6 s) -> 18.8M (2.0e6 s, peak) -> 17.30M (3.67e6 s).
  Late distribution: r < R_c 52k; R_c-1e13 0.23M; 1e13-5e13 2.1M; 5e13-5.87e13 4.7M; 5.87e13-6.5e13 7.8M;
  6.5e13-3e14 2.3M; beyond 3e14 0.08M (0.02-dex bin counts).

Refinement (all inherited from cons12; exact semantics in modules_idort.md "Star refinement bands"):

| criterion | rule |
|---|---|
| mass | split if m > 2 m_ref, merge if m < m_ref/2. m_ref = 1.3434e28 g (6.75e-6 Msun) (`ReferenceGasPartMass`, `TargetGasMassFactor 1`, `RefinementCriterion 1`, `DerefinementCriterion 1`) |
| volume | MaxVolume 1.988e42 (dx 1.26e14); MaxVolumeDiff 5 (neighbour volume ratio); MinVolume 0 |
| around the point mass | `REFINEMENT_AROUND_DM`: 10 cells per softening length (`RefinementCellsPerSoftening`); forced derefinement inside the softening disabled (as in the old code) |

Volume bands (`REFINEMENT_VOLUME_LIMIT`): inside each band the cell volume limit is lowered to V_max; split if
V > 2 V_max, merge only if V < V_max/2; the smallest applicable limit wins.

| band | condition | V_max | dx = V^(1/3) |
|---|---|---|---|
| CSM | 1.174e13 < r < 1.761e14 | 6.779e36 | 1.89e12 |
| photosphere radius | 5.870e12 < r < 6.457e13 | 1.763e35 | 5.61e11 |
| density | 1e-13 < rho < 2e-9 | 2.204e34 | 2.80e11 |
| temperature | 4000 < T < 20000 K and rho > 1e-13 | 2.204e34 | 2.80e11 |

Resulting resolution:
- Envelope cell mass: 1.4-4.8e28 g. Envelope dx: 5e11-1e12 cm at 1e13-3e13 cm, i.e. about 40 cells per radial decade
  at 1e12-1e13 and about 100 at 3e13.
- Core dx: 6.5e10 cm at r ≈ 7e11, 1.1e11 at R_c (IC); so R_c is resolved by roughly 20 cells.
- Photosphere dx: 1.7e12 cm at t = 0 (0.5 cells per density scale height H_rho = 9.2e11 at R_ph). As the band
  refinement acts it drops to 3.0e11 (3.4 cells per H_rho, t >= 1.6e6 s).
- tau of the spherically averaged rho chi falls from about 1e5 to 1e-2 within one 0.02-dex bin: the photospheric
  density cliff (1e-9 -> 1e-17 within about 1e12 cm) is only marginally resolved.
- 1.2-3 R: dx ≈ 2.2e12 cm. Beyond 3.3 R: dx 0.15-0.5 R.

Gravity:
- Gas self-gravity: Arepo tree, hierarchical, non-periodic. ErrTolTheta 0.5, relative opening criterion with
  ErrTolForceAcc 0.001, direct summation below 1000 particles, gravity at cell centres.
- Gas softening: adaptive, 2 x cell radius, minimum 8.26e9 cm.
- Point mass (type 1): exact N^2 gravity on the gas, epsilon 6.288e11 cm. It feels no force, gets the reaction kick,
  and its velocity is zeroed before each drift (`DMFIXED`), so it stays at the centre (offset 0, |v| <= 14 cm/s).
- `GRAVITY_PM_ENERGY_CONSERVING`: point-mass energy correction ported from the old code. Without it the core spun up
  to 17 km/s (run long_fhyb).

Hydro and time steps:
- Second-order finite-volume moving-mesh hydro, Riemann solver HLLD (MHD with B = 0, needed for local steps on
  VORONOI_CELL), CourantFac 0.3.
- Local (hierarchical) time steps (`TIMESTEPS_WAKEUP`, `TREE_BASED_TIMESTEPS`). Top step 7275.58 s, finest bin 227.4 s
  (top/32). 512 top steps and 16,538 radiation solves to 3.725e6 s.
- The signal speed includes radiation pressure support: c^2 = c_gas^2 + (4/9) E_rad/rho (1 - e^-tau), scaled by (1 - f)
  where f > 0. The Riemann solver sees P_gas + f P_rad.

## 4. Microphysics

EOS (`IdortEosMode 2`, the cons12 gas EOS):
- Source table: OPAL EOS2005 `EOS5_data` (Z = 0.02; 15 MB). Path `/orion/u/jinma/AREPO/Arepo-snap-util/eostable/EOS5_data`.
- Outside the OPAL region: an ideal gas in Saha/H2 equilibrium (H2, H, H+, e-, He ladder, C/N/O/Ne ionisation).
- Blend: C2 weight between the two; OPAL for rho > 2e-14 and T above about 2057 K at rho < 2.7e-4
  (old OPAL_EOS_BLEND_IDEAL + OPAL_IDEAL_STITCH_FIX).
- **Ionisation is included; radiation is NOT in the gas EOS** (EOS_NO_RAD). The radiation energy is the IDORT
  intensity field.
- Exception (f-hybrid, point 5): in optically thick cells a share f of a T^4 (energy) and a T^4/3 (pressure) is added
  to the EOS.
- Runtime: a fast table built at start-up, ln T 20 K - 1e8 K (step 0.005), ln rho 1e-26 - 1e4 (step 0.05), one X node.
  Accuracy against the reference EOS: e rms 8e-5, max 3e-3.

Opacities: grey, Rosseland and Planck means.
- Tables: MESA-format `gs98_z0.02_x0.7_ma_{ross,planck}.data`, copied into `rsg-g25-orion/`.
  - Header: "Ferguson, Alexander et al. 2005 with GS98 photospheric metals", X = 0.70, Z = 0.02.
  - Grid: logT 2.7-8.7 (180 points), logR -8 to 1 (37 points), logR = log rho - 3 log T + 18.
  - The Ferguson low-T molecules and dust and the OPAL high-T part are in the tables themselves (MESA's merged
    gs98 tables); we do no splicing at run time.
  - md5: ross 352f7b44fc25b53b518cafc665501e6e, planck 8c5453f63f02fe933261d92a0dd6cece.
- Interpolation: bilinear (`IdortOpacityInterp 1`, = cons12). Off-table points use constant continuation, clamped 1e-6
  of the mean grid spacing inside the edge; the CSM is at logR < -8.
- Scattering kS = 0. The cell X is not used for kappa.
- Old-physics modifiers (`IdortOpacityModifiers 1`; they reproduce the cons12 `set_opacity()` to 3e-14). Applied in order:
  1. floor: kappa >= 1e-5 cm^2/g on both raw lookups;
  2. T floor: for T < 3200 K, kR = max(kR(T), kR(3200 K)) and kP likewise, at the same rho. **Kappa never drops below
     its 3200 K value, so the molecular opacity minimum below 3200 K is cut off;**
  3. dust switch, for r > 1.174e14 cm (= 2 R*): kR -> fd kR + (1 - fd) kP, with fd = 0 below 1200 K and 1 above
     1500 K (linear between);
  4. photosphere switch, rho_s = 1e-10: kP -> fp kP + (1 - fp) kR, with fp = 0 for rho <= 1e-11, 1 for rho > 1e-9,
     linear in rho between. **In the thin atmosphere and CSM the absorption/emission opacity equals the Rosseland mean;**
  5. core: for r < R_c, kR = kappa_eff(r) from the c1s table (log-linear in r) and kP = min(kP, kappa_eff).

Floors:
- no gas T floor;
- no GasTmaxCap;
- no velocity caps;
- MinEgySpec = 1e10 erg/g (Arepo hydro u floor; total EgyInjection from all resets 4.2e39 erg = 1.7e-6 L_c t);
- MinimumDensityOnStartUp 1e-20;
- opacity floors as above.

The CSM gas is free to cool: 170-500 K, 280 K at 3e15 cm. The old code had a 531 K floor (point 8).

## 5. Radiation (IDORT)

| item | value |
|---|---|
| method | implicit discrete ordinates on the Voronoi mesh (Ma, Pakmor, Justham & de Mink 2025, arXiv:2503.16627), module rewrite |
| angles | 24 directions (`IDORT_NUM_ANGLES=24`, BRULS set as cons12 `IdortBrulsNum 2`) |
| frequency | grey: emission and absorption kP, extinction kR + kS |
| speed of light | physical c, no reduced c (`IdortClight -1`) |
| time discretisation | one backward-Euler radiation step per hydro step of each active cell (local steps), dt_i = the cell's hydro step |
| frame | mixed frame, O(v/c) terms with the cell velocity at the start of its hydro step (`IdortUseVelocity 3`) |
| face closure | `IdortTau0 0`, `IdortAlpha 5`. Upwind weight w = 1/(1 + tau_f), tau_f = 5 min(chi_L, chi_R) d (linear form). cons12 used the quadratic w = 1/(1 + tau_f + tau_f^2/30) with mean-kappa face tau |

Gas coupling:
- `IdortCoupleGas 1`, conservative. T_gas is an unknown of each cell block, solved by bracketed Newton.
- The gas receives exactly the lab-frame energy dE and momentum dp that the radiation loses. Opacities are evaluated at
  the start-of-step state.
- The emission term per cell is kP B(T) + (kR + kS - kP) J0; at the photosphere T_rad and T_gas are coupled only
  through this implicit block, with no extra relaxation.
- With the photosphere switch kP = kR for rho < 1e-11.
- Measured (shell averages): T_gas = T_rad to ~1e-3 in the interior, and 2470 vs 2418 K at the photosphere
  (5.89e13 cm, late). In the f-ramp shell at 5.62e13 cm the averages differ (9760 vs 16350 K) because of dense, cold
  cells at the corrugated photosphere; those cells are in equilibrium with the radiation cell by cell. In the CSM
  T_gas and T_rad decouple: at t = 0, T_gas 630 K vs T_rad 1500-2400 K at 7e13-2e14 cm.

Solver settings (production):
- Newton-FGMRES with block Gauss-Seidel preconditioning (`IdortPrecSweepMode 4`: task-local symmetric sweeps, halo after
  each pass); DSA compiled but off (`IdortDsa 0`).
- Tolerances: per-cell content norm `IdortTol 1e-3`, energy-weighted global residual `IdortTolGlobal 1e-4`,
  `IdortKrylovRtol 0.1`, `IdortKrylovMaxIter 100`. Also `IdortGmresOrth 1` (CGS2), `IdortKrylovBasisF 1`.
- Positivity: (cell, angle) pairs with I < -0.1 IdortTol s_i are flagged and their outflow faces fully upwinded; at most
  2 re-solves (`IdortPosRounds 2`, `IdortPosLag 0`, `IdortPosNoise 0.1`). The flags persist and are cleared at full
  syncs every 32 top steps. Nothing is clamped.
- Health over the whole run: max |res_E| 1.2e-16, max |res_p| 6e-18, created energy 0, no FCT, no negative content,
  every solve converged. The 920 "skips" are empty active sets (benign).

Thick-core treatment, f-hybrid (`IdortFhybMode 1`; examples/idort/README.md "Conservative f-hybrid"):
- Weight: in optically thick cells a share f of the radiation energy and pressure is carried in the EOS.
  f = smootherstep(ln(tau_cell/20)/ln(600/20)) with tau_cell = chi V^(1/3). f = 0 below tau_cell 20 and 1 above 600;
  f may grow by at most 0.5 per cell from an f = 0 neighbour.
- Coupling: the DO still solves the full field. The hydro carries the face flux of the equilibrium field f B(T0)
  (central, through the DO faces), and the DO source is compensated by exactly that share. The scheme is exactly
  conservative, with the ledger G = f a T^4 V.
- Where f = 1 in g25: f = 1 out to 5.37e13 cm and f > 0.01 to 5.62e13 cm (0.92-0.96 R), the same at all checked times
  (k = 10-72).
  cons12 instead had a radius ramp: f = 1 for r < 4.99e13, ramping to 0 at 5.40e13.
- What it fixes: without f (DO everywhere, "option b") the thick core develops grid-scale non-radial motions; with
  f-hybrid and the point-mass energy fix, core v_rms equals cons12 (point 7).

Boundary conditions for the radiation:
- Outer: WALL_SOLID cells at the box faces supply incoming intensity 0 (vacuum); outgoing radiation leaves there.
- Inner: none (the core is part of the domain). The heating enters as gas energy (point 6).

L(r) diagnostics: face-flux shell log through 48 log-spaced spheres (8.8e10 - 2.33e14 cm, plus R_c) at every top step
(`idort_shell_flux.txt`), and snapshot fields ERad, FRad (lab), PRad, IdortTemperature, IdortOpacity, IdortFhyb, IV.

## 6. Boundaries, outflow and inner boundary

Outer edge:
- No sponge and no damping.
- The box-face layer of the outer 32^3 background lattice (cell size 5.47e14 cm) is WALL_SOLID: frozen gas,
  rho 1e-20, T ≈ 500 K, radiation vacuum. The next layer is WALL_FLUID (slip wall).
- Gas between the star and the walls is the tenuous IC CSM/background (rho 1e-17 -> 1e-20). It is free to move, and it
  falls in coherently at 2-16 km/s over the run.
- Light-travel delay: the box faces are 8.76e15 cm away, so the escaping L in the budget lags the photospheric L by
  about 3e5 s. Compare L(r <= 3e15) rather than L_esc.
- The old code used IDORT_OUTFLOW at the box boundary (mirror ghost faces) instead.

Inner region (no boundary):
- Gas extends to r = 0 around the fixed point mass.
- Heating: `SPHERICAL_SOURCE` adds L_c = 6.854917e38 erg/s to the gas energy with W(r; h) = cubic spline, h = 3.522e11 cm
  (0.2 R_c), normalised by the discrete sum over the cells. Injected total = L_c t to 3e-6.
- No external gravity table.
- Opacity for r < R_c: kappa_eff (point 4).

## 7. Run status, cost and diagnostics

Status:
- **long_v1 completed** t = 0 -> 3.725e6 s (43.1 d) on 560 MPI ranks: 5 nodes x 112 cores, Intel Xeon Platinum 8480+.
  That is about 0.24 sound-crossing times (the IC builder's acoustic time is 1.53e7 s = 177 d) and about 0.4
  dynamical times sqrt(R^3/GM) = 8.6e6 s.
- 38 snapshots (every 1.0186e5 s), 10 GB each late; `long_v1/output` is 366 GB.
- **This is still the relaxation transient from an unrelaxed IC, not a steady convective state** (point 8).
- Restart: job 1 died at 1.73e6 s (base-code NGBNL "GrDir exhausted", now fixed). Job 2 restarted from 1.31e6 s with the
  classic neighbour gather. The restart left no imprint (v1/m8 ratios 1.000 ± 0.5 %).
- A twin run, long_m8 (earlier build d05fae0cd8 + point-mass energy fix, tighter IdortTol 1e-4), agrees with long_v1 within 1 %.

Cost (wall per 7275.58-s top step from cpu.txt; 1 day = 11.875 top steps; core-h = 560 x wall):

| build / phase | wall per top step | wall per sim day | core-h per sim day |
|---|---|---|---|
| long_v1 whole run (36f0bdb1e0, -O2), 9.30 h total (+0.9 h lost to the crash) | mean 65 s | 0.216 h (12.9 min) | 121 |
| long_v1 early, top steps 1-180 (t < 1.31e6 s, 6-14M cells) | 39.6 s | 0.131 h | 73 |
| long_v1 late, top steps 400-511 (17-18M cells) | 89.4 s | 0.295 h (17.7 min) | 165 |
| 0ab6a8294a (pushed; bitwise speed-up, -24 % per top step on the late bed), late, scaled | ~67 s | ~0.22 h | ~124 |
| same code built -O3 (local commit 07d4413e0f makes -O3 the production build; -29 %), late, scaled | ~63 s | ~0.21 h | ~116 |
| estimate for a full run with 0ab6a8294a / -O3 (1.31x / 1.23x cons12) | | ~0.16 / ~0.15 h | ~91 / ~85 |
| old code cons12, Orion rerun, same 560 ranks, 5.32 h to 3.73e6 s | 41-43 s late | 0.123 h | 69 |

Notes on cost:
- Radiation is 82 % of the wall; about 55-60 % of a late full sync is positivity re-solves.
- The -24 % / -29 % come from a 3-full-sync late-state test bed (`STATUS_speed2.md`); the rows using them are scalings,
  not a measured full run. Run-to-run (node placement) noise is about ±5-15 %.
- Tables: `cost_table.txt`, `seg_cost.txt`.

Diagnostics we track:
- L(r,t): shell-flux log every top step and snapshot flux (lab and comoving).
- R_ph(t): tau = 2/3 on the shell-averaged rho (kR + kS).
- T_eff from L(3e14) and R_ph.
- Mass-weighted v_rms in the core (r < R_c, r < 1e13) and in envelope shells; E_kin.
- Radial profiles of rho, T_gas, T_rad, kappa, f.
- Energy ledgers: radiation closure res_E, global gas + radiation + gravity, L_esc.
- Solver statistics.
- Pulsations: **not analysed** (snapshot cadence 1.0e5 s, run length 0.24 acoustic times).

Files with these:
- `timeseries_v1_m8_old.txt`: per snapshot, long_v1 / long_m8 / cons12; column definitions in the header.
- `vrms_timeseries.txt`: shell v_rms, long_v1.
- `shellflux_Lt_v1_m8_old.txt`: L(t) at 8 radii, every top step.
- `stats_Lr_v1_m8_old.txt`: L(r) lab and comoving, and profile ratios at k = 10, 20, 40, 60, 72.

### Numbers to compare against (long_v1; cons12 in brackets where available)

Luminosity at given radii, from the snapshot flux (lab frame; comoving where it differs). L_c = 6.855e38.

| t [s] | L(1e13) com. | L(5e13) | L(1e14) | L(3e14) | L(3e15) | L(8e15) |
|---|---|---|---|---|---|---|
| 0 | 1.30e38 | 5.33e37 | 4.53e38 | 4.55e38 | 4.32e38 | 4.47e38 |
| 1.02e5 | 7.96e37 | 5.31e37 | 5.43e38 | 5.52e38 | 5.23e38 | 4.51e38 |
| 5.09e5 | 7.92e37 | 5.28e37 | 2.90e38 | 2.79e38 | 3.17e38 | 4.18e38 |
| 1.02e6 | 1.05e38 | 5.29e37 | 2.39e38 | 2.27e38 | 2.38e38 | 2.52e38 |
| 2.04e6 | 1.03e38 | 5.15e37 | 2.73e38 | 2.59e38 | 2.65e38 | 2.63e38 |
| 3.06e6 | 1.09e38 | 5.25e37 | 2.60e38 | 2.46e38 | 2.65e38 | 2.80e38 |
| 3.67e6 | 1.04e38 (1.28e38) | 6.42e37 (5.79e37) | 2.10e38 (2.33e38) | 1.99e38 | 2.14e38 (2.44e38) | 2.32e38 (4.08e38, old T floor) |

Columns L(1e14) and beyond are the surface luminosity. Snapshot L(1e14) and the shell log agree.

Late-time L(r) at t = 3.667e6 s, comoving, in units of L_c:

| r [cm] | 3e11 | 1e12 | R_c | 3e12 | 1e13 | 2e13 | 3e13 | 5e13 | 5.6e13 | 6.5e13 | >= 1e14 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| L / L_c | 0.97 | 0.92 | 0.33 | 0.66 | 0.15 | 0.17 | 0.23 | 0.08 | 0.18 | 0.31 | 0.29-0.34 |

Notes:
- cons12: 1.00 at R_c, 0.90 at 3e12, 0.19 at 1e13; equal to long_v1 from 2e13 outward.
- The rest of L_c across R_c-1e13 is carried by the gas (not extracted directly).
- Inside 1e12 cm the lab-frame L is ±(1-10) L_c with alternating sign: v_r E noise of the ~1 km/s core motions.
  Use the comoving L there.

Escaping luminosity at the box walls, L_esc / L_c (window means). The old code had 0.73-1.24, from its CSM T floor.

| window [s] | L_esc / L_c |
|---|---|
| 0-2.5e5 | 0.68 |
| 2.5-5e5 | 0.75 |
| 5-7.5e5 | 0.55 |
| 1-1.25e6 | 0.37 |
| 2-2.25e6 | 0.39 |
| 3-3.25e6 | 0.41 |
| 3.5-3.75e6 | 0.35 |

Photosphere radius R_ph(t) [cm] and T_eff(t) [K]. cons12 at 3.67e6 s: 5.828e13, 3100 K.

| t [s] | R_ph | T_eff |
|---|---|---|
| 0 | 6.055e13 | 3634 |
| 1.02e5 | 6.046e13 | 3815 |
| 4.07e5 | 6.000e13 | 3321 |
| 8.15e5 | 6.006e13 | 3085 |
| 1.22e6 | 5.993e13 | 3086 |
| 1.63e6 | 5.946e13 | 3142 |
| 2.04e6 | 5.902e13 | 3196 |
| 2.45e6 | 5.871e13 | 3246 |
| 2.85e6 | 5.854e13 | 3212 |
| 3.26e6 | 5.843e13 | 3126 |
| 3.67e6 | 5.833e13 | 3010 |

R_ph shrinks slowly. The mass-weighted radius <r>_m = 2.162e13 cm is constant to < 0.1 %: no envelope expansion.

Velocities (mass-weighted rms |v|, from `vrms_timeseries.txt`) and E_kin:

| t [s] | r < R_c | r < 1e13 | 0.45-0.55 R | 0.4-0.6 R | 0.9-1.0 R | E_kin [erg] |
|---|---|---|---|---|---|---|
| 2.04e5 | 0.33 km/s | 0.28 | 0.33 | 0.33 | 0.35 | 1.0e43 |
| 1.02e6 | 0.46 | 0.29 | 0.079 | 0.078 | 1.3 (radial) | 3.4e43 |
| 2.04e6 | 0.76 | 0.45 | 0.067 | 0.072 | 3.6 (radial) | 8.4e43 |
| 3.67e6 | 1.00 (cons12 1.00) | 0.69 (cons12 0.65) | 0.105 | 0.106 | 5.4 (radial) | 1.13e44 (cons12 1.15e44) |

Notes on the velocities:
- At 0.5 R the motion is mostly radial (v_r rms 0.097 of 0.105 km/s) and weak: convection has not developed yet in 43 d.
- In 0.9-1.0 R the photospheric layer falls in at 10-16 km/s (about 0.2 Msun/yr inward, an IC relaxation transient)
  and decelerates in a standing shock at 0.9-1.0 R.
- The shell-averaged mean v_r at 5.6-6.8e13 cm late is -6 to -17 km/s (profile_late, column vr_mean).

Energy budget at 3.70e6 s, in L_c t:
- injected 1.000, escaped 0.444;
- dE_int +0.29, dE_kin +0.045, dE_pot +0.62, dE_rad +0.60, dG +0.63;
- global non-conservation **+0.374** (cons12 +0.350, of which 0.49 is its f > 0 nonclosure term). See point 8.

## 8. Known issues

1. **Not a relaxed star yet.**
   - The IC is mapped from MESA without relaxation. Within the first ~1e6 s the surface luminosity falls from
     0.68-0.8 L_c to about 0.35 L_c and T_eff from 3600-3800 K to 3000-3250 K (MESA: 4088 K).
   - The radiative L inside the envelope stays at 0.08-0.2 L_c (MESA carries the rest by convection, which has not
     developed: v_rms at 0.5 R ~ 0.1 km/s after 43 d).
   - The photospheric layer falls in (Mdot ≈ 0.2 Msun/yr inward) onto a standing shock at 0.9-1.0 R.
   - The IC CSM (r^-2, rho 1e-17 at R) falls in coherently.
   - 3.7e6 s is only 0.24 acoustic times. The old-code runs on Fritz go to 2.6-3e7 s (results requested, not received).
2. **Photosphere marginally resolved.**
   - 0.5 (t = 0) -> 3.4 cells per H_rho; the density cliff 1e-9 -> 1e-17 lies within one 0.02-dex shell.
   - In the 5.8-6.5e13 band, rho differs from cons12 by +13 % (median, up to 24 %), although both use the same mesh
     rules.
   - The positivity flags (inward rays) all sit in R* < r < R_ph. The luminosity there carries a resolution error of
     a few %, which shrinks under refinement (modules_idort.md, photosphere resolution test).
   - 1.2-3 R is under-resolved for any extended atmosphere (dx/H about 3-10).
3. **Core flux structure.**
   - The comoving radiative L dips to 0.33 L_c at the heating edge R_c (cons12: 1.00); gas carries the rest.
   - Lab-frame L in the core is ±10 L_c noise from v_r E.
   - The structure and v_rms equal cons12 to < 1 %. The difference is the f-hybrid vs the old kappa_eff +
     radiation-in-EOS treatment, not a failure.
   - Before the point-mass energy fix the core spun up to 17 km/s (long_fhyb); fixed with GRAVITY_PM_ENERGY_CONSERVING.
4. **Global energy non-conservation of +0.37 L_c t over the run** (cons12 +0.35). The radiation ledger closes to 1e-16,
   so the residual is on the gas/gravity side (dE_pot +0.62 L_c t). **[uncertain: cause not identified]**
5. **Differences to the old code cons12** (decisions D1-D7 in `PARAMS_cons12.md`):
   - CSM gas T floor (531 K, old) vs none (new): old L(8e15) ≈ 4e38, L_esc ≈ L_c; new CSM 280 K, flat L. Compare
     L(r <= 3e15).
   - Face closure linear (new) vs quadratic (old).
   - f ramp at 5.37-5.62e13 (new) vs 4.99-5.40e13 (old), so L(5e13) is +11 % (new).
   - Opacity evaluation time: start of the radiation step (new) vs one hydro step lagged (old).
   - No dv caps or limiters in the new code.
6. **Opacity physics is simplified.**
   - Grey; a single X = 0.7 table against EOS X = 0.654.
   - kappa is floored at its 3200 K value, so the molecular opacity minimum and the MOLsphere are absent.
   - kP = kR in thin gas (rho < 1e-11).
   - No line cooling and no non-LTE.
7. **Cost.** The new code is 1.7x the old code cumulatively at the same rank count (-O2, before the speed-up); about
   1.3x / 1.2x expected with 0ab6a8294a / -O3. Positivity re-solves are the main remaining cost.
8. **Base-code issues** (Arepo2 node-local neighbour graph, VORONOI_CELL + refinement; `ARE2_BASE_BUGS_NOTE.md`):
   - six fixes (GrBPatQ, edge pool, lone-gravity domain work, GrDir, entry headroom, graph-on restart segfault), all
     on the pushed branch (1-3 with basefix 09479fc441, 4-6 with c68c89994d);
   - long_v1 used build 36f0bdb1e0, which needed `AREPO_NGBNL_GRAPH=0` after the crash;
   - graph on and graph off are not bitwise identical (round-off at refinement thresholds; not investigated).
9. **Small items.**
   - Total mass 0.2 % above MESA (not investigated).
   - IC names bgr2/coreB/opal/c1s inferred from the files (point 2).
   - 920 empty-active-set skips (benign).

### If you build the same star in AthenaK

The minimal set that defines this model:
- **Profile:** the MESA profile (`mesa_g25_p438_selected.txt`).
- **Core:** for r < 2 R_c use the c1s core (`core_c1s_profile.txt`: rho, T, P, M_gas, kap_eff), with M_pm 2.11248e34 g,
  softening kernel h ≈ R_c, and heating L_c 6.855e38 inside 0.2 R_c with a cubic-spline kernel.
- **Composition:** X = 0.6538, Z = 0.0198 uniform.
- **EOS:** OPAL EOS5 Z = 0.02 gas-only + Saha blend, radiation separate.
- **Opacities:** the two gs98 tables with the five modifiers of point 4.
- **CSM:** rho 1e-17 (r/5.87e13)^-2 to a 1e-20 floor, T ≈ 630 K.
- **Outer boundary:** vacuum at 8.76e15 cm.

The IC gas state itself (6.09M cells, point values at the generators) is in the 731-MB `g25_full.hdf5` if you prefer to
map that. Snapshots of long_v1 (rho, u, v, T, ERad, FRad, IV, ...) can be shared on request (10 GB each).

---

## Files

In `docs/handover/rsg-g25-orion/`:

| file | size | content |
|---|---|---|
| `param_long_v1.txt` | 10 KB | the exact param.txt of long_v1 |
| `Config.sh` | 3.7 KB | compile options (identical to the options printed in long_v1/run.log) |
| `PARAMS_cons12.md` | 14 KB | every cons12 option/parameter -> module equivalent, decisions D1-D7 (written 2026-10-14; since then f-hybrid is on (D2) and the solver tolerances are the production ones (D5), see param_long_v1.txt) |
| `convert_ic.py` | 4.5 KB | cons12 IC -> module IC (boundary roles) |
| `create-ics-golden.py` | 9.6 KB | base IC generator (stellar_ics, Orion copy) |
| `mesa_g25_p438_selected.txt` | 414 KB | MESA profile 438, 15 selected columns + header values |
| `core_c1s_profile.txt` | 474 KB | c1s core model (attributes + every 2nd of 4280 rows, 17 columns) |
| `gs98_z0.02_x0.7_ma_ross.data`, `..._planck.data` | 53 + 50 KB | opacity tables |
| `profile_ic_t0.txt` | 46 KB | shell profile at t = 0 (after the steady radiation start) |
| `profile_late_t3.667e6.txt` | 54 KB | shell profile at t = 3.667e6 s |
| `vrms_timeseries.txt` | 6 KB | v_rms (total, radial) in 6 shells + E_kin, every snapshot |
| `timeseries_v1_m8_old.txt` | 38 KB | per-snapshot L(r), R_ph, T_eff, v_rms, E_kin, H/dx: long_v1, long_m8, cons12 |
| `shellflux_Lt_v1_m8_old.txt` | 19 KB | L(t) at 8 radii from the shell-flux logs |
| `stats_Lr_v1_m8_old.txt` | 17 KB | L(r) lab/comoving and profile band ratios at 5 times |
| `cost_table.txt`, `seg_cost.txt` | 2 + 2 KB | wall per top step by phase and component |

On Orion (large):

| item | path |
|---|---|
| MESA profile (full) | `/orion/u/jinma/AREPO/RSGg_g25_IDORT/ics/m22_p438.data` (3.9 MB); all runs in `/orion/ptmp/jinma/Arepo/golden.zip` (2.6 GB) |
| production IC (cons12) | `/orion/u/jinma/AREPO/RSGg_g25_IDORT/ics/rhoin-RSGg_golden25_p438_hr_bgr2_coreB_opal_c1s.ic.hdf5` (755 MB) |
| module IC | `/orion/ptmp/jinma/idort_dev/runs_g25/ics/g25_full.hdf5` (731 MB); 1.19M-cell smoke IC `g25_red.hdf5` |
| base IC | `/orion/ptmp/jinma/Arepo/ICs_golden/ic/RSGg_golden25_p438_hr/` |
| core table | `/orion/u/jinma/AREPO/RSGg_g25_IDORT/ics/coreB_side_g25_c1s.hdf5` (3 MB) |
| OPAL EOS | `/orion/u/jinma/AREPO/Arepo-snap-util/eostable/EOS5_data` (15 MB) |
| run | `/orion/ptmp/jinma/idort_dev/runs_g25/runs/long_v1/` (output 366 GB, snapshots 0-37, logs idort_budget/solver/shell_flux) |
| analysis | `/orion/ptmp/jinma/idort_dev/phys_cmp/v1/` (prof/*.npz shell sums, budget/, figures), `STATUS_phys_cmp.md` |
| old-code rerun | `/orion/ptmp/jinma/idort_dev/cons12_orion/` |
| code | `/orion/u/jinma/AREPO/Arepo2-IDORT-dev` (Bitbucket, branch Arepo2-IDORT-dev, 0ab6a8294a) |
