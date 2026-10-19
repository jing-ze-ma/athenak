# Mass-gaining star in an Algol binary: 2-D r-phi design (RY Per first, u Her second)

Status: design only. No code, no runs. Written 2026-10-06.

Files in this directory:

- `lit/`: the source papers and their text dumps (made with pypdf):
  - `barai2004.txt`: Barai et al. 2004, astro-ph/0309734.
  - `1604.07589.txt`: Van Rensbergen & De Greve 2016.
  - `1408.2681.txt`: Kolbas et al. 2014.
  Line numbers below refer to these dumps.
- `scripts/roche_stream.py`: the ballistic L1 stream integrator. It works in the restricted 3-body problem in the co-rotating frame. Every number marked ESTIMATE in sections 1 and 3 comes from it.

  Usage: `python3 roche_stream.py Macc Mdon a_Rsun P_day Racc_Rsun Tdon`

Conventions:

- **QUOTED** means the number is a verbatim quote from the paper given.
- **ESTIMATE** means I computed it; the formula or script is given.
- Azimuth phi is measured at the accretor:
  - phi = 0 points to the donor.
  - phi > 0 points to +y. With Omega along +z, +y is the **trailing** hemisphere of the accretor, because the accretor moves along -y.

---

## 1. Literature parameters

### 1.1 RY Per (HD 17034)

| Quantity | Value | Status | Verbatim quote and source |
|---|---|---|---|
| Period | 6.863569 d | QUOTED | "we adopted in our radial velocity analysis the orbital period assumed by Olson & Plavec (1997), P = 6.863569 d" (Barai+04 §3; Table 4 "P (d) ... 6.863569") |
| Eccentricity | 0.036 +- 0.005 | QUOTED | Table 4: "e ... 0.0 b 0.036 ± 0.005" (Barai+04; Popper 1989 value fixed at 0). Treated as circular here. |
| M_gainer (spectroscopic) | 6.24 +- 0.24 Msun | QUOTED | Table 4: "MP (M⊙) ... 6.24 ± 0.24" (Barai+04) |
| M_donor (spectroscopic) | 1.69 +- 0.23 Msun | QUOTED | Table 4: "MS (M⊙) ... 1.69 ± 0.23" (Barai+04) |
| Separation a | 30.3 +- 0.6 Rsun | QUOTED | Table 4: "a (R⊙) ... 30.3 ± 0.6" (Barai+04) |
| Inclination | 83.0 +- 0.3 deg | QUOTED | "the orbital inclination found by Olson & Plavec (1997), i = 83.◦0±0.◦3" (Barai+04 §6) |
| Masses and radii (photometric + spectroscopic) | M_P = 6.25 +- 0.16, R_P = 4.06 +- 0.14; M_S = 1.60 +- 0.10, R_S = 8.10 +- 0.17 (Msun, Rsun) | QUOTED | "MP/M⊙ = 6.25 ± 0.16 and RP/R⊙ = 4.06 ± 0.14 for the B4: V gainer and MS/M⊙ = 1.60 ± 0.10 and RS/R⊙ = 8.10 ± 0.17 for the F7: II-III donor star" (Barai+04 §1, citing Olson & Plavec 1997, AJ 113, 425) |
| Teff and log g | 18000 K / 4.02 (gainer); 6250 K / 2.83 (donor) | QUOTED | "assumed effective temperatures of Teff = 18000 and 6250 K and gravities of log g = 4.02 and 2.83 for the primary and secondary, respectively; Olson & Plavec (1997)" (Barai+04 §5) |
| Teff, alternative | log Teff = 4.259 (gainer), 3.802 (donor) | QUOTED | Table 1 row "RY Per 6.68356 1.6 6.25 8.1 4.06 3.802 4.259", columns "P Md Mg Rd Rg Log Teff,d Log Teff,g". The text says "temperatures for RY Per are from Peters & Polidan (2004)" (Van Rensbergen & De Greve 2016, arXiv:1604.07589, Table 1). Their P = 6.68356 d disagrees with 6.8636 d and is probably a typo. |
| Luminosity | M_bol about -3.3 (gainer), -0.2 (donor) | QUOTED (see caveat) | "more luminous (M_bol about -3.3 and -0.2 mag, respectively) than a typical Algol system" (Olson & Plavec 1997, AJ 113, 425, ADS abstract). Caveat: this text came from a search-engine copy of the ADS abstract, because ADS refused the direct fetch. Check it before citing. |
| Luminosity, derived | L_P ~ 1600 Lsun | ESTIMATE | From the M_bol above: 10^(-0.4(-3.3-4.74)) = 1640. From R = 4.06 and Teff = 18000: 4.06^2 (18000/5772)^4 = 1560. |
| v sin i (gainer) | 213 +- 10 km/s (UV); 212 +- 7 km/s (optical) | QUOTED | "best fit ... occurred for V sin i = 213 ± 10 km s−1, which is in good agreement with the value of V sin i = 212 ± 7 km s−1 found by Etzel & Olson (1993) from optical lines" (Barai+04 §4) |
| Synchronicity | 7.2 +- 0.3 | QUOTED | "This implies that the primary is spinning (7.2 ± 0.3)× faster than the synchronous rate" (Barai+04 §4) |
| Synchronicity, older value | ~10 | QUOTED | "the gainer is rotating approximately 10 × faster than the orbital synchronous rate" (Barai+04 §1, paraphrasing Olson & Plavec 1997) |
| v_eq; v_eq/v_crit | 214.6 km/s; 0.40 | ESTIMATE | v_eq = 213 / sin 83 deg. v_crit = sqrt(G M_P / R_P) = 542 km/s. v_sync(R_P) = Omega R_P = 29.9 km/s, which gives the 7.2x consistently. |
| Mass-transfer rate | **no measurement found** | QUOTED (absence) | "No trustworthy determination of the mass transfer rate was found in this case." (Van Rensbergen & De Greve 2016, §3). Their best-fit binary-evolution model gives "RY Per 5.25+2.6 ... 1.4E-5" Msun/yr (Table 3, last column dM/dt). That is a **model** value, not a measurement. Typical Algol range in the same literature is 1e-11 to 1e-7 Msun/yr (search summary only, unquoted). Treat Mdot as unknown, 1e-7 to 1e-5 Msun/yr. |
| Roche radius of gainer | 15.0 Rsun | QUOTED | "the disk is assumed to be relatively thin and to extend from the photosphere to the Roche radius of the primary, RRoche = 15.0R⊙" (Barai+04 §6) |
| Disc / circumstellar evidence | persistent, variable, Keplerian disc; elongated along the line of centres | QUOTED | "the mass gaining primary is surrounded by a persistent but time variable accretion disk" (abstract). "The overall best fit was obtained with Keplerian motion and m = 5.5 ± 1.5" (§6). "The base density in this case was n0 = 9 × 10−4 cm−3, which corresponds to an electron density of Ne ≈ 2 × 108 cm−3" (§6). "The resulting surface density distribution is elongated along the axis joining the stars" (abstract). |
| Disc size | 0.85 of the gainer's Roche radius | QUOTED | Table 2 row "RY Per 0.616 7969 12456 0.85 0.80"; Obs.Size = 0.85. "The size is the radius of the disk divided by the Roche radius of the gainer"; sizes "from ... Sudar et al. (2011) for RY Per" (Van Rensbergen & De Greve 2016). I did not read Sudar et al. 2011 itself. |
| Stream geometry | direct, near-tangential hit on the trailing edge | QUOTED | "the gas stream strikes the star very close to its trailing edge (Fig. 8), and this is clearly a favorable situation to impart angular momentum" (Barai+04 §7). "RY Per falls in the region between the curves [ϖmin and ϖd], and other systems in this same region have variable disks (Peters 2001)" (§7). "RY Per straddles the border between permanent and transient accretion disks" (Olson & Plavec 1997 abstract; same caveat as the M_bol row). |
| Hot impact plasma | O VI (300 kK) emission in totality | QUOTED (secondary) | FUSE spectrum "shows only circumstellar emission lines from ionized species ranging from O VI (formed in a 300 kK plasma) to N II ..., probably in the region of the gas stream impact". Source: ResearchGate abstract of "FUSE Observations of the Active Interacting Binary RY Persei" (Peters & Polidan), seen in a search summary only and not read. |
| Distance | ~840 pc | QUOTED | "Because of its distance (about 840 pc)" (Olson & Plavec 1997 abstract, same caveat) |

Not found or not read: Olson & Plavec 1997 full text; Plavec et al. papers on RY Per; Richards & Albright 1999 RY Per tomogram; Sudar et al. 2011; any measured Mdot or dP/dt. The Barai+04 numbers above already carry the Olson & Plavec values.

#### Derived Roche and stream quantities for RY Per

All values below are ESTIMATES from `scripts/roche_stream.py 6.24 1.69 30.3 6.863569 4.06 6250`.

Inputs and frame:

- q = M_d/M_a = 0.271.
- Omega = 2 pi / P = 1.0595e-5 s^-1.
- a Omega = 223 km/s.

Roche geometry:

- L1 lies at 0.631 a = **19.1 Rsun** from the accretor centre. It is found from the root of the x-acceleration on the line of centres.
- The Eggleton (1983) Roche-lobe radius of the accretor is 0.494 a = **14.98 Rsun**. Barai+04 quote 15.0, which checks.
- The donor's Eggleton lobe is 8.29 Rsun; R_S = 8.10 Rsun, so the donor fills it.

Ballistic stream:

- Donor sound speed: c_s = sqrt(k T / (1.27 m_H)) = 6.4 km/s at 6250 K, mu = 1.27 (neutral photosphere).
- The ballistic start is at L1 with speed epsilon = c_s/(a Omega) = 0.029 toward the accretor.
- Stream-width scale: c_s/Omega = **0.86 Rsun**. This is the Lubow & Shu (1975) scaling; the O(1) prefactor is not computed.
- Closest approach with no star: r_min = **0.0889 a = 2.69 Rsun**.
  - Lubow & Shu fit: 0.0488 q^-0.464 = 0.0895 a.
  - Circularisation radius fit: r_circ ~ 0.0859 q^-0.426 a = **4.54 Rsun**.
  - Both fit formulas are as commonly quoted from Lubow & Shu 1975 and Ulrich & Burger 1976; the fit coefficients were not re-quoted from those papers.
- Since R_P = 4.06 > r_min = 2.69 and R_P < r_circ = 4.54, the stream gives a **direct but grazing impact**. This is the "between the curves" regime of Barai+04.

Impact on the surface:

- Impact position: phi = **+75 deg**, on the trailing side near the limb.
- Velocities:
  - Rotating frame: v_r = -349 km/s, v_phi = +551 km/s.
  - Inertial frame, centred on the accretor: v_phi = 581 km/s.
- The velocity is **58 deg from the inward normal**.
- Specific angular momentum: j_impact / j_Kep(R_P) = **1.07**.

Stream at candidate outer radii (rotating-frame velocities):

| r_out | r (Rsun) | phi (deg) | v_r (km/s) | v_phi (km/s) |
|---|---|---|---|---|
| 0.9 d_L1 | 17.2 | 2.6 | -49 | +23 |
| 0.8 d_L1 | 15.3 | 5.9 | -96 | +50 |
| 0.7 d_L1 | 13.4 | 10.1 | -143 | +82 |

#### Other ESTIMATES

Photospheric scale height of the accretor:

- g = 1.04e4 cm s^-2.
- H_p = k T / (mu m_H g) = 2.3e8 cm = **3.3e-3 Rsun** (mu = 0.62, T = 18 kK).
- H_p/R = 8e-4.

Spin-up time:

- Assume the stream brings j_Kep(R_P) = sqrt(G M R) = 1.5e19 cm^2/s per gram.
- Assume J_* = k^2 M R^2 Omega_*, with k^2 = 0.06 (assumed). This gives J_* = 4.5e51.
- Then t_spin = J_* / (Mdot j_Kep), which is:
  - 1.5e5 yr at 1e-6 Msun/yr;
  - 1.5e4 yr at 1e-5 Msun/yr.
- The simulation therefore measures the **torque** (spin-up rate). It cannot follow the spin evolution.

Accretion luminosity: G M Mdot / R = 48 Lsun x (Mdot / 1e-6 Msun/yr). This is about 3 % of L_P at 1e-6.

### 1.2 u Her (contrast case)

Source for everything quoted: Kolbas et al. 2014, MNRAS 444, 3118, arXiv:1408.2681.

| Quantity | Value | Status | Quote |
|---|---|---|---|
| Period | 2.05102685 d | QUOTED | Table 1: "Orbital period P d 2.05102685" |
| Masses, radii, Teff | M_A 7.88 +- 0.26, R_A 4.93 +- 0.15, Teff 21600 +- 220 K; M_B 2.79 +- 0.12, R_B 4.26 +- 0.06, Teff 12600 +- 550 K | QUOTED | Abstract: "(the mass-gaining star) we find MA = 7.88 ± 0.26 M⊙, RA = 4.93 ± 0.15 R⊙ and Teff,A = 21 600 ± 220 K. For the secondary (the mass-losing star) we find MB = 2.79 ± 0.12 M⊙, RB = 4.26 ± 0.06 R⊙ and Teff,B = 12 600 ± 550 K" |
| Separation | 14.95 +- 0.17 Rsun | QUOTED | Table 2: "Semimajor axis R⊙ 14.95 ± 0.17" |
| q, inclination | 0.354; 78.9 +- 0.4 deg | QUOTED | Table 1: "Mass ratio q 0.354"; "Orbital inclination deg 78.9 ± 0.4" |
| Luminosity | log L = 3.68 +- 0.03 (A), 2.63 +- 0.08 (B) | QUOTED | Table 2: "log L L⊙ 3.68 ± 0.03 2.63 ± 0.08" |
| Rotation | v sin i 124.2 +- 1.8 vs V_synch 121.7 +- 3.5 km/s, so **synchronous** | QUOTED | Table 2: "Veq sin i km s−1 124.2 ± 1.8 107.0 ± 2.0"; "Vsynch km s−1 121.7 ± 3.5 105.0 ± 1.5" |
| Mass-transfer rate | not found | none | Not given in Kolbas+14 (grep for "transfer rate" and "yr−1" found nothing). |

u Her stream ESTIMATES, from `roche_stream.py 7.88 2.79 14.95 2.05102685 4.93 12600`:

- Distances:
  - d_L1 = 9.04 Rsun.
  - Accretor lobe R_L = 7.04 Rsun, so R_A/R_L = 0.70.
  - r_min = 1.16 Rsun.
  - r_circ = 2.0 Rsun.
- Impact: **deep and direct**, at phi = +18 deg.
  - The velocity is 33 deg from the normal.
  - j_impact / j_Kep = **0.65**.
- Check: v_sync(R_A) = 121.6 km/s, which reproduces the quoted 121.7.

Contrast with RY Per: u Her has a deep hit with sub-Keplerian j and a synchronous gainer. RY Per has a grazing hit with j about Keplerian and a gainer at 7.2x synchronous.

---

## 2. Reuse map in AthenaK (rt-integration, read only)

### Coordinates

- Spherical-polar is a real `Coordinates` geometry, switched on with `mesh/use_spherical_polar` (src/mesh/mesh.cpp:233; flag at src/mesh/mesh.hpp:310).
- The axes are fixed: x1 = r, x2 = theta, x3 = phi (src/coordinates/coordinates.cpp:1611 `CoordSphericalPolar`).
- Optional radial stretching is available (mesh.cpp:234-240).
- Geometric sources: coordinates.hpp:99 and :162 (`SrcTermsCurvilinearWB`).
- **A true 2-D r-phi run is impossible.** An x1-x3 plane is fatal (mesh.cpp:664), and every active dimension needs at least 4 cells (mesh.cpp:652-654, :712).
- **Route: 3-D sp with nx2 = 4, a thin theta band symmetric about pi/2.** Existing examples of a symmetric periodic theta wedge:
  - inputs/radiation/he_presn_m1_wedge.athinput:39-48
  - inputs/radiation/bsg_m1_wedge.athinput:20-29
  - inputs/hydro/he4_presn_sp.athinput:72-81
- Cylindrical coordinates do not exist. The alternative is 2-D Cartesian x-y, but the star surface would be staircased, which is poor for a torque measurement.

### Rotating frame

- There is no generic module for it; it lives in the pgens.
- `SourceFunc` in src/pgen/deep_hot_jupiter_rt.cpp:4345. The rotating-frame block is at :4555-4590 (centrifugal Omega^2 r sin theta and Coriolis 2 Omega v_phi, on IM1, IM2, IM3 and IEN).
- `problem/rot_potential` (dhj_rt.cpp:2117) moves the centrifugal term into the potential through `TotPotAt` (:612).
- **The 09-28 double-count fix** is the guard `if (!(rotpot && use_etotgrav))` at dhj_rt.cpp:4574-4577 (commits d0e36fed, merge 8fa76784).
- There is a restart guard on rot_potential and omega in main.cpp:71-80.
- Note: dhj rotates about the planet's own axis. Here the rotation axis passes through the **binary centre of mass**, which is offset from the grid centre. The centrifugal potential must therefore be -1/2 Omega^2 [(x - x_cm)^2 + y^2] with x_cm = a q/(1+q) = 6.46 Rsun (ESTIMATE). Put it in the potential and keep Coriolis explicit.
- No indirect term is needed: the accretor is at rest in the rotating frame, and the rotating frame about the centre of mass is non-inertial only through its rotation.

### Potential and gravity

- src/srcterms has no point mass (srcterms.cpp:40-44), so gravity is handled per pgen through `user_srcs_func`.
- The `etotgrav` machinery takes an arbitrary cell potential `phicc0` plus face potentials `phi0.x1f/x2f/x3f`:
  - flags: hydro.cpp:182, hydro.hpp:179-180
  - allocation: hydro.cpp:249-257
  - `AddGravFlux`: hydro_tasks.cpp:466
  - `AddGravEtot`: :921
  - It does not support SMR/AMR (hydro.cpp:187).
- Fill pattern for the potentials: src/pgen/he_star_m1.cpp:663-692.
- Point-mass helpers:
  - red_giant.cpp:777-781
  - dhj `GravPotAt`/`TotPotAt` (dhj_rt.cpp:612)
- Well-balancing:
  - `wellbalance_dynamic`: hydro.cpp:194
  - `wb_x1`: :195
  - `wb_x2`: :225
  - effective-potential WB: hydro.hpp:333-341
  - There is no `wb_x3`. The Roche potential is phi-dependent, so WB can only take the radial part near the star.

### Boundary conditions

- `user_bcs_func` is declared in pgen.hpp:108-109 and called at hydro_tasks.cpp:724, after the built-in BCs (:719-720).
- The sp template is `HeStarBC` (he_star_m1.cpp:1559): a par_for over (m,k,j) with a test on `mb_bcs.d_view(m, outer_x1) == BoundaryFlag::user`.
  - Its hydrostatic top is at about :1604.
  - Its inner reflecting wall is at about :1693-1705.
- An alternative is `RedGiantBC` (red_giant.cpp:4564), with `problem/inner_bc` wall or open (:74-110).
- The built-in `inflow` sets one constant state per face (bvals.hpp:290; hydro_bcs.cpp:84ff).
- **No existing BC applies inflow on part of a face.** New code: a phi window on `x3v(m,k)` in the outer-x1 ghosts, with outflow and no-inflow elsewhere.

### Floors

src/eos/eos.cpp:38 (`dfloor`), :146-147 (`pfloor`, `tfloor`), :162 (`sfloor`), :39 (`dfloor_keep_velocity`).

### Stellar-envelope pgens

- he_star_m1 requires etotgrav + wellbalance_dynamic + wb_x1 (he_star_m1.cpp:319-322).
  - Gravity: `HeStarGravity` (:1429).
  - Envelope column and WB pair: :663-692.
- red_giant: `RedGiantGravity` (red_giant.cpp:3128), sp via `curv` (:1367).

### Binary / Roche / L1 / stream code

None exists anywhere in src (grep for roche, companion, lagrange, mass_transfer).

### Verdict

Reuse:

- the sp coordinates (thin 4-cell theta band) and the stretched radial grid;
- `etotgrav` + phi0 to carry the full Roche potential: two point masses plus the centre-of-mass centrifugal term, on all three faces;
- the Coriolis block and the double-count guard from dhj_rt;
- the `HeStarBC` skeleton for the user BCs;
- the floors, the `user_*_func` enrolment, and the optional `wb_x1` near the surface.

Write new:

- a pgen `src/pgen/algol_accretor.cpp` (user pgen, `-D PROBLEM=algol_accretor`). Because it is a user pgen, no src/CMakeLists change is needed. It contains:
  - the Roche potential fill;
  - the stream phi-window BC;
  - the rotating inner-surface BC;
  - the mass and angular-momentum flux diagnostics (`user_hist_func`);
- the input file.

---

## 3. Proposed setup (RY Per)

### Frame and geometry

- Co-rotating frame with Omega_orb about the binary centre of mass.
- Grid centred on the accretor; the donor point mass sits at (r = a, phi = 0) and is outside the grid.
- The orbit is circular: the quoted e = 0.036 is ignored in stage 1.

### Domain

Radius:

- r_in = R_P = 4.06 Rsun.
- r_out = 15.3 Rsun = 0.8 d_L1, about the Roche-lobe radius.
  - At r_out the stream crosses at phi = 5.9 deg with v_r = -96 and v_phi = +50 km/s in the rotating frame (ESTIMATE).
  - Option: r_out = 17.2 Rsun (0.9 d_L1), where the stream is slower (-49 / +23 km/s) but sits closer to the nozzle.

Phi: 0 to 2 pi, periodic.

Theta:

- 4 cells about pi/2, each with the theta spacing equal to the phi spacing, so the band is ±0.006 rad wide.
- Use reflecting theta faces. Periodic faces are the precedent from he_presn, but see the open decisions.
- Evaluate the potential and the centrifugal term **at theta = pi/2**: use the in-plane x, y, and give no theta gravity component. This keeps the band a pure 2-D model.

### Grid

- n_phi = 2048, so dphi = 3.07e-3.
- Log-r spacing dr/r = dphi gives n_r ~ 432, rounded to **448**.
- **Total: 448 x 4 x 2048 = 3.7e6 cells** (ESTIMATE).
- Resolution of the stream width (c_s/Omega = 0.86 Rsun):
  - at r_out it is 3.2 deg, about 18 cells;
  - near the impact, dr = 0.0125 Rsun.
- The surface scale height H_p = 3.3e-3 Rsun is **not** resolved: dr = 3.8 H_p at r_in. Resolving it with dr = H_p/4 would cut dt by about 15x, so stage 1 has no envelope (see below).

### Inner boundary (stage 1)

The star is a rigid rotating surface without an atmosphere:

- Gas may only leave the grid: a diode with v_r <= 0, so mass is absorbed into the star.
- Ghost density and pressure are hydrostatic extrapolations.
- Tangential ghost velocity:
  - No-slip: v_phi = (Omega_* - Omega_orb) r. In the rotating frame that is 185 km/s for Omega_* = 7.2 Omega_orb.
  - The alternative is zero-gradient (free-slip). See the open decisions.

### Inner boundary (stage 2, option)

- Add a thin envelope of about 20 H_p, with H_p inflated or genuine, on a stretched radial grid.
- Use the he_star_m1 pattern: etotgrav + wellbalance_dynamic + wb_x1, with an effective potential that includes the stellar rotation.

### EOS and thermodynamics (ESTIMATE)

Post-shock temperature, T_s = 3/16 mu m_H v_perp^2 / k:

- 1.7e6 K for v_r = 349 km/s;
- 6e6 K for the full 657 km/s.

Stream density, taking A ~ (c_s/Omega)^2 = 3.6e21 cm^2 and v = 650 km/s:

- rho ~ 2.7e-10 g/cc x (Mdot / 1e-6 Msun/yr).

Impact-region timescales (post-shock, compression x4):

- Thomson depth τ_es ~ 40 x (Mdot/1e-6).
- Radiative cooling time, assuming Λ ~ 1e-22 n^2: ~1e-2 s, but the region is optically thick.
- Photon diffusion time τ L/c ~ 80 s x (Mdot/1e-6).
- Flow time across the stream width, L/v_perp: ~1700 s.

**Verdict (ESTIMATE):**

- At Mdot <= 1e-6 the impact radiates its heat away about 20x faster than the flow crosses it. It is close to isothermal, consistent with the 300 kK O VI plasma.
- At Mdot = 1e-5 the two timescales are within a factor of 2.
- An adiabatic gamma = 5/3 run would therefore over-inflate the hot spot.

**Proposal:**

- Stage 1: an isothermal EOS (c_iso about 15-20 km/s, set by the accretor's Teff of 18 kK) or gamma = 5/3 with a fast thermal relaxation to T_eq. Radiation transport is not needed in stage 1.
- Later: M1 radiation or FLD, for the Mdot = 1e-5 end and for luminosity predictions.

**Consequence:** isothermal or adiabatic hydro without cooling is scale-free in density. Mdot then only sets the ratio of the stream density to the floor and ambient density. Results scale as torque / Mdot, which matters because Mdot is unknown (section 1).

### Stream inflow BC (outer x1 ghosts)

Inside the window |phi - phi_s| < 3 sigma_phi:

- rho = rho_s exp(-(phi - phi_s)^2 / (2 sigma_phi^2)).
- The values below are for r_out = 0.8 d_L1, all from the ballistic integration:
  - phi_s = 5.9 deg;
  - sigma_phi = (c_s/Omega)/r_out = 0.056 rad;
  - v_r = -96 km/s, v_phi = +50 km/s (rotating frame);
  - T = 6250 K donor, or isothermal c_s.
- Normalisation: rho_s = Mdot / (v_in · 2 pi sigma_w h). Here sigma_w = c_s/Omega, and h is the vertical stream thickness, about c_s/Omega in 3-D. The thin band has no real vertical extent, so state h explicitly.

Elsewhere on the face: outflow with no inflow (v_r >= 0, copy).

Ambient medium and floor:

- Ambient: rho_amb = 1e-5 rho_s, isothermal, at rest in the rotating frame.
- Floor: dfloor = 1e-6 rho_s.

### Time step and cost (ESTIMATE)

Time step:

- dt = 0.3 dr_min / (|v| + c_s), with dr_min = 8.7e8 cm and |v| + c_s ~ 680 km/s.
- That gives dt = 3.8 s and **1.55e5 steps per orbit** (P = 5.93e5 s).
- Work: 3.7e6 cells x 1.55e5 = 5.8e11 cell updates per orbit.

GPU time per orbit, assuming 3e8 to 1e9 zone-cycles/s per MI300A for plain sp hydro (an assumption, not measured; the smoke run measures it):

- 10-30 min per orbit on one GPU;
- 10 orbits: about 2-5 h on one apu node (2 GPUs).

Stage-2 envelope with dr = H_p/4: about 15x more steps.

### Measured quantities

Mass budget:

- Mdot_in through the stream window, Mdot_acc through r_in, and Mdot_out through the rest of r_out.

Angular momentum budget:

- Angular-momentum flux through r_in, in the inertial frame, about the accretor centre:
  - Jdot = ∮ rho v_r (r v_phi,in) r dphi (· band area).
  - Plus the viscous or no-slip surface stress. With a no-slip wall this appears as the momentum flux into the ghost cells, so measure the x3-momentum flux at r_in face by face from the Riemann fluxes.
- Specific angular momentum accreted, j_acc = Jdot / Mdot_acc, compared with:
  - j_Kep(R_P) = 1.5e19 cm^2/s;
  - j_* = R_P v_eq = 6.1e18 cm^2/s (ESTIMATE).
- Spin-up time, from J_* / Jdot.

Disc versus direct impact:

- the phi-distribution of Mdot and Jdot at r_in, split into the impact window around phi ~ 75 deg and the rest;
- the azimuthally averaged Sigma(r) and v_phi(r) / v_Kep, i.e. whether a ring or disc forms;
- elongation along the line of centres, compared with Barai+04.

Diagnostics go through `user_hist_func`. Two-dimensional r-phi maps go to the binary output.

### Staged plan (each stage gated by an apudev smoke at nlim about 10, same binary, input and keys)

0. **Build.** Incremental user pgen. CPU login-node correctness checks:
   - (a) Roche potential and Coriolis only, no stream, uniform gas initially at rest in the rotating frame. It should stay at rest to round-off where the potential is flat. Check the hydrostatic part with etotgrav.
   - (b) A ballistic test particle, i.e. a cold, low-c_s stream with no star, reproduces r_min = 2.69 Rsun from the script.
1. **Smoke:** apudev, nlim 10. This also measures zone-cycles/s, which replaces the cost ESTIMATE.
2. **Steady stream, non-rotating star** (Omega_* = Omega_orb, synchronous, i.e. v_phi = 0 at r_in in the rotating frame). Run 5-10 orbits until Mdot_acc and Jdot are steady. Gives the efficiency j_acc / j_Kep for a synchronous star.
3. **Same with Omega_* = 7.2 Omega_orb** (no-slip). Does the stream still spin the star up at the observed rate, or has it reached the torque-free state? This is the RY Per science result.
4. Resolution check: n_phi 1024 / 2048 / 4096. Then the u Her contrast (deep impact, synchronous).
5. Later: the thin envelope (stage-2 BC), radiation, and 3-D (a theta range of ±0.3 rad with the stream's vertical structure).

---

## 4. Open decisions for the user

1. **Mdot**: no measurement exists. Run scale-free (isothermal or no cooling, results per unit Mdot), or fix 1e-6 (or the model value 1.4e-5) once cooling or radiation enters?
2. **EOS for stage 1**: isothermal (simplest; the ESTIMATE says the impact is near-isothermal for Mdot <= 1e-6) or gamma = 5/3 with a thermal-relaxation source?
3. **Inner boundary**:
   - absorbing rigid surface, no-slip with Omega_*, versus free-slip;
   - or a resolved thin envelope with WB, which costs about 15x the time steps.
4. **Outer radius and stream injection**:
   - ballistic window at 0.8 d_L1 (as briefed), or 0.9 d_L1;
   - alternatively, put r_out just beyond L1 and hold a donor-atmosphere state near phi = 0, so the L1 nozzle forms by itself. This is closer to Lubow & Shu, but the circle then cuts the donor lobe.
5. **Theta band**: reflecting or periodic theta faces, and the in-plane (theta-independent) potential. Alternatively, 2-D Cartesian x-y with a staircased star.
6. **Masses**: Barai+04 spectroscopic (6.24 / 1.69 Msun, a = 30.3) as used here, or Olson & Plavec (6.25 / 1.60 Msun)? The difference is about 5 % in q.
