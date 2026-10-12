# NOTE shd-1012: why the spherical shadow (T7) looks weird, and the vet_gd twin

Worker SH, 2026-10-12, Raven. Work dir /raven/ptmp/jinma/shd_1012. No code changes, nothing committed.
- **Binary:** prodfinal d3bf42d7 (/raven/ptmp/jinma/prodfinal_1011/bin/athena_cpu_d3bf42d7). FIN deck = suiteF_1012 sp runs_js/shd/HFfZ_32 (copied to base.athinput).
- **Runs:** runs/<name>, with cmd.txt (the exact keys) and bin_md5.txt in each.
- **Analysis:** ana/shd.py (metrics -> ana/<run>.json/.txt), ana/plot.py (figure), ana/contrast.py (twin shell-mean diagnostic).
- **Figure:** fig/07_shd_sp.png. The page's fig is NOT overwritten.

## 0. Setup facts (the question assumed some of these differently)

- **Grid is nx1 x nx2 x nx3 = 32 x 48 x 48, not 48 x 48 x 32.** r covers 1-6 uniformly (dr = 0.156). theta covers 70-110 deg (40 deg). phi covers 0-60 deg.
- **Clump (r = 3, a = 0.4) is about 5.1 radial x 18 theta x 12 phi cells.** Cell aspect at r = 3 is dr : r dtheta = 3.6 : 1, so the radial direction is the coarse one. The blocky edge in the maps is these flat-shaded cells.
- **Lateral BCs are PERIODIC in theta and phi, not reflecting** (`ix2/ox2/ix3/ox3_bc = periodic`, sp_x2_periodic_image).
- **Periodic images do not change the exact answer here.** The phi period (60 deg) is an exact rotation symmetry. The theta images sit at +-40 deg. Their penumbrae reach only down to |lat| >= 23 deg (the penumbra half-width at r = 6 is about 17 deg), so the infinite-domain exact solution is the right reference in the 40-deg wedge.
- **vet_gd wall handling:**
  - Wall intensities are re-indexed by the local frame and read from the nearest set direction (rad_m1_vetgd.cpp:32-37).
  - Near-tangent reads beyond the band are clamped. The log shows "band 48 cells (needed 79)": about 3e3 clamped reads per build in the base wedge, none in the wide wedge.
  - The wide-wedge-without-twin control (section 2) shows that neither the walls nor the clamping causes the features.
- **The pile-up is gone:** front E/E_free = 0.92 / 1.00 (mean / max) in every d3bf42d7 run.
- **Equivalence:**
  - d3bf42d7 = JS's j1 for T7: every metric is identical to runs_js.
  - Old 78b2c961 differs only in the front (1.56 / 1.84).

## 1. Run matrix and metrics (exact = visible-star solid angle, s3_eval.shd_exact)

**Metric definitions:**
- **lit L1:** mean |E/E_free - 1| over lit cells (exact f > 0.98) within 0.36 rad of the clump axis, on the shells r = 4.05 / 4.99 / 5.77.
- **max3d:** max E/E_free anywhere outside the clump. Exact is <= 1 everywhere.
- **axis:** E/E_free on the star-clump line (mean of the 2 cells nearest lat 0). Exact is 0.06 at r 4.99 and 0.31 at r 5.77.
- **umbra:** mean E/E_free at r 4.05 where the exact value is < 0.02.

| run | what | lit L1 4.05 / 4.99 / 5.77 | lit max | max3d | axis 4.99 / 5.77 | umbra |
|---|---|---|---|---|---|---|
| old120 (78b2c961, page fig) | FIN 120 | 0.329 / 0.482 / 0.506 | 2.20 / 2.97 / 3.59 | 5.20 | 0.83 / 1.06 | 0.144 |
| **b120** | FIN 120 | 0.329 / 0.482 / 0.506 | 2.20 / 2.97 / 3.59 | 5.21 | **0.80 / 1.02** | 0.138 |
| **b432** | FIN, 18x24 = 432 dirs | 0.077 / 0.085 / 0.187 | 1.20 / 1.26 / 1.59 | 1.65 | 0.14 / 0.26 | 0.038 |
| r0_120 | rotation off (rotate_dt = 0, rotate_every = 0) | 0.379 / 0.548 / 0.592 | 2.31 / 3.91 / 5.20 | 5.31 | 0.03 / 0.07 | 0.035 |
| r0_432 | rotation off, 432 | 0.076 / 0.085 / 0.179 | 1.20 / 1.26 / 1.56 | 1.65 | 0.14 / 0.26 | 0.039 |
| snone_120 | vet_gd_symref = none | = b120 (to 3 digits) | | 5.21 | 0.80 / 1.02 | 0.138 |
| hr120 | 64x96x96, 120, t = 299 (rotation block 1) | 0.355 / 0.478 / 0.466 | 2.79 / 5.60 / 6.35 | 6.44 | 0.07 / 0.13 | 0.080 |
| **hr120b** | 64x96x96, 120, nlim 80, t = 598 (same rotation block as b120) | 0.347 / 0.519 / 0.557 | 2.87 / 3.77 / 5.07 | 9.10 | **1.19 / 1.60** | 0.138 |
| **hr432** | 64x96x96, 432 | 0.077 / 0.114 / 0.227 | 1.20 / 1.39 / 1.74 | 1.78 | 0.04 / 0.18 | 0.016 |
| **wd120** | wide wedge 80 x 120 deg, 96x96 lateral (same cells), 120 | 0.112 / 0.132 / 0.127 | 1.22 / 1.35 / 1.41 | 1.99 | 0.32 / 0.43 | 0.068 |
| **wd432** | wide wedge, 432 | **0.038 / 0.023 / 0.035** | 1.10 / 1.05 / 1.09 | **1.11** | 0.12 / 0.24 | 0.037 |
| **tw0_120** | twin off (vet_gd_twin, _fuse = false, symref none), 120 | 0.179 / 0.201 / 0.297 | 1.34 / 1.49 / 1.63 | 1.95 | 0.52 / 0.65 | 0.107 |
| tw0_wd120 | twin off, wide wedge, 120 | 0.233 / 0.202 / 0.417 | 1.55 / 1.49 / 1.69 | 2.00 | 0.20 / 0.17 | 0.059 |
| tw0_432 | twin off, 432 | 0.074 / 0.095 / 0.142 | 1.21 / 1.28 / 1.44 | 1.44 | 0.15 / 0.29 | 0.038 |
| nc120 | no clump (clump_rho = background), 120 | 0.000 everywhere | 1.000 | 1.000 | (1.00 = free) | - |

## 2. Verdict per feature

### 1. Blocky clump edge: TEST (plotting and resolution)
- **Cause:** the clump spans 5 radial cells. Flat shading of a 3.6 : 1 cell shows the staircase. The exact panel shows the same staircase in its masked interior.
- **Resolution:** at 64x96x96 the edge is half as coarse. The lit metrics do not improve: lit L1 at 432 is 0.077 / 0.114 / 0.227 vs 0.077 / 0.085 / 0.187.
- **Conclusion:** the "weird" look is not spatial resolution.

### 2. Lobes E/E_free 1.4-2.2 (up to 3.6 at r 5.77; 5.2 anywhere) near the lateral walls: METHOD (the vet_gd twin), amplifying a 120-direction ray effect. Not the BC, not the resolution.
- **No clump: exact.** At 120 directions without the clump, E/E_free = 1.000 everywhere. The twin removes the set's ray noise exactly on a laterally uniform state.
- **Twin off is better.** With the clump, at 120 directions:
  - lit L1 drops from 0.33 / 0.48 / 0.51 to 0.18 / 0.20 / 0.30.
  - max3d drops from 5.2 to 1.95.
  - The twin makes the clump case 2-3x worse.
- **A wider wedge cures the twin, not the walls.**
  - Twin on: the wide wedge brings lit L1 to 0.11 / 0.13 / 0.13.
  - Twin off: the wide wedge does not help (0.23 / 0.20 / 0.42 vs 0.18 / 0.20 / 0.30).
  - So the wedge width matters only through the twin. It does not act through the periodic walls or the clamped band reads: tw0_wd120 has no clamping, yet it is no better than tw0_120.
- **Mechanism:** the twin's laterally uniform state is the arithmetic SHELL MEAN of chi and S (sect. 3). In the base wedge the clump covers 7.6 % of the shell area (1.9 % in the wide wedge). The clump shells then get <chi> dr = 12 / 11 / 5.7 against the true mean transmission -ln<exp(-chi dr)> = 0.08 / 0.07 / 0.04, a factor about 150 per shell (ana/contrast_t7.txt).
  - The twin's column optical depth through the clump shells is about 40 (wide: 10), while the real mean transmission is 0.19 (wide: 0.05).
  - So the twin sweeps a star inside an OPAQUE shell. Outside r = 3.4 its field is not the star's ray pattern.
  - Subtracting "twin - shell mean" (LAT slots) and rescaling the half-range J by <J_twin>/J_twin (rad_m1_vetgd.cpp:6459) then injects a noise pattern that does not belong to the main field.
  - The 4x thinner twin shell in the wide wedge keeps more of the star's pattern, so the twin helps there (0.11-0.13 vs 0.20-0.42 without it).
- **What remains is a ray effect.** Without the twin, 120 directions still give lit max 1.3-1.6 and 432 give 1.2-1.4: the usual ray effect of a 0.1 sr direction weight against a 0.09-0.13 sr star (test angular resolution).
- **Best case:** wd432 is close to exact everywhere: lit L1 <= 0.04, max 1.11.

### 3. The axis spike behind the clump (0.80 at r 4.99, 1.02 at r 5.77; exact 0.06 / 0.31): ray effect of the ONE rotated 120 set that is active at the end, amplified by the twin
- **Rotation off removes it:** 0.03 / 0.07.
- **Same orientation at 2x resolution makes it larger:** hr120b, the same rotation block (alpha = 21.2 deg), gives 1.19 / 1.60. With another orientation (hr120, alpha = 55.6 deg) it is gone (0.07 / 0.13).
- **Not vet_col and not the symref col.**
  - vet_gd_symref = none changes nothing (0.80 / 1.02).
  - vet_col does not depend on the set rotation, yet the spike follows the rotation.
- **432 directions:** 0.14 / 0.26 (exact 0.06 / 0.31).
- **Twin off:** it drops to 0.52 / 0.65.
- **Consequence for the page:** "last set wins" makes such single-set artefacts visible in a steady test. AV's rotation averaging (rotavg-1012) is the relevant general cure. Not duplicated here.
- **Not a polar-axis effect.** The line is in the equatorial plane. The product set's z axis is the coordinate polar axis, so the line lies between the n_z = +-0.149 rings, 8.6 deg off the plane.

### 4. Penumbra and umbra noise: angular resolution, plus the twin at 120
- **432 directions:** penumbra L1 0.06-0.17, umbra 0.04 (hr432: 0.016).
- **120 directions:** penumbra L1 0.29-0.50, umbra 0.14.
- **Wide wedge at 432:** penumbra L1 0.057 / 0.061 / 0.071.
- **Breakdown:** the excess at 120 is the twin (sect. 2) plus the set's ray effects. Spatial resolution does not reduce it.

### Test setup recommendation (no test-specific fix)
- **Keep the test, with a wider wedge.** The 40 x 60 deg periodic wedge is a valid test of the method; the exact answer is the same as in the infinite domain. The same test in an 80 x 120 deg wedge, or with the clump covering a few % of the shell, is the physically fairer check of a shell-mean method.
- **Keep the 40-deg wedge as a stress test of the twin.** It is the case where the shell mean is not representative.
- **Use 432 directions for the page,** or 120 once rotation averaging exists. With 120 and "last set wins", single-set ray effects dominate.

## 3. What the twin does (code as of d3bf42d7, wt-prodfinal/src)

**Defaults (rad_m1.cpp):**
- **twin:** vet_gd_twin is ON by default on a fresh vet_gd run (rad_m1.cpp:1149-1158). twin_fuse, twin_det and vet_source_noesrc follow it (defaults-1009).
- **symref:** vet_gd_symref = col is the fresh default where the twin is on (rad_m1.cpp:1172-1205).
- **twin_lowmem:** read only when named (rad_m1.cpp:1229-1236).

**When it is used:**
- **Every vet_gd build.**
  - Stage 0 runs before the real sweep's moments (rad_m1_vetgd.cpp:4941). When fused, it rides in the same sweep.
  - Stage 1 runs after the moments (rad_m1_vetgd.cpp:4969).
- **Only where vet_gd is active:** shells i >= ilo = is + vlat_icut (rad_m1_vetgd.cpp:5908, 6250).
  - vlat_icut is the global-minimum index of the first shell, counted from the top, whose column depth tau_top + dtau/2 reaches vet_col_lat_taucut (30). It is computed in VetLatCsFinish, rad_m1_vetlat.cpp:507-531.
  - Below the cut, D comes from vet_col alone (the radial column f_K): no sweep and no twin.
  - Above the cut, slot 0 holds D_rr from the sweep (LAT0 = D_rr - f_K) and LAT1..5 the off-diagonal and tangential terms (rad_m1_vetgd.cpp:45-55).
  - In production that region is the photosphere and the atmosphere above tau about 30. In T7 nothing reaches tau 30, so icut = 0 and the whole wedge is twin-active.

**Stage 0 (VetGdTwin, rad_m1_vetgd.cpp:5888-6060; fused at 6027-6058):**
- **Shell means:** the band's ln chi and ln S are replaced by ln of the shell's ARITHMETIC means <chi>, <S> (shell_mean with ex = true, rad_m1_vetgd.cpp:6038-6043, 6066).
- **Twin sweep:** one sweep with the same directions and angle gives the twin's LAT slots (vgd_twl) and the twin's half-range ratios, C1 and J (vgd_hrt).
- **Symref:** with symref col, a 1-D p-ray solution of the same shell means is solved on the host (VetGdSymRef, :5761-5850).

**Stage 1 (rad_m1_vetgd.cpp:6095-6123):**
- **Tensor:** LAT0 -= (twin LAT0 - ref), with ref = <twin LAT0> (or the symref f, weighted by VgdSrW = min(J_twin/J, 1), :104-115). LAT1..5 -= twin LAT1..5.
- **Half-range moments:** VetGdHalfRange(2) (:6228-6475) does the same for the half-range h+-:
  - r -= (r_twin - <r_twin>), then clamped.
  - C1 -= its pattern.
  - J *= <J_twin>/J_twin (:6459). This is the factor that becomes noise when J_twin is an attenuated, unrepresentative field.

**Why a shell-mean reference:**
- The exact laterally uniform D is diagonal and laterally constant. The twin's lateral pattern on that state is therefore the set's pure sampling noise, and subtracting it makes uniform states exact (nc120: 1.000).
- **Assumption:** this is right only if the main field's ray pattern equals the twin's, i.e. the shell is laterally near-uniform in transmission.

**Cost and memory:**
- **CPU (T7, this work):** the fused twin is 17 % of the run time at 120 directions (b120: 19.3 of 112.8 s; wd120: 27.4 of 161 s) and 24 % at 432 (b432: 67.4 of 284 s).
- **GPU (bsgaccel notes):** fusing saves 15 % per cycle vs the unfused twin. The twin's own intensity array is the memory item: twin_lowmem = false is -5.6 % time where it fits (bsgaccel C2/lowmem rows).
- **Twin off saves the twin sweep and its intensity array.** Picard went 2.94 -> 3.63 at 120 (T7), 3.01 -> 3.15 at 432.

## 4. Production exposure (ana/contrast.py; AG Car A = rampA_1011 full-size 3-D dumps, 480 x 128 x 128)

**Method:** chi = kappa_R(rho, T) rho from the dump (EOS and TOPS table as rampA reduceA.py). The cut is recomputed as in VetLatCsFinish.

| snapshot | icut r [cm] | max per-shell <chi>dr / (-ln<e^{-chi dr}>) | max tau_twin / tau_eff (column, tau_eff > 0.01) | tau = 1 radius twin / eff / median column | max sigma(B)/<B> |
|---|---|---|---|---|---|
| T7 base 40 x 60 | (all) | **150** | **~500** | twin 3.27 vs never | 0 |
| T7 wide 80 x 120 | (all) | 155 (thinner shell) | ~520 | twin 3.11 vs never | 0 |
| AG Car A t = 2.0e6 s | 2.643e13 | 1.04 | 1.2 | 2.708 / 2.705 / 2.708 e13 | 0.37 |
| AG Car A t = 4.0e6 s | 2.769e13 | 1.14 | 2.06 | 2.919 / 2.910 / 2.913 e13 | 0.69 |
| AG Car A t = 5.2e6 s | 2.913e13 | 1.39 | 2.05 | 3.031 / 3.022 / 3.026 e13 | 0.87 |

**AG Car A, opacity:**
- Per shell, the twin's opacity is within 4-39 % of the true mean transmission. T7 is off by a factor 150.
- The twin's photosphere is 0.1-0.3 % above the effective one.
- The column tau ratio peaks at about 2 in the tau 1-30 layer just below the photosphere, e.g. tau_twin 11.4 vs tau_eff 5.8. That is where columns of very different depth coexist (chi_std/mean up to 2-5 in the H recombination zone).

**AG Car A, source contrast:** the S (B = a T^4) contrast is large, sigma/mean up to 0.87 near the photosphere. T7 has no S contrast, so T7 does not test this part.

**BSG:**
- The unseeded R (bsgrec_1011) stays 1-D (lateral KE <= 4e-16), so the twin is exact there.
- No seeded run with the final method exists yet (bsgprod P is prepared, not run).
- As a structural proxy I used the seeded (1 % T) 256 x 256 BSG repro of 10-06, /raven/ptmp/jinma/bsg_raven_1002/repro_sf_cont (older RT method, convecting, lateral v_rms 14 km/s). Ideal EOS T, Ma+2026 Rosseland table.

| snapshot | icut r [cm] | max per-shell ratio | max column tau_twin / tau_eff | tau = 1 radius twin / eff / median | max sigma(B)/<B> |
|---|---|---|---|---|---|
| BSG seeded t = 4.54e6 s | 3.355e12 | 1.9 | 5.2 | 3.568 / 3.547 / 3.558 e12 | 0.85 |
| BSG seeded t = 4.96e6 s | 3.366e12 | 1.75 | 5.15 | 3.590 / 3.558 / 3.568 e12 | 0.76 |

- **Thick layer:** in the tau 5-150 layer the twin's column depth follows the MEDIAN column (e.g. r 3.45e12: twin 62, median 59), while the mean transmission is set by the thinnest columns (tau_eff 13.6). Ray noise is negligible there (diffusive field).
- **Photosphere:** the twin's tau = 1 radius is 0.6-0.9 % above the effective one (0.3 % above the median).
- **Exposure ranking:** about 2-3x that of AG Car A, and 2 orders of magnitude below T7.

**AG Car A A/B, reduced (worker AC, via coordinator):**
- Setup: E0r32 restart at t = 3.0e6 s, 32 x 32 lateral, convective KE 1.1e44 erg, 240 cycles, Viper apudev, ed32f0fa. OFF = twin false + symref none + noesrc false + fuse/det false.
- **L_top: max relative difference 3.1e-5. L_mid: 1e-6. Picard 7.03 vs 7.01.**
- Files: /viper/ptmp2/jinma/agcmodel_1012/twin_{ON,OFF}.

**Why AG Car A is insensitive:** the twin IS active there (icut sits at tau about 30 below the photosphere: r 2.64-2.91e13, R_ph 2.7-3.0e13). But its shell-mean state is representative:
- per-shell <chi>dr is within 1.04-1.39 of the mean transmission;
- the twin photosphere is within 0.3 %;
- in the twin-active layer, the chi contrast is in the column-depth distribution, not in a few opaque cells covering the shell (T7).
- The S contrast (sigma(B)/<B> up to 0.87) is large. Its effect on the result is nil in AC's A/B: the twin only subtracts the twin's own ray pattern, which on a near-uniform-transmission shell matches the main field's.

**Twin on/off A/B at full size (AG Car A, 480 x 128 x 128):**
- Setup: rampA rst10 (t = 5.0e6 s), binary 1fbe8621 (the rampA binary), Raven gpudev 1 node x 4 A100, about 11 min each, output dt 2e4 s.
  - ON = /raven/ptmp/jinma/shd_1012/agc/on (job 31109912).
  - OFF = agc/off (job 31109728): vet_gd_twin, _fuse = false, vet_source_noesrc = false (this binary predates symref).
  - Comparison: ana/agcab.py, ana/abstat.py.
- **Reach:** ON reached t = 5.415e6 s (cycle 4947), OFF t = 5.471e6 s (cycle 5005). Cost: 9.89e6 vs 1.05e7 zone-cycles/cpu-s, so the twin costs about 6 % per cycle at full size on A100.
- **Picard:** mean 14.1 vs 13.9 (max 28 both, 0 non-converged).
- **Luminosity over the common window (48 hst rows, 5.0-5.47e6 s):**
  - L_top: ON - OFF max 1.1 %, mean 0.11 %. The window means differ by -0.09 %.
  - L_mid: window means differ by -0.05 %. KE_int: -0.07 %.
- **Noise floor:** the ON rerun itself differs from the original rampA ON link (same binary, same restart) by up to 2.9 % in L_top (2.2 % L_mid, 3.0 % KE: GPU/convective chaos). The temporal std of L_top in the window is 2.3 %.
- **So the twin has no measurable effect in AG Car A at full resolution either.** Dumps: shell-mean E, F1 differ <= 0.3 %, at the outer boundary region r > 7e13.
- **AC's reduced A/B** (E0r32, 32 x 32, Viper, ed32f0fa) agrees: L_top 3.1e-5.

## 5. General fix candidates (proposal, nothing implemented)

The one rule must hold on all tests and production, both closures (vet_sc has no twin; this is vet_gd only), same keys.

**A. Twin off.**
- **Pros:** simplest and saves its cost.
- **T7:** 120: lit L1 0.18-0.30 (vs 0.33-0.51 on). 432: unchanged or better (0.07-0.14).
- **Cons:** it loses the exact uniform-state noise cancellation (nc120 1.000 -> the bare 120-set ray noise). Every laterally uniform suite case would show set noise again at 120: thin_L16, core, pulse, relax... Those cases are where the twin was introduced (vgdlat/defaults-1009 records).
- **Needs:** the sp suite rerun at 120 before it could become a default.

**B. A transmission-consistent shell state (smallest change, local in code).**
- **Rule:** keep the twin, but build its uniform state from the shell's mean TRANSMISSION, not its mean opacity: chi_twin = -ln<exp(-chi dr)>/dr. The source likewise is the transmission-weighted mean <S (1 - e^{-chi dr})>/<1 - e^{-chi dr}>.
- **Uniform shells:** identical to now (bitwise when chi is laterally constant).
- **T7:** the clump shell's twin dtau goes from 12 to 0.08. The twin then sees the star pattern, as in nc120.
- **AG Car:** <= 40 % per shell change, at the photosphere 0.3 %.
- **Cost:** two extra exp/log per cell per build. The same gather; no memory.
- **Caveat:** it still assumes one pattern per shell. In a channel (bright beam through an opaque shell) the twin then describes the bright part and mis-corrects the opaque part. That is the case VgdSrW already guards for the symref shift.

**C. A contrast-triggered local weight (complements B).**
- **Rule:** apply the twin's pattern subtraction with a per-cell weight w = min(J_twin/J, J/J_twin), or the VgdSrW form, on LAT0..5 and on the h+- / C1 / J subtraction (today only the symref shift is weighted). Where the main and twin fields disagree (shadow, beam, channel), the cell falls back to the bare sweep. Near-uniform shells keep the full cancellation.
- **Cost:** the J of both sweeps is already stored (vgd_hr/vgd_hrt M1_VGD_HRJ). No memory.
- **Risk:** w is itself built from noisy J. It needs a floor and a check on the uniform suite cases (should stay bitwise there, since w = 1 when J = J_twin).

**D. A local (not shell-mean) reference: per-cell or per-patch lateral smoothing of chi and S, swept as the twin.**
- **Most general,** but it needs either a per-patch twin (memory x patches) or a smoothed full field: the same memory as the twin, one more pass for the smoothing.
- **No clear gain over B + C** for the T7 / AG Car contrasts measured above.

**Recommendation:** B + C. A only if the 120-direction uniform cases are shown not to need the twin.
- **Validation:** sp suite (FIN keys), T7 base and wide at 120 / 432, nc120 (must stay 1.000), and the AG Car / BSG A/B. Same keys in all.

**Interaction with UN's cut change (vet_col_lat_cut_thin, unify-sc-gd-1012):**
- **The twin-active region gets deeper.** The cut never sits above a thin cell, so in a cavity or sshell the shells inside the opaque layer that contain the thin cavity now join the vet_gd/twin region.
- **Exactly the T7 failure.** Those shells are the most laterally inhomogeneous (thin cavity next to thick material). The arithmetic <chi> makes the twin opaque while the main field streams through the cavity.
- **So UN's change raises the twin's exposure.** B (transmission mean) or C (weight) should be in place, or checked on UN's sshell case with twin on/off.
- **The cut logic itself is untouched here.**

## 6. twinfix-1012: implementation and validation (B = vet_gd_twin_state = trans; C = vet_gd_twin_weight = j; A = twin off)

**Branch:** twinfix-1012 (local, not pushed), worktree /raven/u/jinma/ATHENAK/wt-twinfix, from prodfinal-1011 2b84ffc6.
- 2fc2edae: the keys, default mean/none.
- 182153fe: fresh default `vet_gd_twin_state = trans`.

**Binaries:**
- shd_1012/bin/athena_cpu_t1 / athena_he_cpu_t1 (2fc2edae, md5 9214f149 / 7ce6310c). All validation named the keys explicitly.
- athena_cpu_t2 (182153fe, fresh-default check).
- athena_cpu_r1 (twrot-test-1012 = 2fc2edae merged with AV's rotavg e584cfe1, test only).

**Code (all in src/rad_m1):**
- **Keys** (rad_m1.cpp, after the vet_gd_twin key): `vet_gd_twin_state = mean | trans`, `vet_gd_twin_weight = none | j`, `vet_gd_twin_wq` (0.5). They are read when vet_gd_twin is on; a non-default value without the twin is FATAL.
- **B:** `VetGdTwinTrans` (rad_m1_vetgd.cpp).
  - It is called in VetGdTwin stage 0 on both the fused and the unfused/lowmem path, after the shell means and before VetGdSymRef, so the symref uses the same state.
  - Per shell it sums t = e^{-chi dr}, w = 1 - t, S w, and min/max ln chi in a fixed order, rank-ordered over ranks (min/max by Allreduce).
  - Then chi_tw = -ln<t>/dr (log1p(-<w>) where <t> >= 0.5), and S_tw = <S w>/<w>.
  - A shell with min ln chi == max keeps the arithmetic means: bitwise on uniform shells.
- **C:** `VgdTwW` = min(1, min(J/J_tw, J_tw/J)/wq). It multiplies the twin subtraction of LAT0..5 (VetGdTwin stage 1) and of r+-, C1 (additive) and J ((mo/ht)^w) in VetGdHalfRange(2), on both the det-MPI and serial paths. Where w = 1 the old expressions are used.
- **Cost / memory of B:** one extra fixed-order pass over the band per build (2 exp + 1 expm1 per cell), plus 5*c1 Reals. No new big array.

**Key sets:** all new runs use FIN2 (stretch 0.02, slope correct, js_face hcface). The control is the same binary with mean/none.

### Bitwise status
- **FIN keys, mean/none (t1) vs d3bf42d7:** T7 E, F1 bitwise; F2/F3 differ at 1e-14 because the runs used different OMP thread counts (8 vs 14).
- **No-clump T7, trans + j vs control:** all of m1 bitwise. Uniform shells are untouched.
- **Laterally symmetric suite cases, B+C vs control:** E, F1, dens, eint, velx L1 and orders identical to the printed digits in every xord case (rsw, rswt, rswv, rswtv, marsh, rlat*, patmT). The patm, thin, core, dif, relax and pulse exact metrics are identical too.
  - Only the lateral noise components of the symmetric cases differ. These are F2/F3, vely/velz of relax1_T and rswtv, with L1/mean ~ 1, i.e. round-off noise relative to a zero mean: chi is not exactly laterally constant there, so trans is active at round-off.
- **Reduced AG Car A and bsgR, B+C vs control:**
  - AG Car: L_top/L_in 1.00022 (min 0.99965); Picard 3.234 vs 3.232; E/J-1 at R_ph -1.60e-4 in both.
  - BSG: L_top identical to 9 digits; Picard 2.817 both.

### T7 (spherical shadow), FIN2

| arm | dirs | lit L1 4.05/4.99/5.77 | lit max | max3d | penumbra L1 | umbra | axis 4.99/5.77 (ex 0.06/0.31) | Picard |
|---|---|---|---|---|---|---|---|---|
| control (mean) | 120 | 0.329/0.482/0.506 | 3.59 | 5.21 | 0.29-0.50 | 0.138 | 0.80/1.02 | 2.95 |
| **B trans** | 120 | **0.006/0.008/0.009** | 1.02 | 1.19 | 0.12-0.18 | 0.060 | 0.32/0.49 | 2.96 |
| C weight j | 120 | 0.181/0.214/0.301 | 1.77 | 1.97 | 0.16-0.45 | 0.112 | 0.55/0.69 | 5.20 |
| B+C | 120 | 0.006/0.008/0.009 | 1.02 | 1.20 | 0.12-0.18 | 0.068 | 0.35/0.50 | 4.39 |
| A twin off | 120 | 0.179/0.201/0.297 | 1.63 | 1.95 | 0.16-0.45 | 0.107 | 0.52/0.65 | 2.86 |
| control | 432 | 0.077/0.085/0.187 | 1.59 | 1.65 | 0.06-0.17 | 0.039 | 0.15/0.27 | 3.00 |
| **B trans** | 432 | **0.003/0.004/0.006** | 1.00 | 1.01 | 0.05-0.07 | 0.037 | 0.14/0.24 | 2.99 |
| B+C | 432 | 0.003/0.004/0.006 | 1.00 | 1.01 | 0.05-0.07 | 0.037 | 0.14/0.24 | 2.90 |
| A | 432 | 0.074/0.095/0.142 | 1.44 | 1.44 | 0.06-0.13 | 0.039 | 0.15/0.29 | 3.15 |
| wide control | 120 | 0.112/0.132/0.127 | 1.41 | 1.98 | 0.12-0.21 | 0.069 | 0.33/0.44 | 3.12 |
| wide B+C | 120 | 0.004/0.004/0.005 | 1.03 | 1.47 | 0.13-0.18 | 0.067 | 0.25/0.31 | 2.72 |
| wide A | 120 | 0.233/0.202/0.417 | 1.69 | 2.00 | 0.17-0.53 | 0.059 | 0.20/0.17 | 2.62 |
| wide control | 432 | 0.038/0.023/0.035 | 1.10 | 1.11 | 0.06-0.07 | 0.037 | 0.12/0.24 | 3.01 |
| wide B+C | 432 | 0.002/0.004/0.006 | 1.00 | 1.00 | 0.06-0.07 | 0.037 | 0.13/0.24 | 3.00 |
| wide A | 432 | 0.085/0.094/0.147 | 1.41 | 1.42 | 0.07-0.14 | 0.039 | 0.15/0.29 | 3.14 |

- **B alone does the job:**
  - Base and wide wedges now agree (domain independence restored).
  - The 120-direction >1 lobes are gone (lit max 1.02).
  - The remaining 120-direction error is penumbra and axis ray noise of one set.
- **C** gives no accuracy beyond B and raises Picard (2.95 -> 4.4-5.2 at 120). It stays off.
- **Separating the lobes.** The twin-off lobes (A: lit max 1.4-1.7, L1 0.18-0.42) are the bare ray noise of the set on the star field. A representative twin cancels them, so they are what the twin buys. The control's lobes (up to 3.6) were that ray noise plus the opaque-twin error.
- **Figure:** fig/07_shd_sp_twinfix.png (control vs trans, 120 and 432).

### What the twin buys (A = twin off + symref none + noesrc false), FIN2
- **sp suite, P1 patm atmosphere (patmT_L16_s):**
  - tau < 0.1 E error at N = 256: 4.5e-5 (C) vs 7.9e-3 (A), 175x worse; it does not converge with A (order -0.5).
  - E/J-1 at R_ph: -2.9e-4 vs +2.2e-3.
  - Picard: 1.9 vs 2.9.
- **Other symmetric cases:**
  - marsh_r1 F1 L1 is 35 % worse with A.
  - patmT F1 order drops from 2.00 to 1.42.
  - core_u is mixed: A has the better F r^2 (1.002 vs 0.977), the worse E (+3.4e-2 vs -2.2e-2).
- **AG Car A reduced (300 cycles):**
  - L_top/L_in min/max: 0.99286 / 1.00141 (A) vs 0.99965 / 1.00022 (C).
  - E/J-1 above 1.1 R_ph: max 3.95e-2 vs 3.1e-3.
  - Picard: 4.10 (max 18) vs 3.23 (max 6). Wall time is about equal (1633 vs 1651 s).
- **bsgR reduced (30 cycles):** L_top differs by 7e-6; Picard 2.98 vs 2.82; cost -11 % (1242 vs 1398 s CPU).
- **AG Car full size (A100, sect. 4):** twin cost about 6 % per cycle; no measurable L_top effect.
- **Conclusion:** A loses the optically thin accuracy and the Picard robustness of the uniform-state noise cancellation, so A is rejected. B keeps all of it and fixes the shell-mean failure.

### Cost
- **B vs control:**
  - T7 CPU: Picard unchanged.
  - AG Car reduced: wall 1657 vs 1651 s (B+C). bsgR: 1413 vs 1398 s CPU (B+C, +1 %).
  - The B pass alone is a small part of that; B-only timings are in tw12.
  - Memory +5*c1 Reals.

### P2 wide channel (pinc_ha20_sr30, 64; metric p2eval M1 L / L_ex, stagger true / false; FIN arm strings as AV)
Binary r1 = twinfix 2fc2edae + AV rotavg e584cfe1 (test merge twrot-test-1012, not for integration). Runs in shd_1012/runs_p2.

| arm | stagger true | stagger false | Picard |
|---|---|---|---|
| R0, control (mean) | 1.096 | 0.775 | 2.18 / 2.09 (AV: 1.104 / 0.777) |
| R0, B trans | 0.994 | 0.897 | 2.97 / 2.18 |
| R0, C weight j | 0.994 | 0.993 | 3.14 / 2.26 |
| R0, B+C | 0.995 | 0.981 | 3.10 / 2.18 |
| rotate_avg 4, control | 1.071 | 1.306 | 2.11 / 2.12 (AV: 1.07 / 1.31) |
| **rotate_avg 4, B+C** | **0.998** | **1.000** | 2.13 / 2.14 |
| rotate_avg 4, twin off (AV) | 1.022 | 1.021 | |
| rotate_avg 4, B only | pending (tw12, job 31110760) | | |

- **B makes the twin average-able in P2.** With B+C, rotate_avg 4 is within 0.2 % for both sets, against 7 % / 31 % for the mean twin.
- **The C weight does matter in a channel.** Where the beam and the twin differ locally, B alone with a single set (R0, stagger false) is 0.897 and B+C is 0.981.

### T7 + vet_rotate_avg 4 (r1 binary, FIN2, 120 directions)

| arm | lit L1 4.05/4.99/5.77 | max3d | penumbra L1 | umbra | axis 4.99/5.77 (ex 0.06/0.31) | Picard |
|---|---|---|---|---|---|---|
| control (mean) | 0.094/0.175/0.223 | 2.18 | 0.12-0.27 | 0.062 | 0.24/0.29 | 2.26 |
| **B+C** | **0.004/0.004/0.005** | 1.09 | 0.08-0.12 | 0.044 | 0.16/0.18 | 2.26 |
| A (twin off) | 0.093/0.089/0.081 | 1.24 | 0.08-0.15 | 0.079 | 0.38/0.59 | 2.27 |

- **Rotation averaging removes the single-set axis spike and halves the 120-direction penumbra error.** It does so on top of B.
- **With averaging, C costs no Picard here.**

### Fresh default and restart
- **Fresh default:** 182153fe makes `vet_gd_twin_state` default `trans` on a FRESH vet_gd run with the twin. Checked: binary t2 on the T7 FIN deck without the key reproduces the trans results (lit L1 0.006/0.008/0.009, Picard 2.94).
- **Restart rule:** a restart whose file lacks the key keeps `mean` (bitwise older runs). This is by construction, not run-tested.
- **Weight:** `vet_gd_twin_weight` stays `none` by default.

### Not done / pending (scope freeze)
- **tw12 (B-only arm), pending in the queue:** suite, reduced AG Car / bsgR, P2 + rotate_avg 4, T7 + rotate_avg 4.
- **UN's sshell deck with/without UN's cut change:** not run (no hash received).
- **future:** a production-scale C weight with a smaller Picard cost. Candidates are the w floor/q and applying the weight only to LAT0 and h+-, not to J.
