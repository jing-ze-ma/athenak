# NOTE thinfix (thinfix-1011): the thin top of the atmosphere and the collimated-beam energy of the implicit VET

Branch `thinfix-1011` (local, NOT pushed) in the worktree /raven/u/jinma/ATHENAK/wt-thinfix, on thinsw-1010 61cab6c0.
Work dir: /raven/ptmp/jinma/thinfix_1011.

**Commits**
- ea3ba71e: `implicit_thin_bq = lin` and `implicit_thin_recon = js`.
- cfb3258f: `implicit_thin_recon = jse`, which is js weighted by the thermal fraction. This is the final commit.

**Binaries**
- bin/athena_cpu_cfb3258f, md5 0897daaf... It is identical to athena_cpu_jse1, so every jse1 run is valid for cfb3258f.
- bin/athena_cpu_ea3ba71e, md5 2fb64cf3... It is identical to athena_cpu_js1.
- bin/athena_he_a100_cfb3258f, md5 5c43fec7...

## Verdict

**Final variant (arm `tJe`).** It is tEs plus two new keys:
```
implicit_flux = blend   implicit_flux_faces = all   implicit_flux_beam = halfrange   implicit_hr_recon = dc
implicit_blend_xthin = 0   implicit_thin_switch = fse   implicit_thin_corr = hcs
implicit_thin_recon = jse   implicit_thin_bq = lin                 (new; implicit_thin_jslim = 2 default)
+ hesdirk2 2nd order in time: time2_vet_sc = predict, implicit_vimp = true, time2_vstage = true
```

| criterion | tEs (thinsw best) | **tJe (final)** | verdict |
|---|---|---|---|
| atm thin top, tau < 0.1, self-conv. orders 32..2048 | 0.74 0.60 1.52 0.95 0.62 | **2.04 2.02 2.09 2.01 1.99** | **FIXED** |
| atm tau 0.1-1 / 1-10 | 2.26 2.41 1.82 **0.58 0.13** / 2.0-2.4 | **1.92 2.03 2.01 2.00 2.00** / 2.00 | **FIXED** |
| atm tau > 10 (thick bottom) | 2.01 2.05 1.78 1.50 1.48 | 2.01 2.04 1.77 1.49 1.47 | = cen (1.57 1.38), not this fix (see D3) |
| Hopf L1 at 2048 (angular-quadrature floor) | 3.454e-4 | 3.455e-4 | = |
| 1-D X/T/C ladders (pulses k0.128/k12.8/k1280, rw_t10, rw_t1000; T with predict) | | **identical to tEs, digit for digit** (122 lines) | PASS |
| rw_t10 X, last order (E) | 1.70 (cen 1.96) | 1.70 raw; **2.02 / 1.94 with the mean removed (cen 2.03 / 1.95)** | measurement artifact (D2) |
| gates G1 cfl 1e2 / 1e4 | L1 7.08e-4; 1e4 NON-CONV 34 | L1 7.08e-4; 1e4 NON-CONV **16** | PASS |
| gate G3 c1 / c10 | 0.0143 / 0.0148 | **0.0139 / 0.0142** | PASS |
| gate G5 k1280 / k12800 (k12.8 / k128 fail for all arms) | 1.0167 / 0.99996 (0.508 / 0.534) | same (0.508 / 0.534) | PASS |
| beams: stable 5/5 | yes | yes, 0 NON-CONVERGED | PASS |
| beam E/J_exact median, xb20 / ba0 / ba20 (exact 1) | 0.209 / 0.683 / 0.110 | **0.754 / 1.083 / 0.485** | improved 1.6-4.4x; limited by the formal solution (D4) |
| beam luminosity F/H_exact at x = 0.5 | 0.20 / 0.80 / 0.135 | 0.67 / 1.01 / 0.62 | |
| cyl <\|E/J-1\|> / <\|F\|/\|cH\|> | 0.346 / 0.707 | 0.453 / 0.592 | worse (D4) |
| shd3b edge at x = 0.65 / 0.85 (exact 0.184 / 0.295) | 0.236 / 0.343 | **0.195 / 0.316** | better |
| beam cost, Picard / BiCGStab inner (xb20, ba0, ba20, cyl, shd3b) | 5.49/12.2, 4.00/18.0, 5.51/16.0, 4.99/10.0, 2.02/6.4 | 4.86/13.0, 4.19/17.1, 4.88/15.6, 4.02/12.4, **3.00/16.6** | OK; shd3b is +30 % wall (355 s vs 274 s) |
| GPU AG Car A, 4 A100, 30 cycles: s/cycle and Picard (all: 0.409, 4.33) | 0.407, 4.43 | 0.408, 4.43, **bitwise = fses** | PASS (vet_gd: the new keys are inactive) |
| default behaviour | | bitwise: 6 tests | PASS |

**Problem 1 is fixed.**
- The cause was the outer Marshak face of the half-range solve, not the switch.
- The fix gives clean 2nd order in every atmosphere region above tau = 10, at levels up to 2048.

**rw_t10 1.70 is not a scheme error.**
- It is a uniform drift of the background E between levels, 1.4e-10 relative, which the self-convergence norm divides by the 1e-3 to 1e-5 wave amplitude.
- With the level mean removed, every arm gives 2.02 / 1.94, the same as cen.

**Problem 2 is improved but not closed.**
- The beam E/J error of tEs is mostly the vet_sc formal solution itself, not the face flux.
  - The SC beam luminosity H_sc/H_exact is 0.23 / 0.82 / 0.24 with nphi = 32, and 0.23 / 0.90 / 0.25 with nphi = 128.
  - tEs followed it: E/J_sc is about 1.
- jse makes the thermal half-range flux 2nd order where E is proportional to J. That makes the moment beam 3-4x brighter than tEs, and on xb20 / ba20 brighter than the SC beam itself.
- Collimated beams are now 0.49-1.08 of exact (they were 0.11-0.68).
- cyl moves away from exact (0.35 -> 0.45 <|E/J-1|>). It is a broad isotropic emitter, where the SC J is the reference that E follows (far E/J_sc goes from 1.06 to 0.90).
- An exact beam E/J would need a less diffusive formal solution: bilinear short characteristics lose 75 % of a 20-degree beam in the 5-cell channel. That is outside this fix.

## D. Diagnosis

### D1. Thin top of the atmosphere (problem 1)

**Atm ladder at 32..2048 cells.**
- Self-convergence is computed per tau band (ana_atm_self.py, last-level orders). With 2048 cells, the band tau < 0.1 holds 479 cells.
- RESULTS: cpu/run/RESULTS/atm_self_regions.txt.

| arm | what it isolates | tau < 0.1 | 0.1-1 | 1-10 |
|---|---|---|---|---|
| cen | central, Marshak q = 1/2 | 1.06 1.02 1.01 1.01 1.00 | ~0.95 -> 0.99 | 2.39 1.79 0.81 0.58 0.83 |
| tZ | half-range path, w == 0 (central interior), top q = h+_N (cell) | 0.85 0.91 0.95 0.98 0.99 | 0.56 -> 0.92 | 2.1-2.4 -> 0.99 |
| tE / tA | fse / hj without correction | 0.85 1.06 0.69 0.99 1.15 | ~1 | |
| **tEs / tAs** | fse + hcs | 0.74 0.60 1.52 0.95 0.62 | 2.26 2.41 1.82 0.58 0.13 | 2.0-2.4 |
| tQh | fsq + hc (w ~ 1 through the top, from C1) | 1.91 1.96 1.98 1.99 2.00 | 2.0 | 2.0 |
| hr | production (xthin, dt-dependent) | 0.26 0.79 1.70 1.96 -10.4 | | |
| tW / tWs | w == 1 everywhere, hj h0 = -1 (with hcs) | O(1) errors, no steady state / Picard DIVERGED at 1024 | | |
| **tZb** | **w == 0 + bq = lin** | **2.05 2.04 2.02 2.01 2.01** | **2.0** | **2.0** |
| tEsb | tEs + bq = lin | 2.04 2.02 2.09 2.01 1.99 | 2.0 | 2.0 |
| tJb / tJsb / tJe | fse + js (no hcs) / + hcs / jse + hcs, all with bq = lin | 1.98-2.09 | 2.0 | 2.0 |

**The switch does not flicker.**
- atm_w.py reconstructs w = smoothstep(\|H\|/J) from the m1_vsc output.
- In tau < 0.1 the weight is 0.54..0.98, smooth in z, and the same at every level: 0.541..0.983 at 512, 0.540..0.984 at 2048.

**The cause is the half-range outer Marshak q = H+_N/J_N.**
- The ratio is taken at the top CELL, which is O(dx) off the face value H+_f.
- This first-order boundary error sets E_N. Every arm whose top-band flux is not fully the formal solution's own carries it into the whole thin region: cen, tZ, tE, and tEs through its (1 - w) central part.
- tQh hid it because w ~ 1 there. The thinsw "erratic" tEs orders are a first-order term with a small (1 - w) coefficient beating against second-order terms. Asymptotically they would be 1st order.
- **Fix:** `implicit_thin_bq = lin`, with q = H+_f/J_N and H+_f = (3 H+_N - H+_{N-1})/2 limited to [1/2, 2] H+_N.
  - The vacuum face carries H+ only.
  - With it even w == 0 (tZb) is 2nd order.
  - The half-range reconstruction (dc vs plm) and the hcs correction are not the cause.
- **Measurement.**
  - The band mask (tau from the analytic profile at the coarse-cell centres) and the number of cells (7 to 479) are adequate.
  - The Hopf reference cannot measure order: its floor is the angular quadrature, 3.45e-4 at 2048 for every arm.

### D2. rw_t10 X order 1.70 vs cen 1.96

- tZ (w == 0) and cen differ only at the 256 level: E(tZ) - E(cen) = -4.1e-8 uniformly, on a background E = 300 and a 1.7e-4 wave.
- That is a 1.4e-10 relative drift of the total energy, i.e. round-off over 4200 steps. The two linear-solver paths differ: BiCGStab breaks down on the y-uniform 2-D problem and falls back to line_jacobi in every step, for both arms.
- The self-convergence norm divides by the wave amplitude, so the drift shows up as a 2.5e-4 "error".
- With the level mean removed (ana/od_dm.py), the last orders are:

| | be X E | hd2 X E | rw_t10a (amp 1e-3) be X |
|---|---|---|---|
| cen | 2.03 | 1.95 | 2.03 |
| tZ | 2.02 | 1.94 | 2.03 |
| tEs | 2.02 | 1.94 | 2.03 |
| tJe | 2.02 | 1.94 | |

- At amplitude 1e-3 the raw tZ is 1.79. So the drift is real, but it is a mean offset, not the wave.

### D3. atm tau > 10

- The last orders fall to about 1.5 for every arm, cen included (1.91 1.57 1.38).
- The cause is not relaxation: a 10x longer run (case atmL, T = 20000) gives identical numbers.
- It comes from the thick bottom of the vet_sc / central part, is common to all arms, and is unchanged here.

### D4. Collimated-beam energy (problem 2)

**New diagnostics.**
- ana/beamej.py: beam E against the exact J and against the code's own formal-solution J_sc (m1_vsc).
- ana/beamlum.py: column luminosity at x = 0.5 / 0.7 / 0.9 for the E-solve face flux, the SC H1 and the exact Hx, plus emitter-cell E/J_sc.
- ana/emj.py.

**1. The formal solution itself misses the tilted beams.**

| test | J_sc/J_ex median (beam) | H_sc/H_ex luminosity, x = 0.5 | same with vet_nphi = 128 |
|---|---|---|---|
| xb20 | 0.209 | 0.226 | 0.234 |
| ba0 | 0.655 | 0.818 | 0.901 |
| ba20 | 0.151 | 0.242 | 0.250 |

- Angular resolution is not the cause.
- The bilinear (first-order) short-characteristics interpolation smears a 20-degree beam laterally into the absorbing channel walls.
- tEs's moment E followed the SC: E/J_sc = 0.99 / 1.03 / 0.75.

**2. The emitter.**
- Emitter-cell E/J_sc: tEs 0.51-0.55, cen 0.71-0.82, tJe 0.61-0.64.
- The emitter is under-resolved: tau_cell = 1, 5 cells wide.
- **Suggested cap.** I tested implicit_thin_lte = 0.8: a thermal cell with J_sc close to S (min/max > 0.8) has its thinness scaled to 0 (central).
  - Beams, cyl and shd3b came out bitwise unchanged (arms tJl8, tEl8).
  - The interior emitter faces are already central under fse (C1 small, J ~ S). The surface face takes w = 1 from the vacuum cell.
  - So capping w in thick thermal cells cannot change E/J. The code was removed.
- The deficit is the central flux at tau_cell = 1 plus the surface face. A cap on the surface face itself would give the cen-like surface flux: cen's emitter E/J_sc is only 0.71.

**3. Better surface/thin treatment: `implicit_thin_recon = js`.**
- The half-range face ratios use the face mean of the formal-solution half-range fluxes: h+ = (H+_L + H+_R)/(2 J_L), h- = (H-_L + H-_R)/(2 J_R).
  - The factor is limited to [1/2, 2] and |h| <= 1.
- Where E is proportional to J, the upwind face flux is the SC's own central flux c H_f, which is 2nd order. The operator stays homogeneous in E: no fixed source, M-matrix, |F| <= cE donor.
- It is dt-free and tau_cell-free. The face factor is 1 + O(dx) in smooth fields, so the weight does not move toward upwind as dx -> 0.
- **Effect:**
  - E/J_ex goes from 0.21 / 0.68 / 0.11 to 0.75 / 1.08 / 0.49. F/H_ex at x = 0.5 goes to 0.67 / 1.01 / 0.62.
  - No superluminal faces (xb20 0 %; ba0 / ba20 1.2-1.4 % as before).
  - Direction error unchanged (4.4 / 4.7 / 11.7 deg).
  - shd3b edges sharper.
- **But js slowed the time-dependent scattering pulse:** G5 k12.8 ratio 0.51 -> 0.28. The lagged J (S = E^n) does not describe a time-resolved E, as for hc.
- **Final `jse`** blends js with dc by the thermal fraction eps of the face.
  - Thermal media (beams, Marshak) get js. Scattering media get dc, and in steady state hcs + bq give 2nd order there.
  - Result: G5 back to 0.508, every 1-D ladder identical to tEs, beams identical to js.

**4. What remains.**
- The tilted-beam E/J (0.75 / 0.49) and cyl (<|E/J-1|> 0.45) errors are set by the formal solution's spatial interpolation.
- The cure is a higher-order (e.g. monotone quadratic) SC interpolation in rad_m1_vet.cpp. That is new work and is not attempted here.

## Runs, scripts and results (all under /raven/ptmp/jinma/thinfix_1011)

**cpu/ (copy of thinsw cpu/ with additions)**
- od.py: OD_NX env; new arms tW tWs tAs tZs tZb tEsb tJ tJb tJsb tWb tWJb tJe, plus the P/Pvs forms; case atmL.
- battery_raven.sh: BTAG / BEXTRA / BNLIM for beams; arms tW tWs tEsb tJ tJb tJsb tJl8 tEl8 tJe.
- New scripts: cmd.sbatch, ana/atm_w.py, ana/od_dm.py, ana/beamej.py, ana/beamlum.py, ana/emj.py.
- The atm input has an added m1_vsc tab output.

**cpu/run/RESULTS**
- atm_self_regions.txt, atm_hopf.txt
- order_eval_new.txt (tEsb, tJsb, all ladders), order_eval_tJe.txt (tJe, with the demeaned rw)
- gates.txt (tEs, tEsb, tJsb, tJe)
- beams.txt

**Beam runs:** cpu/run/beams/<test>/{tJ, tJb, tJsb, tJe, tEsb, tJl8, tEl8, *_p128}.

**Default-bitwise checks**
- verify/ (ea3ba71e) and verify_cfb/ (cfb3258f) against the base binary 9d8fd8c9 = 61cab6c0 code.
- Tests: ba0 hr / cen / tEs (bin), G1 atm2d hr (tab + hst), atm tEs 128 (tab), pulse_k0.128 hr hd2 X128 (tab).
- Result: 0 differing data files.

**GPU**
- Script: gpu/gpu_arms_raven.sh with arm fsesje; input gpu/inp/agcar_A_thf.athinput, which has the new keys at their defaults.
- Job 31061760, results in gpu/RESULTS_gpu.txt.
- all 0.409 s/cycle, Picard 4.33. fses 0.407, Picard 4.43. fsesje 0.408, Picard 4.43, bitwise = fses.
- Diffs vs all are the same as thinsw's fses: m1_e < 3e-5 inside 0.95 R_ph.

**Rules followed**
- Held jobs and other dirs/worktrees were untouched.
- Commits are local only, not pushed.
- `git submodule update --init kokkos` was needed in wt-thinfix to build.
