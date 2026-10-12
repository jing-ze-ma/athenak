# NOTE-bubbles (worker BB, 2026-10-12): IDORT multi-source bubbles (Ma+2026 Sect. 4.2.2)

Branch **bubbles-1012**, worktree /raven/u/jinma/ATHENAK/wt-bubbles (kokkos symlinked, not committed). Local only, not pushed.
- Base: bsgcand-1012 d3bf42d7.
- 0aedee55: vet_sc rotation + vacuum transverse faces.
- fe39c749: merge of defaults2-1012 3cd2a2fb (fresh default `vet_gd_rotate_slope correct`) + one beam-pinning rule.

Work dir /raven/ptmp/jinma/bubbles_1012.
- Binaries: bin/athena_cpu_{bc3 = 0aedee55, m4 = fe39c749}, bin/athena_he_cpu_{bc3, m4}.
- Inputs inp/, runs runs/, analysis ana/, figure fig/bub_outflow.png.

## 1. Verdict: TEST (boundary condition). The "front lag" was the periodic box.

**Single disc** (T 10, r 0.03, ρκ 1e5, FIN keys, ring 100). Runs: runs/one (periodic unit box, 256²) and runs/oneo (2x2 box, outflow).

- **Periodic box (as in the suite deck).** The steady vet_sc formal solution wraps rays through x1 (vet_x1_periodic, npass 4) and x2.
  - In vacuum, J_sc is flat at about 2600: 12x the exact T⁴ asin(R/d)/π at d = 0.47. H/J is only 0.03-0.17.
  - This is the lossless periodic equilibrium, not the retarded field. So the Eddington tensor is near-isotropic far from the source.
  - The moments therefore diffuse instead of streaming: F/cE 0.2-0.3. E piles up near the disc (1.4-1.9x exact). The front is at about 0.12 when ct = 0.25. E outside the disc is 0.46 of exact.
  - BE and hesdirk2 are identical; c dt/dx 0.26 and 2.56 are identical.
- **x1 opened** (outflow + marshak; vet_x1_npass 1 gives the same result). Along the x2 axis: H/J 0.99, F/cE 0.95-1.0, front at about ct.
  - Along the x1 axis it is still wrong (J/Jex up to 3.7). Steep rays still wrap in x2.
- **All faces open** (this branch, runs/oneo, 2x2 box): the front is at ct, F/cE 0.98-1.0, and E/Eex is flat in radius at 0.77-0.78.
  - The ratio is constant in time, so energy is no longer lost; what remains is a luminosity deficit, section 4.
- **IDORT setup.** The user (IDORT author) checked: the radiation boundaries are OUTFLOW on all faces (paper Sect. 4.2.2). PeriodicBoundariesOn 1 in param.txt is the Arepo mesh flag only.
  - So the suite deck (periodic) was not the IDORT test. With a steady formal-solution closure, a lossless periodic vacuum box is a different problem: the steady periodic field is the 4π equilibrium, not the retarded one.
- **The rest of the setup matches IDORT** (create.py):
  - 20 discs (logistic map), r ≤ 0.03, T ≤ 10, ρκ = 1e5 inside, vacuum κ = 0;
  - c = 10, a = 1, dt_max 1e-3, TimeMax 0.1;
  - scored at 0.05 / 0.25 / 1 t_c = snapshots 1 / 5 / 20.
  - Differences: IDORT starts with E = T⁴ inside the discs (we start from the floor; the discs equilibrate in one step at ρκc dt = 1e4). IDORT uses 1024² and 200 angles in 2-D; we use 512² and ring 100 / 200.

## 2. Changes (general, not bubble-specific)

1. **Vacuum transverse faces** (0aedee55).
   - New keys `<rad_m1>/implicit_bc_x2min|x2max|x3min|x3max` = auto (default: periodic on a periodic mesh, else reflect = zero flux; unchanged and bitwise) | reflect | marshak.
   - marshak puts the free surface F = c marshak_q E (outward) in the implicit transverse row: Cartesian, sp and cs rows (rad_m1_implicit.cpp, m1_impl_tcell).
   - vet_sc on a non-periodic x2 (x3) mesh: the SC sweep treats the face as vacuum (a foot point outside the box gets I = 0, with the edge cell's χ and S), and the tensor ghosts are edge-copied. This needs implicit_bc_x2*/x3* = marshak (FATAL otherwise) and one MeshBlock.
   - vet_gd FATALs on transverse marshak (its sweep has no vacuum transverse face).
2. **vet_sc angle-set rotation on the vet_gd clock.**
   - VetRotInit (rad_m1_vetgd.cpp) is the one parse of vet_gd_rotate_dt|t0|align|stretch|slope for both closures, with the same fresh defaults: auto = 20 x first dt, align true, stretch 0.02, slope correct (since the defaults2 merge). Restarts lacking the keys keep the old behaviour (vet_sc: no rotation).
   - The set changes in one place, RadiationM1::VetRotUpdate → VetBuildAngles(beta), using VetGdBlock / VetGdAngleBlk (golden-ratio block sequence):
     - double-Gauss: azimuth about x1 turned by (2π/nphi) q;
     - ring: 0.8 (q − ½) (2π/N), so a ray never reaches ±x2, which an x1 sweep cannot carry.
   - VetGdRotStep / Time2RotCorrect now act for vet_sc too.
   - Banded multi-block sweep: the local ray copy is refreshed, and the band covers every rotated set.
   - Hook for AV's averaging: VetRotUpdate (set-change event) and vsc_beta (state).
3. **One pinning rule for both closures.** Boundary beams (vet_gd_beam; implicit_beam_n), and vet_sc's vet_milne diagnostic, inject along fixed rays of the set. They therefore get no default rotation, and a named rotation with them is FATAL (as vet_gd_beam already was).
   - Affected decks: cb, cbg, cbh, shd, shdg, cloud. They are unchanged by construction.

## 3. Bubbles with outflow (IDORT setup)

Deck inp/bubu.athinput: unit box, outflow + marshak on x1 and x2, vet_x1_periodic false.
- Cross-check: inp/bubo.athinput (2x2 box; open-equivalent to t = 0.11) gives the same numbers to 3 digits.
- Reference: retarded, OPEN (no images): idort_page/figsrc/bubref.py 'open'. Metric: bubeval sampling every 4. Script ana/bubo_eval.py.

| arm (binary) | L1(E) at 0.05 / 0.25 / 1 t_c | E_tot | L1(F) |
|---|---|---|---|
| ring 100, no rotation (bc3) | 0.238 / 0.306 / 0.257 | 0.83 / 0.76 / 0.75 | 0.58 / 0.46 / 0.30 |
| ring 100, rotation (m4, slope correct) | 0.233 / 0.301 / 0.257 | 0.83 / 0.76 / 0.75 | 0.57 / 0.45 / 0.30 |
| ring 200, no rotation (bc3) | 0.233 / 0.299 / 0.252 | 0.83 / 0.77 / 0.75 | 0.57 / 0.45 / 0.29 |
| ring 200, rotation (m4) | 0.235 / 0.298 / 0.251 | 0.83 / 0.77 / 0.75 | 0.57 / 0.45 / 0.29 |
| Gauss 8x32 ρκ 1e3, no rotation (m4; vs E3) | 0.183 / 0.253 / 0.182 | 0.86 / 0.84 / 0.83 | 0.52 / 0.41 / 0.23 |
| Gauss 8x32 ρκ 1e3, rotation (m4; vs E3) | 0.183 / 0.254 / 0.187 | 0.86 / 0.84 / 0.83 | 0.52 / 0.41 / 0.23 |
| old: periodic suite deck FIN r100 k5, vs periodic ref | 0.33 / 0.61 / 0.67 | 0.71 / 0.52 / 0.46 | |
| old: periodic run, vs open ref (Q2 box2) | 0.31 / 0.57 / 0.69 | 0.72 / 0.58 / 0.87 | |

- **Rotation and angle count matter little:** at most 0.005 in L1. The ring lumps are not the dominant error.
- **IDORT's own result.** IDORT's notebook snapshot at 1 t_c (bubble_2d/plot.ipynb output, 'hot' log colour bar 10..1.5e5), colour-inverted (ana/idortcmp.py):
  - The inversion is calibrated on the isolated discs: E = T⁴ inside, ratio 1.00 ± 0.01.
  - IDORT's field is a median 1.24x the open exact field (p10 1.00, p90 1.36), and L1 vs exact is 0.27.
  - Ours (r100) is 0.73x with L1 0.25. **At 1 t_c we are as close to the exact solution as IDORT itself, but on the low side.**
- **Remaining morphology** (fig/bub_outflow.png): fronts at ct but not perfectly circular. There are elongated small bubbles along x1 at 0.05 t_c and a straight diagonal edge on the large bubble at 0.25 t_c: grid anisotropy of the half-range transport.
- **Verdict on the suite criterion.** L1 ≤ 0.20 (set by us) is met by Gauss k3 at 0.05 and 1 t_c (0.18, 0.18), not by the ring k5 arms (0.23-0.31). IDORT itself scores 0.27 at 1 t_c against the same reference. I propose PASS-with-note: "agrees with the exact retarded field as well as IDORT does; residual = disc luminosity (sect. 4)".

## 4. OPEN (routed to UN): luminosity of small opaque sources

Single disc in the 2x2 open box (runs/oneo; ana/lum.py = ∮F_r dl / (2 R c T⁴); ana/lumfs.py = the same with the formal solution's H):

| run | moments L/L_ex (r = 0.04-0.25) | formal solution L_FS/L_ex | J_FS/J_ex |
|---|---|---|---|
| ρκ 1e5, dx 1/256 (b3k5) | 0.75-0.80 | 0.56-0.61 | 0.57-0.62 |
| ρκ 1e7 (b3k7) | 0.78 (same) | | |
| ρκ 1e3 (b3k3) | 0.80-0.85 | 0.77-0.83 | 0.79-0.84 |
| ρκ 1e2 (b3k2) | E_out 0.84 | | |
| ρκ 1e5, dx 1/128 (b3k5n256) | 0.73-0.84 | | |
| ρκ 1e5, dx 1/512 (b3k5n1024) | 0.76-0.81 | 0.61-0.65 | 0.61-0.63 |

- F/cE is 0.98-1.0 beyond r = 0.06, so transport is fine. The deficit is set at the disc surface and is constant in time.
- The moments emit more than their own formal solution: the half-range upwind flux of the disc cell, c h⁺ E, is close to the right surface flux.
- The formal solution under-emits from an optically thick disc by 40 %, and that does not converge with dx (0.57 → 0.62 from dx 1/256 to 1/512).
- The coordinator notes that resolution independence argues against a pure interpolation error. Candidates for UN:
  1. the SC bilinear interpolation of χ and S at the foot point across a sharp opaque edge, which kills grazing vacuum rays and replaces them with a·S. This could be resolution independent, because the staircase geometry is self-similar;
  2. the thin-face hc correction using the FS's own (low) face flux;
  3. the gate q = min(E/J, J/E)².
- **Decks for UN:** inp/oneo.athinput (2x2, 512², disc T 10, r 0.03 at the centre, outflow, marshak on all faces; 40 steps, about 1 min on 12 threads). Arms in jobs_o3.txt; analysis ana/oneo.py, lum.py, lumfs.py.

## 5. Other decks with transverse boundaries (task iii)

- **shd / shdg / cloud: keep periodic x2.** IDORT shadow_2d is periodic in y: the irradiated-cloud halos wrap across y in IDORT's own snapshots (NOTE-check). These are beam decks, so they are pinned and never rotate.
- **pipe: no change needed.** The x2 = ±2 edges lie inside the thick wall material (Gentile geometry mirrored), so the BC does not matter.
- **hcyl.** Its exact reference (ana/hcyl.py) is DEFINED with periodic images in x2 (a row of cylinders), and the steady problem has no light-travel issue. So the periodic deck is self-consistent and stays.
  - As a check of the new BC on a steady test, hcyl with open x2 (inp/hcylo.athinput, marshak x2) was scored against an open reference (ana/hcylo.py, kmax = 0):

  | run | E_L1_out | F_L1_out | rF ratio | fxx(3) code / exact |
  |---|---|---|---|---|
  | open g8x32 | 0.020 | 0.014 | 0.983 | 0.641 / 0.6415 |
  | open n512 | 0.013 | 0.013 | 0.991 | 0.637 / 0.6415 |
  | periodic deck, runs_js g8x32 | 0.018 | 0.012 | 0.988 | 0.294 / 0.297 |
  | periodic deck, runs_js n512 | 0.014 | 0.013 | 0.993 | 0.292 / 0.297 |

  Same quality; the vacuum faces work in a steady test.

## 6. Validation

- **Bitwise, rotation off** (`vet_gd_rotate_dt=0` named on the new binary; ana/bwcmp.py, all data arrays). bc3 vs d3bf42d7:
  - Cartesian eq1, cb r4, shd r12g2, hcyl g8x32: identical. marsh c1, nltek3, sdifk: identical files.
  - Reduced CPU AG Car A, FIN deck, 60 cycles: 12 dumps + both hst identical, Picard 2.722/6 both.
  - Reduced BSG (bsgR_FIN, 4 blocks, 30 cycles): hst written every cycle, identical; Picard 2.797/9 both.
  - Job 31108235; runs in suiteF_1012/runs_bbw_{d3,bc3} and bw/{agc,bsg}_{d3,bc3}.
- **Restart.** d3bf42d7 hcyl g8x32 wrote a restart at cycle 32. The new binary (m4) restarted from it to cycle 64: no rotation (the key is absent from the file), and every dump is identical to the d3bf continuous run (rst/hcyl, job 31109247).
- **Cartesian suite with the new fresh default** (vet_sc rotates unless pinned; m4; runs in suiteF_1012/runs_bb): see sect. 7.
  - Key set: the suite's FIN arm, which names no rotate keys, so the code defaults apply: stretch 0.02, slope correct, js_face hcface. That is FIN2-equivalent for vet_sc.
  - Baseline: runs_js. That code does not rotate vet_sc, so runs_js is also the FIN2 control of the unchanged code.

## 7. Cartesian suite: m4 (rotating default) vs runs_js

All 40 runs: rc 0, 0 NON-CONVERGED.
- 29 rotate (every non-beam vet_sc deck). The 11 beam decks are pinned, and their numbers are identical to runs_js.
- Metrics files: ana/cart_{js,bb}.txt (ana/evalS.sh).

**Changes against runs_js, to the printed digits** (everything not listed is identical):

| test | runs_js | runs_bb (rotating) |
|---|---|---|
| hcyl g8x32: E_L1_out / F_L1_out / rF | 0.0178 / 0.0122 / 0.988 | 0.0172 / 0.0116 / 0.988 |
| hcyl n512 | 0.0137 / 0.0127 / 0.993 | 0.0144 / 0.0104 / 0.996 |
| eq3t, column 3 | 4.462475 | 4.462474 |
| pipe c2200, last probe T | 0.557 | 0.556 |
| bub (OLD periodic deck) r100k5, L1 at 0.05 / 0.25 / 1 t_c | 0.246 / 0.571 / 0.714 | 0.237 / 0.579 / 0.719 |
| bub r100k3, L1 | 0.248 / 0.604 / 0.762 | 0.246 / 0.616 / 0.770 |
| bub g8x32k3, L1 | 0.271 / 0.612 / 0.769 | 0.271 / 0.606 / 0.765 |

- **Picard.** Within ±2 % everywhere, except hcyl: 2.20 → 2.41 (+10 %), from the slope-correct solves at block jumps.
- **Steps and time.** Alignment adds steps where dt is not a divisor of the 20 x dt0 block: ddf 6100 → 6156 (+0.9 %, +9 % wall), eq3t 341 → 344, rs_m2 +1; elsewhere identical. Wall time −11 % to +9 %.

**Default policy: rotation as the vet_sc fresh default is neutral on the Cartesian suite.**
- No failure, no metric worse than noise, small gains on hcyl F.
- It gives no measurable gain on the bubbles (≤ 0.005).
- It costs about 10 % Picard on hcyl and about 1 % more steps on dt-limited decks.
- It is the shared one-method default (same as vet_gd), so I left it ON as committed.
- If the user prefers vet_sc not rotating, the single switch is the `vet_sc && !sc_pin && !rev_named && fresh` branch in VetRotInit (rad_m1_vetgd.cpp).

## 8. Open items / not done

- **Moment Marshak transverse face:** wired into the E row (diagonal + tt). The persistent boundary face flux F0 stays 0, as for a reflecting face, so the cell-centred F of a boundary cell misses the outflow half. This only affects output F and the gas momentum split at that face.
- **vet_sc vacuum faces need one MeshBlock** (FATAL otherwise). The banded multi-block sweep still wraps.
- **GPU (A100) bitwise of the Cartesian kernels** was not run (CPU only). The new branches are separate statements guarded by run-time flags that are false by default.
- **hesdirk2 fallbacks.** In the single-disc runs the hesdirk2 stage is "NOT ADMISSIBLE" on most steps and falls back to BE: min(E, T) ~ 1e-13 in vacuum cells ahead of the front. This happens in the bubbles too (about 18 of 100 steps). It is pre-existing and unrelated to this branch.
- **Page / idort_results.json:** not edited. Suggest that the 12_bub page owner replace the periodic deck with inp/bubu.athinput, score against bubref 'open', and add IDORT's own inverted snapshot (ana/idortcmp.py) as a column. Figure: fig/bub_outflow.png.
