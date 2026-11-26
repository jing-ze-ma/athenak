# NOTE 2026-11-27za orion -> Caltech: HLLC deep band in m35mms: resolution test (HLLRES, 4 fresh runs)

From agent HLLCRES (orion), on the USER's order: "run the resolution test for the HLLC deep band in M35 MMS". Your ring work (26h-27s) found: the band at 0.72-0.76 R_ph,IC is set by the hydro (R6 hlle removes it; rays, WB, dt, face-flux form, FOFC do not), and M35C (128 x 128) put it near the same cells per wavelength (9.1 against 10-11, k ratio 0.56-0.61), not at the same l. Open: (a) only one coarser lateral grid, on another binary than the release; (b) no finer grid; (c) no radial refinement (the band's radial width / nodes in cells or in R?). This test closes those, all on the RELEASE binary with HLLC (the default).

## Choice of grids (why a 10 x 10 deg wedge)
A 2x lateral run of the full 20 x 20 deg wedge (512 x 512) needs 4x the memory of production (about 150 GB per H200 on 4 GPUs: does not fit) or 4 nodes for ~10 h. Instead all arms use a **10 x 10 deg wedge centred on the same equatorial patch** (theta 85-95 deg, phi 0-10 deg): its 256 x 256 grid has the lateral cell of a 512 x 512 full wedge. The 128 x 128 narrow arm (W1, production's lateral cell, M35C's cost) is the control for the wedge change: it must match prod_c1_1124's 256 run (same lateral cell) in band k20, cells per wavelength and growth, otherwise the narrow-wedge ladder is not used and I will ask for the full-wedge 512 run on 4 nodes instead. Radial grid: the stretch map is a function of the uniform index fraction, so nx1 = 672 with the same f_stretch_r splits every production cell exactly in two (old faces are the even new faces).
Restarts cannot change the grid, so all arms are FRESH starts from the IC with the same seed keys (as M35C); the band grows from 1e-8 to 1e-2 km/s in 3e5 s (10+ e-folds, e-fold 1.1e4-5e4 s), so its growth, k and width can be read on all grids at equal times, and against prod_c1_1124 (256, full wedge) and M35C (128, full wedge) dumps you already have.

## Binary
athena_rel6_01deec542 (md5 7b0900a6, CUDA 12.9 + g++13, star_vet, your 27v build in wt_orion1121/rel6). Release tip aa2646b92 differs from 01deec542 in two docs files only (`git diff --name-only 01deec542 aa2646b92`: docs/vet_clean/PERFORMANCE.md, docs/vet_clean/performance/wall_times.csv), so this IS the release code. Never a CUDA 12.2 build.

## Decks
Base: the release template inputs/stars/star_m35mms.athinput at aa2646b92, TABLES/ replaced by your table copies, nothing else changed except the lines below. Fresh starts: no rst fallback applies, but name every line below anyway so the effective values are in the log.

Common to all arms (add / replace):
```
<mesh>
x2min = 1.48352986419518            # 85 deg  (HLLRES: 10 x 10 deg wedge, production 80-100)
x2max = 1.65806278939461            # 95 deg
x3min = 0.0
x3max = 0.174532925199433           # 10 deg  (production 20 deg)
f_stretch_r = -0.455726942889053    # unchanged, named
<time>
tlim = 3.5e5                        # as M35C
cfl_number = 0.3                    # release default for a fresh VET run, named
<hydro>
rsolver = hllc                      # release default (recipe.cpp), named
reconstruct = plm                   # inert on spherical polar (your 27s), named
fofc = true                         # release star value (KEYS.md), named
fofc_retest = true
<rad_vet>
vet_gd_gl_nmu = 7                   # the default set, named
vet_gd_gl_nphi = 18
vet_gd_gl_stagger = true
vet_gd_wall_interp = true           # release defaults for a fresh SP wedge; REQUIRED here because
vet_gd_wall_sym = true              #   10 deg is half the GL-18 azimuthal period (KEYS.md)
vet_gd_theta_seam_sym = true        #   theta wedge centred on the equator, nmu odd
clean_stag_mg_levels = 4            # unchanged; tile 64 -> coarsest 4 cells (as M35C)
<problem>
star_seed_band = n2
star_seed = 1.0e-3
star_seed_lmin = 16                 # = production's printed window, fixed on every grid
star_seed_lmax = 64
star_seed_rng = 1234
<output2>                           # hydro_w every 1800 s (frequency needs it), float32 as production
dt = 1800.0
<output3>
dt = 1.0e30                         # vet off
<output5>
dt = 1.0e30                         # rst: only t = 0 and the final one
<output6>
dt = 1.0e30
<output7>
dt = 1.0e30
```
If vet_gd_wall_sym on the 10 deg wedge FATALs or gives a CLEAN WARNING at start-up, use vet_gd_gl_nphi = 36 instead (whole period; your R1 showed the ray set does not change the band) and say so; do not use any other fallback.

Per arm (grid lines):
| arm | mesh nx1 x nx2 x nx3 | meshblock | layout_tile | ranks / GPUs | lateral cell vs production | radial cell vs production | est. wall to 3.5e5 s |
|---|---|---|---|---|---|---|---|
| W05 | 336 x 64 x 64 | 336 x 64 x 64 (1 block) | 64 | 1 | 2x (= M35C) | 1x | ~2-2.5 h |
| W1 (control) | 336 x 128 x 128 | 336 x 64 x 64 (2 x 2) | 64 | 4 | 1x (= prod 256) | 1x | ~2.1 h (M35C cost) |
| W2 | 336 x 256 x 256 | 336 x 128 x 128 (2 x 2) | 128 | 4 | 0.5x | 1x | ~6-8 h (production cost per cycle) |
| W1R2 | 672 x 128 x 128 | 672 x 64 x 64 (2 x 2) | 64 | 4 | 1x | 0.5x | ~8 h if dt halves (2x cells, 2x cycles) |
| W1E (optional, last) | as W1 | as W1 | 64 | 4 | 1x | 1x | ~2.1 h; ONLY change hydro/rsolver = hlle |
(The seed needs one MeshBlock across x1: every arm has it. W05 on one GPU keeps tile 64 and mg 4; do not split it into 32-wide blocks.)
Nodes: W1R2 and W2 first (critical path ~8 h), W1 and W05 together on a third node (W05 1 GPU, W1 4 GPUs: if they do not fit together run W05 after W1). W1E only if a node is idle after these; it is the fresh-start HLLE counterpart of W1. Directory suggestion /resnick/groups/carnegie_poc/jingze/hllres_1127/<arm>. Only idle test nodes; no production stopped. Dumps: W2 ~86 GB, W1R2 ~43 GB, W1 ~21 GB, W05 ~5.5 GB (~160 GB total).
Guard: stop an arm at its first FATAL / CLEAN WARNING; after 200 cycles post s/cycle, dt and the wall estimate per arm. If W1R2's estimate exceeds 12 h, set its tlim to 2.6e5 s (enough e-folds) and say so.

## Checks at start-up (post in the first NOTE)
- the seed line of every arm (l window, rng, amplitude, radial band r_lo-r_hi; production's was 8.9756e11-9.34835e11, 29 cells); W1R2's radial band in cm must be the same within one production cell (say if not);
- the wall / seam key lines (vet_gd_wall_sym, theta_seam_sym) and that rsolver = hllc (W1E: hlle) is in the log;
- the effective key list of W1 against the stored keys of prod_c1_1124's m35mms rst (keys that differ, by name and value): this says what the release changes against the R0 harness;
- dt and s/cycle after 200 cycles.

## Measurements (definitions as your spec_k.py of a0dc2503; v_r minus the shell mean, 2-D FFT per shell)
Wavenumber in all tables as **k20 = k per 20 deg** (k20 = 2 x k per 10-deg wedge) so that the narrow and full wedges compare directly, plus **cells per wavelength** (= lateral cells per 20 deg / k20) and l ~ 18 k20. Band = the peak-centred band k20_peak / 1.25 .. 1.25 k20_peak (fixed windows do not transfer between grids); also give the A/B windows of 26k converted to k20 for continuity.
1. Band rms at r/R_ph,IC = 0.72, 0.74, 0.76 every 1e4 s from 5e4 s to the end, for every arm, next to prod_c1_1124 (256, full) and M35C (128, full) at equal times (within 900 s).
2. Peak k20 and power-weighted k20 (and cells per wavelength) at 1.8e5, 2.6e5, 3.46e5 s at the three radii; axis share (|kx| <= 1 or |ky| <= 1) with its counting share.
3. Growth e-fold of the band rms in the windows 1-2e5, 2-2.6e5, 2.6-3.5e5 s (as 26k), and ln-residual scatter.
4. Radial structure at 2.6e5 and 3.46e5 s: radius of the peak of the band rms, its FWHM in R_ph,IC AND in local radial cells, the node radii and the node spacing lambda_r (phase jumps of the strongest coefficient), for W1 against W1R2 (and W05/W2).
5. Frequency of the strongest coefficient from its phase (1800 s dumps) against N at 0.74 (1.60e-4 rad/s), as 26h; mark the alias caveat.
6. Run health: dt, s/cycle, wall, CW / FATAL / DT_RETRY, and the hst means of KE_int, KEr_int, L_top/L_in per unit solid angle against prod_c1_1124 over 2-3.5e5 s (the convection zone in a 10 deg wedge).

## Figures (plot style: docs/handover/plot_style, `import idort_style; idort_style.use()` before plotting, as NOTE 24zzzzv)
F1 band rms at 0.74 vs t, log y, all arms + prod_c1 + M35C. F2 shell power spectrum at 0.74, t = 2.6e5 and 3.46e5, x = k20 (top axis: cells per wavelength), one line per grid. F3 band rms vs r at 3.46e5, normalised to its peak, against r/R_ph,IC (left) and against the radial cell index from the wall in production cells (right), W1 and W1R2 (+ W05, W2). F4 v_r maps at 0.74 at 3.46e5 for W05, W1, W2 (same physical patch, same colour range per panel noted). PNG + the numbers as text, in caltech_results/HLLRES/.

## Fixed reading (written before the runs)
- **Lateral grid artefact:** across W05 / W1 / W2 the peak k20 scales with the lateral cell (about x2 per halving, i.e. the same 9-11 cells per wavelength), and/or the growth rate or the amplitude at equal t changes monotonically with resolution by more than 25 % (M35C was 3.3-4.2x larger than the 256 run at 3.46e5 s).
- **Physical:** the peak stays at the same l (k20 22-26) on all three grids and the growth e-fold agrees within 25 % (and W1 = prod_c1 within 25 %).
- **Radial:** the band's FWHM / lambda_r is the same in local cells in W1 and W1R2 (radial grid sets it) or the same in R_ph (physical); a band that vanishes in W1R2 or W2 counts as grid-set.
- W1 must match prod_c1 (same lateral cell) within 25 % in growth and 20 % in k20; if not, the wedge change matters and this ladder is reported as inconclusive.

Report in one NOTE caltech-orion-hllres (start-up checks first, readings when W1/W05 end, final when W2/W1R2 end). Thank you.

-- HLLCRES (orion)
