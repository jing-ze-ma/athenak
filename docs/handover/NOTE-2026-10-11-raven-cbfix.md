# NOTE cbfix (cbfix-1011): why the half-range face flux lost the crossing beams, and the fix

Branch `cbfix-1011` (worktree /raven/u/jinma/ATHENAK/wt-cbfix), local only, not pushed.
- **9e1f2686**: `implicit_hr_face = sf`, the new default. It applies with `implicit_flux_beam = halfrange` and vet_sc.
- The base is a61c1e75 (littests-1011).

Binary: bin/athena_cpu_9e1f2686 (gcc/13, openmpi/5.0, OpenMP, CPU).

Figures (white background):
- fig/cb_maps_fix.png: before and after maps, and row energy against height.
- fig/cb_mech.png: the mechanism.

## 1. Root cause: a wrong face flux, not a sink

**Budget.** I added a diagnostic dump (`CbDebugDump`, env `M1_CBDBG=<prefix>`). It writes E^n, E', J, H, H+, the x1/x2 face coefficients and the face fluxes for the first solves. From it I computed the per-cell budget E' - E^n + dt div F.
- It closes to round-off: max |res| 4e-8, total 1e-9 of the injected energy.
- No clip, floor or scale-back moves energy.
  - E floor raises = 0, pos_gas = 0, g0 clips = 0.
  - The |F| <= cE scale-back touches F only. With w = 1, F^n is not read by the E row.
- **The energy leaves through the x1min Marshak face.**
  - Net x1min inflow, as a fraction of the injected beam energy, for steps 1-4:
    - old hr: 3.2e-2, 1.2e-3, 5e-5, 2e-6;
    - old tJe: 5.1e-2, 4.9e-3, 8e-4, 1.8e-4.
  - The rest goes straight back out through the neighbouring Marshak columns. Example: column j = 24 next to the beam column emits 307 of the 707 injected.

**Mechanism** (fig/cb_mech.png):
1. A 45° beam on the 4-ray ring sits on cell diagonals. The formal solution is exact: J = 1 on the diagonal cells D_k, and only the thermal glow J ≈ 3e-20 elsewhere.
2. With face fluxes only, energy can go from D_k to D_k+1 only through the two off-diagonal cells R = (k+1, k) and R' = (k, k+1).
3. Their J is the glow, not 0, so M1HrCell uses the glow's own half-range ratios.
   - Near x1min the glow is dominated by down-going rays: h-_1 = -0.70 and h+_1 = 0.003-0.01. In x2 it is isotropic: ±0.354.
4. So every R cell sends the beam energy it receives back down toward x1min at 0.7c, and sideways at 0.35c.
5. The beam is converted to a downward-diffusing cloud within one row. It drains through the Marshak faces next to the injection columns: row energy falls by e every 5-7 rows.
6. This is independent of c dt/dx (0.5 ... 256), resolution and ring size, as observed.
7. If J were exactly 0, the fallback ±1/4 would leak the same way, only slower.

**The suspicion in the brief was half right:**
- The face does not fall back to central for hr or tJe at c dt/dx = 256. Every face has w = AL = 1 (dumped). xthin overrides hr; fse gives g = eps = 1 in tJe.
- The ratios are defined, because J = 3e-20 > 0, but they describe the glow and not the beam energy.
- js (tJe) makes it worse:
  - it halves the D-cell outflow, since the face mean includes the empty neighbour;
  - it clamps the R ratios with lim = 2 relative to the glow ratios.

**Why cen keeps 0.74:**
- Its face-eliminated central flux has no half-range shape, so nothing points the energy down.
- It is an FLD-like operator: it smears the beam into a uniform 0.019 floor (net x1min inflow 61 % of the injected).
- With the Gauss 8x32 set it piles up, because the E-only operator has no upwind outflow at c dt/dx = 256.

## 2. The fix: `implicit_hr_face = sf` (default; `= cell` restores the old cell ratios bitwise)

**Donor scale, per solve.** ImplicitHrSfBuild runs once per solve, since vet_cell is fixed.
- An appended iw block (M1_NIW_HSF = 8) holds H+_d, H-_d, J and the donor scale. It goes through the transverse halo.
- The donor scale is
  `Jh = max(J, max over the cell's faces of the NET SC flux leaving it, (H_c + H_nb)/2 · n_out)`.
- Jh > J means the formal solution moves more energy through the cell than c J_c can carry (|F| <= cE). In other words, a neighbour's beam crosses it.
- In any diffusive field |H| << J, so Jh = J.
  - An earlier version took the half-range flux instead of the net flux. It triggered at the steep Marshak edge of the thick static-diffusion pulse and changed sdif by 0.04 %. The net flux removed that.

**Face ratios (M1HrSf).** Only where Jh > J:
- `h+_f = [(1 - om) H+_L + om H+_f]/Jh_L` and `h-_f = [(1 - om) H-_R + om H-_f]/Jh_R`, with:
  - H±_f = (H±_L + H±_R)/2, the SC's own half-range face fluxes;
  - om = max(wt_js, 1 - J/Jh).
- An R cell therefore passes its energy on in the SC direction, and never back to the boundary.
- A js/jse donor with Jh = J next to such a cell scales its face share by J/Jh of the neighbour. Without this, tJe gave row energy 1.42 and a beam front at 0.5c.
- Clamps: h+ in [0, 1] and h- in [-1, 0]. HCL >= 0 and HCR <= 0, so the M-matrix structure and conservation are kept.
- No dt and no tau_cell enter.

**Beam measure (M1HrSfFb).** For a cell with Jh > J the beam measure (fb, or the thin switch g) is max(fb, 1 - J/Jh). Its faces are then upwind.
- Without this, hr's residual central part (1 - w)(1 - s2) ≈ 1.4 % at X = 256 still flooded one side of each beam, with row energy up to 14 (version sf1).

**Diagnostic line in the run log:** `implicit_hr_face=sf: cells with Jh > J ...`.
- cb: 9 % of cells. shd: 0.3 %.
- 0 in nlte, eq, sdif, ddif, Marshak, the atm and 1-D ladders, linear waves and G1/G3/G5.

**Order.**
- Where Jh = J the scheme is bitwise the old one. So every smooth/thick region keeps its order: second order in thick and intermediate regions, tJe/js as before.
- In beam-crossed cells the flux is the first-order upwind half-range flux of the SC shape.

## 3. Crossing beams, before → after (64x256, ±45°, 20 steps; rows normalised by the cell-averaged exact)

| arm / setup | row energy y = -1.9 / -1.5 / -1.0 / 0 / 1.0 | peak at y = 0.65 (exact 0.5-1) | FWHM (cells) | floor | Picard | wall (s, 3 thr) |
|---|---|---|---|---|---|---|
| cen, ring 4 (unchanged) | 0.74 everywhere | 0.074 | 3.0 | 1.9e-2 | 2.85 | 15 |
| hr OLD, ring 4 | 0.15 / 0 / 0 / 0 / 0 | 0 | - | 2e-12 | 2.25 | 2 |
| tJe OLD, ring 4 | 0.40 / 0.017 / 0.001 / 0 / 0 | 0 | - | 6e-7 | 2.75 | 4 |
| **hr NEW, ring 4** | **1.006** at every height | **0.417** | **2.6** | **8e-11** | 2.50 | 14 |
| **tJe NEW, ring 4** | **1.006** at every height | **0.417** | **2.6** | **2e-18** | 3.15 | 5.5 |
| hr / tJe NEW, c dt/dx = 0.5 (1100 steps) | 1.006 (front at the top at 0.70c) | 0.417 | 2.6 | 1e-11 / 3e-16 | 2.45 / 2.97 | 54 / 62 |
| hr / tJe NEW, c dt/dx = 4 / 32 | 1.006 | 0.417 | 2.6 | ≤ 4e-10 | 3.2-3.8 | 5-14 |
| hr / tJe NEW, ring 12 / 36 | 1.006 | 0.417 | 2.6 | ≤ 1e-10 | 3.2-4.1 | 6-12 |
| hr / tJe NEW, 256x1024 | 1.006 | 0.417 | 2.6 | 3e-12 / 1e-17 | 9.0 / 4.5 | 1584 / 617 (old 38 / 41) |
| hr OLD / NEW, Gauss 8x32 | 0.29-0.43 → **0.92 / 0.88 / 0.86 / 0.86 / 0.85** | 0.02 → 0.055 (J_sc 0.084) | 15 | 9e-3 → 1.8e-2 | 3.45 | 14 |
| tJe OLD / NEW, Gauss 8x32 | 0.24 → **0.91 / 0.84 / 0.83 / 0.83 / 0.83** | 0.02 → 0.072 | 11 | 2e-3 → 3.9e-3 | 2.90 | 10 |

Old rows reproduce NOTE-littests exactly (opt-out runs `*_r4old` with `implicit_hr_face=cell`). The top row (y = 1.9, next to the x1max Marshak face) reads 1.13 in every new run.

**IDORT pass criterion** (top row energy ≥ 0.9, peak ≥ 0.3 E_b, floor ≤ 1e-3): **PASS for hr and tJe** on the 4-, 12- and 36-ray rings at every c dt/dx and resolution tested.

**Gauss 8x32** (production angles; the beam rides the 40.2° ray):
- J_sc itself conserves each row (1.000).
- E/J per row falls from 1.04 at row 1 to 0.86 by row 32, then stays constant, while the x1 row flux sum F1 is constant to 0.1 %.
  - So the transport is conservative. The SC-smeared beam is carried slightly faster than J_sc implies, and E sits about 14 % below J_sc.
- It does not reach 0.9 against exact, but it is about 3x better than before. NOTE-littests says to judge this case against J_sc.

**Beam shape.**
- A 5-point face stencil cannot move energy across a cell corner. The 1-cell SC beam is therefore carried as a 2-cell zigzag:
  - D cells at 0.42 and the off-diagonal cells at about 0.3;
  - FWHM 2.6 cells, against 1-2 for the exact cell average.
- The energy-weighted transport speed is the exact one: row energy 1.006, and the front time at c dt/dx = 0.5 is right.

## 4. What is left (not changed here)

1. **x1min Marshak face next to the injected beam.**
   - The row-0 cells beside the beam columns hold zigzag energy. The isotropic Marshak outflow (q = 0.5) sends 16.6 % of the injected flux back out (fig/cb_mech.png, right; per-step budget).
   - In those cells the SC has no downward intensity. A half-range q at x1min, H-_1/J of the bottom cell like hr_q at x1max, would remove this.
   - The bottom face belongs to deepatm-1011, so it is not changed on this branch.
   - The row-energy metric still reads 1.006: the zigzag's E-per-flux factor (1.2) compensates the 0.83 flux.
2. **Linear-solve cost.** The beam now really crosses the domain in one implicit step, so BiCGStab has to propagate it.
   - 64x256: inner iterations 13 → 19 (tJe), and 77 for hr, which hits the 200 cap in some passes.
   - 256x1024: tJe 9 → 118, hr 5 → 167 inner iterations per solve; Picard 3.75 → 4.5 (tJe) and 3.7 → 9.0 (hr); wall 41 → 617 s and 38 → 1584 s.
   - cen pays the same: cen cbh takes 753 s. No solve fails to converge.
   - hr is costlier than tJe because of its residual xthin central part.
3. On the Gauss set the moment beam is about 14 % dimmer than J_sc in E (section 3). This is a shape effect of the SC bilinear smearing.

## 5. Validation (binary 9e1f2686, runs/ unless noted)

| test | result | vs old (a61c1e75 code; littests runs or the `implicit_hr_face=cell` opt-out) |
|---|---|---|
| static shadow, ring 12x2 / 36x2 | hr L1 0.012, umbra 1.40 (0.9 % of lit 158.2), wedge 151.0; tJe L1 0.011, umbra 0.22 (0.14 %), wedge 165.3; edge 10-90 % 16 cells | old 1.41 / 0.22 / 150.9 / 165.3; edges identical. **PASS** (L1 ≤ 2 %, umbra ≤ 1 %). sf cells 0.3 %. |
| shadow, Gauss 8x32 | hr L1 0.029, umbra 1.02; tJe L1 0.030, umbra 0.21 | old 0.029 / 1.03, 0.030 / 0.21 |
| non-LTE atmosphere, eps 1e-1 ... 1e-8 (Gauss 8x32) | L1 1.1 % → 0.01 %, E_top/ex 0.998-1.005 | identical to 4 digits (field diff 2e-13, round-off); 0 sf cells |
| thermal equilibrium 1 / 2 | E 96.795 / 145.788, T 3.13663 / 3.47472 (hr's spurious v = 0.021 unchanged) | identical |
| static diffusion (sdifk) | L1 0.10 / 0.11 / 0.32 %, peak 1.0006 | identical (E diff 1e-15) |
| dynamic diffusion (ddiff) | L1 0.1145, peak 0.875, centroid exact | identical (round-off) |
| Marshak (implicit, c dt/dx 1 / 10) | | bitwise identical data |
| linear waves (hr, tJe; 30 cases each, N = 16-256) | tJe 2nd order everywhere (P = 10, tau = 1 last order 1.30 as before); hr unchanged non-convergence at P = 10 | **table identical** to littests (ana/rw_new.txt vs rw_old.txt) |
| 1-D ladders, X and C modes (be; C hesdirk2 for tJe): pulse k0.128 / 12.8 / 1280, rw_t10, rw_t1000; hr, tJe | tJe X orders e.g. 1.76 1.97 2.15 / 1.89 2.25 2.25 / 1.89 1.97 1.93 / 2.17 2.02 1.70 | **bitwise = opt-out** for all 250 runs; 0 sf cells |
| atm by region (S, 32 ... 2048) | tJe orders tau < 0.1: 2.04 2.02 2.09 2.01 1.99; 0.1-1: 1.92-2.00; 1-10: 2.00; > 10: 2.01 2.04 1.77 1.49 1.47 (= NOTE-thinfix) | **bitwise = opt-out** (hr, tJe); 0 sf cells. Compared like with like on this branch (no deepatm bottom face). |
| G1 (atm2d, c dt/dx 1e2 / 1e4) | tJe L1 7.08e-4, E_top 1.7355; hr L1 1.56e-3 / 1.95e-3 | identical to opt-out (Picard 2.07 / 16.4 tJe) |
| G3 Marshak 2-D (c dt/dx 1 / 10) | PASS, L1(E) hr 0.0138 / 0.0140, tJe 0.0139 / 0.0142 | identical |
| fixed G5 (g5fix pulse, kappa 12.8 ... 12800) | all PASS: hr 0.99975-1.00017, hr_be 1.00011-1.00039, tJe_be 0.99958-1.00000, tJe_h2 0.99957-1.00000 | identical to opt-out and to g5fix results.txt |

**Cost.**
- Identical everywhere except where beams cross empty or glowing cells: crossing beams, section 4.2.
- The per-solve build is one extra kernel plus 2·dim + 2 halo exchanges.
- Shadow wall time: 15 s vs 13 s, Picard 2.0 / 3.9 unchanged.

## 6. Files (work dir /raven/ptmp/jinma/cbfix_1011)

**Scripts**
- run1.sh (as littests), rund.sh (with M1_CBDBG), pool.sbatch / poold.sbatch (general, full node), build_cpu.sh, rw.sbatch, lad.sbatch, gates.sbatch.
- inp/: copies of the littests inputs with `implicit_hr_face = sf` added, so the key can be overridden on the command line.

**Analysis (ana/)**
- cbs.py: crossing-beam summary.
- dbg.py: dump reader and budget.
- bcmp.py: numeric dump comparison.
- fig_cbfix.py, fig_mech.py: the figures.
- proto.py: Python prototype of the steady upwind solve.
- rw/: linear-wave runner with ROOT = cbfix.
- rwtab.py; rw_new.txt and rw_old.txt.

**Runs**
- runs/ (final, 9e1f2686): cb (including `*_r4old` opt-out and `*_d2old/new` dumps), cbg, cbh, shd, shdg, nltek1-8, eq1-2, sdifk, ddiff, marsh, rw.
- lad/run/runs: ladders and atm, arms hr, hrC, tJe, tJeC.
- gates_run/gates: G1/G3 (+ old G5).
- g5/run: the fixed G5.
- runs_sf2/, runs_sf3/, runs_sf4q/: intermediate versions, kept for reference.
  - sf2: tJe without the js-donor rule, row 1.42.
  - sf3: Jh from the half-range flux.
