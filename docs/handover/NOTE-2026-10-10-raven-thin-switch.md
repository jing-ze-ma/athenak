# NOTE thin-switch (thinsw-1010): a dt-free, step-frozen optically-thin switch for the implicit VET face flux

Branch `thinsw-1010` (local, NOT pushed) in worktree /raven/u/jinma/ATHENAK/wt-thinsw, on fork/xthinfix-1009 c696c7d1
(= e50ec7d6 code) + cherry-pick of torder 8ee1e217 (time2_vet_sc). Work dir /raven/ptmp/jinma/thinsw_1010.
Final commit c3ea5bec (code = 9d8fd8c9 + a comment fix). Binaries: bin/athena_cpu_9d8fd8c9 (md5 43d5732e...),
bin/athena_he_a100_9d8fd8c9 (md5 a1cd652b...).

## Verdict

**Best variant: `implicit_thin_switch = fse` + `implicit_thin_corr = hcs`** (arm tEs). It is on top of the hr keys,
with `implicit_blend_xthin = 0` (xthin is ignored in the new mode anyway):
```
implicit_flux = blend   implicit_flux_faces = all   implicit_flux_beam = halfrange   implicit_hr_recon = dc
implicit_thin_switch = fse   implicit_thin_corr = hcs   (defaults: implicit_thin_h0 0.3, implicit_thin_h1 0.6)
+ for hesdirk2 2nd order in time (torder): time2_vet_sc = predict, implicit_vimp = true, time2_vstage = true
```
`fse` alone (tE) is the minimal variant. It gives the same results everywhere except the scattering atmosphere,
where `hcs` turns the thin top from 1st order into roughly 2nd order.

**Status against the user criteria:**
- **PASS:** the 1-D X ladders = cen, T ladders with predict 2, gates, all 5 beams stable and cheaper than hr, GPU
  cost = `all`.
- **Marginal:** the atm thin-top order (tau < 0.1: 0.60 then 1.52, erratic). rw_t10 X last order 1.70 (cen 1.96);
  this comes from the half-range path with w = 0 (control arm tZ), not from the switch.
- **Open accuracy issue:** collimated-beam energy E/J. Every half-range arm misses J here, hr included.

| arm (keys) | thin X (pulse k0.128 hd2; atm tau<0.1, Kn>0.3) | intermediate X (k12.8 hd2; atm tau 0.1-1, 1-10) | thick X (k1280; atm tau>10) | rw_t10 / rw_t1000 X (E) | T hd2 with predict (k0.128/k12.8/k1280/rw10a/rw1000a E) | gates G1 (max dE/E) / G3 / G5 thick | beams xb20 ba0 ba20 cyl shd3b (Picard / BiCGStab inner) | GPU AG Car s/cyc, Picard (all 0.403-0.407, 4.33) | verdict |
|---|---|---|---|---|---|---|---|---|---|
| cen (ref) | 1.94 2.04 / 1.01, 0.94 | 2.25 2.24 / 0.97, 0.81 | 1.93 / 1.91 | 1.96 / 1.92 | 2.00 / 2.02 / 2.00 / 2.00 / 2.09 | FAIL Hopf 0.17 / PASS / PASS | 5/5, 2.2-4.7 / 47-106 | (not hr) | ref |
| hr (prod., xthin 30) | 1.90 **0.69** / 1.70, **-0.01** | 2.24 **1.07** / **0.27, 0.06** | 1.80 / **-0.42** | **-0.79** / 1.34 | 1.06 / 1.06 / **-1.17** / **0.40** / 1.88 | PASS | 5/5, 3-5 / 17-45 | 0.403, 4.33 | FAIL (1) |
| hrx0 (xthin 0) | = cen / 0.68, 0.72 | = cen / 0.77, 0.32 | = cen / 1.98 | 1.70 / 1.92 | (= tA) | PASS | 5/5 but 4.8-9.4 / 58-182, 39-56 % faces superluminal | 0.403, 4.43 | FAIL cost on beams |
| tA `hj` | = cen / 0.69, 0.74 | = cen / 0.79, 0.34 | = cen / 1.98 | 1.70 / 1.92 | 2.00 / 2.02 / 2.00 / 1.91 / 2.09 | PASS (1.68e-2) | 5/5 but 3-6.2 / 56-170 (3-4x hr wall), 39-56 % faces superluminal | 0.406, 4.43 | stable, beams costly |
| tB `kn` (Kn0 1, p 4), tAB `hjkn`, tFB `fskn` | **0.71** | **0.03** | 1.93 | 1.70 | - | - | tB xb20 **DIV c4** | - | FAIL (1) |
| tF `fs` (C1) | **0.33** | 2.24 / - | = cen | 1.70 | - | - | 5/5, 2-5.6 / 6-18 | = hj | FAIL thin pulse |
| tFl `fs+lag` (C4a) | -0.26 | 2.24 | - | 1.70 | - | - | **DIV** xb20 ba0 ba20 shd3b (c1-c7) | - | FAIL |
| tFs `fs+sc` (C4b) | -0.23 | 2.24 | - | 1.70 | - | - | 5/5 but E/J ~ 30-150 | - | FAIL |
| tFh `fs+hc` (ungated) | 0.32 | 2.24 | - | 1.70 | - | - | **DIV c0** xb20 ba0 ba20 | - | FAIL |
| tQ `fsq` (qp 8) | = cen / 1.31 overall | = cen | = cen | 1.70 / 1.92 | = cen | PASS (1.63e-2) | **DIV** xb20 c15, ba0 c42, ba20 c16 | - | FAIL beams |
| tQh `fsq+hc` | = cen / atm 1.98, 1.98 | = cen / 1.98, 2.00 | = cen / 1.77 | 1.70 / 1.92 | = cen | PASS | **DIV** xb20 ba0 ba20 | - | FAIL beams |
| tI `fsi` (+hc tIh, qp 2 tI2) | **0.26** | 2.24 | = cen | 1.70 | - | - | 5/5, 2-4.8 / 6-19 | - | FAIL thin pulse |
| **tE `fse`** | = cen / 0.69, 0.74 | = cen / 0.79, 0.34 | = cen / 1.98 | 1.70 / 1.92 | 2.00 / 2.02 / 2.00 / 1.91 / 2.09 | PASS (1.68e-2) / PASS / PASS | **5/5, 2.0-5.6 / 6-18 (= fs bitwise)** | 0.409/0.421, 4.43 (= hj bitwise) | PASS except atm thin top (0.69) |
| **tEs `fse+hcs`** | = cen / **1.52, 1.81** | = cen / **1.82, 2.06** | = cen / 1.78 | 1.70 / 1.92 | 2.00 / 2.02 / 2.00 / 1.91 / 2.09 | PASS (1.97e-2) / PASS / PASS | **5/5 = tE bitwise** | 0.425, 4.43 (= tE bitwise) | **best** |

Notes on the table:
- **Orders:** the last X order of each ladder (32..512; rw_t1000 32..256). The atm regions are the last order of
  ana_atm_self.py (128 -> 256 -> 512).
- **X ladders = cen** means the L1 numbers are identical to cen's: all pulses be/hd2 XTC, k1280 and rw_t1000.
  pulse_k0.128 hd2 X is 1.75 1.94 2.04 for cen, tA, tQ, tE and tEs.
- **T ladders without predict** equal cen's for every switch arm (hd2 T k0.128 1.57, k12.8 1.02, k1280 1.90,
  rw_t10 1.03, rw_t1000 1.01). hr is worse there.
- **T ladders with predict (torder 8ee1e217):** cen, tA, tQ, tQh, tE, tEh and tEs all give the numbers shown.
  - rw_t10a F1 is 1.87 at the last level (cenPvs 1.97) and dens 2.08.
  - hr stays at about 1 (k1280 -1.17, rw_t10a E 0.40). This confirms the torder finding that the dt-dependent
    xthin weight is what kills hr's time order.
- **atm Hopf L1 at 512:** cen 8.4e-3 (top cell 0.15); hr 4.2e-3; tA/tE 3.72e-4; tEs 3.72e-4 (top 2.8e-5).
  The Hopf floor is the angular quadrature; the region orders above are self-convergence.
- **Gates** (G1 cfl 1e2 / 1e4, G3 c1 / c10, G5 k12.8..12800) for tE / tEs:
  - G1 max|dE/E| 1.68e-2 / 1.97e-2, L1 6.6e-4 / 7.1e-4;
  - G3 L1 0.0143 / 0.0140;
  - G5 k1280 1.0167, k12800 0.99996.
  - k12.8 / k128 FAIL (0.51-0.55) for every arm including cen and hr, as known.
  - G1 cfl 1e4 NON-CONVERGED 34 (implicit_maxit 30 in the gate input): the same for every arm.
- **Beams:** 100 cycles, c dt/dx ~ 300, BiCGStab inner mean per pass.
  - hr: Picard 3-5, inner 17-45.
  - tE / tEs: xb20 5.57/12, ba0 4.00/18, ba20 5.42/16, cyl 4.99/10, shd3b 2.02/6.4. That is fewer inner
    iterations than hr on every beam, and wall 0.6-0.9x hr.
  - The cyl metric matches hr (<|E/J-1|> 0.346 vs hr 0.377), with no face |F| > cE (hj / hrx0: 1.1 %).
  - shd3b edges are as hr (0.236 / 0.343).
  - Collimated beams xb20 / ba0 / ba20: the beam E/J median is 0.21 / 0.68 / 0.11 for tE / tF. For hr it is
    2.9 / 7.1 / 2.4, and for hrx0 1.45 / 7.1 / 1.8; the exact value is 1.
    - Every half-range arm misses J. The idort factor at the thick emitter surface decides the sign of the
      error: hrx0 with alpha = 0, which is not frozen, gives E/J equal to hj's.
    - So this is an emitter-surface flux issue of the half-range blend, not of the new switch. **Open.**
    - The superluminal-face fraction is 0 for fs / fse (hj / hrx0: 39-56 % of all faces, max |F|/cE 2776).
- **GPU AG Car A** (4 A100, 30 cycles, gpudev):
  - On vet_gd (src 2) the formal-solution C1 / eps are not computed, so `fse` and `fse+hcs` reduce to `hj` there.
    The rerun bin data are bitwise equal.
  - s/cycle 0.406-0.425 vs all 0.403-0.407 (run-to-run spread ~4 %). Picard 4.43 vs 4.33; inner 3.21 vs 3.33.
  - Diffs vs all: below 3e-5 inside 0.95 R_ph (m1_e), 0.95-1.05 R_ph velx 0.14 and m1_f2/3 0.1. Above 1.05 R_ph,
    vely/velz 0.6 x max|all| (sq: 130-900x). Full table: gpu/RESULTS_gpu.txt, RESULTS_gpu2.txt, RESULTS_gpu3.txt.

**Key findings:**
1. **`fse` (C1 non-locality weighted by the thermal fraction) + `|H|/J`.** It is dt-free and frozen per solve
   (built at Picard pass 0, so the operator is fixed). It keeps central everywhere the field is diffusive or the
   SC source is the lagged E^n (scattering). It goes upwind where the formal solution's own thermal source shows
   non-local transport (beams, shadows, vacuum). This is the only variant here that is both 2nd order in 1-D and
   cheap and stable on the beams.
2. **Why the plain C1 weight (`fs`) fails the thin pulse k0.128.** That pulse is scattering-dominated
   (S = E^n) and time-resolved (c dt/dx < 1). Any weight that upwinds it gives a pre-asymptotic order of about
   0.3 (fs 0.33, fsi 0.26, kn 0.71, hr 0.69). The quasi-static formal solution of t^n is O(1) wrong there, and
   even the isotropic +-1/4 fallback (fsi) gives 0.26.
3. **J-vs-E^n quasi-steadiness gates (fsq) look dt-free but fail the beams.**
   - From t^n formal-solution data alone, a time-resolved thin transient (J << E^n) and a steady shadow (J ~ 0,
     E^n > 0 by numerical diffusion) look the same locally.
   - Gating them central at c dt/dx = 300 diverges (xb20 / ba0 / ba20).
   - Only the step's light-crossing ratio (dt) separates them physically, and the brief forbids dt. The thermal
     fraction is the dt-free proxy that works for these tests.
4. **The Knudsen number of E^n (kn) is a poor thinness measure.** |grad ln E|/chi is large in the Gaussian
   tails of a diffusive pulse, which makes pulse_k12.8 order 0.03.
5. **C4, the formal-solution flux in the thin part.**
   - **lag** (C4a, F_sc - F_hr(E^n) as a fixed source) diverges on 4/5 beams. Its explicit -F_hr(E^n) is an
     explicit transport term at c dt/dx 300.
   - **sc** (C4b) is stable but wrong (E/J 30-150; atm order 0.13).
   - **hc**: the SC's own central-minus-upwind face flux, gated by q = min(E^n/J, J/E^n)^2. Ungated it diverges
     at cycle 0 (light fronts). Gated, it gives 2nd order in every atm region, but still diverges in the thermal
     beams with fse/fsq.
   - **hcs** (hc times the scattering fraction) keeps the atm benefit and is inactive in the thermal beams.
   - **C4c** (SC flux scaled by the iterate) is the existing half-range flux itself: c (h+ E_L + h- E_R) with
     h = H/J.
6. **Freezing is cheap and harmless.** Faces are built at pass 0 of each solve. hrx0 with alpha = 0, rebuilt
   every pass, gives the same beam metrics as hj. The new arms' Picard counts on the beams are at or below hr's.

**Caveats / to do:**
- fse relies on the thermal fraction. A thin beam in a SCATTERING medium falls back to hj: stable, but with
  hj's beam cost (3-4x hr) and superluminal dark faces. This was not tested; there is no scattering beam test in
  the battery.
- On vet_gd (sph, AG Car) only |H|/J is used (no C1, no eps, no correction). It is cheap there, but the C1
  benefit would need the sum w|I - S| and eps in the vet_gd sweep (rad_m1_vetgd.cpp VetGdHalfRange).
- C1 is accumulated in the single-block VetShortChar and the multi-block VetSweepMB cell kernels only, not in
  the banded VetMomAdd path (vet_mb_mom_fuse). There it stays 0, and fse degrades to hj.
- The atm thin top with tEs is erratic (0.74 0.60 1.52). With tE it is about 0.7. The criterion asks >= ~1, so
  this is borderline.
- rw_t10 X last order 1.70 vs cen 1.96 is common to every half-range arm with w = 0 (tZ: hj with h0 = 10,
  identical). It is a baseline difference of the blend path, not the switch.
- The time order needs torder's time2_vet_sc = predict (cherry-picked here as 10d36189). Without it every
  arm's T ladders equal cen's lagged ~1 values.

## Commits (thinsw-1010, local)
- 4676ba1e implicit_thin_switch (hj | kn | fs | hjkn | fskn), implicit_thin_corr (lag | sc), M1_VET_NL sweep sum (C1)
- 7a106224 fix: ifw2/ifw3 carry the DG slot (was written out of bounds -> heap corruption at exit; found with a
  Debug/bounds-check build)
- 066118ca implicit_thin_corr = hc
- cab4b34a quasi-steadiness gate q on lag / hc
- e9014976 implicit_thin_switch = fsq
- d868d46c implicit_thin_switch = fsi
- 10d36189 cherry-pick of torder 8ee1e217 (time2_vet_sc = lag | predict | rebuild)
- 5b4d0e8f implicit_thin_switch = fse (M1_VET_EPS)
- 9d8fd8c9 implicit_thin_corr = hcs
- c3ea5bec comment line length (cpplint: no new errors in the touched files)

Default bitwise: implicit_thin_switch = none is the default. 4676ba1e and 9d8fd8c9 reproduce
athena_cpu_e50ec7d6 bin and tab data bitwise on pulse_k0.128 (hr arm), ba0 (hr and cen, 20 cycles) and G1 atm2d
(hr, cfl 1e2) (verify/, verify2/). GPU: all_1 vs all_2 bitwise. The new keys are read with GetOrAdd, but a fresh
input that names them on the COMMAND LINE needs them in the file: fixkeys.py does that for the battery, and the
GPU input copy has `implicit_thin_switch = none`, `implicit_thin_corr = none`. Other keys are read only when the
switch is on: implicit_thin_h0/h1, kn0, knp, fsa, fsj, qp.

## Runs / scripts
- CPU battery copy: cpu/ (battery_raven.sh with prep_fast, od.py with the t* arms, fixkeys.py, run.sbatch).
  - Run tree: cpu/run (runs/, beams/, gates/, RESULTS/: order_eval_full.txt, order_eval_fse.txt,
    order_eval_tEs.txt, order_T_predict.txt, atm_self_regions.txt, atm_hopf.txt, collref_*.txt, fast_eval1.txt).
  - cpu/run_bug: the first round, run on the out-of-bounds binary 4676ba1e. It is invalid; do not use it.
  - All runs: Slurm general, 72 cores.
- Experiment exp1/: freeze-vs-idort check (ba20 / xb20, 30 cycles).
- GPU: gpu/gpu_arms_raven.sh (ARMS hj, hjkn, hrx0, fse, fses, all), gpu/inp/agcar_A_ths.athinput. Jobs
  31054015, 31055398, 31055770.

## What was tried, in order (log)
1. hj / kn / fs / hjkn / fskn plus C4 lag / sc, all frozen at pass 0.
   - The first binary wrote DG out of bounds into ifw2/ifw3: rc 139 at exit on the 512 levels. A Debug build
     found it. Fixed in 7a106224 and everything rerun.
   - Results: kn fails k12.8, fs fails the thin pulse, lag diverges, sc is wrong. hj = hrx0: correct orders,
     costly beams.
2. hc (ungated): diverges at cycle 0 on beams (light front, J >> E^n). Adding the q gate made it stable in
   scattering media but not in thermal beams.
3. fsq (C1 gated by J~E^n): 1-D orders = cen, but beams diverge (shadows look like transients).
4. fsi (isotropic h+- fallback where J != E^n and not beamy): beams fine, but the thin pulse is still 0.26 and
   atm worse. This shows the thin-pulse failure is upwinding per se, not the angular shape.
5. Cherry-picked the torder fix: every switch arm gets T 2 with predict; hr does not.
6. fse (C1 times the thermal fraction): everything passes except the atm thin top (~0.7).
7. fse + hcs: atm thin top / intermediate become ~2nd order; beams and GPU unchanged (bitwise).


## Design (written before coding)

### Where w enters (existing code, src/rad_m1/rad_m1_implicit.cpp at c696c7d1)
- Every face of the all-faces half-range blend carries AL = w, HCL = w ch h+_L, HCR = w ch h-_R
  (M1HrFace, l.256; x1 faces in the `m1_impl_aphll` kernel ~l.9988, x2/x3 faces in
  ImplicitLatFaceCoef ~l.2837). The face flux is (1 - w) F_central(E') + c (h+_L E'_L + h-_R E'_R) w,
  F_central the face-eliminated staggered flux (theta = 1 - AL on transverse faces).
- h+- = H+-/J of the vet_sc formal solution (vet_cell M1_VET_HP1+d, M1_VET_H1+d, M1_VET_J;
  M1HrCell src 1) or vet_gd (vgd_hr, src 2). |H|/J = fb.
- The x1 face also has a lagged additive slot M1_IFW_DG (the plm deferred correction): it enters
  only the right-hand side (rr -= nu DG) and the reconstructed face flux. The x2/x3 face arrays
  ifw2/ifw3 have the same slot, unused so far.
- Existing weights: idort_f w = 1/(1 + x + x^2/tau0) * smoothstep(fb; 0.3, 0.6), x = tau_f
  (tau_cell-based; -> smoothstep(fb) as dx -> 0), then the xthin override
  w <- 1 - (1 - w)(1 - X^2/(X^2 + X0^2)), X = (c dt/dx)/(1 + c dt chi_f) (dt dependent).

### New key: <rad_m1>/implicit_thin_switch = none (default, bitwise) | hj | kn | fs | hjkn | fskn
Needs implicit_flux = blend + implicit_flux_beam = halfrange (+ faces = all on multi-D).
When on, implicit_blend / idort / xthin are ignored on the half-range faces (no dt anywhere):
- hj: w_A = S((g - h0)/(h1 - h0)), g = max over the two cells of fb = |H|/J (formal solution),
  S the smoothstep 3s^2 - 2s^3 clamped to [0,1]. Defaults h0 = 0.3, h1 = 0.6
  (implicit_thin_h0/h1). Exactly 0 below h0.
- kn: w_B = K^p/(K^p + Kn0^p), K = |E^n_R - E^n_L|/(tau_f max(E^n_L, E^n_R)) the face Knudsen
  number of E^n (M1_IW_EN, the start-of-step state of the solve), tau_f = chi_f dx_f. tau_f = 0
  (vacuum) or max E^n = 0 (dark) gives w = 1 (nothing diffuses: upwind, realisable).
  Defaults Kn0 = 1 (implicit_thin_kn0), p = 4 (implicit_thin_knp).
- fs (coordinator C1, formal-solution non-locality): g = max(fb, C1/1.5, fsa * an, fsj * dJS) with
  - C1 = sum_rays w |I - S| / J (new vet_cell slot M1_VET_NL, accumulated in the vet_sc sweep; in
    the diffusion limit C1 = 1.5 |H|/J, hence the 1/1.5 so one threshold pair serves both);
    C1 also flags crossing beams (H ~ 0, I far from S) and vacuum (S = 0: C1 = 1);
  - an = |D - I/3|_F / sqrt(2/3) (raw K/J anisotropy, 1 for a single beam), weight implicit_thin_fsa (default 0);
  - dJS = |J - S|/max(J, S), weight implicit_thin_fsj (default 0; a fallback where C1 is not swept);
  w = S((g - h0)/(h1 - h0)).
- hjkn / fskn: w = max(w_A, w_B).
All of fb, C1, an, dJS, E^n are fixed through the step (the VET rays are built from E^n), tau_f
moves only with implicit_opac_update. **Frozen operator**: in the new modes the face coefficients
(x1 and x2/x3) are built at the first Picard pass of each solve and reused, so the Picard loop
sees a fixed linear transport operator (no iterate-dependent switch).

Resolution behaviour: fb, C1, an, dJS and K are continuum quantities of the field at t^n (up to
the angular quadrature), independent of dx and dt. In a thick/intermediate diffusive region
fb ~ Kn/3 < h0 => w = 0 exactly => the compact central flux, 2nd order. Where w > 0 (Kn >~ 1:
transport regime, the "thin" part by the criterion's own definition) w is O(1) and fixed =>
1st-order upwind error there.

### Coordinator C4 (formal-solution flux in the thin part), key implicit_thin_corr = none | lag | sc
Conservative (face flux form):
- lag (C4a): DG = w ch (Hf - h+_L E^n_L - h-_R E^n_R): the face flux becomes
  (1 - w) F_cen(E') + w [F_hr(E') + (F_sc^n - F_hr(E^n))], a FIXED source in the implicit solve.
  Hf = (H_L + H_R)/2 the face mean of the formal-solution flux (E units), clamped to [-E^n_R, E^n_L]
  (donor admissibility). Note F_hr(E^n) = c (H+_L E^n_L/J_L + H-_R E^n_R/J_R): the existing
  half-range flux IS the SC flux scaled by the iterate (coordinator C4c) with an upwind face value,
  so C4c needs no new code; C4a adds the central (2nd-order) face value plus the J/E^n mismatch.
- sc (C4b): HCL = HCR = 0, DG = w ch Hf: pure lagged SC flux where w = 1.
vet_sc only (src 1); with vet_gd the correction is 0.
Time: F_sc^n is the quasi-static formal solution of the t^n sources (lagged one step, and no
light-travel time): O(dt) in time, O(1) wrong in a time-resolved thin transient.
