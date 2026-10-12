# NOTE sporder: the remaining hesdirk2 time-order floors on the spherical wedge (vet_gd)

Worker X, 2026-10-12, MPCDF Raven. Work dir /raven/ptmp/jinma/sporder_1012. Branch sporder-1012 (worktree
/raven/u/jinma/ATHENAK/wt-sporder) from prodfinal-1011 4253af30. No push.

## 0. Setup
- Binary bin/athena_cpu_4253af30 (copy of prodfinal_1011/bin).
- Harness sp/ana (copy of suiteF sp/ana): arm `HFfZ` = the candidate with every key named
  (hesdirk2, time2_vet_col predict, implicit_vimp, time2_vstage, xthin 0, fse / hc / dc / cell, vet_bc_diff face,
  vet_seg emit, time2_pred_cap 0.5, conservative2, cs2_shift false, cs2_hll drag, vimp_bc true, time2_start be).
  Overrides `HFfZ+tag[+tag]` (spl.VARS): r0 = vet_gd_rotate_every 0, f0z = implicit_f0_init zero, sl = slope start,
  nvs / lag / nts / c1 / m / tol / tw0 / nb ...
- tprof.py / xprof.py: per-radius T / X difference profiles.

## 1. Findings so far (fresh T ladders l = 1..16, orders over 1-2 / 2-4 / 4-8 / 8-16 shown for the last three)

### 1a. rswv (radial shell + v_r = 0.01 sin), conservative2
| arm | dens | velx | m1_e | m1_f1 |
|---|---|---|---|---|
| HFfZ (candidate) | 1.27 1.15 1.18 | 1.81 1.58 1.42 | 1.77 1.58 1.44 | 1.78 1.68 1.25 |
| + rotate 0 | same | same | same | same |
| + implicit_f0_init zero | **2.02 2.01 1.97** | 1.92 1.74 1.59 | 1.92 1.74 1.45 | 1.87 1.85 1.34 |
| + f0z + r0 | 2.01 2.01 1.96 | 1.90 1.78 1.56 | 1.91 1.77 1.42 | 1.90 1.78 1.56 |
| + f0z + r0 + tol 1e-14/1e-15 | identical to f0z + r0 | | | |
| + f0z + r0 + lag / thin_switch none | identical (within 0.01) | | | |
| + f0z + r0 + time2_vstage false | 1.63 1.43 1.24 | 1.29 1.16 1.08 | 1.36 1.21 1.11 | 1.00 1.00 1.00 |
| + f0z + r0 + cs2_dbg 1 / matched | 2.0 1.3 0.4 | 1.6 0.6 0.3 | 1.5 0.6 0.4 | 1.6 0.5 0.4 |
| + slope start + r0 | 1.32 1.20 1.20 | 1.88 1.79 1.64 | 1.88 1.80 1.69 | 1.89 1.79 1.63 |
| Cartesian twin rswv_cart (vet_sc) | 1.31 1.16 1.21 | 1.88 1.76 1.45 | 1.85 1.75 1.59 | 1.89 1.85 1.76 |
| rswv_cart + f0z | 2.02 2.03 1.86 | 1.91 1.77 1.26 | 1.93 1.85 1.75 | 1.92 1.86 1.80 |

- rho floor = start-up: implicit_f0_init auto sets F0 = F - aE from the pgen's F = 0 under a moving gas (radiation
  at rest in the lab), a stiff initial layer (1/(c kappa) ~ 1e-3 < dt). `zero` (radiation co-moving, the well-prepared
  IC) gives rho 2.0. The slope start does not cure it.
- Remainder (v, E, F1 1.4-1.8 at l 8-16, errors ~4e-7 rel): not rotation, not tolerance, not vet_gd (the Cartesian twin
  has it too); the post-solve split is far worse. Same family as rlat's 1.8 (in-solve E part of dA; worker T).
- rswv X (dens 0.13, v 0.7): a physical boundary layer at the outer wall (gas piles against the reflecting wall,
  the heated radiation leaks through the Marshak top: width sqrt(c t/(3 kappa)) ~ 0.026 = 3 cells at N 256).
  Cells r < 2.75 only: dens 1.78, v 1.95, E 1.99, F1 1.95 (r < 2.5: 1.80 / 1.97 / 2.04 / 1.97; r < 2.9 already
  takes in the layer: 0.9-1.0).

### 1b. rswt / rswtv (tau 20 shell, top bath)
- rotate_every 10 is most of the floor: F1 1.0/0.0 -> 2.39/0.14 (finest error 3e-5 -> 2.7e-6), E 0.84/0.19 -> 1.61/0.52.
  Rotating at the same physical times on every level (rotate_every 10 l) = rotate 0.
- The rest of the floor is vet_gd-specific: the Cartesian twin rswt_cart_T is 1.8-1.9 in every variable.
- No effect: thin_switch none, vet_gd_twin off.

### 1c. rlatt_om5 (thin lateral wave)
- rotation 0: errors 10-100x lower but a floor remains (E 0.87 0.69 0.74, F1 0.98 0.54 0.69).
- rlatt X F1 0.6: F1 is a tiny quantity (mean |F1| 1.2e-5 vs E 3e-2); its absolute X error is at the T-floor level.

### 1d. Marshak (marsh_r1_T): flat in the Cartesian twin too (marsh_cart_T 0.2-0.5); rotation / thin switch no effect.

## 2. The relaxed VET source of the predicted build (the vet_gd thin floor)

- `time2_vet_col_source` default `relaxed`: the predicted build's source relaxes P (~ t^{n+1}) by a further dt
  (VcolRelaxedSource / VlatRelaxedSource), an O(dt) source error where the exchange is resolved (vcpred-1010).
- With rotation off, `state` makes the thin tests 2.0:
  - rswt_T r0+state: dens 2.41 2.19 1.81, v 2.02 1.93 1.80, E 2.01 2.00 1.97, F1 2.05 2.04 2.07 (relaxed r0: E 1.61 0.52 0.64);
  - rlatt_om5_T r0+state: dens 2.05 2.01 2.01, v 1.92 1.96 1.97, E 1.90 1.97 2.03, F1 1.92 1.96 1.98, F3 1.81 1.91 1.94
    (relaxed r0: E 0.87 0.69 0.74, F1 0.98 0.54 0.69);
  - rswv (thick, f0z r0): state = relaxed (the source is irrelevant there).
- `state` on reduced CPU AG Car A (candidate deck, conservative2, 120 cycles, job 31091608): Picard 8.57 (max 200),
  2 NON-CONVERGED, 2.3x cost, v changes O(1): still unusable (as vcpred-1010, pred_cap does not save it).
- New opt-in keys (sporder-1012, default bitwise: rswv_T HFfZ level 2 tab IDENTICAL; AG Car R = base all bins identical):
  - `time2_vet_col_source = extrap`: 2 S*(dt) - S*(2 dt), clamped to [S*(0), S*(dt)]. AG Car: Picard 3.03 = relaxed,
    max change 7e-8 (v), 2e-12 (E). rswt r0: E 1.73 1.21 1.33 only (gas-T relaxation x = 1.25/l: pre-asymptotic).
  - `time2_vet_col_source = blend` (+ `time2_vet_col_source_x0`, default 1): S*(0) + w (S*(dt) - S*(0)),
    w = x^2/(x0^2 + x^2), x = chat dt rho kappa_E + c dt rho kappa_P 4 a T^3/((1+x_E) rho c_v): the state's source
    (O(x^3) from it) where the exchange is resolved, the relaxed source where stiff. Results in sect. 3.

## 3. blend source and time-based rotation (binary xs3 = blend + extrap + my rotate_dt; the commit 9164fea8 keeps
## ONLY blend (rotate_dt dropped for worker Z's vgdrot2-1012 a78e1977, same key); bin/athena_cpu_bl1 = 9164fea8 is
## bitwise equal to xs3 on rswt_T r0+blend levels 1-2; defaults bitwise: rswt_T / rswv_T level-2 tabs IDENTICAL to
## 4253af30, reduced AG Car R = 4253af30 all bins identical)

| case (orders 1-2 2-4 4-8 [8-16 / 16-32]) | E | F1 | v | dens |
|---|---|---|---|---|
| rswt_T candidate | 0.84 0.19 0.28 | 1.00 0.00 0.69 | 1.79 0.50 0.73 | 2.46 1.52 0.55 |
| rswt_T + rotate 0 | 1.61 0.52 0.64 0.82 | 2.39 0.14 0.44 0.73 | 1.86 0.87 0.61 0.80 | 2.47 1.84 0.86 0.81 |
| rswt_T + rotate 0 + state | 2.01 2.00 1.97 1.92 | 2.05 2.04 2.07 2.13 | 2.02 1.93 1.80 1.68 | 2.41 2.19 1.81 1.45 |
| rswt_T + rotate 0 + extrap | 1.73 1.21 1.33 1.60 | 2.41 1.09 1.16 1.51 | | |
| rswt_T + rotate 0 + blend | 1.88 2.02 2.03 1.96 | 2.23 1.95 1.84 1.92 | 1.92 1.91 1.81 1.68 | 2.45 2.19 1.77 1.42 |
| **rswt_T + blend + rotate_dt 10 dt0** | 1.88 2.00 1.96 | 2.18 1.78 1.44 | 1.92 1.91 1.81 | 2.43 2.16 1.73 |
| rswt_T + blend (rotate 10 cycles) | 0.99 0.19 0.30 | 1.00 0.02 0.65 | | |
| rlatt_om5_T candidate | 0.68 0.10 0.60 | 1.16 0.55 0.41 | 1.36 0.90 0.69 | 1.37 0.73 0.48 |
| rlatt_om5_T + rotate 0 + state | 1.90 1.97 2.03 | 1.92 1.96 1.98 | 1.92 1.96 1.97 | 2.05 2.01 2.01 |
| rlatt_om5_T + rotate 0 + blend | 1.51 2.11 2.42 | 1.54 2.14 2.28 | 2.32 2.01 1.29 | 2.29 2.01 1.87 |
| rsw_T + blend + rotate_dt | 1.99 1.98 1.97 | 2.00 2.00 2.00 | 2.00 2.00 2.00 | 2.04 2.08 2.14 (= candidate) |

- Reduced CPU AG Car A (candidate deck agc/inp/A0c = prod/agcar_A_fresh + conservative2 / cs2_shift false / cs2_hll
  drag / vimp_bc true, 480x32x32, 120 cycles; jobs 31091608, 31092677):
  relaxed Picard 3.029 (max 6); state 8.57 (max 200, 2 NONCONV, 2.3x cost, v O(1)); extrap 3.029, max change 7e-8 (v);
  **blend 3.029, max change 6.9e-8 (v), 3e-14 (E), cost +0.5 %**.

- rotate_dt rows above used my (dropped) vet_gd_rotate_dt = 10 dt0 (block of t^n); Z's vgdrot2 implements the same key
  (+ rotate_t0, rotate_align). With it + blend, rswt_T E 1.88 2.00 1.96, F1 2.18 1.78 1.44 (vs rotate-every-10:
  E 0.99 0.19 0.30, F1 1.00 0.02 0.65; finest F1 error 3.3e-5 -> 4.6e-7).

## 4. Other items

- **rlatt X (phi ladder, radial grid fixed 16 cells, dt ~ dphi)**: F1 0.46 / v 0.62 with rotation off and blend (E 2.02,
  F3 1.93, dens 2.40). Not the twin (tw0 bitwise-equal orders), not the thin switch (nts equal). 432 directions (Q):
  F1 0.70, v 0.82, errors -30 %: the fixed vet_gd direction set / lateral interpolation (V's vgdlat: first order in the
  lateral cell, ray effects) sets F1, a tiny quantity (mean |F1| = 4e-4 of c E). Thin region (tau 10 over the shell).
- **rswv X**: wall pile-up layer at the outer reflecting gas wall under the Marshak radiation top (last cell rho 1.062 vs
  1.021). Order ~0.6 out to N 1024 (rswv_cart the same: 0.5), interior r < 2.75: 1.8-2.0. Not the time scheme.
- **rswv T remainder** (f0z): v/E/F1 1.90 1.78 1.56 1.28 (l to 32), Cartesian twin the same, restart ladder the same;
  not the source (state/blend equal), not cs2_dbg 2, vet_seg, vet_bc_diff, tolerance. matched + implicit_rad_inertia
  false: dens 2.03 2.06 2.10, v 1.97 1.94 1.89. v scan: v 0.001 1.95 1.89 1.75, 0.01 1.90 1.78 1.56, 0.02 1.74 1.49
  1.15; kappa_s 100 / 1e4 similar. -> the radiation-inertia (dA) treatment under conservative2 (worker T, rlat).
  implicit_vimp false / vimp_bc false: 1.0 (they are needed).
- **Marshak T (spherical and Cartesian, hesdirk2 and BE)**: flat (H 0.2-0.4, BE 0.5-0.6 dens, 0.2-0.3 E) to l = 64.
  Not: rotation, thin switch, pred_cap 0, slope start, tol 1e-14, time2 source (state/blend), gas_newton off, hr_pos
  bound, pos_floor off, blend gate (idort without f, fmode mean: bitwise equal), outflow hydro walls, warm start
  (T0 0.3), smooth Gaussian start in equilibrium; the order is already low at t = 0.0125 and the differences decay
  with time. BE itself is below first order, so a non-smooth ingredient of the discrete Marshak problem shared by
  both schemes (candidate: the T = (E/a)^(1/4) equilibrium manifold at the cold foot; not resolved). OPEN.
- Restart ladders (rlad) are not usable for order: the restarted dt sequence differs from dt0/l (marsh: 804 steps for
  768) and vet_gd restarts carry a cold-start term (rswt_T restart 0.2-0.7).
- rswv fresh-start rho: implicit_f0_init zero (well-prepared IC) is the test-side fix; the pgen F = 0 with
  f0_init auto puts the radiation at rest in the lab under a moving gas (stiff initial layer).

## 5. blend regression (bl1 = 9164fea8, candidate HFfZ vs HFfZ + blend)
- relax1/relax2 early T (vs the ODE): both at round-off (1e-12 - 1e-13), identical behaviour.
- patmT L16 poly (P1T, steady coupled): dens identical, E/v differences <= 2e-9 (vely is a zero quantity).
- rlat_T (thick lateral): 1.83 1.82 1.74 both arms; differences <= 1e-7.
- pinc (P2, held gas): bitwise identical (Picard 6.59, the known 1 NONCONV in both).
- rsw_T: identical orders (2.0).
- Reduced CPU AG Car A: Picard 3.029 both, max change 6.9e-8 (v), 3e-14 (E).

## 6. Recommendation
- Production candidate: add `time2_vet_col_source = blend` (x0 1). With Z's `vet_gd_rotate_dt` (vgdrot2-1012) the thin
  top-bath / lateral floors (rswt, rlatt) become ~2 at production-like rotation.
- For the rswv family use `implicit_f0_init = zero` in the test IC (well-prepared); the remaining 1.3-1.8 is the
  in-solve radiation inertia under conservative2 (worker T).
- Marshak T flat and rswv / rlatt X are not the time scheme (see sect. 4); Marshak stays OPEN.

## Files
- sp/ana: spl.py (arm HFfZ + VARS tags incl. bl, st, xs, r0, rl, rd, f0z ...; NEWTAGS -> SPT_BIN2), spt.py (case tokens
  vNN / kNN / wNN / sm / ho, rlad), tsum.py (order summary), tprof.py / xprof*.py / tdump.py (profiles).
- sp/jobs_*.txt, sp/raw_*.txt (job lists), logs/, sp/runs/<case>/<arm>_<lev>; agc/ (reduced AG Car: agc.sh, inp/A0c,
  cmp.py; runs/j.31091608 relaxed/state/extrap, j.31092677 relaxed/blend).
- bin: athena_cpu_4253af30 (base), xs1/xs2/xs3 (scratch builds), **bl1 = 9164fea8** (+ he_cpu).
