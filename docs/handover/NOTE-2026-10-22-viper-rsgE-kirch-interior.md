# NOTE 2026-10-22 viper -> caltech (USER): RSG E N=32: kirchhoff fixes the cliff but NOT the seam imprint; a cube-locked interior flow from t=0 feeds the imprint

Runs: RSG E deck, cs N=32, f55be371 / kirch-1022 (61580891), viper; restart arms from our base rst at 40.0 d. Scripts scr/imprint.py,
scr/radfe.py, scr/interior.py in /viper/ptmp2/jinma/photinst_1022 (logs/imprint_*.txt, logs/interior_base.txt).

## 1. kirchhoff restart 40 -> 59 d vs base
- Cliff fixed: cliff-cell e/rho sd 0.015 at 58 d (base 0.21); 0 dt retries (base: retries from 56.1 d, dt 1000 -> 447 s); dt 1000 s;
  no measurable cost penalty on MI300A so far.
- Panel contrast NOT reduced at equal time (full field, FC |x|,|y|<0.3 vs ED max>0.8):
  | t | L(1.2-1.4 R_ph) fc/ed-1 kir / base | T_eff fc/ed-1 kir / base |
  |---|---|---|
  | 58 d | +3.1e-3 / +3.6e-3 | +3.1e-4 / +2.6e-4 |
  | 59 d | +3.3e-3 / +5.1e-3 | +3.3e-4 / +4.0e-4 |
  Lateral L sd at 1.3 R_ph larger with kir (2.4-2.6 % vs 1.3-1.7 %). With kir the contrast grows smoothly from +7e-4 at 41 d (e-fold
  ~12 d), without the cliff noise. Base after 56 d: L contrast 6e-3 (60 d), 1.3e-2 (67 d), T_eff contrast 1.0e-3 (68 d), retries.
- Metric caveat: removing l <= 8 before fc/ed also removes the cube pattern (cubic harmonics l = 4, 6, 8); use the full-field fc/ed.

## 2. Source: a cube-locked interior flow present from t = 0 (base; identical to 3 digits in the kir arm below ~0.95 R_ph)
Per shell z = (mean_FC - mean_ED)/lateral sd (a purely cell-geometric pattern gives z ~ +0.20-0.21 with these masks; random convection ~0).
- t = 0: Trad and e/rho have z = +0.20 on every shell 0.3-0.95 R_ph (the IC's lateral variation is purely geometric).
- From 3 d: v_r has z = +0.21 and rms(v_r)_FC/rms(v_r)_ED - 1 = +0.30 on every shell: face centres move radially ~30 % faster, mean
  outflow at face centres and inflow at edges (a cube-locked l = 4-like circulation).
- Amplitude grows secularly: v_r rms at 0.9 R_ph 1.7e2 (3 d), 5.6e2 (12), 1.1e3 (21), 2.0e3 (30), 3.2e3 (39 d) cm/s; at 0.3 R_ph
  1.5e3 -> 3.3e4 cm/s by 68 d.
- Where resolved convection develops z drops toward 0, then the cube pattern RE-EMERGES at 0.95 R_ph (z_vr 0.17 at 54 d -> 0.21 at
  60-68 d; rms ratio 0.03 -> 0.33), exactly when the L panel contrast passes 1e-3 (55 d); same in the kir arm.
- Not the IC seed (he_seed: random-sign Y_lm, l 1-8 -> random-sign z), not the cliff.
Reading: a spherically symmetric hydrostatic state is not exactly in equilibrium on cs cells (pressure gradient vs gravity /
geometric source terms / cell volumes and face areas differ between face-centre and edge cells); the residual drives a cube-locked
circulation that later shows up as the seam imprint. Onset ~55 d at N=32 vs ~49 d in your N=128 prodF: similar at 4x lateral
resolution, consistent with a cell-shape-set amplitude.
So kirchhoff (or harmonic) fixes the column-local cliff runaway (noise + retries) but on its own it will not remove this imprint.

## 3. Viper next (USER GO)
- Static test: a spherically symmetric HSE star on the cs mesh, hydro only (no radiation, no seed), N=32 and 64, with and without the
  well-balancing options: measure the cube-locked v_r growth. If it reproduces, the fix is a well-balanced HSE that is exact on cs cells.
- Harmonic arm (d3df1401): key-off gate PASS, key-on C==C2 and R==C bitwise; restart arm 40 -> 60 d starting.

## 4. Ask of Caltech
Please run the same interior metric on your N=128 RSG E prodF (and HG) dumps from t = 0 (copy scr/interior.py from viper
/viper/ptmp2/jinma/photinst_1022/scr, or the definition above: per-shell z of v_r, Trad, e/rho and rms(v_r)_FC/_ED - 1, FC/ED masks
as above) to see whether the cube-locked v_r pattern is present at N=128 and how its amplitude scales with N.
