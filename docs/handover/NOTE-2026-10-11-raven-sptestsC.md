# NOTE sptestsC: spherical-wedge (vet_col + vet_gd) suite, matched vs conservative gas source

Worker C, 2026-10-11. Work dir /raven/ptmp/jinma/sptestsC_1011.

- Binary: `bin/athena_cpu_530c38e2`, a copy of prodfinal_1011/bin (md5 656e4177...). It is prodfinal-1011 530c38e2:
  consrc merged opt-in, and `vet_seg` defaults to `emit`.
- Harness: a copy of sptests_1011/ana (spt.py, spl.py, s3_eval.py, patm_eval.py, pin_eval.py, lr_eval.py) plus new
  code:
  - ana/cmpmk.py: bitwise and difference comparison, X/T self-convergence, totals of E, momentum and L_z.
  - ana/evalall.sh: runs every evaluation. Outputs are in ana/out/*.txt.
  - Run tree: runs/.
- Arms (all fse+hc with the production keys):
  - `HFf` = GD + idort_f blend + `implicit_blend_xthin 0`, `implicit_thin_switch fse`, `implicit_thin_corr hc`, plus
    `time_scheme hesdirk2`, `time2_vet_col predict`, `implicit_vimp`, `time2_vstage`, with `vet_bc_diff face`.
    vimp and vstage are dropped where gas_feedback = false, as in sptests.
  - Suffix `M` = `implicit_gas_src = matched`, `K` = `conservative`. Both are named in the deck.
  - `Ff` = the be variant (T ladders only). `HCFf` = the Cartesian vet_sc analogue, used for the large-radius
    controls.
- Every input sets `ic_rad = none`.
- No hr or cen arms (the old numbers are in sptests_1011/NOTE-sptests.md).
- Jobs: sCsmk 31078197, sCA 31078204, sCB 31078255, sCC1 31078444, sCC2 31078592, sCC3 31079072. All ran on
  `general -c 72`, none waited more than a few minutes, so no Caltech TASK was needed.

## 1. Verdict

1. **Held gas (gas_feedback = false): bitwise.** These are P1 patm, P2 pinc, L (pulse r0/w 1 and 1000, nlte2 r1000,
   and the Cartesian twins), T4 dif plus its T ladder, T6a thin (4x4 and 16x16), T6b core (4x4 and 16x16), and T7
   shd. M and K give byte-identical tab/bin data at every level (ana/out/bit.txt), as designed. The exact-solution
   numbers below therefore hold for both arms.
2. **Coupled gas in static equilibrium: identical to round-off.**
   - patmT (photosphere, radiative equilibrium): M vs K differ by <= 2e-12 in E, rho and eint.
   - T1 relax: identical to the printed digits.
3. **Coupled gas with radiation-driven motion: K differs from M by O(E/rho c^2), and K is the physically right one.**
   - rsw is the radiation-dominated shell with E/(rho c^2) = 0.03. There M and K differ by 3.1 % L1 in v and 0.05 % in
     E and rho, at every N (64-256).
   - The difference is the same in the Cartesian vet_sc twin and in sp r0/w = 1000, so it is not a vet_gd effect.
   - The wave-front shift is M - K = +0.0049 (N 128, t = 2). The prediction of the radiation inertia
     (E + P_r)/(rho c^2) = 0.04 in the sound speed is +0.0047: c_s 0.12180 -> 0.11944.
   - So `conservative` carries the radiation inertia, and `matched` does not.
4. **Momentum conservation, new lateral test rlat_om5** (section 4):
   - Total L_z (gas + radiation) drifts by **-4.1 to -5.2 % of the radiation L_z with matched**. The amount is
     independent of N (32-128), of be vs hesdirk2, and of thick vs thin. It is set at start-up, when the gas drags the
     radiation up to F = 4/3 E v and the gas over-pays.
   - **With conservative the drift is +5e-6** (thick, closed in L_z) and **-1.3e-4** (thin, where L_z can leave
     through the Marshak top).
5. **Orders: X unchanged, but K lowers the hesdirk2 time order.**
   - X: K is equal to M or slightly better: rsw interior 2.0, lateral rlat 1.85-1.9, both arms. K's L1 error is 2-12 %
     smaller.
   - T, be: K = M, order 1.0 everywhere.
   - **T, hesdirk2 on rsw:**
     - M gives dens 2.02/1.90, v and F1 1.73/1.23.
     - **K gives dens 1.66/1.47, v and F1 1.39/1.18**, with a 2-3x larger dt error (v: 1.8e-5 vs 6.7e-6 at the finest
       pair).
     - This is the consrc rw_t10a F1 drop, reproduced on the production path, and it affects every variable.
   - On rlat and relax, H is the same for M and K.
6. **Pass/fail (fse+hc production keys, both arms):**

| stage | test | result (both arms unless noted) | verdict |
|---|---|---|---|
| P1 | patm L16 poly, photosphere | E/J-1 at R_ph -0.47 / -0.19 / -0.014 % (N 64/128/256); band 1-10 order 2.3/1.2; max\|F r^2-1\| 3.8e-4 (N 256); env +0.45/+0.88/+1.14 % (lateral 16x16 floor, known) | PASS (M = K bitwise) |
| P1T | patmT L16 poly (coupled gas) | E/J-1 at R_ph -0.52 / -0.16 / +0.04 %; dens/eint X order 2.0; M-K <= 2e-12 | PASS |
| P2 | pinc thermal channel beam, c dt/dx 300 | L/L_ex 1.50 / 1.53 / 1.56 at r 4.3/4.6/4.9 (be was 2.25-2.35); L_fs/L_ex 0.53; 1 NONCONVERGED step (Picard 200 at cycle 82), mean 6.6 | known excess, PASS stability; NONCONV is new with hesdirk2 |
| L | pulse r0/w 1 / 1000 vs Cartesian | L1 8.8e-2 (true focusing) / 1.4e-3 (solver floor, as before) | PASS (bitwise M = K) |
| L | nlte2 r0/w 1000 (vet_col_source gas) vs exact | per band -0.26 / -0.03 / +0.07 % | PASS |
| L | rsw r0/w 1000 vs Cartesian rsw | L1 2-5e-4 of the perturbation, both arms | PASS |
| T1 | relax1/relax2 early T ladder | error to round-off (2e-12) by dt = tau/8; M = K | PASS |
| T4 | dif r0 = 3 (exact) | L1 7.2e-3 / 1.8e-3 / 4.6e-4, X order 2.01 / 2.00; T ladder flat (space floor) | PASS |
| T5 | spherical Marshak wave (new, c 100, rho kappa_a 10) | X order 1.0 (dens, v), 0.5-0.75 (E, eint), as expected at a front; r0/w 1000 vs Cartesian 1.3-1.9 %; M-K 2.5e-5 in v | PASS (first order at the front; not consrc-sensitive, E/rho c^2 ~ 1e-6) |
| T6a | thin black-body sphere, 4x4 / 16x16 | F r^2 second order (1.9e-2 -> 1.3e-3); E L1 4x4 1.2/2.4/3.2 %, 16x16 0.18/0.62/0.95 % (grows with N_r, the known lateral-interpolation error) | flux PASS; E open (known) |
| T6b | IDORT sphere, 4x4 / 16x16 | interior and edge converge (r<0.9 -0.12 %; edge -1.4 / -2.4 % at N 256); exterior r>3 +7.6 % (4x4), **+4.2 % (16x16)**, not converging; L 1.03-1.04 | interior PASS; thin exterior FAIL (lateral error, 16x16 does not cure it) |
| T7 | shadow, 48x48, 120 directions | umbra E/E_free 0.049; lit L1 0.39-0.52; penumbra 0.50-0.68 (ray effects, as in sptests; 432 dirs gave 0.08-0.14) | PASS umbra, ray-limited elsewhere |
| T9 | rsw (rho kappa_s 1000) | interior X order 2.0 (both); T: be 1.0, H M 1.9-2.0 (dens) / 1.2-1.7 (v), **H K 1.5-1.7 / 1.2-1.4** | PASS X; T-order drop with K |
| T9v | rswv: rsw + radial flow v_r = 0.01 sin (moving shell, v tau/c ~ 20) | X order 1.8-2.0 both; T be 1.0, H ~1.0-1.1 (both arms) | X PASS; H T order 1 in both arms (open) |
| T9t | rswt / rswtv (rho kappa_s 10, tau 20) | **no convergence** in X or T, either arm: see section 5 (vet_gd ignores the top bath) | test-design artefact, excluded |
| new | rlat lateral radiation-acoustic wave (thick), +rotation | X order 1.85-1.9 both; T be 1.0; H 1.5/1.3 both, with an H dt-error 9x that of be | X PASS; H lateral T order open |
| new | dynamic diffusion with rotation (drot) | smoke only: hydro itself does not conserve E/L_z for a super-sonically rotating cold heavy gas (centrifugal collapse onto the walls); replaced by rlat_om5 | dropped |
| T8 | injected beams | skipped (as instructed) | - |

## 2. Held-gas stages (P1, P2, L, T4, T6a, T6b, T7): bitwise

ana/out/bit.txt has 32 lines of `IDENTICAL`. Picard counts are identical too. The numbers in the table come from
ana/out/exact.txt and ana/out/lr.txt. Notes on the final code and hesdirk2, compared with the sptests be numbers:

- P1 patm. At R_ph: be -0.73 / -0.26 / -0.05 % (sptests); here hesdirk2 + vet_seg emit gives -0.47 / -0.19 /
  -0.014 %. The thin-envelope plateau is +1.14 % at N 256, unchanged (lateral 16x16).
- P2 pinc. L/L_ex is 1.50-1.56 under hesdirk2 (be 2.25-2.35). There is one NONCONVERGED Picard (200 passes, resid
  1.8e-6, cycle 82) that be did not have.
- T6b at 16x16. The thin exterior is +3.3 % (1.1-3) and +4.2 % (r > 3) at N 256, against +8.8 / +7.6 % at 4x4. It
  does not converge with N_r. 16x16 roughly halves the error but does not remove it. This is the per-shell lateral
  interpolation error of sptests P1b.

## 3. Coupled-gas stages: M vs K

L1 |M-K| / mean at the finest level (ana/out/diff.txt, ana/out/lrc.txt):

| case | N | v (velx) | E | rho | comment |
|---|---|---|---|---|---|
| rsw (thick, E/rho c^2 0.03) | 64 / 128 / 256 | 3.1 / 3.1 / 3.1 % | 4.6e-4 | 3.5e-4 | N-independent: different continuum limit (radiation inertia) |
| rsw r0/w 1000 | 128 | 3.0 % | 4.4e-4 | 3.3e-4 | same in the plane limit |
| rsw Cartesian vet_sc | 128 | 3.0 % | 4.5e-4 | 3.3e-4 | same on vet_sc: not a vet_gd effect |
| rswv (moving shell) | 256 | 0.79 % | 5.1e-4 | 3.9e-4 | |
| rlat / rlat_om5 (lateral, thick) | 128 | 9.3 % | 1.7e-3 | 1.4e-3 | longer t = 4, lateral phase |
| rlatt_om5 (lateral, thin) | 128 | 3.1 % | 9e-4 | 1.5e-3 | |
| marsh (c 100) | 64-256 | 2.5e-5 | 5e-7 | 5e-7 | E/rho c^2 ~ 1e-6 |
| patmT | 256 | (noise) | 2e-12 | 0 | static |
| relax1/2 | all | 0 to print | | | |

**Which arm is right? rsw wave speed** (Cartesian twin and r0/w 1000; the right-moving v peak at t = 2, with
parabolic sub-cell peak location):

| | N 64 | N 128 |
|---|---|---|
| M peak | 2.25155 | 2.25499 |
| K peak | 2.24675 | 2.25007 |
| M - K | 0.0048 | 0.0049 |
| predicted shift from the radiation inertia: 2 t (c_s0 - c_s0 / sqrt(1 + (E+P_r)/rho c^2)) | 0.0047 | 0.0047 |

- The LTE mixture gives Gamma1 = 1.349 and c_s0 = 0.1218.
- K's v amplitude is 2 % lower (3.94e-3 vs 4.02e-3). The radiation now carries its share of the momentum.
- `matched` therefore leaves out the radiation inertia (E + P_r) v / c^2, which is an O(E/rho c^2) error in the
  dynamics. `conservative` contains it.
- For AG Car, E/(rho c^2) ~ (v/c)^2-scale ~ 1e-6, so the two coincide there, consistent with the consrc AG Car result
  (dumps within 1e-7).

## 4. Conservation diagnostics

**L_z, rlat_om5 / rlatt_om5** (ana/out/cons.txt):

- Setup: rsw medium, a lateral mode 1 + 0.1 cos(6 phi), and solid-body rotation Omega = 0.005. phi is periodic, the
  theta and radial walls reflect, so L_z = Sum r sin(theta) (rho v_phi + F_phi/c^2) dV is conserved.
- Totals are from the bin dumps over the full 3-D wedge.

| arm | N_phi | dL_tot / max\|L_rad\| at t = 0.5 | at t = 4 |
|---|---|---|---|
| HFfM thick | 32 / 64 / 128 | -4.17e-2 (64) | -4.10 / -4.15 / -4.16e-2 |
| HFfK thick | 32 / 64 / 128 | +9.6e-7 (64) | +6.5e-6 / +5.4e-6 / +4.7e-6 |
| FfM / FfK thick (be) | 64 | -4.18e-2 / +8.4e-7 | -4.16e-2 / +5.6e-6 |
| HFfM thin (tau 10) | 32 / 64 / 128 | -4.22e-2 | -5.22 / -5.23 / -5.23e-2 |
| HFfK thin | 32 / 64 / 128 | -1.3e-5 | -1.1e-4 / -1.3e-4 / -1.5e-4 (the open Marshak top) |

- matched: the gas gives 9.66e-4 of L_z while the radiation receives 9.27e-4 (64, t = 0.5). The excess, 4 % of the
  radiation momentum, disappears at start-up. The drift then stays flat, so it is the uncharged radiation inertia
  dF/dt, not a per-step leak.
- conservative closes L_z to 1e-6 to 5e-6 of L_rad, i.e. 2e-7 of the total L_z (Picard-level, tol 1e-12).
- Total energy is not a closed check in these tests: every case has a Marshak top or bath. M and K agree to <= 2e-5
  relative everywhere (rsw -3.2e-6 for both; relax1 -1.796e-4 for both, the bath; marsh identical).

## 5. Findings that are not consrc (present in both arms)

1. **vet_gd ignores `implicit_ebath_x1max`: its formal solution always takes a vacuum top**
   (rad_m1_vetgd.cpp:3754 only has an inner bath). Where the top is not thick, a uniform equilibrium with a top bath
   is not stationary:
   - rswt_q (no pulse, rho kappa_s 10, tau 20, bath = E): E falls by 7 / 10 / 17 % near r = 2.8 at N 64/128/256. The
     gas is pushed into the top wall (rho 1.14 -> 1.53 there). The error grows with N.
   - This holds with be and hesdirk2, and with rotate 0 and 10.
   - rsw_q (tau 2000) keeps the equilibrium to 2e-4.
   - The Cartesian vet_sc twin keeps it exactly.
   - So rswt, rswtv and the thin rlatt runs (Marshak bath 0.03 at the top) are not valid convergence tests. Their M/K
     comparison is still like-for-like.
   - Production uses a vacuum top, so this does not affect AG Car. Supporting a top bath in vet_gd would only matter
     for test problems.
2. **hesdirk2 on the lateral radiation-acoustic wave (rlat_T, thick, vet_col with D = I/3):**
   - The order is 1.5/1.3, not 2.
   - The dt error is 9-10x that of be at c dt/dx_phi ~ 1 (v: 3.3e-2 vs 3.5e-3 at dt0).
   - It is the same for M and K, and the same with rotation.
   - The radial rsw is 1.9-2.0 (dens) under hesdirk2, so the problem is specific to the lateral (offdiag / operator)
     path.
   - **Open; this is worth a look before relying on hesdirk2 time accuracy for lateral dynamics.**
3. rswv (moving shell) under hesdirk2 has T order 1.0-1.1 in both arms (be 1.0). The second order seen on the static
   rsw does not survive a background radial flow v_r = 0.01 c. Open.
4. pinc under hesdirk2 has one NONCONVERGED step (Picard 200); be had none.
5. T5 Marshak (new): r0/w = 1000 vs Cartesian differs by 1.3-1.9 % at N 128, against the 1.1e-3 floor of the pulse.
   This is the vet_gd vs vet_sc difference at a thin Marshak front (tau 10 total). Not chased.

## 6. Recommendation

- **Keep `implicit_gas_src = conservative` opt-in for now. Do not switch the default yet.** The arguments for it are
  strong:
  - it is bitwise on every held-gas test and round-off on static coupled ones;
  - it conserves L_z to 5e-6 of L_rad where matched loses 4-5 % of the radiation momentum;
  - it gets the radiation-inertia sound speed right to 4 % of the shift (0.0049 vs 0.0047), where matched is 2 % fast
    at E/rho c^2 = 0.03;
  - X orders are unchanged.
- But it **lowers the hesdirk2 time order on rsw** from ~2 to ~1.5 in dens and to 1.2-1.4 in v and F, with 2-3x the
  dt error. This is the rw_t10a F1 drop, now seen on the production path and in every variable.
- The likely place is the stage F^old / `csinc` handling of the inertia term in Time2FormStage (rad_m1_time2.cpp):
  the gas pays d(aE)/dt per stage. A ladder with `implicit_rad_inertia` off, or a check of the stage-increment
  algebra, is the next step.
- For AG Car the choice is immaterial: E/(rho c^2) ~ 1e-6, and the consrc AG Car dumps agree to 1e-7.
- Once the hesdirk2 order is restored, make conservative the default for gas_feedback = true.
- Independent of consrc:
  - look at the hesdirk2 lateral time order (finding 2);
  - optionally add a top bath to the vet_gd formal solution so that test problems with a bath are consistent
    (finding 1).

## 7. Files

- ana/spt.py, new cases:
  - `rswt` (rho kappa_s 10), `rswv` / `rswtv` (gas_vr 0.01), `_r<k>` / `_cart` for rsw, `_q` (no pulse);
  - `marsh_r<k>` / `marsh_cart` (T5);
  - `rlat` / `rlatt` (lateral wave), with `_om<k>` (rotation) and `_rt` (reflect top; refused by vet_gd);
  - `drot` (dropped);
  - token `h<k>` (k outputs).
- ana/spl.py: arm suffixes M / K.
- ana/cmpmk.py: bit / diff / xord / tord / cons / lrc. ana/evalall.sh runs all of them.
- Outputs: ana/out/{bit,exact,lr,xord,tord,cons,diff,lrc}.txt.
- Job lists: jobs_smoke.txt, jobs_A.txt, jobs_B.txt, jobs_C1.txt, jobs_C2.txt, jobs_C3.txt. Logs are in logs/.
