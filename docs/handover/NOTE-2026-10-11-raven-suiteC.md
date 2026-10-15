# NOTE suiteC (suiteC-1011): the Cartesian IDORT/Jiang suite with the conservative gas source

Work dir: /raven/ptmp/jinma/suiteC_1011. MPCDF Raven, CPU (`general`, 72 cores, OpenMP, one MeshBlock).

**Question.** Is `rad_m1/implicit_gas_src = conservative` (consrc-1011, merged opt-in in
prodfinal-1011) safe on every test of /raven/ptmp/jinma/suite_1011/NOTE-suite.md? Why does the
rw_t10a F1 time order drop from 2.0 to 1.7 with it? Can it become the default?

**Short answer.**
- **Do not make `conservative` the default.** It is a *model* change, not a re-bookkeeping of the
  same model: it gives the gas the extra inertia (4/3) E/c^2 everywhere, whatever the optical depth.
  - On every test with gas_feedback = false or static gas it is **identical** to matched (cb, cbg,
    cbh, shd, shdg, nLTE, sdif, ddn, hcyl, bubbles, pipe, G5: 0 difference or bitwise).
  - On the gas-coupled tests with small E/(rho c^2) it **agrees with matched to <= 1e-7 ... 1e-9**
    (shocks, Marshak, dynamic diffusion, cloud), and it fixes eq3t (P_tot 3.000000, E 1.13260 = the
    analytic 1.1326; matched 2.99938 / 1.13407).
  - It **FAILS the linear waves at P = 10** (c = 10, E/(rho c^2) = 0.1): L1(rho) stalls at 0.11-0.21
    for N = 16 ... 256 in all three tau, omega off by 7 % at tau = 0.01 and 1, BE and hesdirk2.
    Matched passes (2nd order, omega to 1e-4).
  - With hesdirk2 it is **first order in time**, in every variable, with an error coefficient
    proportional to E/(rho c^2) (P = 10 T-ladder: order 0.96-1.07; matched 1.9-2.0 in rho).
- **The rw_t10a F1 drop is the small-epsilon end of the same defect.** At c = 1000, P_rad/P_gas = 100,
  E/(rho c^2) = 3e-4: the conservative solution differs from matched by a dx- and dt-independent
  frequency shift that is exactly the dispersion shift of a gas inertia rho -> rho + 4E/(3c^2)
  (measured -1.0595e-3 + 8.090e-4 i, predicted -1.059e-3 + 8.09e-4 i), plus a first-order-in-dt term
  of size ~E/(rho c^2): F1 order 1.65 at l = 8-16 and **0.96 at l = 16-32**. Not start-up, not the
  Picard tolerance (tol 1e-10/1e-12/1e-14 identical), not implicit_rad_inertia (bitwise with it off).
- **Where it comes from.** The conservative deposit charges the gas -(F'_c - F^old_c)/(c chat), and
  the cell flux F_c = F0 + a(w) E carries a(w) = (4/3) w of the *kicked gas velocity*. The face
  equation of F0 never pays -d(aE)/dt (the f0 model; matched omits the momentum, which is the moveq
  leak). So the velocity part of d(aE)/dt goes 100 % to the gas, at every optical depth. That is the
  physical radiation-enthalpy inertia only in the diffusion limit: at P = 10, tau = 100 RECc gives
  omega = 14.31 against the exact O(v/c) transport root 14.385 (matched 15.31, the f0-model root
  15.315), but at tau = 0.01 and 1 the transport root sits at the matched value (7.991, 6.288) and
  RECc is 7 % off (7.44, 5.84). And because d(aE) contains the Heun (explicit) hydro velocity change
  of the stage, the reaction to the explicit force is integrated with the ESDIRK stage weights:
  first order in time with coefficient E/(rho c^2).
- **Recommendation.** Keep matched as the default (as the user decided). Before conservative can be
  default, the velocity part of dA = a(w') E' - a(w^n) E^n must get the same Theta split that
  implicit_rad_inertia gives the E part (thick: gas, thin: F0 keeps it), and the explicit hydro
  Delta v must be taken out of F^old (charge only the implicit-stage change). Neither is in the code.
  For AG Car the effect is negligible (E/(rho c^2) ~ (c_s/c)^2 ~ 1e-7), as consrc's AG Car test found.

## 0. Binary, arms, runs

- **Binary** bin/athena_cpu_530c38e2 (md5 656e41770164cd24a592a97d97290968), a copy of
  /raven/ptmp/jinma/prodfinal_1011/bin (prodfinal-1011 530c38e2 = prodcand2 767d3313 + consrc merge
  1cb04a19 + vet_seg default emit). Used for every run. Its `m1_test = cyl` already allows 24
  objects, so the bubbles need no test build.
- **Arms** (arms.sh; NOTE-suite definitions, cen/hr not rerun):
  - REC = tJe keys (`implicit_thin_corr=hcs`, jse, lin, fse, xthin 0) + hesdirk2 + time2_vet_sc
    predict + time2_vstage + implicit_vimp (on where gas_feedback and dynamic gas).
  - RECbe = the same with time_scheme = be.
  - RECc / RECcbe = REC / RECbe + `rad_m1/implicit_gas_src=conservative`.
  - REC and RECbe were rerun on the new binary so that every matched/conservative pair is on the same
    code. Differences against NOTE-suite numbers come from vet_seg = emit and predfix (sect. 3).
- **Jobs** (all `sbatch -p general -c 72`, <= 3 at a time): pool jobs 31078139 (jobs_p1.txt, 104
  runs), 31078140 (jobs_p2.txt, 37 runs), misc 31078141 (ladders, linear waves, G5), diagnosis
  31079095 / 31079254 / 31079384. Nothing went to Caltech (no job waited > 20 min).

## 1. Matched vs conservative, per test

"identical" = the last dumps differ by 0 (the gas momentum source is not used: gas_feedback = false
or kinematic gas). Numbers are REC / RECc unless noted; pass criteria as in NOTE-suite.

| test | criterion | REC (matched) | RECc (conservative) | RECc vs REC | status RECc |
|---|---|---|---|---|---|
| crossing beams ring 4, c dt/dx 256 | row 1.0, peak >= 0.3, floor <= 1e-3 | 1.006 / 0.417 / 4.5e-18 (now completes: predfix) | same | identical | PASS |
| crossing beams c dt/dx 0.5, ring 12, 256x1024, Gauss 8x32 | " | 1.006 / 0.417 | same | identical | PASS |
| static shadow ring 12x2 | L1 <= 2 %, umbra <= 1 % | L1 0.011, umbra 0.22, edge 16 | same | identical | PASS |
| static shadow Gauss 8x32 | robustness | L1 0.030 | same | identical | = REC |
| thermal eq. 1 | E, T -> 96.795, 3.13663 | 96.79505, 3.13663 | same | E bitwise, v noise 1e-15 | PASS |
| thermal eq. 2 | 145.788, 3.4747, v = 0 | 145.78781, 3.47480, max v 1.7e-11 (predfix: was v 1.7) | same, 1.7e-11 | E bitwise | PASS |
| moving eq. eq3t, v 3, c 10 | E 1.1326, v 2.9554, P 3 | 1.13407 / 2.95471 / 2.99938 | **1.13260 / 2.95538 / 3.000000** | conservative better | PASS (better) |
| moving eq., BE | " | 1.13411 / 2.95469 / 2.99937 | 1.13262 / 2.95537 / 3.000000 | " | PASS |
| non-LTE eps 1e-1 ... 1e-8 | L1 <= 2 %, E_top 1 % | L1 1.11 % ... 0.01 %, E_top 0.9994-1.0000 | same | identical (kinematic) | PASS |
| static diffusion | L1 <= 1 % | 0.28 % at t_end, peak 1.0001 | same | identical | PASS |
| dynamic diffusion ddf | centroid at vt, peak 5 % | 0.9700 / -1e-4 / L1 0.0237 / D_eff 0.984 | identical to 4 digits | <= 1e-4 | PASS |
| ddiff / ddn | " | 0.9700 / 0.9983 | same | ddn bitwise | PASS |
| BE ddf / ddn | " | 0.8748 / 0.9296 (BE time error, as before) | same | identical | as REC-be |
| linear waves P = 0.01, tau 0.01/1/100 | L1 ~ N^-2, omega | orders 1.99 / 2.05 / 2.12 at N 256 | 2.42 / 2.53 / 2.71 (closer to the reference) | omega shift 5e-4 | PASS |
| **linear waves P = 10, tau 0.01/1/100** | " | orders 2.17 / 2.06 / 2.02, L1(256) 2.9e-4 / 1.5e-4 / 2.3e-4 | **L1 0.178 / 0.115 / 0.213, no convergence**; omega 7.44 / 5.84 / 14.31 vs ref 7.988 / 6.234 / 15.315 | 7 % in omega | **FAIL** |
| linear waves P = 10, BE | " | 1st order (BE) | no convergence | | FAIL |
| Marshak c dt/dx 1 / 10 | L1(E) <= 0.02 | (NOTE-suite 0.0139 / 0.0141) | E to 1e-13, v to 3e-9 rel. | ~ identical | PASS |
| irradiated cloud | qualitative | x_c 0.0088 / 0.0157 / 0.0276 | same | 1e-7 | reproduced |
| G5 (static, 4 kappa) | ratio within 2 % | 0.99957 / 0.99987 / 1 / 1 | same (bitwise numbers) | identical | PASS |
| radiation shocks M 1.2 / 2 / 3, 512 | L1 <= 2 % vs IDORT transport | <= 0.0004 / 0.0014 / 0.0012 | same | <= 1e-9 | PASS |
| shocks M 2, 3, 1024 | " | <= 0.0007 / 0.0017 | same | identical digits | PASS |
| homogeneous cylinder 256 / 512 | E, rF out 10 % | E out 4.4 % / 3.9 %, F 5.7 % / 6.0 %, f_in 0.3344 / 0.3331 | same | identical | PASS (256 now passes: vet_seg emit) |
| crooked pipe c dt/dx 2200 | dt self-consistency | probes 0.469 / 0.428 / 0.387 / 0.345 / 0.278 keV at t 10 | same | identical (static gas) | reproduced |
| crooked pipe c dt/dx 220 | " | probes 0.469 / 0.428 / 0.387 / 0.345 / 0.277 keV at t 10 (NOTE-suite 0.468 / 0.427 / 0.385 / 0.343 / 0.276 at 9.4) | same | identical (REC and BE pairs) | reproduced |
| multi-source bubbles, ring 100, rho kappa 1e5 / 1e3, Gauss 1e3 | L1 <= 20 % | L1 0.67 / 0.76 / 0.76 at 1 t_c, E_tot 0.52 / 0.68 / 0.57 | same | identical | FAIL (as before; better than 0.21 E_tot: vet_seg emit) |

Old cen/hr numbers: see NOTE-suite sect. 1 (not rerun).

**Acceptance criteria** (2nd order in space and time, 1st order in space allowed only in the optically
thin part): REC passes them on the waves, diffusion and shocks. RECc fails the time criterion
(hesdirk2 1st order, sect. 2.3) and the consistency criterion (P = 10 waves do not converge to the
model's own dispersion root).

## 2. The F1 time-order drop (rw_t10a) and what it is

### 2.1 Ladders on the new binary (ladder/, od.py, rw_t10a and rw_t1000a, X/T/C, mt vs cs)

- All X and C orders, and all BE orders, are identical between matched and conservative (lad_mt.txt,
  lad_cs.txt). rw_t1000a: identical orders in every mode.
- rw_t10a hesdirk2 (tJePvs = REC) T mode, levels l = 1, 2, 4, 8, 16, 32 (dt = T/(1100 l), N = 128):

  | var | matched L1 (successive) | orders | conservative L1 | orders |
  |---|---|---|---|---|
  | E | 9.5e-4 1.5e-4 3.7e-5 9.9e-6 2.4e-6 | 2.67 2.00 1.91 2.05 | 9.5e-4 1.5e-4 3.7e-5 9.9e-6 2.4e-6 | 2.67 2.00 1.92 2.07 |
  | F1 | 4.0e-4 1.0e-4 2.8e-5 7.0e-6 1.6e-6 | 1.99 1.84 2.01 2.15 | 3.7e-4 8.8e-5 2.4e-5 7.5e-6 **3.8e-6** | 2.04 1.91 1.65 **0.96** |
  | rho | 5.3e-5 ... 2.0e-7 | 1.98 1.96 2.08 2.05 | 5.2e-5 ... 1.6e-7 | 2.00 1.99 2.16 2.24 |
  | v | 3.6e-5 ... 2.1e-7 | 1.91 1.87 1.89 1.76 | 3.6e-5 ... 1.9e-7 | 1.92 1.88 1.92 1.86 |

  The extra level shows that the conservative F1 has a genuine first-order term (it flattens at
  ~4e-6), not a fluctuation. Matched at l = 32 continues at 2.15.
- **Not the Picard tolerance:** implicit_tol 1e-10 (tJePvsL) gives the same digits as 1e-12;
  1e-14 / lin 1e-15 (tJePvsQ) gives F1 1.79 (mt) / 1.21 (cs) at 16-32.
- **Not start-up:** the cs-mt difference grows steadily over the period (ana: csdiff.py; F1 3e-4 at
  1/16 period, 1.3e-3 at the end), and is the same at every dt.
- The fine-level F1 orders of *both* paths also oscillate over the wave phase (1.4-2.3 per dump at
  l = 8-16, csdiff.py); the period-end value sampled 2.01 (mt) vs 1.65 (cs) in consrc. The l = 32 level
  separates the two.

### 2.2 A dt- and dx-independent model offset

- cs - mt in the last dump (rw_t10a, all modes and levels) is a pure k = 1 mode of relative size
  1.4e-3, the same at N = 32 ... 512 (X and C) and l = 1 ... 32 (T), in BE as well as hesdirk2, and in
  the consrc cen ladders too.
- As an effective frequency (k = 1 coefficient, IC -> one period): matched X512 differs from the code
  model's root vetqsc_d (5.961672 - 2.017107 i) by +1.1e-4 + 7e-5 i (discretisation); conservative
  by -9.5e-4 + 8.8e-4 i. The cs - mt shift is -1.0595e-3 + 8.090e-4 i at every T level.
- The dispersion relation of the code model (radwave_disp.py, 'f0' form) with the gas inertia
  rho -> rho + 4 E0/(3 c^2) gives a shift of **-1.059e-3 + 8.09e-4 i**: an exact match.
- The same at P = 10, c = 10 (E/(rho c^2) = 0.1): f0 model + inertia predicts 7.489 - 0.528 i /
  5.862 - 0.769 i / 14.384 - 0.447 i for tau = 0.01 / 1 / 100; RECc measured 7.441 - 0.522 i /
  5.839 - 0.751 i / 14.312 - 0.457 i (N = 256).
- implicit_rad_inertia = false does not change it (RECc P = 10 omega 5.84 / 14.31 again; runs/rwD),
  and in the T ladders the rad_inertia-off rows are bitwise equal: in a linear wave the E-change part
  of dA (which rad_inertia handles) is second order; the velocity part is the one at work.

### 2.3 The first-order time error is proportional to E/(rho c^2)

T ladders of the P = 10 wave (inp/rw_p10t{1,100}, N = 128, dt = T/(256 l), l = 1 ... 16):

| case | arm | rho orders | v orders | E orders | F1 orders | rho L1 at 8-16 |
|---|---|---|---|---|---|---|
| tau 1 | matched hesdirk2 | 1.97 1.92 1.55 | 1.74 1.63 1.41 | 2.00 0.63 0.66 | 1.35 1.11 1.08 | 1.6e-5 |
| tau 1 | conservative hesdirk2 | 1.24 1.00 0.96 | 1.19 1.11 1.05 | 1.17 1.12 1.07 | 1.03 0.99 1.01 | 8.8e-5 |
| tau 100 | matched hesdirk2 | 2.02 2.00 1.89 | 1.08 1.17 1.02 | 2.32 2.27 1.04 | 1.10 1.18 1.03 | 3.3e-6 |
| tau 100 | conservative hesdirk2 | 0.92 0.95 0.98 | 1.09 1.05 1.03 | 1.15 1.08 1.04 | 1.09 1.05 1.02 | 1.4e-4 |
| both | either, BE | 1.00 | 1.00 | 1.00 | - | |

- Conservative + hesdirk2 is first order in every variable, with errors 5-40x those of matched.
- The first-order coefficient scales with E/(rho c^2): rw_t10a (3e-4) shows it only in F1 at
  l >= 16, P = 10 (0.1) everywhere. Worker C's spherical rsw (sptestsC) is in between.
- Side finding: **matched hesdirk2 is not clean either at P = 10** in this fixed-N T ladder (v and F1
  ~1.0-1.2 at tau = 100, E 0.6 at tau = 1, small coefficients). The C ladders and the suite waves
  (dt ~ dx) still give 2nd order, so it is a small dt-only term (predict / vstage path?), not
  investigated further here. Without vimp (tJeP) both paths are first order (as known from
  torder-1010: hesdirk2 needs vimp + vstage), so that variant does not separate them.
- **Mechanism (code reading, rad_m1_impl_wb_kernel.hpp ~l. 289-420, rad_m1_time2.cpp ~l. 483):**
  F^old_c is the stored cell flux at the stage start + csinc, i.e. a(w_old) E_old; F'_c uses the
  kicked velocity w of the gas, whose momentum already holds the Heun stage's explicit hydro update.
  The gas therefore pays (4/3) E Delta v_hydro / c^2 inside the implicit solve: the inertia reaction
  to the explicit force enters through the ESDIRK slopes (csk1/csk2) while the force itself enters
  with the Heun weights. The two quadratures differ at O(dt), so the error is (E/(rho c^2)) O(dt).

## 3. Changes against NOTE-suite on the new binary (matched)

- predfix (time2_pred_cap 0.5 etc.): **thermal equilibrium 2 now passes under REC** (v 1.7e-11,
  E 145.78781) and **crossing beams at c dt/dx 256 complete** (was FATAL).
- vet_seg = emit (new default): **homogeneous cylinder 256^2 E outside 12 % -> 4.4 %** (now PASS),
  bubbles at rho kappa 1e5 E_tot 0.21 -> 0.52 at 1 t_c (still FAIL), the opaque-emitter defect of
  NOTE-suite sect. 4 is what emit addresses.
- eq3t matched: E 1.13407 (NOTE-suite 1.13329), P_tot leak 6.2e-4.

## 4. Recommendation

1. Keep `implicit_gas_src = matched` as the default. Conservative fails the P = 10 linear waves
   (no convergence, 7 % omega) and is first order in time with hesdirk2 wherever E/(rho c^2) is not
   tiny. It must stay opt-in, labelled "diffusion-limit inertia; not for radiation-dominated thin gas".
2. It is harmless for AG Car (E/(rho c^2) ~ 1e-7) and better on eq3t, so production may keep either.
3. A fix to make it default-capable (not implemented):
   - split the velocity part of dA = a(w') E' - a(w^n) E^n with the face Theta exactly as
     ImplicitRadInertia does for the E part (thick -> gas, thin -> F0), so the inertia appears only
     where the radiation is coupled (this is what the transport root does: inertia at tau = 100, none
     at tau <= 1);
   - form F^old_c with the stage's post-hydro velocity, a(w_stage) E_old, so the explicit hydro
     Delta v is not charged at the implicit stage (removes the O(E/(rho c^2) dt) term).
   Test: the P = 10 waves (omega -> transport root at tau 100, -> f0 root at tau 0.01), the P = 10
   T ladder (order 2), eq3t (P_tot), rw_t10a T to l = 32.
4. The matched P = 10 T-ladder irregularity (v ~1.0 at tau = 100) deserves its own look.

## 5. Files

- Scripts: arms.sh, run1.sh, pool.sbatch, misc.sbatch, jobs_p1.txt, jobs_p2.txt; g5/run_g5.sh (+
  tJe_h2c/tJe_bec arms); ladder/{od.py (+ tJePvsQ/L, tJePvsNI, tJeP cases rw_p10t1/100),
  lad.sbatch, diag.sbatch, diag2.sbatch, jobs_*}; ana/rw/cases_C.txt, cases_D.txt.
- Analysis: ana/pair.py (last-dump comparison), ana/eqbin.py, ladder/ana/csdiff.py (cs-mt per level and
  dump), ladder/ana/csdx.py; tables ladder/lad_mt.txt, lad_cs.txt, ana/rw_table_suiteC.txt;
  fig/suiteC_cloud.png.
- Runs: runs/<test>/<arm>_<tag>, runs/rw, runs/rwD, ladder/run_{mt,cs}/runs, g5/run.
- old_scaffold/: the unused job lists of the earlier attempt.
