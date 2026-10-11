# NOTE consfix (consfix-1011): the conservative gas source made consistent — `implicit_gas_src = conservative2`

Worker G, 2026-10-11. Work dir /raven/ptmp/jinma/consfix_1011. Branch consfix-1011 from prodfinal-1011 8279c0a6,
worktree /raven/u/jinma/ATHENAK/wt-consfix, code commit **27c35877** (not pushed).
Binaries: bin/athena_cpu_27c35877 (built-in pgens), bin/athena_he_cpu_27c35877, bin/athena_he_a100_27c35877
(md5 a1b0c047 / 81b61d4b / d6d1f80a).

## 1. Result

- New opt-in value `<rad_m1>/implicit_gas_src = conservative2`. `matched` stays the code default.
  `matched`, `conservative` and `conservative_full` are **bitwise unchanged**: 704/704 eq3f tab files against
  the parent binary 18ff91c8, and 8/8 P = 10 wave runs against suiteC's 530c38e2 runs.
- conservative2 fixes both defects that suiteC found in `conservative`:
  - The P = 10 linear waves converge at second order, to the lab-F moment-model root (`vetqs_d`). At tau = 100
    that root is the transport root, 14.385. The old conservative did not converge, and its omega was 7 % off.
  - hesdirk2 is second order in time again. The F1 drop on rw_t10a is gone: 2.26 at l = 16-32, against
    0.96 for conservative and 2.15 for matched. On worker C's rsw (vet_gd, the production path) every
    variable has order 2.0 with hesdirk2, against M 1.2-1.9 and K 1.2-1.5.
- It conserves total momentum like conservative:
  - eq3f P_tot to round-off (8.9e-16), E 1.13260 (analytic 1.1326).
  - rlat_om5 L_z drift +5.3e-6 of L_rad, against -4.15e-2 for matched.
- Reduced CPU AG Car A (production FRESH deck, 120 cycles):
  - stable; same dt, budget, L_top and E_rad as matched;
  - Picard 2.98 vs 2.94;
  - cost +2.8 %.
- **Recommendation:** usable for production now (named in the deck). Make it the code default for gas_feedback
  runs after H's hesdirk2 fix is merged and the T-ladder is re-checked. A ~1e-6 first-order component is left in
  the P = 10 T-ladder. matched has a 10-50x larger one, so it is part of the common hesdirk2 issue, not this source.

## 2. Diagnosis (what was really wrong with `conservative`)

1. **First order in time = vimp inconsistency, not F^old.**
   - The vimp E rows predict the gas velocity change from the drag rule (dt k F0/c).
   - The conservative write-back gives a different one: it also includes -(dA)/(c chat).
   - So the enthalpy flux that moved E used a different velocity from the one stored in F and in the gas, an
     O(E/rho c^2 dt) mismatch every step.
   - Proof with diagnostic bits on rw_p10t100, hesdirk2 T, L1(rho) at l = 8-16:
     - conservative: 1.4e-4, order 0.98 (bit 7 is bitwise equal to it);
     - conservative + the vimp prediction alone (bit 5): 3.5e-6, order 1.8-2.1.
2. **The model error (100 % of the velocity part of d(aE)/dt to the gas) needs the lab-F face equation.**
   - The face equation must be F0' = theta (F0^n - dA_f - ...).
   - The split must sit **inside** the solve.
   - A first attempt applied the split after the solve, as ImplicitRadInertia does for the E part. The E rows
     then transport F0 without -theta dA while the stored faces have it. That gave order ~0.7, with errors
     3x those of conservative.
3. B's proposal (b), "F^old at the post-hydro velocity", breaks momentum conservation as stated.
   - The conservative equivalent is to keep the lab F and give F0 the hydro's explicit tendency -dA/dt.
   - It changes only the ESDIRK slope bookkeeping: one stage solve is identical.
   - It is implemented as (c) below.
   - It matters little: with and without it, the P = 10 combined waves and rsw are equal to 3 digits. The
     fine-dt T-ladders split as follows:
     - tau 100: without it, order 2.0 to l = 32;
     - tau 1: with it, 1.86 vs 1.43 at l = 32-64.
   - It is kept as the IMEX-consistent choice. Diagnostic bit 4 turns it off.

## 3. What conservative2 does (rad_m1_implicit.cpp, rad_m1.hpp)

- **(a) ImplicitCs2Pass**, once per Picard pass after the vimp build:
  - dA_c = A'_c - A^old_c:
    - A' = the enthalpy flux the write-back stores: a(v + dv_pred) E' (hesdirk2 + vimp + vstage), else a E';
    - A^old = F^old - F0^old face mean, with F^old exactly the conservative write-back's old cell flux.
  - Its face mean goes into F0^n of every face that has a face equation (x1, x2, x3; tracked in cs2f*).
  - So the rows (d), the face kernels (g) and the transverse terms all use the same F0.
  - Z_f adds back the true F0^n, so conservation by construction is unchanged.
  - F0^n is restored before the write-back, so the slopes and the boundary terms see the true old value.
  - ImplicitRadInertia is skipped (its E part is inside dA).
  - The drag then hands the gas (1 - theta) of dA:
    - thick: all of it (radiation inertia);
    - thin: the lab flux is kept.
  - AP-HLL faces: only the (1 - AL) part has a face equation. The HLL part keeps its drag-only rule, as in
    conservative.
- **(b) ImplicitCs2Vimp**, inside ImplicitVimpBuild:
  - dv^k becomes the conservative write-back's velocity change:
    rho dv = rho dv_drag - [(1 - c) dA(dv) - sum_f (1/2) om_f w_f dA_nbr]/(c chat), with om_f = 1 - AL.
  - The own dv is implicit (eta = (1 - c)(1 + D_dd) E'/(rho c chat)). The neighbours' dA and the D_de cross
    terms are lagged by one pass.
  - The P rows are scaled by 1/(1 + eta).
- **(c) ImplicitCs2Shift**, at the solve start, before the F0^n copy (not at the first conservative solve):
  - per cell, s = a(v_start) E_start - (F_stored - F0 face mean);
  - F0 faces -= the face mean of s;
  - the lab F (u0) is untouched.
- **Cost:**
  - 6 iw components (A^old, dA) and the cs2f face arrays, allocated only with conservative2;
  - per pass, 1-3 scalar halo exchanges (dA).
- **Not supported:** the cubed sphere (fatal).
- **Diagnostic key** `implicit_cs2_dbg` (bitmask): 1 no split (= conservative + vimp prediction), 2 no vimp
  prediction, 4 no hydro shift.

## 4. Tests (final binary 27c35877 unless noted)

### 4.1 P = 10 linear waves (c = 10, E/rho c^2 = 0.1)

Combined mode, nt = 2N. Conservative2 runs start from the vetqs_d eigenmode and are measured against it
(`ref=vetqs_d`). Matched and conservative are from suiteC, against vetqsc_d.

| tau | arm | L1(rho) N 16 / 64 / 256 | orders | omega(256) | model root |
|---|---|---|---|---|---|
| 0.01 | matched h2 | 6.3e-2 / 5.5e-3 / 2.9e-4 | 1.54 1.97 2.08 2.17 | 7.98859 - 0.53315i | vetqsc_d 7.98801 - 0.53314i |
| 0.01 | conservative h2 | 0.158 / 0.175 / 0.178 | none | 7.441 - 0.522i | |
| 0.01 | **conservative2 h2** | 6.4e-2 / 5.6e-3 / 2.9e-4 | 1.54 1.96 2.09 2.17 | 7.99244 - 0.51941i | vetqs_d 7.99183 - 0.51937i |
| 1 | matched h2 | 3.8e-2 / 2.8e-3 / 1.5e-4 | 1.69 2.08 2.14 2.06 | 6.23459 - 0.87552i | vetqsc_d 6.23399 - 0.87560i |
| 1 | conservative h2 | 0.088 / 0.112 / 0.115 | none | 5.839 - 0.751i | |
| 1 | **conservative2 h2** | 3.9e-2 / 2.9e-3 / 1.7e-4 | 1.67 2.08 2.14 1.94 | 6.28412 - 0.85814i | vetqs_d 6.28356 - 0.85844i |
| 100 | matched h2 | 6.1e-2 / 3.7e-3 / 2.3e-4 | 2.01 2.03 2.02 2.02 | 15.31398 - 0.43696i | vetqsc_d 15.31514 - 0.43666i |
| 100 | conservative h2 | 0.266 / 0.216 / 0.213 | none | 14.312 - 0.457i | |
| 100 | **conservative2 h2** | 6.2e-2 / 3.8e-3 / 2.3e-4 | 2.00 2.03 2.02 2.01 | 14.38446 - 0.37946i | vetqs_d 14.38553 - 0.37925i (transport 14.385 - 0.441i) |

- With BE, conservative2 is first order (as matched BE).
- At P = 0.01, conservative2 converges at order 1.99-2.11.
- vetqs_d is the moment model with the lab-F flux equation. It is the model that conserves momentum, and its
  real part is the transport root at every tau. Its damping at tau = 100 differs from transport by the moment
  closure, not by the gas source.

### 4.2 T and X ladders (ladder/, od.py of suiteC, N = 128)

Space (X): rw_p10t100, rw_p10t1 and rw_t10a all have order 2.0-2.3 in every variable.

**rw_p10t100, T, hesdirk2**, L1 at l = 1-2 ... 32-64:

| arm | rho | v | F1 |
|---|---|---|---|
| matched | 1.9e-4 ... 1.3e-7, order 2.0-2.3 | 4.3e-4 ... 1.2e-5, **order 1.0** | 4.4e-4 ... 1.1e-5, **order 1.0** |
| conservative (suiteC) | order 0.92-0.98, 1.4e-4 at l = 8-16 | 1.0 | 1.0 |
| **conservative2** | 1.7e-4, 3.7e-5, 9.3e-6, then ~1e-6 (irregular) | 1.5e-4 ... 8.5e-7, order 1.8 2.2 1.2 1.2 1.1 | order 1.8 2.2 1.1 1.1 1.1, 9.2e-7 at the end |
| conservative2 without the shift (bit 4) | order 2.04 2.05 2.11 1.98 1.26 | 2.0 to l = 32 | 2.0 to l = 32, 1.5e-7 at the end |

**rw_p10t1, T, hesdirk2**, L1(rho) at l = 1-2 ... 32-64:

| arm | rho |
|---|---|
| matched | order 1.97 1.92 1.55 1.18 1.03 (3.4e-6 at the end) |
| **conservative2** | order 2.00 1.99 1.95 1.91 1.86 (8.3e-7 at the end); v 1.97 ... 1.74 |

**BE, T:** conservative2 has order 1.00 (rho) at tau 1 and tau 100.

**rw_t10a, T, hesdirk2**, l = 1 ... 32:

| arm | F1 orders | other variables |
|---|---|---|
| matched | 1.99 1.84 2.01 2.15 | |
| conservative | 2.04 1.91 1.65 0.96 | |
| **conservative2** | 1.90 2.07 0.99 2.26 (the l = 8-16 dip is the phase sampling seen in both arms) | E 2.0, rho 2.0, v 1.7-1.9 (as matched) |

**rw_t10a, X and C:** order 2.0.

### 4.3 eq3f moving equilibrium (consrc inputs, tJe, c = 10, t = 1)

| arm | E / v | max \|dP_tot\| | max \|dE_tot\| |
|---|---|---|---|
| matched h2v | 1.13406 / 2.95471 | 6.2e-4 | 9.4e-11 |
| conservative h2v / be | 1.13260 / 2.95538 ; 1.13262 / 2.95537 | 8.9e-16 | 9.4e-11 / 1.4e-9 |
| **conservative2 h2v / be** | 1.13260 / 2.95538 (both) | 8.9e-16 | 3.0e-11 / 1.4e-10 |
| conservative2 h2v, tol 1e-12 | 1.13260 / 2.95538 | 1.3e-15 | 2.0e-14 |

The target is E 1.1326 and v 2.9554.

### 4.4 Spherical wedge (vet_gd, fse + hc, worker C's sp harness in sp/, runs_fin; M and K are C's runs)

**rsw_T, hesdirk2**, orders at l = 2-4 / 4-8:

| arm | rho | v | E | F1 | v error at l = 4-8 |
|---|---|---|---|---|---|
| M | 2.02 / 1.90 | 1.73 / 1.23 | | | 6.7e-6 |
| K | 1.66 / 1.47 | 1.39 / 1.19 | | | 1.8e-5 |
| **K2** | 2.04 / 2.07 | 1.98 / 1.97 | 1.98 / 1.96 | 1.99 / 1.97 | 3.5e-6 |

- BE (FfK2): order 1.00.
- rsw: K2 - K is 0.07 % in v, against M - K 3.1 %. K2 carries the radiation inertia in the thick shell, as K
  does.
- rlat_om5_T, hesdirk2: 1.55 / 1.36 for all three arms (C's open lateral issue, not the gas source).
- patmT (static coupled): K2 - M ≤ 8e-10.

**L_z, rlat_om5**, dL_tot / max|L_rad| at t = 4:

| arm | thick, hesdirk2 | thick, BE | thin (rlatt_om5) |
|---|---|---|---|
| M | -4.15e-2 | -4.16e-2 | -5.2e-2 |
| K | +5.4e-6 | +5.6e-6 | -1.3e-4 |
| **K2** | +5.3e-6 | +5.5e-6 | -1.3e-4 (open Marshak top) |

### 4.5 Reduced CPU AG Car A

- Setup: athena_he_cpu_27c35877; a copy of prodfinal_1011/prod/agcar_A_fresh.athinput; 480x32x32, 4 ranks x 6
  threads per arm; 120 cycles to t = 1.503e5 s; only implicit_gas_src varied; runs in agcar/{mt,cs,c2}.

| arm | dt | Picard mean / max | zone-cycles/cpu-s | budget residual (rel. Etot / rel. int L_top) |
|---|---|---|---|---|
| matched | 1252.3 | 2.94 / 5 | 1.787e5 | -1.8e-9 / -1.2e-5 |
| conservative | 1252.3 | 2.94 / 5 | 1.775e5 | identical |
| **conservative2** | 1252.3 | 2.98 / 5 | 1.738e5 (-2.8 %) | identical |

- Last dumps, c2 vs mt:
  - dens, eint and E within 7e-8 (max-relative);
  - F within 7e-7;
  - v1 1.4e-3 of max|v1| (cs vs mt: 1.6e-3);
  - v2, v3 4e-5.
- No failures or NONCONVERGED steps.

### 4.6 Cartesian suite subset (run1.sh arm RECc2, against suiteC RECc)

Every run completed (rc 0). The table gives the last-dump max|a-b|/max|a| of RECc2 against RECc (ana/pair.py).
Velocity components that are noise around zero are left out.

| test | RECc2 vs RECc |
|---|---|
| eq3t (+ RECc2be) | E 3.8e-6, F1 2.5e-5, eint 6e-7: the eq3f totals above |
| eq1, eq2 (static coupled) | dens, eint, E bitwise |
| ddn (gas_feedback false) | hydro bitwise, E 6e-15 |
| nltek4 (kinematic) | E 1e-12 |
| shd (static) | E 4e-11 |
| ddf (dynamic diffusion) | dens 7e-8, E 1.8e-6 |
| rs_m2 (radiation shock, 512) | 1e-8 |
| cloud | 4e-7 |
| Marshak c dt/dx 1 / 10 | E 7e-6 / 8e-7, F1 7e-4 / 1.4e-4, eint 4.5e-4 / 4.7e-5 (at the front) |

Marshak is the only visible change. Where the realisability clip binds, the clip remainder now enters dA, and
the thin faces keep it as lab flux. In conservative it was charged back to the gas every step. The change is far
below the test's L1(E) criterion (0.02).

## 5. Files

- Code: src/rad_m1/rad_m1.hpp, src/rad_m1/rad_m1_implicit.cpp (27c35877).
- Scripts:
  - build*.sh, smoke.sbatch, rwd4.sbatch, agcar.sbatch, pool.sbatch / run1.sh / arms.sh (RECc2, RECc2be);
  - reg/run_cs.sh (eq3f);
  - ladder/{lad.sbatch, ana/od.py (+ tJePvsD<k>, tJeD<k> arms)};
  - ana/rw/run_radwave.py (+ `ref=` key), ana/rwtab.py;
  - sp/ana/spl.py (+ arm suffixes K2, Z<k>).
- Runs:
  - runs/rw_27c35877, runs/rwbit_27c35877, reg/runs_cs/eq3f/*_fin;
  - ladder/run_fin, ladder/run_w4..w7 (diagnostic bits), ladder/run_mt;
  - sp/runs_fin, agcar/.

## 7. Follow-up (consfix2): the shift is off by default; rlat

Commits on top of prodfinal ed29da1a, which already merges 27c35877:

- **ad8bd227**: ImplicitCs2Shift becomes opt-in, `<rad_m1>/implicit_cs2_shift` (default false).
- **1fbe8621**: diagnostic `implicit_cs2_dbg` bits 8 (no x1 split) and 16 (no x2/x3 split). The default is bitwise
  equal to ad8bd227.

Binaries in bin/ (`*_cf2a` is the ad8bd227 build):

- athena_cpu_1fbe8621 (a copy of cf2b);
- athena_he_cpu_1fbe8621;
- athena_he_a100_1fbe8621.

Runs are in h2/ (h2order FIN decks: time2_start = slope + implicit_vimp_bc; h.py tord). Orders are over the level
pairs 1-2, 2-4, 4-8 (and 8-16).

### 7.1 Shift off vs on

| case | matched FIN | c2, shift off (now the default) | c2, shift on |
|---|---|---|---|
| rw_p10t100, tol 1e-14 | 2.00 (N) | E 2.00 2.00 2.00; F1 1.99 1.96 1.92; v 2.01 2.01 1.95 | E 1.88 1.70 1.48; rho 2.02 1.39 0.68 |
| rw_p10t1 | E 1.79 1.73 1.29; rho 1.99 2.06 2.00 | E 1.81 1.92 2.00; rho 1.91 1.84 1.70; F1 1.86 1.67 1.56 | rho 1.99 1.98 1.94; F1 1.79 1.58 1.08 |
| rsw_T (vet_gd) | E 1.82 1.47; v 1.78 1.29 | **2.00 in every variable** | |
| rw_t10a | 2.0 | E, rho 2.0; F1 1.90 1.98 1.29 | |
| rlat_T | 1.94 / 2.16 | 1.81 / 1.73 | 1.81-1.84 / 1.69-1.77 (N) |

Conservation and AG Car with the shift off:

- eq3f: P_tot 8.9e-16 (hesdirk2 and BE), E 1.13260.
- rlat_om5 L_z: +5.2e-6 (hesdirk2) and +5.5e-6 (BE).
- Reduced AG Car A, 120 cycles (agcar/cf2_*, a deck copy with implicit_cs2_shift named):
  - c2 with the shift off: Picard 2.94, the same as matched (2.98 with the shift on).
  - Budget, L_top and E_rad are identical to matched.
  - Dumps vs matched: v1 1.6e-3 of max, the rest ≤ 7e-7.
  - The shift-on run is bitwise equal to 27c35877.

### 7.2 rlat_T at about 1.8

rlat_T: rho kappa_s 1000, c dt rho kappa about 8-60, theta 0.02-0.1. All variants below have the shift off.

| variant | orders | Picard mean |
|---|---|---|
| c2 | 1.81 / 1.73 | 4.1 |
| tol 1e-14 | identical to c2 | |
| bit 8 (no x1 split, so not the Marshak end face) | 1.81 / 1.72 | |
| bit 16 (no x2/x3 split) | 1.81 / 1.73 | 3.4 |
| bits 8+16 (no split at all) | 1.81 / 1.73 | |
| bit 1 + implicit_rad_inertia false | 1.81 / 1.73 | |
| **bit 1** (no velocity split; the E part split after the solve, as matched / conservative) | **1.96 / 2.12** | 3.17 |

- So the deficit is how the E part a (E' - E^n) of dA is treated: inside the solve, or not at all, gives 1.8; the
  post-solve ImplicitRadInertia rule gives 2.0.
- It is not the shift, the end face or the Picard tolerance.
- Its size is small: the finest-pair error is 2.2e-6 against matched 1.65e-6.
- Possible follow-up, not done: keep the E part on the post-solve rule and put only the velocity part into the
  solve's dA. This needs the vimp prediction adapted to the drag share of the E part.

## 8. A100 smoke (1fbe8621, gpudev 1 node x 4 A100, full-size AG Car A, 60 cycles; gpu/runs/j.31087519)

| run | deck | rc / FATAL / NON-CONV | Picard mean / max | s/cycle (cpu time / 60) | L_top | vs reference |
|---|---|---|---|---|---|---|
| c2P | prod FRESH copy + he_sponge_dmax 2e-16 + conservative2, implicit_cs2_shift = false named | 0 / 0 / 0 | 2.63 / 6 | 0.458 | 2.8430e38 | production smoke |
| mtA | pfmerge agcar_A_dc60 (+ the shift key), matched | 0 / 0 / 0 | 2.80 / 6 | 0.417 | 2.8430e38 | **BITWISE** vs N's ed29da1a mtA (dump data, hst) |
| c2A | same, conservative2 (shift off) | 0 / 0 / 0 | 2.80 / 6 | 0.419 | 2.8430e38 | vs N's c2A (shift on): v1 1.0e-3 of max, F1 1.3e-6, rest ≤ 1e-7 |
| c2As | same, conservative2 with the shift on | 0 / 0 / 0 | 2.94 / 6 | 0.417 | 2.8430e38 | **BITWISE** vs N's c2A |

Notes:

- The bin files differ by bytes only in the parameter header (the extra shift key).
- The energy budget residual is identical in all runs: -4.9e-10 of Etot.
- c2P differs from the others also by he_sponge_dmax.
