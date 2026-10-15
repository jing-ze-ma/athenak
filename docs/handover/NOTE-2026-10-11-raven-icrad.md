# icrad-1011: a VET-consistent radiation start for AG Car A

Branch `icrad-1011` in worktree `/raven/u/jinma/ATHENAK/wt-icrad`, based on vetacc-1011 bb1ebff8.
| commit | change | binaries in `bin/` |
|---|---|---|
| **7045bc14** | the opt-in `ic_rad` keys | `athena_he_a100_7045bc14`, `athena_he_cpu_7045bc14` |
| **99fb1407** | `ic_rad = fs_relax` becomes the default when eligible; relax_n 40, relax_tol 1e-4 | `athena_he_a100_99fb1407`, `athena_he_cpu_99fb1407`, `athena_cpu_99fb1407` (no pgen) |
| **32546788** | vet_source_noesrc in every backward-Euler solve under hesdirk2 (fixes the hesdirk2 dip); pskip is off during the IC passes | `athena_he_a100_32546788`, `athena_he_cpu_32546788`, `athena_cpu_32546788` |

None of these is pushed.

## 1. Problem (coordinator diagnosis, re-checked)

The IC column sets the envelope E with an Eddington/diffusion-like relation. Taking shell means
from vetacc `j.31064135/be`, E(0)/E(50) is 1.08 at 0.99 R_ph, 2.4 at 1.009, 2.6 at 1.05–1.1,
2.0 at 1.6 and 1.35 at 2.9 R_ph. F1 is already right. The gas also starts too hot in the
envelope: eint(0)/eint(50) is 1.27 at R_ph, 1.65 at 1.05–1.1, 2.3 at 1.6 and 35 at 2.9 R_ph.
The relaxed run has S_fs = E to 1e-3, so the gas follows the radiation down. The surplus leaves
within about one light-crossing time, which produces the +54 % L_top spike and a transient of
about 15 cycles.

The steady state of the discrete M1-VET system is **not** E = J. At cycles 50 and 200,
<E>/<J_fs> is 1.000 below 0.99 R_ph, about 1.021 at 1.02–1.1 R_ph and about 1.010 at 2 R_ph.

## 2. Implementation (opt-in; keys read only when named; default bitwise)

`src/rad_m1/rad_m1_icrad.cpp` adds a hook in `Driver::Initialize` (Step 2b). It runs only when
`res_flag` is false, so **never on a restart**. It runs after the first dt and before the
cycle-0 outputs. Afterwards Steps 1 and 2 are redone (boundaries, primitives, dt).

| key `<rad_m1>/` | meaning |
|---|---|
| `ic_rad` | `none` / `fs` / `relax` / `fs_relax`. Default since 99fb1407: `fs_relax` when transport = implicit, closure = vet_col, vet_gd and a fluid are all present; `none` otherwise |
| `ic_rad_fs_iter` (3) | (a) E := J passes; the tensor is rebuilt after each pass |
| `ic_rad_relax_n` (40), `ic_rad_relax_tol` (1e-4) | (b) max steps; stop at \|dL_top/L_top\| < tol |
| `ic_rad_relax_dtfac` (1) | relaxation dt / run dt |
| `ic_rad_diag` (0) | N: append shell means of E, J_tw and F1 to `icrad_ej.txt` at every IC pass and at the first N builds of the run (read-only) |

* **(a) fs.** This is one ordinary m1 task-list pass. Closure limits, opacity, the
  start-of-step vet_col/vet_gd build, then `IcRadHook`: E := J (u0, iw EN/EP, halos), rebuild,
  repeated `fs_iter` times. The solve is then skipped, and the fluid refresh and ghost exchange
  of the list run as usual. J is twin-corrected under the fused twin, as vgd_hr's J:
  J *= <J_tw>/J_tw. The shell means use a fixed order plus a rank-ordered Allgather sum, so the
  result is run-to-run reproducible on GPU and does not depend on the decomposition. The IC flux
  and the gas are kept. The first gd shell keeps its E.
* **(b) relax.** Up to N implicit radiation steps (ordinary BE m1 task lists,
  `t2_solve = none`) with no hydro stages and no advance of time or cycle. For these steps only,
  `dbg_gas_force(_trans) = false` and `implicit_vimp = false`, so the gas momentum and velocity
  stay the IC's. The gas–radiation **energy** exchange stays on, so the envelope gas relaxes
  together with E. The esrc (frozen MLT deposit) is a function of t only, so it stays the IC's.
  The stopping test is L_top = Σ f0x1(ie+1) A, a reproducible sum. Afterwards the predictor
  history is dropped (`pred_ok = pred2_ok = false`) and the flags are restored.

**Choice.** Both are implemented. The default is fs_relax, (a) then (b), by the user's decision after the restart check. (b) carries the result:

* It converges to the discrete solver's own steady state for the frozen hydro, including the
  2 % E/J offset and the gas temperature.
* (a) alone leaves the gas hot: E/J is about 1.01 right after the passes, but the first step heats E back and L_top still overshoots by +8.8 % (see sect. 3).

(a) in front costs 0.6 s and gives the same result as relax alone, to 1 %.

## 3. Results (AG Car A, fixed dt 1239.996 s)

Sources:
* **DeltaAI** (GH200, binary 7045bc14), jobs 3357130 and 3357234. Notes are on origin/bsg-files-1009 under
  `docs/handover/NOTE-2026-10-11-deltaai-agcarA-{face-icrad,icrad-restart}.md` and
  `deltaai-agcarA-1011/analysis.txt`.
* **Raven** GPU job 31067115 (A100, 7045bc14), arms `befrr` (fs_relax with restart files), `berx`, `hpfr`, `rs0` and
  `rs10`. Analysis is in `ana/ana_31067115.txt` and `ana/gas_befrr.txt`.
* The default-IC references are vetacc `j.31064135/{be,beH,hp}` and `j.31064136/hpH`.

An earlier Raven attempt, 31066834, failed at startup: `rad_m1/ic_rad` cannot be set on the command line unless the
key already exists in the deck. The deck copy `inp/agcar_A_icr.athinput` now names it.

### 3.1 L_top per cycle (`fig/ltop_60.png`)

The ratio is L_top/L_in, from Raven be runs.

| cycle | 0 | 1 | 2 | 3 | 4 | 5 | 10 | 30 | 60 |
|---|---|---|---|---|---|---|---|---|---|
| default | 1.0000 | **1.5356** | 1.4449 | 1.2491 | 1.1083 | 1.0403 | 1.0014 | 1.0004 | 1.0003 |
| relax | 1.0006 | 1.0006 | 1.0006 | 1.0005 | 1.0005 | 1.0005 | 1.0004 | 1.0003 | 1.0003 |
| fs_relax | 1.0006 | 1.0006 | 1.0005 | 1.0005 | 1.0005 | 1.0005 | 1.0004 | 1.0003 | 1.0003 |

DeltaAI gives the same picture. fs alone cuts the spike only to +8.8 %, because it leaves the gas hot and E/J
about 1.15 off. relax and fs_relax are flat to 6e-4 from cycle 0. All variants agree within 7e-6 by cycle 60.

### 3.2 E(0)/E(50) and flux factor (`fig/e0e50_ff.png`)

| r/R_ph | 0.99 | 1.0 | 1.009 | 1.05 | 1.1 | 1.6 | 2.0 | 2.9 |
|---|---|---|---|---|---|---|---|---|
| E(0)/E(50), default | 1.080 | 1.087 | 2.395 | 2.611 | 2.615 | 2.004 | 1.687 | 1.355 |
| E(0)/E(50), fs_relax | 1.015 | 0.974 | 1.001 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| F/cE at cycle 0, default | 0.034 | 0.24 | 0.25 | 0.26 | 0.28 | 0.45 | 0.56 | 0.73 |
| F/cE at cycle 0, fs_relax | 0.036 | 0.27–0.54 | 0.59 | 0.69 | 0.74 | 0.91 | 0.95 | 0.99 |

The fs_relax flux factor at cycle 0 equals the relaxed default's at cycle 50 at every radius in this table.
F(0)/F(50) = 1.000.

### 3.3 E/J_tw (shell means) at the first build of cycle 0 (`fig/ej_c0.png`, `icrad_ej.txt`)

| band (R_ph) | 1.0–1.05 | 1.05–1.5 | 1.5–2 | 2–3 |
|---|---|---|---|---|
| default IC (the first build = fs0 / relax0) | 1.966 | 2.045 | 1.520 | 1.196 |
| after fs (3 passes) | 1.011 | 1.010 | 1.009 | 1.009 |
| after relax / fs_relax: cycle 0+ build | 1.0206 | 1.0181 | 1.0095 | 1.0125 |
| relaxed default (cycle 50/200, m1_fs) | 1.020 | 1.019–1.021 | 1.011 | 1.013 |

The relaxed IC reproduces the solver's own steady E/J offset of about 2 %. E := J alone gives 1.01 and leaves the
gas hot.

### 3.4 Structure at cycle 50 vs the relaxed default

Metric: shell-mean relative difference (L1 in brackets), fs_relax vs default be.

| band | < 0.8 R_ph | 0.8–0.95 | 0.95–1.05 | > 1.05 |
|---|---|---|---|---|
| dens | 2e-7 (1e-5) | 4e-6 (3e-5) | 9e-5 (6e-5) | 1.3e-2 (8.5e-4) |
| eint | 1.5e-7 | 3.7e-6 | 3.3e-5 | 1.3e-2 |
| m1_e | 1.8e-7 | 5.8e-6 | 1.3e-5 | 3.0e-5 |
| m1_f1 | 1.8e-4 | 7.3e-5 | 1.3e-5 | 1.2e-5 |

The interior and the photosphere are the relaxed default's to ≤ 1e-4. The ~1 % in dens and eint above 1.05 R_ph is
floor-density gas (P_rad/P_tot = 1.00 there).

### 3.5 Envelope gas and hydrostatics (`ana/gas_befrr.txt`, DeltaAI addendum)

At cycle 0, rho and v are bitwise those of the default. eint changes by up to −39 % at 0.95–1.05 R_ph, −47 % at
1.05–1.5 and −97 % at 2.9 R_ph. This is the same state the default IC reaches by itself by cycle 50–60 (ratio
0.997–1.009).

E falls 62 % (P_rad). Above 1.05 R_ph, P_rad/P_tot = 1.00, so the gas-pressure change is ≤ 1e-3 of P_tot there. In
0.95–1.05 R_ph it is ≤ 0.22 P_tot, an upper bound from (2/3)eint. F1 changes by < 0.6 % below 0.95 and by 4 % at
0.994 R_ph. The radiative force follows F and does not change (rho κ F/c).

There is **no velocity transient**. KE_int is 9 % below the default at cycle 1 and 0.6 % below at cycle 60, and both
grow smoothly from the seed. max|v_r| at cycle 50 is 4.76 vs 4.78 km/s (0.95–1.05) and 7.58 vs 7.32 km/s (>1.5).

### 3.6 dt convergence of L_top over cycles 1–30 (`fig/ltop_dtconv.png`, DeltaAI dt/2 arms)

Metric: max |L(dt) − L(dt/2)|/L. Default **8.2e-2**, fs 9.4e-3, relax 2.65e-4, **fs_relax 2.63e-4** (rms 8.3e-5).
The Raven default pairs give the same: be 1.07e-1 and hp 1.77e-1 over cycles 1–30, falling to 7e-6 and 4e-5 after
cycle 30. With the relaxed start the first 30 cycles become a meaningful dt test.

### 3.7 Cost

| arm | IC step | relax steps | breakdown |
|---|---|---|---|
| fs_relax, A100 (Raven) | 6.1 s | 14 | fs 0.6 s + relax 4.5 s |
| relax, A100 | 5.4 s | 14 | — |
| fs_relax, GH200 | 4.3 s (dt), 5.5 s (dt/2) | — | — |

That is a one-off cost of about 15–20 cycles of run time.

### 3.8 Restart

* **Raven 31067115.** rs0 restarts from the post-IC cycle-0 rst of befrr and rs10 from its cycle-10 rst. Neither
  runs the IC passes (the log says `[RESTART: the IC passes are skipped]`). All 10/10 restarted hst rows are bitwise.
  The only differing row is the end-of-run duplicate, whose Picard column is 0, as on every run.
* **DeltaAI 3357234.** Bitwise, including the cycle-20 bins and the rst.
* **Old restart files (no key) under the new default (99fb1407).** Reduced CPU AG Car A, job 31067568: an old-binary
  (7045bc14) run wrote rst files, and the new binary restarted from cycle 2. It does not run the IC passes, and cycles
  3–4 are bitwise vs the straight old run. A fresh start of the production deck with no key runs fs_relax.
* **CPU bitwise of non-eligible problems.** No-pgen builds 99fb1407 and 32546788 vs 42aa6bbc, jobs 31067321 and
  31067999: atm_tJe, ba0_tJe, g1_hr, pulse_tJe (hesdirk2, vet_sc, noesrc off), pulsebe_tJe and linear_wave_hydro are
  bitwise in all bin, hst and tab data. The hook is only armed for implicit + vet_col + vet_gd + fluid runs.

## 4. hesdirk2 + relax dip: diagnosis and fix (commit 32546788)

**Observation** (DeltaAI hprx; Raven hpfr with 7045bc14). Under hesdirk2 + predict + vimp + vstage, the relaxed
start dips by −11 % in L_top at cycle 4 and rings over cycles 2–8.

* The relax takes 23–25 steps and ends at 1.0019 L_IC (under be: 14 steps, 1.0006).
* The IC E/J profile shows the cause: E/J_tw is **0.12** in the whole envelope after the relax. The formal J is 8 S
  even though E = S.

**Cause.** fixbundle F1 (`vet_source_noesrc`: the VET source reads the physical E^n = u0, not the solve's old vector
E^n + dt·esrc) was applied only when `time_scheme == be`. Under hesdirk2, every backward-Euler solve still built S
from E^n + dt·esrc, where esrc is the frozen-MLT deposit, large at the CZ top (0.98–0.995 R_ph). That makes S ≫ E
there, so the envelope is lit by a bright layer. The affected solves are the BESTORE start step, the fallbacks and
the IC passes.

Reduced CPU AG Car A (480×32×32, job 31067406) confirms it at the very first build. be with
`vet_source_noesrc=false` gives exactly the hesdirk2 fs0 numbers: max|E/J−1| 0.893, top 0.658, vs 1.053 / 1.026.

The hesdirk2 STAGE solves are clean, because Time2VetColAt puts the predicted state, without esrc, into iw EN.

**Fix.** `VsrcNoes()` applies noes to every solve except the hesdirk2 STAGE solves. The opt-out is
`vet_source_noesrc_h2be = false`. In addition, pskip is off during the IC passes.

Reduced CPU (job 31067567, `fig/ltop_h2_noesfix_cpu.png`), L_top/L_in for cycles 0–6:

| arm | 0 | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|---|
| be + fs_relax | 1.0006 | 1.0006 | 1.0005 | 1.0005 | 1.0005 | 1.0005 | 1.0005 |
| hesdirk2 + fs_relax, **fixed** | 1.0006 | 1.0006 | 1.0005 | 1.0005 | 1.0005 | 1.0005 | 1.0005 |
| hesdirk2 + fs_relax, old | 1.0020 | 1.0019 | 0.9877 | 0.9276 | **0.8906** | 0.9803 | 0.9894 |
| hesdirk2, default IC, fixed | 1.0000 | 1.5283 | 1.5499 | 1.0979 | 0.9221 | 1.0389 | 0.9803 |
| hesdirk2, default IC, old | 1.0000 | 1.5626 | 1.5891 | 1.0857 | 0.8549 | 1.0406 | 0.9830 |

* With the fix, the hesdirk2 relax is the be relax: 14 steps, 1.000625 L_IC, identical numbers. The fixed hesdirk2
  run stays flat to 2e-4 over 12 cycles.
* The old reduced run reproduces the GPU dip (0.891 vs 0.891).
* With the default IC, hesdirk2 still rings after its +55 % spike. This is the start-up overshoot of the spike
  itself (DeltaAI/coordinator point 1); the noes fix reduces the dip from −14.5 % to −7.8 %.
* The 2nd-order ladders are not touched: STAGE solves are unchanged, and only BE solves of noesrc-on runs change.
* A GPU confirmation of hesdirk2 + fs_relax with 32546788 at full resolution is still open.

## 5. Files

| what | where |
|---|---|
| figures | `fig/ltop_60.png`, `fig/e0e50_ff.png`, `fig/ej_c0.png`, `fig/ltop_dtconv.png`, `fig/ltop_h2_noesfix_cpu.png` |
| analysis | `ana/ana_icr.py`, `ana/gas_icr.py`, `ana/ana_31067115.txt`, `ana/gas_befrr.txt` |
| GPU runs | `runs/j.31067115/` |
| CPU runs | `cpu/runs/j.31067284` and `j.31067406` (noes isolation), `j.31067567` (dip fix), `rst.31067568` (old restart) |
| CPU bitwise checks | `cpu/verify/` (logs `cpuv.31067321`, `cpuv.31067999`) |
| scripts | `gpu_icr.sh`, `cpu/dbg1.sh`, `cpu/rstold.sh`, `cpu/verify.sh` |
