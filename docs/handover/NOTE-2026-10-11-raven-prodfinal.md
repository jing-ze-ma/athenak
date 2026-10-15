# NOTE prodfinal (prodfinal-1011): integration of the final candidate, new defaults, AG Car production package

Branch `prodfinal-1011`, worktree /raven/u/jinma/ATHENAK/wt-prodfinal, MPCDF Raven. Work dir /raven/ptmp/jinma/prodfinal_1011.
Base: prodcand2-1011 9b6d572a (code 767d3313, NOTE-2026-10-11-raven-prodcand2.md).

**Status.**
- The regression passes: with the old values named, both new binaries are **bitwise** identical to 767d3313 on all 33
  literature runs and on the reduced CPU AG Car A.
- The conservative gas source conserves total momentum to round-off on eq3f (8.9e-16).
- The FRESH and RESTART packages are ready. Nothing was submitted.
- **Two open items need a user decision:**
  - `he_ic_rad_hse = true` **fails on AG Car B** (§5.3).
  - The A100 package test is still in the gpudev queue (§6).

## 1. Commits

| commit | content |
|---|---|
| 1cb04a19 | merge consrc-1011 26fa1610: opt-in `rad_m1/implicit_gas_src = matched (default) \| conservative \| conservative_full` |
| 530c38e2 | `rad_m1/vet_seg` default **emit** (fresh runs). A restart whose file lacks the key keeps mean. |
| **18ff91c8** | fresh-run defaults (user request, approved): `implicit_thin_switch = fse`, `implicit_thin_corr = hc`, and he_star_m1 `he_ic_rad = sph1d`, `he_ic_rad_gas = re_atm`, `he_ic_rad_hse = true` (§3) |
| (this NOTE) | docs only |

None of these commits is pushed; the coordinator pushes.

**Binaries** (bin/; READY lists the md5s; nothing is overwritten):

| file | md5 | built from |
|---|---|---|
| athena_cpu_530c38e2 | 656e4177 | 530c38e2, CPU, built-in pgens (workers B and C use this one) |
| athena_he_cpu_530c38e2 | b3357e79 | 530c38e2, CPU, PROBLEM=he_star_m1 |
| athena_he_a100_530c38e2 | 0c1a6c14 | 530c38e2, A100, he_star_m1 |
| athena_cpu_18ff91c8 | 737a6282 | 18ff91c8, CPU, built-in pgens |
| athena_he_cpu_18ff91c8 | 7c516866 | 18ff91c8, CPU, he_star_m1 |
| **athena_he_a100_18ff91c8** | **1a55649f** | 18ff91c8, A100, he_star_m1. **The production binary.** |
| athena_he_cpu_98835d99 | ef08e07f | the held chains' code (hrdet2), CPU. Used only to make a held-equivalent rst (§4.3). |

Build jobs: 31078032 (530c38e2) and 31078644 (18ff91c8).

## 2. Defaults (new in prodfinal)

The rule is the same for every row:
- A **restart whose file lacks the key keeps the old value**.
- The resolved value is echoed through GetOrAdd, so later restarts keep it.
- Explicit input overrides the default.

| key | fresh-run default (was) | scope | opt-out |
|---|---|---|---|
| `rad_m1/vet_seg` | **emit** (mean) | vet_sc, vet_gd, vet_col_lat | `vet_seg = mean` |
| `rad_m1/implicit_thin_switch` | **fse** (none) | Only where the switch is accepted and has ray data: `implicit_flux = blend` with `implicit_flux_beam = halfrange`, `implicit_hr_recon` dc or not named, and vet_sc, vet_gd or `implicit_hr_model`. Elsewhere it stays none, silently (central, ap_hll, berthon, a blend without halfrange). | `implicit_thin_switch = none` |
| `rad_m1/implicit_thin_corr` | **hc** when the switch is fse (none) | as above | `implicit_thin_corr = none` |
| `problem/he_ic_rad` (he_star_m1) | **sph1d** (file) | Every fresh he_star_m1 start with `opacity = table`, no `he_wind_ic` and no `he_ic_rad_thin`; file otherwise. The column is rebuilt on every start, restarts included; a restart keeps file. | `he_ic_rad = file` |
| `problem/he_ic_rad_gas` / `he_ic_rad_hse` | **re_atm / true** (column / false) | sph1d | name them |
| `rad_m1/implicit_gas_src` | matched (unchanged; conservative is opt-in) | implicit M1 | - |

The defaults carried over from prodcand2 are unchanged:
- `time2_pred_cap = 0.5`;
- `vet_x1_periodic_inflow = source`;
- `vet_bc_diff = face`;
- `ic_rad = fs_relax`;
- `vet_source_noesrc_h2be`;
- `implicit_hr_face = sf`;
- the moveq trio;
- `time_scheme` hesdirk2 where it is accepted.

**Command-line keys still have to exist in the deck.**

### 2.1 BSG true reproduction (bsg3d_truerepro2*)

Neither Raven deck (truerepro2_1009/inputs/bsg3d_truerepro2_raven{,_hr}.athinput) names `he_ic_rad*` or `ic_rad`. Both
use `opacity = table`. With 18ff91c8 a fresh start would therefore take the sph1d IC.

To keep the original reproduction IC, add to **both** decks (fork bsg-files-1009 docs/handover/bsg-hrdet-1009/*.in;
the coordinator adds them):

```
<rad_m1>
ic_rad = none
<problem>
he_ic_rad = file
```

`ic_rad` (default fs_relax since prodcand) also changes the radiation start of a fresh vet_col + vet_gd run.

Additionally, `bsg3d_truerepro2_hr_lm` (and the `_raven_hr` deck) is blend + halfrange + vet_gd. It would take the new
thin switch, so add `implicit_thin_switch = none` under `<rad_m1>`.

If the decks are ever run with a prodcand / prodfinal binary, the older fresh-run defaults also differ from the
reproduction. The decks name none of these keys:
- `vet_seg = mean`;
- `vet_bc_diff = cell`;
- `implicit_hr_face = cell` (hr deck);
- `implicit_f0_init = zero`, `implicit_rad_inertia = false`, `implicit_nofb_work = false`.

`time_scheme = be` is already named in the decks.

## 3. Defaults commit 18ff91c8: the checks

All CPU runs are on `-p general`.

- **(a) A deck that names the old values (none / file) is bitwise identical to 530c38e2.**
  - The literature suite (reg/jobs_nd, same 4-thread pool): all 46 runs with every key named, plus the hr arms with
    `implicit_thin_switch = none` named. All the new-vs-530c38e2 pairs are **bitwise** (reg/cmp2.txt: 53 of 54
    lines; the one non-bitwise line is an emit-vs-mean Marshak pair, 5e-15).
  - eq3t cen without any key: bitwise.
  - Reduced AG Car A: the deck agcar_A_va_z names thin none / he_ic_rad file.
    - Its own keys (be, hr xthin 30), 30 cycles: **bitwise**.
    - With the production keys, 60 cycles: **bitwise** vs 530c38e2.
- **(b) A fresh AG Car A deck without the five keys reproduces the explicit run bitwise.** The test uses CPU reduced
  480x32x32, 30 cycles, production dt.
  - fresh0_nokey (18ff) = fresh0 (18ff) = fresh0 (530c38e2) = prod/agcar_A_fresh (every key named, 18ff): all
    **bitwise**.
  - The A100 version of this check (frA0 vs frA) is in the pending job 31079211 (§6).
- **(c) A restart from a held-chain rst without the keys keeps the old behaviour.**
  - The held chains have **no rst yet**: runA_rc / runB_rc have never run, and their run.cfg has `RST=` empty.
  - A held-equivalent rst was therefore made from the held chains' code 98835d99 (CPU build) and their deck
    agcar_rcxA_..._mg_hr: reduced, 3 cycles.
  - Restarting it without keys to cycle 8 on 530c38e2 and on 18ff91c8 gives **identical** rst data and hst files
    (job 31079412).
  - The new rst echoes `implicit_thin_switch = none`, `vet_seg = mean`, `vet_bc_diff = cell` and
    `he_ic_rad = file`.
- **Effect of the new thin default** (hr deck without the keys, xthin 30): fse + hc vs none.
  - eq3t: 1.4e-3. Crossing beams r4: 9e-8. eq2: VET noise only.
  - A no-key run is **bitwise** with an explicit `fse`/`hc` run (eq2, cb, eq3t).

## 4. Regression on 530c38e2

### 4.1 Explicit `vet_seg = mean` + `implicit_gas_src = matched` vs 767d3313

The literature suite was rerun with prodcand2's arms and inputs, the same pool and 4 threads (reg/, job 31078169;
Marshak rerun with `implicit_vimp = false` as in prodcand2's c2b). Every pair is **bitwise**:
- crossing beams r4 / r12 / r4 c0.5 / r4 c32 / cbg;
- shadows shd r12g2 / r36g2 / shdg;
- eq1, eq2 (REC and lag), eq2c VC / VCbe, eq3t;
- nltek1-8, sdifk, ddiff / ddn / ddf;
- Marshak c1 / c10;
- radiation shocks M1.2 / M2 / M3 n512, and M2 / M3 n1024.

Reduced CPU AG Car A, agcar_A_va_y = prodcand2's agcar_A_va_x plus `implicit_gas_src = matched`, production keys, fixed
dt, 60 cycles: **bitwise** vs prodcand2 new_candx (767d3313), 14 files.

### 4.2 `vet_seg = emit` effects (default vs mean, same binary)

| test | effect |
|---|---|
| eq1, eq2 | E mean / min / max identical to 6 digits (96.795052, 145.787811). Only VET noise in F. |
| eq3t, nltek4, Marshak c1 | ≤ 4e-12 (eq3t 1e-15) |
| shadow r12g2 | metrics identical: L1 0.011, umbra 20.5 / 165.3 |
| **crossing beams r4** (c dt/dx 256) | final fields agree to 6e-15. **Picard mean 33.3 → 3.1, NON-CONVERGED 6 → 0, 81 s → 6 s.** Under mean, stage 1 fails 6 times (residual 1e25) and falls back to BE; emit removes this. |
| AG Car A reduced (e_o: `vet_seg = emit`, 530c38e2) | **bitwise** vs prodcand2 new_emit. The effect is therefore as in NOTE-prodcand2 §6: envelope E -1.4e-4, E/J closer to 1, Picard +2 %. |

Two cautions:
- **agcar_A_va_x and va_y name `vet_seg = mean`.** The pf_mean / pf_emit pair there is mean in both arms. The emit arm
  is e_o.
- Restarts from files without the key keep mean (§3c).

### 4.3 `implicit_gas_src = conservative`

The eq3f moving equilibrium (consrc inputs: 32 x2-row slices give exact domain totals, c = 10, t = 1) on 530c38e2,
reg/runs_cs:

| arm | E / v | max \|dP_tot\| | max \|dE_tot\| |
|---|---|---|---|
| cen be matched | 1.13410 / 2.95469 | 6.3e-4 | 2.2e-9 |
| cen be conservative | 1.13262 / 2.95537 | **8.9e-16** | 2.0e-9 |
| cen be conservative, tol 1e-12 | 1.13262 / 2.95537 | 8.9e-16 | 3.0e-14 |
| cen hesdirk2 conservative | 1.13262 / 2.95537 | 1.3e-15 | 8.9e-9 |
| tJe be conservative | 1.13262 / 2.95537 | 8.9e-16 | 1.4e-9 |
| tJe h2v (REC: predict + vstage + vimp) matched | 1.13406 / 2.95471 | 6.2e-4 | 9.4e-11 |
| **tJe h2v conservative** (mean / emit identical) | 1.13260 / 2.95538 | **8.9e-16** | 9.4e-11 |
| hr be conservative | 1.12966 / 2.95672 | 1.3e-3 | 1.3e-9 |

What the table shows:
- With the conservative source, eq3 gas + radiation momentum is conserved to **round-off** on every central / tJe arm,
  BE and hesdirk2.
- The energy defect scales with the Picard tolerance, as in NOTE-consrc.
- hr keeps its AP-HLL leak, as in consrc §2.
- On the suite eq3t (single-slice tabs), conservative gives E 1.13260 against 1.13406 matched.

Reduced AG Car A, conservative vs matched (60 cycles):
- max relative difference 3.5e-4, in velx;
- L_top identical to 5 digits (0.99980 at cycle 60);
- Picard 3.35 / 7 in both arms;
- with emit: 1.8e-4 (velx).

## 5. Production package (prod/; NOT submitted)

### 5.1 FRESH decks

prod/agcar_A_fresh.athinput and prod/agcar_B_fresh.athinput start from the held decks
(rc_1009/inputs/agcar_rcx{A,B}_ge_raven_accel_st_mg_hr.athinput) with:
- the production keys: `time_scheme = hesdirk2`, `time2_vet_col = predict`, `implicit_vimp = true`,
  `time2_vstage = true`, `implicit_blend_xthin = 0`, `implicit_thin_switch = fse`, `implicit_thin_corr = hc`,
  `vet_bc_diff = face`, `ic_rad = fs_relax`;
- `vet_seg = emit`, `implicit_gas_src = matched`, `time2_pred_cap = 0.5`, `implicit_hr_excess = none`;
- the new IC: `problem/he_ic_rad = sph1d`, `he_ic_rad_gas = re_atm`, `he_ic_rad_hse = true` (**B: false, §5.3**);
- **every remaining `<rad_m1>` and `<problem>` key** named at the binary default, from the PAR_DUMP of a fresh run:
  122 rad_m1 and 19-20 problem keys (prod/mkdeck.py).
  - `<hydro>` / `<time>` were left as they are. Naming their dumped defaults would include derived values such as
    `tfloor = FLT_MIN` next to `tfloor_kelvin`.

Production settings unchanged from the held decks:
- adaptive dt (`implicit_cfl = -1`, cfl 0.3);
- tlim (A 8.064e6 s, B 1.3e6 s);
- outputs (hst 1e4 / 1e3 s, bins 1e5 / 1e4 s, rst 5e5 / 5e4 s).

The run configs are prod/run.cfg.fresh_{A,B}. They set BIN = athena_he_a100_18ff91c8 with its md5, IN = the deck, and
`RST=`. They are for rc_raven.sh in a NEW run dir.

CPU checks (reduced 480x32x32):

| deck | binary | cycles | rc | Picard mean / max | NON-CONV | result |
|---|---|---|---|---|---|---|
| A fresh0 (fresh-only keys) | 530c38e2 | 100 | 0 | 3.03 / 6 | 0 | dt 1252 s |
| A, full deck | 18ff91c8 | 30 | 0 | - | - | **bitwise** vs fresh0 |
| B (hse false) | 18ff91c8 | 100 | 0 | 2.94 / 8 | 0 | dt 242 s; **bitwise** vs 530c38e2 |

### 5.2 RESTART

- The held chains have **no rst** (§3c), so there is no newest rst to switch.
- prod/mkrestart.py `<rst> <fresh deck>` prints the KEYS that switch a restart to the production set: every
  rad_m1 / hydro / problem key that differs from the rst's PAR_DUMP or is missing there.
  - IC-only keys are skipped.
  - The restart-aware defaults come out explicit: vet_seg emit, vet_bc_diff face, pred_cap 0.5, thin fse / hc.
- Templates: prod/run.cfg.restart_{A,B}.template.

Test on a COPY of the held-equivalent rst (CPU reduced, 18ff91c8, job 31079412, runs/agc/j.31079412/rsH_p):
- 27 keys (restart_keys_A_heldequiv.txt), from cycle 3 to cycle 13;
- rc 0, Picard 5.21 / 7, 0 NON-CONV (the first step after the restart is BE);
- the rst echoes hesdirk2 / emit / face / fse / hc / pred_cap 0.5 "added on the command line of a restart", and keeps
  `he_ic_rad = file`.

The B evolved state:
- prod/rstcopy/runB.00010.rst is a copy of agcar_1008/runB at t = 4.624e5 s (central flux, older keys).
- prod/restart_keys_B_runB.txt has 61 keys, including new opacity tables.
- Its A100 restart test is in pending job 31079211 (rsB).

### 5.3 Finding: `he_ic_rad_hse = true` fails on AG Car B

- On AG Car B (hot state, log L 6.17), the sph1d IC with `he_ic_rad_hse = true` stops at init with "sph1d: the VET
  iteration failed" (he_star_m1.cpp:948, after the first HSE re-integration). This happens for every taub (300 / 3000 /
  1e4) and for dtau 0.02.
- The same IC with **hse false** works: re_atm, re and column each take 12 VET passes, 0.5 s.
- A plausible cause (not proven): the column's radiative force exceeds gravity in the Fe-opacity zone
  (Gamma_es ≈ 0.7). The HSE march at fixed T then gives a density inversion.
- **Consequence:** with 18ff91c8, a fresh he_star_m1 start of AG Car B **without the IC keys stops at init**.
- The B FRESH deck therefore uses `he_ic_rad_hse = false` (prod/agcar_B_fresh.athinput, marked in the deck).
- **User decision:** accept B with hse false, or have the HSE march fixed for super-Eddington layers.
- A: hse true works (29 outer iterations).

## 6. Pending (gpudev queue: the next start was about 2 h out when they were submitted; the gpu partition is a day out)

| job | content |
|---|---|
| **31079210** (pf_g1, gpudev 15 min, A100) | mkH (held binary hrdet2 A100 + held deck, 3 cycles → rst); rsH_o / rsH_n (530c38e2 vs 18ff91c8 restart without keys, bitwise expected); rsH_p (RESTART package keys, to cycle 13) |
| **31079211** (pf_g2, gpudev 15 min, A100) | frA (prod/agcar_A_fresh) vs frA0 (no-key deck), bitwise expected; frB (prod/agcar_B_fresh); rsB (runB rst copy + restart_keys_B_runB, +2e3 s) |

Outputs go to gpu/runs/j.<id>. Compare with:
- ana/cmp_pairs.sh (bin/hst);
- ana/rstcmp.py (rst data after the PAR_DUMP).

## 7. Files (/raven/ptmp/jinma/prodfinal_1011)

Scripts and inputs:
- build.sh (env T);
- reg/ (run1b.sh, rb.sh, rbi.sh, run_cs.sh, pool*.sbatch, jobs_*.txt, inp_nokey/);
- agc/ (agc_cpu.sh, agc_nd.sh, agc_x.sh, agc_e.sh, rs2_cpu.sh, pairs_agc.txt);
- icb/icb.sh (B IC variants);
- gpu/gpu_g1.sh, gpu/gpu_g2.sh;
- inp/ (agcar_A_va_y / va_z, agcar_{A,B}_fresh0, fresh0_nokey, fresh0h, B_icb);
- prod/ (decks, run.cfg.*, mkdeck.py, mkrestart.py, restart_keys_B_runB.txt, rstcopy/).

Analysis:
- ana/cmp_pairs.sh, ana/eqmom.py, ana/rstcmp.py;
- reg/cmp2.txt.

Runs:
- reg/runs, reg/runs_cs;
- runs/agc/j.31078170, j.31078671, j.31078941, j.31078942, j.31079251, j.31079412;
- icb/runs/j.31078864.
