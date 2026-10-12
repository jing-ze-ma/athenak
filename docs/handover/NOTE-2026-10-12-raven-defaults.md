# NOTE defaults-1012: the signed-off RT method (method-1012) as FRESH-RUN defaults

Worker D2, 2026-10-12, MPCDF Raven. Work dir /raven/ptmp/jinma/defaults_1012. Code: method-1012
(merge of vimpbcdef-1011 2dd7e434 + commit 22f7021d "rad_m1 defaults-1012"), merged into prodfinal-1011.

## 1. Rule (every key below)
- FRESH run, key not named -> the method (FIN) value, only where the key acts.
- RESTART whose file lacks the key -> the old value (bitwise older runs). 4253af30 rst files already
  carry every key 4253af30 GetOrAdd-echoes (vet_seg emit, implicit_gas_src, thin_recon/bq, ...), so they keep it.
- The resolved value is echoed (GetOrAdd / SetReal), so later restarts keep it. Explicit input always wins.

## 2. Default table

| key | old default (4253af30) | new fresh default | restart, file lacks key | applies to | code |
|---|---|---|---|---|---|
| time_scheme | hesdirk2 (m1-defaults2) | hesdirk2 (unchanged) | be | transport implicit, Cartesian / sp wedge, not tau closure, rk2 | rad_m1_implicit.cpp ImplicitInit |
| time2_vet_sc | lag (read when named) | **predict** (unify-1012, on method-1012) | lag | vet_sc, no time2_vet_extrap, vet_sc_every 1 | rad_m1_time2.cpp |
| time2_vet_col | predict on sp wedge | predict (unchanged) | lag | vet_col on the sp wedge | rad_m1_time2.cpp |
| time2_vstage | true on the sp wedge only | **true on every mesh** (unify-1012) | false | every mesh (acts under vimp) | rad_m1_time2.cpp |
| implicit_vimp | true where valid | **true where valid AND `<time>/evolution = dynamic`** (new: kinematic/static -> false) | false | hesdirk2, multi-D, bicgstab, ng>=2, hydro + coupling + gas_feedback + dbg_gas_force | rad_m1_implicit.cpp:2311 |
| time2_pred_cap | 0.5 | 0.5 (unchanged) | 0 | predict/rebuild builds | rad_m1_time2.cpp |
| time2_start | be (read when named) | be (unchanged) | be | hesdirk2 | rad_m1_time2.cpp |
| implicit_gas_src | matched | **conservative2** | matched | hydro + coupling + gas_feedback (else matched) | rad_m1_implicit.cpp:2509 |
| implicit_cs2_hll | drag | drag (unchanged) | drag | conservative2 | |
| implicit_cs2_shift | false | false (unchanged) | false | conservative2 | |
| implicit_vimp_bc | false (read when named) | **true** (vimpbcdef-1011 2dd7e434, merged) | false | where implicit_vimp is on | rad_m1_implicit.cpp:2593 |
| implicit_blend_xthin | 0 | 0 (unchanged) | 0 | blend flux | |
| implicit_thin_switch | fse | fse (unchanged) | none | blend + halfrange beams + (vet_sc / vet_gd / hr model), hr_recon dc | |
| implicit_thin_corr | hc | hc (unchanged) | none | thin_switch fse | |
| implicit_thin_hcsf | (new key) j | j (unchanged) | j | thin switch on | |
| implicit_thin_recon | dc | **jse** (unify-1012) | dc | thin switch on | |
| implicit_thin_js_gd | (new key; js inert on vet_gd) | **x1** (unify-1012) | none | js/jse on vet_gd | |
| implicit_thin_bq | cell | **lin** (unify-1012) | cell | thin switch on | |
| implicit_thin_bq_gd | (new key) | cell (unchanged) | cell | bq lin on vet_gd | |
| implicit_thin_hcface | (new key, read when named) | **tau** | mean | thin_corr hc / hcs (else read when named) | rad_m1_implicit.cpp:2091 |
| implicit_thin_hcface_gd | (new key) | x1 (unchanged) | x1 | hcface tau on the sp wedge | |
| vet_seg | emit | **jump** | mean (files since prodfinal-1011 echo emit) | VET closures (vet_sc, vet_col/vet_gd) | rad_m1.cpp, rad_m1_vet.cpp |
| vet_bc_diff | face | face (unchanged) | cell | VET closures | |
| vet_gd_symref | (new key, read when named) | **col** | none | vet_gd with vet_gd_twin, no twin_full, no beam | rad_m1.cpp:1161 |
| time2_vet_col_source | relaxed (read when named) | **blend** | relaxed | vet_col (incl. vet_gd) with time2_vet_col predict / rebuild | rad_m1_time2.cpp:264 |
| vet_gd_rotate_dt | (new key) none: no rotation unless vet_gd_rotate_every named | **auto = 20 x the first step's dt** (recorded as the number) | no rotation (old) | fresh vet_gd naming neither rotate_dt nor rotate_every, no beam / async | rad_m1.cpp:1004, rad_m1_vetgd.cpp:787 |
| vet_gd_rotate_align | true (with rotate_dt) | true (unchanged) | true | physical-clock rotation | |
| vet_gd_rotate_stretch | 0 | **0.02** | 0 | physical-clock rotation | rad_m1.cpp |
| vet_gd_rotate_slope | keep | keep (unchanged) | keep | physical-clock rotation | |

Notes on the table:
- **implicit_vimp** (behaviour change beyond the FIN list): the default now also needs `<time>/evolution = dynamic`
  ("where gas live"). Before, a fresh kinematic deck with gas_feedback (only the crooked pipe in the suite) took
  vimp = true; the FIN harness named false there. Restarts keep the echoed value.
- **vet_gd_rotate_dt = auto** is implemented (clean and deterministic): resolved in `VetGdAlignDt` from the dt of the
  first `Mesh::NewTimeStep` (after its global MPI min, before any alignment / tlim clamp; same on every rank),
  `vet_gd_rotate_dt = 20 x dt`, written into the parameters with `SetReal` (round-trip digits) and printed
  (`<rad_m1> vet_gd_rotate_dt = auto: 20 x the first dt ... = <value> (recorded ...)`). The first restart file
  already carries the number, so restarts keep it. `vet_gd_rotate_dt = auto` may also be named. With ic_rad the
  value comes from the pre-IC dt (the first NewTimeStep). Until resolved, `VgdRotOn()` is true and the block is 0
  (exactly as with a named value at t = 0), so the run is bitwise the named-value run.
  - Naming `vet_gd_rotate_every` (any value, incl. 0) is explicit input: the cycle clock (or no rotation) is kept,
    no auto. Production decks still name `vet_gd_rotate_every = 10` -> they must DROP it or name `vet_gd_rotate_dt`.
  - Ladders that need the same rotation times on every level (T ladders: F used 10 x dt0) must keep naming
    `vet_gd_rotate_dt`; auto = 20 x each run's own first dt.
  - Steady sp tests whose first dt differs from the later dt (core_u: first dt 52.1, F's 20 x radiation dt = 800)
    get a different rotation interval than F's harness; production (AG Car A: first dt 1252.56 -> 25051.19 s,
    F named 2.5e4) is essentially F's choice.
- Two-stage restarts with a changed cfl (sp patm _tb/_s stage B) keep stage A's recorded rotate_dt unless named.
- Existing decks that name only SOME old values (the prodfinal production deck run/agcA/agcar_A_fresh.athinput:
  names matched, emit, dc, cell, rotate_every 10, but not hcface / symref / time2_vet_col_source / js_gd) give a
  hybrid on a fresh start: the unnamed keys take the new defaults. Production decks should be rewritten with the
  FIN values (or with the keys removed).

## 3. Verification (binaries `*_22f7021d` = method-1012 source with the defaults; CPU sbatch -p general -c 72)

### (a) fresh deck, no method key named (U) == FIN keys named: BITWISE everywhere
U = F's FIN deck / command line with every FIN key removed (and on vet_gd decks also vet_gd_rotate_every / _t0);
setup keys (closure, faces, beam, blend idort_f, vet_ring, ...) stay. N = F's FIN deck with rotate_dt = U's auto
value and stretch 0.02 (= every key named at the new default). FX = U + F's rotate_dt and stretch 0 = F's FIN arm
literally, compared with F's own outputs (78b2c961).

| test | comparison | result |
|---|---|---|
| Cartesian cb r4 (vet_sc, static gas) | U vs F runs_fin/cb/FIN_r4 | BITWISE (4 dumps) |
| Cartesian eq3t (moving gas) | U vs F FIN_t | BITWISE (402 files) |
| cbg Gauss 8x32 | U vs F FIN_g8x32 | completes, BITWISE (4 dumps) |
| crooked pipe c220 (kinematic + gas_feedback: vimp default off) | U vs F FIN_c220 | BITWISE on the 42 files to t = 1.35 of 10 (run cut: F's took 1.8 h) |
| sp patm_L16_s 64 (two-stage, stage B a restart of the new binary's own file) | U vs N; FX vs F | BITWISE; BITWISE |
| sp rsw_T l 1 | U vs N; FX vs F | BITWISE; BITWISE |
| sp T6b core_u 64 | U vs N; FX vs F | BITWISE; BITWISE |
| reduced CPU AG Car A, 480x32x32, 4x18, 60 cyc | U vs N (14 files) | BITWISE; Picard 2.722 / 6, 0 NC; auto rotate_dt 25051.19 |
| same | FX vs F's FIN (dumps 0-3, hst rows) | BITWISE on 11 of 12; the 12th is hydro_w #1 = an output-schedule artifact (my output2 dt 1e30 -> dump at t = 7.5e4, F's at 1e5) |

### (b) restart from an rst written by 4253af30 (no new key in the file)
- Decks: rsw_T PRODd (suiteF sp/runs_d HFfW, + rst at t = 1, m1 bin every 0.5); AG Car A reduced R deck =
  F's agcA_OLD minus the keys 4253af30 does not know (symref, js_gd, hcsf) and minus time2_vet_col_source; rst at
  3.75e4 (cycle 30).
- **New binary restart (RN) == 4253af30 restart (R4): BITWISE** on both (rsw_T 5 files, AG Car 9 files). The
  restarted new binary takes every old value (vet_seg emit from the file, hcface mean, symref none, relaxed source,
  cycle-clock rotation, no stretch).
- **4253af30 itself is not a bitwise continuation** (R4 vs C4 = 4253af30 continuing): AG Car A m1_e 2e-12 (t 5.0e4)
  -> 1.7e-8, F 6-8e-8 relative to max (t 7.5e4); Picard 3.08 vs 3.02 over cycles 30-60; rsw_T E 2.4e-6, F2/F3
  (lateral, ~0) O(1). Pre-existing in 4253af30 (restart = backward-Euler first step without the stored slope;
  the vet_gd stored intensities are re-converged at the first build), not caused by defaults-1012. The criterion
  "new binary restart == 4253af30 continuing" can therefore only be met as "== 4253af30 restarting".

### (c) A100 full AG Car A production deck, keys named vs unnamed
(see sect. 4)

### (d) old explicit-key decks reproduce 4253af30
- cb OLD r4, eq3t OLD t (F's OLD arm): BITWISE vs F's 4253af30 sanity runs (4 / 402 files).
- F's sp rlat_T HFfWO l 1, 2 and AG Car A agcA_OLD as they are: NOT bitwise (dump 0 already differs on AG Car,
  through ic_rad). Cause: these decks did not name `implicit_thin_hcface` (and the sp decks not `vet_gd_symref`,
  `time2_vet_col_source`) because those were opt-in; they now default to tau / col / blend. With those keys named at
  their old values (mean / none / relaxed): **rlat_T l 1, 2 BITWISE (27 files each), AG Car A OLDx BITWISE (14 files,
  Picard 3.016 / 6 = 4253af30)**. A deck reproduces 4253af30 only if it names every changed key (table).

## 4. A100 (c)
At the merge the A100 check was still queued: Raven gpudev job 31104217 (gpu.sbatch) and the DeltaAI race
(TASK-2026-10-12-deltaai-defaults, bsg-files-1009 docs/handover/defaults-1012/). Both run the full production deck
(prodfinal run/agcA/agcar_A_fresh.athinput copy, 480x128x128, 4 GPUs, 60 cycles): GU = every method key and
vet_gd_rotate_every removed, GN = every key named (rotate_dt = GU's auto value, stretch 0.02); bitwise m1 / hydro_w
dumps every 2.5e4 s + hst, Picard, s/cycle. Result: logs/gpu.31104217.out (and the DeltaAI NOTE), appended by the
coordinator.

## 5. Files
- /raven/ptmp/jinma/defaults_1012: build.sh, small.sbatch, small2.sbatch, agc.sbatch, old2.sbatch, gpu.sbatch;
  cart_a.sh, sp_a.sh, old_d.sh, rst_b.sh; ana/deck.py (strip / set the FIN keys), ana/cmpout.py (bitwise run-dir
  compare: bin data, tab, hst rows by time); agc/inp (U, N, FX, R, OLDx), gpu/inp (GU, GN), old/inp (rlat_T OLDx);
  logs/ (a_*, b_rsw, d_small, old2, agc, gpu); deltaai/ (TASK, job.sub, decks).
- Binaries: bin/athena_{cpu,he_cpu,he_a100}_22f7021d (method-1012); prodfinal binaries in
  /raven/ptmp/jinma/prodfinal_1011/bin with the prodfinal merge tag (same src tree).
