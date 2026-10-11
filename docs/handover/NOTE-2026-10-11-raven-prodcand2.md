# NOTE prodcand2 (prodcand2-1011): the final production candidate for rad_m1 (implicit VET)

Branch `prodcand2-1011`, worktree /raven/u/jinma/ATHENAK/wt-prodcand2, MPCDF Raven. Work dir: /raven/ptmp/jinma/prodcand2_1011.
Code commit **767d3313**. b89c137e adds only the sptests 929805c7 docs, and this NOTE commit is docs only.

**Status.** Every merged branch's headline test reproduces on the new binary with all new keys at their defaults. Where a result is not
bitwise, the difference is explained, and in every such case the explanation was checked with the cause switched off (§4).
- The two literature failures of NOTE-suite are fixed:
  - eq2 now passes with REC;
  - the crossing beams r4 / r12 at c dt/dx 256 now complete with REC.
- On the reduced CPU AG Car A, the candidate is **bitwise** identical to prodcand 558be1c2.
- On 4 A100 the candidate (cand1) is **bitwise** identical to prodcand e66ee669 cand1, so the production run is unchanged.
- **Not merged: consrc-1011.** The user's instruction says it is still being worked on by another agent. The user will decide on it.

## 1. Merge log

Base: origin/prodcand-1011 **686f6a24**. Its NOTE is NOTE-2026-10-11-raven-prodcand.md.

| commit | merge | conflicts and resolution |
|---|---|---|
| b8837405 | predfix-1011 f2d7c75c (fix eb75fa96): `time2_pred_cap` 0.5 default, `vet_x1_periodic_inflow = source` default, no fatal on an already-inadmissible hesdirk2 stage | clean |
| 55431efc | segfix-1011 0c061fff: `vet_seg = mean \| split \| emit` in vet_sc, vet_gd and vet_col_lat; rad_m1_vetseg.hpp; cyl pgen 24 objects | `rad_m1_vet.cpp` VetShortChar locals: kept predfix's `x1src`/`dbgc`/`dcl_` and segfix's `spl`/`spe` side by side. The segment code merged cleanly: predfix touches only the boundary half-cell entry intensity, segfix only the interior segment. Default stays **mean**. |
| 9851f232 | icgen-1011 709921a3: opt-in `problem/he_ic_rad = sph1d` (+ `_gas`, `_hse`, `_dump`), he_star_m1.cpp only | clean |
| 767d3313 | sptests-1011 9c800be8: opt-in `implicit_hr_excess = iso` (vet_gd, gated by `implicit_hr_excess_cfl` 10) + NOTE stages C/X/3a | clean; `M1_VGD_NMAX` stays **432** |
| b89c137e | sptests-1011 929805c7 (NOTE stage 3b, docs only) | clean |

Commits are by Jingze Ma; the coordinator pushed b89c137e.

**Binaries** (bin/, built from 767d3313):

| file | md5 | build |
|---|---|---|
| athena_cpu_767d3313 | 196a4b3d | CPU, built-in pgens (includes cyl24) |
| athena_he_cpu_767d3313 | 616ef614 | CPU, PROBLEM=he_star_m1 |
| athena_he_a100_767d3313 | 3bdad89d | A100, PROBLEM=he_star_m1 |

## 2. Defaults: new in prodcand2, and the ones carried over

| default | from | opt-out (old behaviour, bitwise) | acts on |
|---|---|---|---|
| **`time2_pred_cap = 0.5`** (fresh runs; a restart without the key gets 0) | predfix | `rad_m1/time2_pred_cap=0` | hesdirk2 predict (vet_sc and vet_col / vet_col_lat / vet_gd) where one step changes E or eint by more than half. The limited builds are counted in the run.log line `time2_pred_cap=... limited cell-builds=`. |
| **`vet_x1_periodic_inflow = source`** (fresh runs; a restart keeps `boundary`) | predfix | `rad_m1/vet_x1_periodic_inflow=boundary` | `vet_x1_periodic` (single-MeshBlock tests only: eq*, rw, pulse ladders) |
| **no fatal on an already-failed hesdirk2 stage** | predfix | none (the BE redo now runs instead of exit) | `implicit_resid_fatal(_masked)` when the stage is already NOT ADMISSIBLE |
| `vet_seg = mean` (unchanged behaviour) | segfix | - | - |
| `vet_bc_diff=face` | deepatm (prodcand) | `rad_m1/vet_bc_diff=cell` | vet_sc, vet_gd, vet_col_lat, vet_col core |
| `ic_rad=fs_relax`, he_star_m1 builds, fresh implicit vet_col + vet_gd + fluid starts | icrad (prodcand 558be1c2) | `rad_m1/ic_rad=none` | stellar envelopes |
| `vet_source_noesrc_h2be=true` | icrad | `rad_m1/vet_source_noesrc_h2be=false` | hesdirk2 BE solves |
| `implicit_hr_face=sf` | cbfix | `rad_m1/implicit_hr_face=cell` | halfrange + vet_sc |
| `implicit_f0_init=auto`, `implicit_rad_inertia=true`, `implicit_nofb_work=true` | moveq | `=zero`, `=false`, `=false` | moving gas |

**Command-line keys.** AthenaK refuses a command-line key that the deck does not name: "Parameter ... on command line not found".
- To use `vet_seg`, `implicit_hr_excess`, `ic_rad` or `problem/he_ic_rad*`, name them in the deck at their defaults. The run tools' fixkeys.py does that.
- inp/agcar_A_va_x.athinput is agcar_A_va with `vet_seg = mean` and `implicit_hr_excess = none` added. On CPU it is bitwise identical to agcar_A_va (new_candx vs new_cand).

## 3. Production key set (AG Car A)

Deck: vetacc agcar_A_va.athinput, with fixed dt `implicit_cfl=1251.3` (dt = 1239.996 s). Add:

```
rad_m1/time_scheme=hesdirk2 rad_m1/time2_vet_col=predict rad_m1/implicit_vimp=true rad_m1/time2_vstage=true
rad_m1/implicit_blend_xthin=0 rad_m1/implicit_thin_switch=fse rad_m1/implicit_thin_corr=hc
```

Everything else stays at the defaults of §2:
- face, fs_relax, h2be, hr_face sf, moveq;
- pred_cap 0.5;
- vet_seg mean;
- vgdfuse and vgdrot.

**Opt-in options, with recommendations** (not defaults; the user has not approved switching any of them):

| option | keys | recommendation | evidence |
|---|---|---|---|
| segment at opaque -> thin faces | `rad_m1/vet_seg=emit` | **recommended** for new vet_sc / vet_gd / vet_col_lat runs (segfix). Never `split` in production (it blackens heated walls). | Hot emitters fixed (slab J_sc 0.05 -> 0.4996 at tau_cell 10; P2 beam 2.3x -> 1.5x); heated walls unchanged. AG Car A: only a slight change, Picard +2 % (§6) |
| transport-consistent IC | `problem/he_ic_rad=sph1d problem/he_ic_rad_gas=re_atm problem/he_ic_rad_hse=true` (+ `ic_rad=fs_relax`, the default) | **recommended** for fresh AG Car starts (icgen) | Static photosphere: max\|v\| at cycle 60 is 0.29 km/s against 2.67; rho drift 2.6 % against 40 %. Picard 2.79 against 3.25. Solve 12.2 s on the host (§6) |
| P2 consistency limiter | `rad_m1/implicit_hr_excess=iso` | **optional**, for vet_gd thermal-channel / beam problems (sptests P2: L_M1/L_ex 2.3 -> 0.31-0.69). Harmless on AG Car A. | AG Car A changes are ≤ 2e-4; there is one +3e-4 L_top blip at cycle 2 (§6). The 120-direction pinc run has one Picard-198 step (sptests). |

## 4. CPU re-check (compute nodes, `-p general`)

Data arrays were compared with ana/cmpdir.py: bin, tab and hst rows, as max|a-b|/max|a|. Bitwise pairs were run with the same OMP thread count.

### 4.1 Literature suite, REC arm (lit/run1b.sh = suite run1b.sh with BIN = athena_cpu_767d3313; pool 18 x 4 threads)

The reference is predfix eb75fa96 (`*_v2` runs, same thread count). It equals the suite's 4c033dcf except where predfix changed the result.

| test | vs predfix eb75fa96 | result on 767d3313 |
|---|---|---|
| **eq2** REC (and REC lag) | **bitwise** | **PASS**: E 145.7878, uniform; \|v\| ≤ 3e-12 (suite REC: 144.347, v 1.706) |
| eq2c (vet_col, hesdirk2 predict / BE) | E bitwise uniform 145.787811; only F/v round-off noise (\|v\| ≤ 1e-12) | PASS (cap limits as in predfix) |
| **crossing beams** cb REC r4 / r12 at c dt/dx 256; r4 c0.5, r4 c32, cbg g8x32 | **bitwise** | **r4 and r12 complete**. Final outputs = suite / RECbe: row energy 1.006, peak 0.417 |
| static shadow shd r12g2 / r36g2, shdg | bitwise | PASS: L1 0.011, umbra 0.14 % |
| eq1 | bitwise | 96.79505 / 3.13663 PASS |
| eq3t | bitwise | E 1.13407, spatial sd 4e-12 (uniform), PASS |
| non-LTE nltek1-8 | bitwise | L1 1.1 % ... 0.01 %, E_top 0.9994-1.0000 PASS |
| sdifk, ddiff, ddn, ddf | bitwise | as suite / predfix (sdifk L1 ≤ 0.28 %, ddf peak 0.970) PASS |
| Marshak c1 / c10 (`implicit_vimp=false`, as suite) | bitwise | L1(E) 0.0139 / 0.0141 PASS |
| radiation shocks M 1.2 / 2 / 3, n512, and M2 / M3 n1024 | **bitwise vs suite 4c033dcf** at the suite's 9 threads (4-thread runs: 1e-9 to 4e-5, OpenMP round-off) | L1 ≤ 0.0014 (n512), ≤ 0.0017 (n1024), PASS. pred_cap: 0 limited builds |

Linear waves (rw ladder, REC = tJe hesdirk2 predict vstage vimp, N 16-256; ana/rw_767d.txt). Last orders:

| P, tau | 767d3313 | suite 4c033dcf |
|---|---|---|
| 0.01, 0.01 / 1 / 100 | 1.99 / 2.05 / 2.12 | 1.99 / 2.03 / 2.12 |
| 10, 0.01 / 1 / 100 | 2.17 / **2.06** / 2.02 | 2.16 / **1.43** / 2.02 |

- The change is predfix's `vet_x1_periodic_inflow = source`: rw is x1-periodic.
- P = 10, tau = 1: L1 at N = 256 drops from 2.63e-4 to 1.51e-4, and the order goes from 1.43 to 2.06. omega is unchanged to 4e-5.

**Fixed G5** (tJe be / h2, kappa 12.8 / 128 / 1280 / 12800): be 0.999581 / 0.999972 / 1 / 1; h2 0.999565 / 0.999874 / 1 / 1. **All PASS**, identical to predfix and suite.

**Grey atmosphere by region, face** (deepatm bundle, tJe and cen, 32-2048): all 14 runs are **bitwise** with prodcand dv_new. Self-convergence (ana/atm_self_767d.txt), tJe:

| band | orders 32 -> 2048 |
|---|---|
| tau 0-0.1 | 2.04 2.02 2.09 2.01 1.99 |
| tau 0.1-1 | 1.92 2.03 2.01 2.00 2.00 |
| tau 1-10 | 2.00 1.99 2.00 2.00 2.00 |
| tau > 10 | 1.99 2.00 2.00 2.00 2.00 |
| Kn < 0.03 / 0.03-0.3 / > 0.3 | ~2.00 in all |

cen tau > 10: 2.04 2.08 2.20 2.36 2.18, as before.

### 4.2 vet_seg = emit (segfix inputs; vs segfix wip2 a6a97649 runs)

| run | vs segfix | note |
|---|---|---|
| slab RECbe emit, tau_cell 1 / 100 / 1e4 | **bitwise** | J_sc = 0.5 on the thin side, as segfix |
| hcyl RECbe c1000 emit | **bitwise** | |
| slab RECbe tau 100, default (no key) vs segfix mean | E, F1 bitwise; F2 (~0) 3e-4 round-off | default = mean |
| **crooked pipe** REC c220 emit, t = 10 (16 threads) | not bitwise: T_end probes within 6e-4. 0.25 keV arrival 0.1 0.1 0.3 0.5 3.9 (segfix 0.1 0.1 0.4 0.7 4.0); T(t=1) probe 5 0.216 keV (segfix 0.211) | **pred cap**: 59 830 of 1.3e8 builds limited at the fronts. **With `time2_pred_cap=0` it is bitwise with segfix emit.** |
| wall1d REC f1 emit | not bitwise: E ≤ 0.4 % early, 0.14 % later | **pred cap** (388 of 84 672 limited builds). With cap 0 (9 threads), E and dens are bitwise and F1 differs by 1e-7. |

Merged emit therefore behaves exactly as in segfix. On REC runs with steep fronts, the predfix cap adds a sub-percent change on top: the arrival in the pipe is a little earlier.

### 4.3 sptests (spt.py, SPT_BIN = athena_cpu_767d3313)

- pinc FfX 64 (`implicit_hr_excess=iso`) vs xs2 runs_xs2: **bitwise**, so L_M1/L_ex 0.31 / 0.32 stands.
- pinc Ff 64 vs runs_pc: **bitwise**.
- pulse_r1 FfX / Ff 128 vs runs_xs2: **bitwise**, so the gate leaves the pulse unchanged.

### 4.4 Reduced CPU AG Car A (480x32x32, 4 ranks x 3 threads, fixed dt 1251.3, 60 cycles; agc/agc_cpu.sh, agc2.sh)

| arm | vs reference | Picard mean / max, NON-CONV | L_top/L0: cycle 1 / 10 / 60 |
|---|---|---|---|
| **new_cand** (production keys, defaults) | **bitwise** vs old_cand (558be1c2): 14 files incl. hst | 3.35 / 7, 0 | 0.99996 / 0.99983 / 0.99980 |
| new_candx (agcar_A_va_x deck) | bitwise vs new_cand | same | same |
| new_emit | 1.2e-5 (phot.) / 1.3e-4 (env.) shell-mean E; F1 ≤ 1e-5 | 3.44 / 8, 0 | 1.00002 / 0.99984 / 0.99981 |
| new_iso | ≤ 1.3e-5 E, ≤ 4e-7 dens | 3.35 / 8, 0 | identical to 5 digits |
| new_ic (sph1d re_atm hse + fs_relax) | - | 2.95 / 6, 0 | L/L_IC after fs_relax 0.99962 (7 steps); flat to 3e-4 |

- pred_cap limited 0 of 2.9e7 cell-builds, and there were no stage fallbacks.
- The ic_rad start of new_cand took 61.5 s on CPU: 14 relax steps, L/L_IC 1.00054.

icgen headline (deck agcar_A_icg, its own keys, 30 cycles), compared with NOTE-icgen 4.3:
- L_top/L0 at cycle 1 / 2 / 3 / 5 / 10: default 1.52832 / 1.44194 / 1.24905 / 1.04075 / 1.00138, against NOTE 1.5283 / 1.4419 / 1.2490 / 1.0407 / 1.0014.
- sph1d re+hse: 1.00577 / 1.00222 / 1.00071 / 1.00004 / 0.99981, against NOTE 1.0058 / 1.0022 / 1.0007 / 1.0000 / 0.9997.
- Picard 4.23 / 9 and 4.53 / 7, identical to the NOTE.

## 5. Literature-suite status (REC) after prodcand2

| test | status |
|---|---|
| eq1, eq2, eq3t (moving), nLTE, sdifk, ddiff / ddn / ddf, Marshak, shadow, crossing beams (all c dt/dx incl. 256), linear waves, radiation shocks, fixed G5, atmosphere | **PASS** |
| crossing beams in vacuum | stage 1 is still inadmissible on 16/20 steps (offdiag positivity fallback), so REC runs as BE there. Harmless; not fixed (predfix). |
| multi-source bubbles | FAIL (moment front, NOTE-suite 3.4; emit improves the total energy, segfix) |

## 6. GPU: full AG Car A, 4 x A100 (gpudev jobs 31071238 and 31071688, binary athena_he_a100_767d3313)

Setup:
- fixed dt: `implicit_cfl=1251.3`, dt = 1239.996 s;
- 60 cycles, hst every cycle, dumps at cycles 0 and 60;
- runs: gpu/runs/j.31071238 (cand1, cand2, ic) and j.31071688 (emit, iso, with the deck agcar_A_va_x).

The emit and iso arms of 31071238 stopped at the parser (keys not named in agcar_A_va) and were rerun in 31071688.

| arm | rc / FATAL / NON-CONV | Picard mean / max | s/cycle steady / rotation / avg(2-59) | cost s per 1e4 s simulated (steady / avg) | init (ic) |
|---|---|---|---|---|---|
| (1) cand1 | 0 / 0 / 0 | 3.25 / 7 | **0.347 / 0.433 / 0.358** | **2.80 / 2.89** | 31.4 s; fs_relax 4.32 s (13 steps), L/L_IC 1.000612 |
| (2) cand2 | 0 / 0 / 0 | 3.25 / 7 | 0.347 / 0.434 / 0.355 | 2.80 / 2.87 | the same |
| (3) emit | 0 / 0 / 0 | 3.32 / 7 | 0.353 / 0.445 / 0.362 | 2.85 / 2.92 | 31.1 s |
| (4) ic | 0 / 0 / 0 | **2.79 / 6** | 0.330 / 0.436 / 0.340 | 2.66 / 2.74 | 42.3 s: sph1d 12.2 s (331 VET passes, 29 outer), fs_relax 3.08 s (8 steps), L/L_IC 0.99967 |
| (5) iso | 0 / 0 / 0 | 3.30 / 7 | 0.347 / 0.433 / 0.356 | 2.80 / 2.87 | 32.0 s |

Notes on the table:
- **cand2 is bitwise identical to cand1** (14 files: bins hydro_w, m1, m1_thsw, m1_fs at cycles 0 and 60, plus the hst files).
- **cand1 is bitwise identical to prodcand e66ee669 cand1** (j.31068627). The production run is therefore unchanged by predfix / segfix / icgen / sptests at the defaults.
- The A100 timing for the final candidate is 0.347 s per steady cycle (prodcand: 0.345).
- The last cycle (with the dumps) is 4.1-6.1 s and is excluded from the averages.

**L_top per cycle (/ L_top(cand1, cycle 0)):**

| cycle | 0 | 1 | 2 | 3 | 5 | 10 | 20 | 30 | 40 | 50 | 60 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| cand1 (= cand2) | 1.00000 | 0.99993 | 0.99990 | 0.99988 | 0.99983 | 0.99977 | 0.99973 | 0.99972 | 0.99972 | 0.99972 | 0.99973 |
| emit | 1.00000 | 0.99994 | 0.99990 | 0.99988 | 0.99984 | 0.99977 | 0.99973 | 0.99972 | 0.99972 | 0.99972 | 0.99973 |
| iso | 1.00000 | 0.99994 | **1.00022** | 0.99981 | 0.99983 | 0.99977 | 0.99973 | 0.99972 | 0.99972 | 0.99972 | 0.99973 |
| ic | 0.99906 | 0.99908 | 0.99910 | 0.99911 | 0.99913 | 0.99917 | 0.99923 | 0.99927 | 0.99931 | 0.99933 | 0.99935 |

ic starts 0.09 % lower: this is its own transport luminosity. It then rises by 3e-4 over 60 cycles.

**Envelope E/J at cycle 60**, median (p5-p95):

| arm | 1.0-1.05 R_ph E/J_fs | 1.05-1.5 E/J_fs | 1.05-1.5 E/J_tw | 1.5-2 E/J_tw | 2-3 E/J_tw |
|---|---|---|---|---|---|
| cand1 | 1.0024 (0.9862-1.0179) | **1.0008** (0.9715-1.0300) | 1.0005 (0.9999-1.0015) | 1.0011 | 1.0100 |
| emit | 1.0019 (0.9854-1.0172) | **1.0001** (0.9708-1.0292) | 0.9998 (0.9992-1.0007) | 1.0004 | 1.0093 |
| ic | 1.0026 (0.9863-1.0178) | 1.0008 (0.9715-1.0299) | 1.0005 | 1.0011 | 1.0100 |
| iso | 1.0024 (0.9862-1.0179) | 1.0008 (0.9715-1.0299) | 1.0005 | 1.0011 | 1.0100 |

**Differences vs cand1 at cycle 60.** Values are the shell-mean relative difference, max over the shells in the band (ana_va.py vsame "Sh"):

| arm | quantity | < 0.8 R_ph | 0.8-0.95 | photosphere 0.95-1.05 | envelope > 1.05 |
|---|---|---|---|---|---|
| emit | dens / E / F1 / velx | 4e-11 / 1e-11 / 1e-8 / 8e-8 | 2e-9 / 2e-8 / 4e-7 / 3.6e-4 | 2.2e-5 / 1.2e-5 / 7.9e-6 / 2.8e-3 | 1.9e-6 / **1.4e-4** / 8.0e-6 / 5.9e-5 |
| ic | dens / E / F1 / velx | 8e-8 / 7e-8 / 5e-5 / 5.5e-3 | 1e-5 / 1.1e-5 / 3.2e-4 / 0.32 | **7.5e-3 / 2.9e-3** / 4.1e-4 / 0.45 | 8e-5 / 1.9e-3 / 3.9e-4 / 3.4e-4 |
| iso | dens / E / F1 / velx | 4e-12 / 8e-12 / 3e-9 / 1e-8 | 2e-10 / 7e-10 / 5e-9 / 1.2e-5 | 3.7e-7 / 2.2e-7 / 1.3e-7 / 5.1e-5 | 1.2e-5 / 3.4e-5 / 1.4e-7 / 2.1e-4 |

What the differences mean:
- **emit.** The change sits in the thin envelope: E -1.4e-4, and E/J 0.9998-1.0001 against 1.0005-1.0008, so a little closer to 1. Everything else is ≤ 3e-5.
- **ic.** It changes the photosphere, as icgen intends.
  - Photosphere (0.97-1.028 R_ph): max|v| at cycle 60 is **0.29 km/s against 2.67 km/s** (cand1), and max|rho60/rho0 - 1| is **2.6 % against 40 %**.
  - The interior is unchanged to 1e-4, and Picard is 14 % lower.
- **iso.** No change, apart from the cycle-2 L_top blip.

## 7. Files (/raven/ptmp/jinma/prodcand2_1011)

Scripts:
- build.sh (tag T);
- lit/run1b.sh (copy of suite run1b.sh: BIN, INDIR, predfix arms.sh), lit/pool.sbatch, lit/jobs_*.txt, lit/cap0.sbatch;
- cpu/j2.sbatch (atm, G5, rw), cpu/j4.sbatch (emit pipe / slab / hcyl / wall1d, sptests);
- agc/agc_cpu.sh, agc/agc2.sh;
- gpu/gpu_pc.sh, gpu/gpu_pc_b.sh;
- inp/agcar_A_va_x.athinput.

Analysis:
- ana/cmpdir.py, ana/cmpv.py, ana/ltop.py, ana/photv.py, ana/eqstat.py;
- ana/atm_self_767d.txt, ana/rw_767d.txt.

Runs:
- lit/runs/<test>/REC_*_c2 (and _c2b, _o9, _cap0);
- runs/dv_767d (atm, g5), rw/runs_767d, sp/runs;
- runs/agc/j.31071127, runs/agc/j.31072285;
- gpu/runs/j.31071238, gpu/runs/j.31071688 (gpu/runs/cmp links all five arms).

DeltaAI: TASK-2026-10-11-deltaai-prodcand2-agcarA.md (bsg-files-1009, posted by the coordinator) asks for the same five arms on GH200.
