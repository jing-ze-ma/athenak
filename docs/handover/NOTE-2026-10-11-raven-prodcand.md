# NOTE prodcand (prodcand-1011): the production-candidate branch for rad_m1 (implicit VET)

Branch `prodcand-1011` (worktree /raven/u/jinma/ATHENAK/wt-prodcand). Work dir: /raven/ptmp/jinma/prodcand_1011. MPCDF Raven.

**Status.** The merges are done and validated. Every feature is bitwise against its own branch once the other branches' new defaults are opted out. The combined defaults pass every criterion.

**Before production, two re-checks are still needed:**
- **Predictor fix.** The literature suite (NOTE-2026-10-11-raven-suite.md) found that `time2_vet_sc=predict` (hesdirk2) gives a spurious flux in IDORT thermal-equilibrium case 2, and crashes the crossing beams at c dt/dx = 256. Both bugs were already in the per-branch binaries. The predfix-1011 agent is fixing the predictor on the vet_sc and vet_col paths.
- **GPU re-check.** The final candidate will include predfix, segfix, icgen and the ic_rad-scope change below. It needs one more short CPU and GPU re-check. Take the A100 timing in that run.

## 1. Merge log

Base: origin/vgdrot2-1011 **3e790a6c**, which is deepatm-1011 with:
- the thin switch fse/hcs/jse/bq;
- hesdirk2 predict;
- gdthin;
- `vet_bc_diff=face` as the default;
- vgdfuse;
- the vgdrot speed-up.

| commit | merge / change | conflicts and resolution |
|---|---|---|
| 19de78fc | icrad-1011 32546788: `ic_rad=fs_relax` default, `vet_source_noesrc` for every BE solve under hesdirk2, the G5 gate fix bb1ebff8 | Clean. All four noes sites (vet_col, vet_lat, vet_sc x2) already go through `VsrcNoes()`. vet_gd takes its sources from these. |
| c71ce4ce | cbfix-1011 8c05d344: `implicit_hr_face=sf` default. It brings littests-1011: implicit_beam_n, vet_ring, pgen options. | `rad_m1_vet.cpp` VetShortChar locals: kept both deepatm's `bfc` (face) and littests' boundary-beam locals. The beam intensity is added after the face-corrected diffusion intensity, as on each branch. The beams are single-MeshBlock vet_sc only, so the MB sweeps need no change. |
| 4c033dcf | moveq-1011 544e1f5b: `implicit_f0_init`, `implicit_rad_inertia`, `implicit_nofb_work` | Clean. All moveq hooks sit on the one write-back path shared with predict, vgdfuse and vgdrot: F0 init at it == 0, the nofb row, rad inertia before the write-back, and the wb kernel included in both `m1_impl_wb` variants. |
| 080dc504 | sptests-1011 1294e272: opt-in `vet_gd_inner=emit`, `vet_gd_beam`; sph_shell/sph_atm options | `rad_m1.hpp`: kept vgdrot's VgdTabHost/VgdTabDev members and added `VetGdBeamInit`. **Semantic fix:** vgdrot moved the direction vector `d` out of VetGdTables into VgdTabCompute. The first-call `VetGdBeamInit(d)` now rebuilds `d` from the host table `H.dir`. The values are the same, and `vet_gd_beam` requires rotate_every = 0. |
| 710b43b5 | torder-1010 9bc090c6 (docs) | Its code, 8ee1e217, is already in as the cherry-pick 10d36189 via thinsw. The `rad_m1_time2.cpp` conflict was resolved to ours. |
| b9289f7a, f8c9c2e4, 875d1228 | vetacc-1011 05f0f72f, vgdrot-1011 0d8eef6e, icrad-1011 31de33b9 | docs only |
| 51d68fbe | sptests-1011 9b8552ff: the `m1_gdhr` output, and M1_VGD_NMAX 432 -> 2048 | clean |
| e66ee669 | **M1_VGD_NMAX kept at 432** | NMAX sizes `Real buf[M1_VGD_NMAX]` (VetGdWall) per thread. 2048 would spill on A100. The 4x-GL angular study needs a local build with 2048. |
| 558be1c2 | **ic_rad default narrowed to he_star_m1 builds** | See §4. |
| 2e95be82 | docs: NOTE-2026-10-11-raven-suite.md (literature suite on 4c033dcf) | |

Not merged: icgen-1011 (in progress). Commits were made by Jingze Ma. Nothing was pushed by me; the coordinator pushed e66ee669.

Binaries are in bin/:
- **558be1c2**:
  - athena_cpu (built-in pgens), md5 0c911578;
  - athena_he_cpu, md5 8c0d7823, byte-identical to e66ee669's;
  - athena_he_a100, md5 23e40b85.
- **e66ee669**: athena_cpu 39ab667b, athena_he_cpu 8c0d7823, athena_he_a100 d7d1de92.
- **4c033dcf**: athena_cpu 1f797264, athena_he_cpu 503d072b, athena_he_a100 7352ae9e.

## 2. New defaults and their opt-outs

| default (new) | from | opt-out (old behaviour, bitwise) | acts on |
|---|---|---|---|
| `vet_bc_diff=face` | deepatm/gdface | `rad_m1/vet_bc_diff=cell` (a restart without the key keeps cell) | vet_sc, vet_gd, vet_col_lat, vet_col core |
| `ic_rad=fs_relax`, he_star_m1 builds only, fresh starts of implicit vet_col + vet_gd + fluid | icrad, narrowed in 558be1c2 | `rad_m1/ic_rad=none` | stellar-envelope runs |
| `vet_source_noesrc_h2be=true` | icrad | `rad_m1/vet_source_noesrc_h2be=false` | the BE solves under hesdirk2 (BESTORE start, fallbacks, ic_rad passes) |
| `implicit_hr_face=sf` | cbfix | `rad_m1/implicit_hr_face=cell` | `implicit_flux_beam=halfrange` with vet_sc |
| `implicit_f0_init=auto` | moveq | `rad_m1/implicit_f0_init=zero` | the cold-start face flux (moving gas) |
| `implicit_rad_inertia=true` | moveq | `rad_m1/implicit_rad_inertia=false` | the radiation momentum of the enthalpy flux (moving gas) |
| `implicit_nofb_work=true` | moveq | `rad_m1/implicit_nofb_work=false` | `gas_feedback=false` with moving gas |

With all seven opt-outs, the merged binary is bitwise identical to the base littests binary a61c1e75 and to cbfix 9e1f2686 with `hr_face=cell`. This was checked on eq3t cen be and nofb (lit/, runs `*_nall`, `*_cbc`, `*_lit`).

## 3. Production key set (AG Car A candidate)

Input: vetacc agcar_A_va.athinput. It holds the production blend, the be/hr keys and `implicit_blend_xthin=30`. Add:

```
rad_m1/time_scheme=hesdirk2 rad_m1/time2_vet_col=predict rad_m1/implicit_vimp=true rad_m1/time2_vstage=true
rad_m1/implicit_blend_xthin=0 rad_m1/implicit_thin_switch=fse rad_m1/implicit_thin_corr=hc
```

Everything else stays at the new defaults:
- face, fs_relax, h2be, hr_face sf;
- the moveq fixes;
- vgdfuse (halo_exact 2) and the vgdrot speed-up. vet_gd_halo_ipc stays off.

The run used fixed dt (`implicit_cfl=1251.3`, dt = 1239.996 s).

## 4. ic_rad default scope (558be1c2)

The spherical suite found that the `fs_relax` default fired on every fresh implicit vet_col + vet_gd + fluid start. That included a time-dependent pulse test, which lost 62 % of its energy before cycle 0.

The default is now compiled in only for PROBLEM=he_star_m1. CMake sets `M1_ICRAD_DEFAULT_FS_RELAX` on `rad_m1_icrad.cpp`. Every other pgen defaults to `none`, and an explicit `ic_rad` key wins everywhere.

Check (icchk/, job 31069504):

| run | 558be1c2 | 4c033dcf |
|---|---|---|
| he_star AG Car A reduced, fresh start, no ic_rad key | `ic_rad = fs_relax`, done 21.1 s, L/L_IC = 1.000624 | the same, identical L_top |
| built-in pulse (sptests pulse_r1000 deck, `pgen_name=rad_m1_beam`, vet_col + vet_gd) | no ic_rad at all | `ic_rad = fs_relax` (the bug) |
| the same pulse with `rad_m1/ic_rad=fs_relax` | `ic_rad = fs_relax` (explicit key wins) | - |

## 5. CPU validation (4c033dcf; e66ee669 / 558be1c2 checked bitwise against it)

Configurations:
- **ref**: the branch's own binary on its own keys.
- **new**: the merged binary on the same keys with the new defaults.
- **newx**: new plus the opt-outs of the *other* branches' defaults. This should reproduce ref bitwise.

Data arrays were compared with ana/cmpdir.py: bin data, tab and hst rows, as max |a-b|/max|a|.

**Determinism caveat.** CPU OpenMP runs with 4 threads can differ from 1 thread at round-off (1e-15). The cbfix shadow differed from the cbfix agent's own run by 1e-5 in F2 (near zero). So the bitwise pairs were run under identical settings in this work dir.

### 5.1 Crossing beams and static shadow (cbfix; ref = cbfix_1011/bin/athena_cpu_9e1f2686)

| run | new vs ref | newx vs ref | criterion (new) |
|---|---|---|---|
| cb hr / tJe / cen, vet_ring 4 | **bitwise** | bitwise | Row energy 1.006 at every height for hr and tJe (≥ 0.9 PASS); peak 0.417, floor 8.5e-11 / 4.2e-19. cen 0.743, as before. |
| shd hr / tJe, ring 12x2 | **bitwise** | bitwise | as cbfix: hr L1 0.012, umbra 0.9 %; tJe 0.011 / 0.14 % (PASS) |

New is bitwise even with face, because these tests have no diffusion bottom boundary. The same holds at OMP=1 for 558be1c2/e66ee669 vs 4c033dcf.

### 5.2 Moving equilibrium eq3t and ddn (moveq; ref = moveq_1011/bin/athena_cpu_mq3)

| run | new vs ref | newx vs ref | eq3t at t = 1 (new): E / v / P_tot |
|---|---|---|---|
| eq3t cen be, bev, h2, h2v; hr be; tJe be; tJe h2v; cen ovc | 1.1e-5 to 1.3e-5 (face) | **bitwise** | cen be 1.13333 / 2.95510 / 2.99984; cen h2 1.13334 / 2.95510; tJe h2v 1.13329 / 2.95512 (target 1.1326 / 2.9554: **PASS**); hr be 1.1261 / 2.9586 (known hr 1e-3 leak, as on moveq) |
| eq3t cen nofb (gas_feedback=false) | 5.6e-3 | 5.6e-3 | **see below** |
| ddn cen be, cen h2, tJe be, tJe h2 (gas_feedback=false) | 1.4e-5 / 5.4e-6 | the same | see below |

**nofb / ddn are not a merge effect.** The mq3 reference binary does not match the committed moveq-1011 source on the `gas_feedback=false` path. It was built at 04:15, one minute after commit 7a9af8bc, which added the nofb work row.

A fresh build of origin/moveq-1011 544e1f5b (bin/athena_cpu_moveq544e, from git archive):
- is **bitwise** to newx on eq3t nofb, ddn cen be and ddn tJe h2;
- is bitwise to mq3 on eq3t cen be and tJe h2v.

The merged code therefore reproduces moveq exactly. The new-vs-ref differences on the feedback runs (1e-5) come from `vet_bc_diff=face` alone, because newx is bitwise.

### 5.3 deepatm bundle: atm, 1-D ladders, G1 / G3, fixed G5 (ref = deepatm_1011/gdface/bin/athena_cpu_wipB = e84d2513)

94 runs per configuration:
- atm tJe and cen, 32 to 2048;
- pulse_k0.128 and rw_t10, tJe, be and hesdirk2, X/T/C;
- G1 1e2 / 1e4, G3 c1 / c10, the old G5;
- the fixed G5 for tJe_be, tJe_h2 and cen at 4 kappa.

Results:
- **newx: 94/94 bitwise.**
- **new: 62/94 bitwise.** All atm, pulse, G1 and G5 runs are bitwise. 30 rw_t10 runs differ by 4e-13 to 1.2e-10, and G3 by 1e-15. The cause is `implicit_rad_inertia` acting on the small gas velocities, as in NOTE-moveq §4.5; newx removes it.

**atm self-convergence by region, new (face), tJe** (`ana/atm_self_new.txt`; identical to deepatm's face table):

| band | orders 32→2048 |
|---|---|
| tau 0-0.1 | 2.04 2.02 2.09 2.01 1.99 |
| tau 0.1-1 | 1.92 2.03 2.01 2.00 2.00 |
| tau 1-10 | 2.00 1.99 2.00 2.00 2.00 |
| **tau > 10** | **1.99 2.00 2.00 2.00 2.00** |
| Kn < 0.03 | 1.99 2.00 2.00 2.01 2.00 |

cen tau > 10: 2.04 2.08 2.20 2.36 2.18.

**Ladders:**
- pulse_k0.128 eval tables are identical to ref.
- rw_t10 orders: tJe be X E 2.17 2.03 1.73 becomes 2.17 2.02 1.66 (finest pair, L1 3.7e-4 to 3.9e-4). hesdirk2 X/T/C are unchanged to two digits.

**Gates:**
- G1 Hopf L1 6.99e-4 (cfl 1e2 and 1e4).
- G3 Marshak L1(E) 0.0139 / 0.0142 (≤ 0.02, PASS).
- **Fixed G5: all 12 PASS.** The ratio ranges from 0.99957 to 1.00000 (within 2 %).
- The old in-battery G5 still fails at k 12.8 and 128, as on every branch. It is superseded by the fixed G5.

### 5.4 Reduced CPU AG Car A (icrad; ref = icrad_1011/bin/athena_he_cpu_32546788)

Setup:
- mesh 480x32x32, 4 ranks x 3 threads;
- `ic_rad=fs_relax`, fixed dt (cfl 1251.3), 12 cycles;
- run with hesdirk2 + predict + vimp + vstage, without and with the fse/hc thin keys.

Results:
- **newx (vet_bc_diff, hr_face and the moveq opt-outs) vs ref: bitwise** for both arms, on 14 files including the bins and both hst files. vgdfuse and vgdrot are therefore bitwise on CPU too.
- new vs ref: up to 3.8e-2 in the m1_fs formal-solution H at cycle 0. This comes from face and from hr_face sf, which acts with `flux_beam=halfrange`.
- L_top / L_top0 is the same to 5 digits: 1.00000 → 0.99978 (cycle 10) → 0.99993 (cycle 12) for both. The cand arm goes 0.99992 → 0.99990.
- Picard: 4.27 for new vs 4.30 for ref, with NON-CONVERGED 0.

### 5.5 Post-merge bitwise check of the defaults

e66ee669 vs 4c033dcf, both at OMP=1:
- cb hr / tJe / cen r4 and shd hr / tJe: bitwise;
- eq3t, all 9 arms: bitwise.

558be1c2 differs from e66ee669 only in the `ic_rad` default scope (§4). Its athena_he_cpu is byte-identical.

## 6. GPU: full AG Car A, 4 x A100 (Raven gpudev, job 31068627, binary he_a100_e66ee669)

Setup:
- fixed dt: `implicit_cfl=1251.3`, dt = 1239.996 s;
- 60 cycles, hst every cycle, dumps at cycles 0 and 60;
- runs: gpu/runs/j.31068627.

The ref arm failed at startup on the deck's `m1_thsw` output, which d28dacd8 does not have. Its rerun (31069219, with a deck without output6 and output7) was cancelled because the gpudev wait was about 2 h. **The ref numbers below are DeltaAI's**: NOTE-2026-10-11-deltaai-prodcand-agcarA.md, 4 GH200, same candidate commit.

| run | hw | rc / FATAL / NON-CONV | Picard mean / max | s/cycle steady / rotation / avg | cost (s per 1e4 s simulated, steady / avg) | ic_rad start |
|---|---|---|---|---|---|---|
| cand1 | A100 | 0 / 0 / 0 | 3.25 / 7 | **0.345 / 0.431 / 0.352** | 2.78 / 2.84 | 4.29 s (fs 0.44, relax 13 steps 3.17), L/L_IC 1.000612 |
| cand2 | A100 | 0 / 0 / 0 | 3.25 / 7 | 0.344 / 0.430 / 0.352 | 2.77 / 2.84 | the same |
| cand1 | GH200 | 0 / 0 / 0 | 3.25 / 7 (132 solves) | 0.226 / 0.223 / 0.248 | 1.82 / 2.00 | 2.98 s |
| ref (d28dacd8, be / hr xthin 30) | GH200 | 0 / 0 / 0 | 3.78 / 9 (60 solves) | 0.170 / 0.169 / 0.232 | 1.37 / 1.87 | - |

Notes on the table:
- The A100 average excludes the last cycle (4.2 s with the dumps) and starts at cycle 2. Init is 30.9 s.
- **cand2 is bitwise identical to cand1** on A100: 14 files (bins hydro_w, m1, m1_thsw, m1_fs at cycles 0 and 60, plus both hst files). The same holds on GH200.
- **Cost.** cand costs +33 % per steady cycle against ref on GH200 (two hesdirk2 stages, at lower Picard) and +7 % on the average.
- **A100 rotation cycles** cost about 0.43 s against 0.345 s for steady cycles. The A100 cand steady cycle (0.345 s) is faster than vgdrot's 0.422 s for the same keys on the deepatm base.

**L_top per cycle (A100 cand1, L/L_top(0); cand2 identical).** L_top is flat from cycle 0 and falls by at most 2.8e-4:

| cycle | 0 | 5 | 10 | 15 | 20 | 30 | 40 | 50 | 60 |
|---|---|---|---|---|---|---|---|---|---|
| L_top / L_top(0) | 1.00000 | 0.99983 | 0.99977 | 0.99974 | 0.99973 | 0.99972 | 0.99972 | 0.99972 | 0.99973 |

The ref arm (GH200, default IC) has a spike: ×1.536 at cycle 1, 1.0010 at cycle 10, 1.00035 at cycle 60. At cycle 60, cand / ref = 0.99999.

**E/J in the envelope at cycle 60, cand1** (A100; GH200 identical to 4 digits). Values are median (p5-p95):

| | 1.0-1.05 R_ph | 1.05-1.5 | 1.5-2 | 2-3 |
|---|---|---|---|---|
| E/J_fs | 1.0024 (0.9862-1.0179) | **1.0008** (0.9715-1.0300) | 0.9991 | 1.0064 |
| E/J_tw | 1.0023 (1.0018-1.0028) | **1.0005** (0.9999-1.0015) | 1.0011 | 1.0100 |
| ref E/J_fs (GH200) | 1.0197 | **1.0187** | 1.0078 | 1.0101 |

**cand vs ref, shell means at cycle 60 (GH200):** max over the shells in each band of |Δ<q>| / |<q>_ref|.

| quantity | < 0.8 R_ph | photosphere 0.95-1.05 | envelope > 1.05 |
|---|---|---|---|
| dens | 1.4e-5 | 2.95e-2 | 1.54e-2 |
| E | 1.2e-5 | 5.90e-2 | 1.94e-2 |
| F1 | 6.4e-3 (at the base) | 3.6e-5 | 1.1e-5 |
| velx (max\|Δ\| / max\|ref\| in band) | 2.4e-2 | 4.9e-2 | 3.2e-2 |

The interior is unchanged to 1e-5. The photosphere carries the change: ref is still relaxing from its +54 % step-1 spike, and the envelope E/J goes from 1.019 to 1.001.

## 7. Files

Scripts:
- build.sh: the three builds, tag T.
- build_mq.sh: the moveq source rebuild check.
- cpu/: copy of the deepatm bundle with an `EXTRA_KEYS` hook in one.sh, battery_raven.sh and run_g5.sh, plus cpu_val.sh.
- lit/: run1.sh, wrap.sh (configs ref, new, newx, e66, nall, cbc, lit, mqo, mqb) and the pool*.sbatch job files.
- agc/agc_cpu.sh, gpu/gpu_pc.sh, gpu/gpu_pc2.sh, icchk/icchk.sh.

Analysis:
- ana/cmpdir.py (data compare), ana/ltop.py;
- ana/dv_cmp_{new,newx}.txt, ana/atm_self_new.txt, ana/lad_*.txt.

Runs:
- runs/lit/{cbfix,moveq}/…, runs/dv_{ref,new,newx}/, runs/agc/j.31068235/;
- gpu/runs/j.31068627/, icchk/j.31069504/.

DeltaAI TASK: bsg-files-1009 docs/handover/TASK-2026-10-11-deltaai-prodcand-agcarA.md (pushed by the coordinator, 58ee4a14).
