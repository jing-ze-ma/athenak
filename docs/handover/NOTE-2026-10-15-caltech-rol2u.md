# NOTE caltech -> raven: ROL2U (fin51-r 1bd0c2ce) gates on the newest He giant + AG Car E1 rst (TASK-2026-10-15-caltech-rol2u)

Work dir `/resnick/groups/carnegie_poc/jingze/rol2u_1015` (bundle at `caltech_rol2u/`, JOBS.txt). Production dirs were only read;
the seed rst were copied once.
USER 10-15 (via coordinator): section 4 (window B/RW and the optional H arm) is SKIPPED, so there is no window table and no HSTCMP.
The chains were submitted held and then released by the coordinator (user decision) without waiting for the gates. The He giant
chain runs on 2 H200 (user change).

## 1. Build
- fin51-r @ `1bd0c2cefbd178ef94bc7de5b4e691f8e90ea974`, incremental in builds_fin51 (43 objects recompiled; fin51 recipe, nofma,
  gcc 13.2 + cuda 12.9 + hpcx 2.17.1, Kokkos 4.6.2).
- `/resnick/home/jingze/ATHENAK/builds/bin/athena_gpu_he_star_m1_1bd0c2ce_nofma`, md5 **dd5f68261ad6a2d3d8050cfefe3cb2db**.

- Header fix (user 10-15, coordinator build): branch **fin51-r-hdr @ 5ab8c87e1e3f953752268d0ff89bf0636bc1fe70** = 1bd0c2ce +
  src/parameter_input.cpp only (`<par_end>` search limit 10 x 4096 -> 256 x 4096 bytes, ~1 MB; no effect on results).
  Binary `builds/bin/athena_gpu_he_star_m1_5ab8c87e_nofma`, md5 **0677cd3be292391bc64225823bbf6bb4**. hgR2w run.cfg now uses
  it (BIN/BIN_MD5), so the first hgR2w link runs 5ab8c87e. The gates below ran 1bd0c2ce.

## 2. Seeds (copies)
| case | source rst | t | cycle | md5 | `<par_end>` |
|---|---|---|---|---|---|
| hg | hegiant_1014/wide4/run/rst/hegiant.00060.rst | 2.592e6 (30 d) | 147411 | 9557d36afc0cd734ba0696ec511c7924 | 41554 |
| e1 | e1s_1014/prod/rst/agcar3d.00052.rst | 2.6e7 | 38823 | a99fc21c632965769f0acf22ee2f6fff | 39049 |

## 3. Gates (debug QoS, hpc-sm-02-10)
s/cycle = `cpu time used` / cycles, so it includes the rst read and setup; the L arms run only half the cycles and therefore look
slower. GPU peak = nvidia-smi peak per GPU.

### He giant, job 4552182 (4 H200, 30 cycles 147411 -> 147441): **GATE PASS**
| arm | rc | wall s | s/cyc | Picard mean/max | NON-CONV | lagged | e_floor clips | DT_RETRY | GPU peak GiB |
|---|---|---|---|---|---|---|---|---|---|
| P d505d21f | 0 | 79 | 1.03 | 7.23 / 9 | 0 | - | 0 | 0 | 43.1 |
| B base | 0 | 76 | 1.01 | 7.23 / 9 | 0 | - | 0 | 0 | 29.7 |
| C ROL2U | 0 | 66 | 0.80 | 3.06 / 4 | 0 | 62 | 0 | 0 | 29.7 |
| L1 | 0 | 60 | 1.24 | 3.10 / 4 | 0 | 31 | 0 | 0 | 29.0 |
| L2 | 0 | 59 | 1.27 | 3.03 / 4 | 0 | 31 | 0 | 0 | 29.7 |
| F forced | 0 | 67 | 0.87 | 3.06 / 4 | 0 | 62 | 0 | 1 (forced) | 29.7 |

- GATE-BASE P=B: rst BITWISE, hydro.hst + user.hst BITWISE (30 rows), bins CMPLAST BITWISE.
- GATE-RESTART C=L2: rst BITWISE, both hst BITWISE (15 rows), bins BITWISE.
- GATE-RETRY C=F: rst BITWISE, both hst BITWISE (30 rows), bins BITWISE; exactly 1 `DT_RETRY ... reason=forced`.
- FATAL 0 and NaN 0 in all arms. `<par_end>` C 42512, L2 42512 (< 44000). The absent dt_retry key was added on restart as expected.

### AG Car E1, job 4552183 (2 H200, 20 cycles 38823 -> 38843): **GATE FAIL** (GATE-BASE and GATE-RETRY)
| arm | rc | wall s | s/cyc | Picard mean/max | NON-CONV | lagged | e_floor clips (solved / eint<=0) | DT_RETRY | GPU peak GiB |
|---|---|---|---|---|---|---|---|---|---|
| P 5b7e1bc5 | 0 | 110 | 3.68 | 39.4 / 134 | 0 | - | 5 / 17 | 2 | 58.2 |
| B base | 0 | 109 | 3.66 | 39.4 / 134 | 0 | - | 5 / 17 | 2 | 35.1 |
| C ROL2U | 0 | 58 | 1.15 | 5.27 / 8 | 0 | 41 | 22 / 4 | 0 | 34.6 |
| L1 | 0 | 51 | 1.53 | 4.90 / 7 | 0 | 21 | 18 / 3 | 0 | 34.8 |
| L2 | 0 | 52 | 1.72 | 5.65 / 8 | 0 | 20 | 4 / 1 | 0 | 34.6 |
| F forced | 0 | 61 | 1.21 | 5.27 / 8 | 0 | 41 | 22 / 4 | 1 (forced) | 34.6 |

- GATE-BASE P=B: **DIFF** in the rst, both hst and the bins (CMPLAST 9 differ).
- GATE-RESTART C=L2: BITWISE (rst, both hst over 10 rows, bins).
- GATE-RETRY C=F: **DIFF** in the rst, both hst and the bins (8 differ); exactly 1 forced DT_RETRY line.
- rc 0, FATAL 0, NaN 0 in all six arms. `<par_end>` C 39580, L2 39382.

**Both failures start right after a DT_RETRY at cycle 38828 (= seed + 5).**
- P and B take the same natural retries:
  - cycle 38828: dt 340.426 -> 172.264, resid 0.0211119; P prints `reason=DIVERGED`, B `reason=resid>implicit_resid_fatal`.
  - cycle 38841: NON-CONVERGED, dt 331.063 -> 167.863.
- First differing hst row is row 6 (t = 2.6001874391181756e7, the first row after the 38828 retry), at round-off level:
  - P=B hydro.hst col 6: 8.7875502689635318e33 vs 8.7875502689641221e33.
  - P=B user.hst cols 4-5: 1.2642173950256634e38 vs ...638e38, and 1.2632912460411525e38 vs ...629e38.
  - C=F user.hst first diff row 6 (col 6: 1.2592999301793033e38 vs ...035e38); hydro.hst first diff row 9.
- Max rel diff in the final bins:

| field | P vs B (t 2.60044e7) | C vs F (t 2.60111e7) |
|---|---|---|
| dens | 1.2e-7 | BITWISE |
| velx | 1.2e-7 | 1.1e-7 |
| vely | 1.4e-7 | 1.9e-7 |
| velz | 3.5e-7 | 2.7e-7 |
| eint | 4.2e-7 | 1.2e-7 |
| m1_e | 7.6e-7 | 1.2e-7 |
| m1_f1 | 5.0e-6 | 1.2e-7 |
| m1_f2 | 1.8e-3 | 1.4e-3 |
| m1_f3 | 7.4e-4 | 7.7e-4 |

- Reading: 1bd0c2ce with default keys follows the production binary exactly until a DT_RETRY. After a retry (natural in B, forced
  fac 1.0 in F) the state differs at round-off, which chaos then grows. The He giant (no natural retry; forced retry BITWISE) does
  not show it. The suspect is the retry restore path changed in 1bd0c2ce (rad_m1_dtretry_visit.hpp, or the opacity / gas-newton
  state that is not saved/restored for the redone step). It needs a Raven check.
- Note: the RST format and the restart itself are fine (C=L2 BITWISE). The E1 ROL2U arm C itself ran cleanly: Picard 39.4 -> 5.3,
  0 retries (BASE 2), 3.7 -> 1.15 s/cycle.

## 5. Chains
- He giant **hgR2w** `rol2u_1015/hg/chain`:
  - Setup: 1 node x 2 H200, 2 ranks x 2 MeshBlocks; link_hgR2w.sh = link_hg4w_ext.sh with -n 2, gres 2, mem 240G.
  - run.cfg: BIN 5ab8c87e (hdr fix; was 1bd0c2ce at submission), XKEYS tlim 1e9 + dt_retry + ROL2U keys, prune keep 3 on this dir.
  - Seed: hard link of OUR seed copy. Submitted with `-t 7-00:00:00 --time-min=12:00:00`.
  - Jobs: 4552184 -> 4552185 -> 4552186 (afterany), released (coordinator, user decision), pending at the time of this note.
  - rmax guard: rmax_guard_r2.sh, RLIM 150, jobs 4552218 + afterany 4552229. Own guardwork/ cache, merged page run =
    He giant page run up to the seed + this chain.
  - 2-GPU fit: 4-GPU gate peak 29.7 GiB/GPU; FIN5 2-H200 smoke 84 GB/GPU; expected ~60-85 GiB on 141 GB.
- AG Car E1 **agcE1r** `rol2u_1015/e1/chain`:
  - Jobs 4552187-89 were submitted and released by the coordinator; **removed from the queue after GATE FAIL, before any started.**
  - Held until Raven/user decide.
  - run.cfg: BIN 1bd0c2ce, deck 680dadb0, KEYS = ROL2U, launch2.sh, E1C_MARGIN 1200; currently E1C_GATE=0, set by the coordinator.

## 6. Analysis
Window analysis: skipped together with section 4 (user). Long-term side-by-side tables follow in a later NOTE once hgR2w has
output.
