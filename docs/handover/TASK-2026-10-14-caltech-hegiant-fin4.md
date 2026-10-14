# TASK raven -> Caltech: He giant on FIN4 (fin4-1014 @ d505d21f): smoke + restart-bitwise + 15-min background run NOW, then a HELD production chain

From the Raven coordinator (worker HG). USER 10-14: the He giant is now FIRST priority; Caltech runs the FIN4 validation (Raven queues
are too slow). Keep the BSG / Linhao work running as it is. Your FIN3 He giant run (fin3, link 1 4454171) stays stopped; keep its data.
Background: NOTE-2026-10-14-raven-caltech-hegiant-heating-reply.md (the base heating is MESA's own L(r_in); USER keeps the setup).

## 1. Build (report md5)
- Fork branch **`fin4-1014` @ `d505d21f`** (= FIN3 9972a56d + OT time2_enth_be / implicit_cs2_vimp_npass + OD mbband + DI dic-prod).
  Same toolchain and flags as your 9972a56d build (`-D PROBLEM=he_star_m1 ... -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`).
  Name `athena_gpu_he_star_m1_d505d21f_nofma`. (If you already built d505d21f for the FIN4 suite TASK c5787e33, reuse it.)
- Raven / Viper use the same hash (A100 cd049ef5, MI300A a94677b7).

## 2. Deck
- `docs/handover/hegiant-1014/hegiant_fin4.athinput.in` (md5 1b4d0592) -> `sed "s#@FILES@#<your He giant files dir>#g"`. Change no key.
  Same IC/tables as before (ic_giant_own.txt ddc72692, rosseland c1ef063c, planck b699c60d).
- = the FIN3 deck you ran (hegiant-1013/hegiant_fin3.athinput.in) + 3 named FIN4 keys: `time2_enth_be = true`, `implicit_cs2_vimp = true`,
  `implicit_cs2_vimp_npass = 3`. he_ic_equil = false, vet_gd_halo_ipc = false, lowmem not set (= false).
- Layout as your FIN3 run: 1 node, 2 ranks x 2 MeshBlocks (897x64x64), 1 rank per H200, your CUDA_VISIBLE_DEVICES wrapper, --cpu-bind=cores.

## 3. Validation (one job is fine; report each item)
A. Smoke C: fresh, `time/nlim=24 time/ndiag=1 output1/dt=1e-6 output2/dt=1e30 output3/dt=1e30 output4/dt=1e30 output6/dt=1e30 output5/dt=240`
   (rst 00001 at t = 240 s = cycle 12). Report rc, FATAL / NaN / DIVERGED / NON-CONVERGED counts, Picard mean/max, dt, s/cycle (cycles 3-24),
   GPU peak memory, and the startup line `time2_enth_be=true implicit_cs2_vimp_npass=3`.
B. Restart R: NEW srun `-r <C>/rst/hegiant.00001.rst time/nlim=24` into a fresh dir. The user.hst + hydro.hst rows after cycle 12 and the
   final bin dumps (every variable, data after the header) must be BITWISE vs C. Not bitwise -> STOP and report.
C. Background run BG (~15 min wall): fresh, deck as is plus `problem/he_seed=0.0 output2/dt=1000 output3/dt=1e30 output6/dt=1e30` on the
   command line, `-t 00:15:00` (it writes a final rst/dump). Then:
   `BASE=hegiant LOGF='run*.log' VISPY=<repo>/vis/python python3 docs/handover/hegiant-1014/bg.py 4.30368e12 1.25226e13 <BG dir> 1`
   Report:
   - the energy budget from hegiant.user.hst (wedge values; x 12 = full sphere): L_in + P_esrc vs L_MESA(r_in) = 6.74e6 Lsun full sphere
     (= 2.15e39 erg/s wedge); dEtot/dt; L_top;
   - the band table at the last dump, compared with the Raven FIN3 reference E0 at EQUAL model time
     (`hegiant-1014/bg_E0_ref_raven.txt`: same deck physics, no seed, FIN3 9972a56d-equivalent dad00933, Raven A100). E0 at t 1.0e4 s:
     0-0.5 R_ph +0.29 km/s / 0.002 | 0.99-1.01 +2.69 / 0.510 | 1.01-1.03 -0.11 / 0.008 | 1.03-1.10 -0.17 / 0.035 | 1.1-1.5 -0.15 / 0;
     the photosphere band expansion (+2.7 km/s at 1e4 s, +5.2 at 2.3e4 s in E0; MESA surface v_r +6.6 km/s);
   - fofc report line (last), Picard, NC, s/cycle.

## 4. Production chain: prepare it HELD, release ONLY on a GO file
- New run dir e.g. `/resnick/groups/carnegie_poc/jingze/hegiant_1014/fin4`; fresh t = 0 to tlim 2.592e6 s; your link script pattern
  (run.cfg read at job start; FRESH if no rst else the newest valid rst + rst_info.py last_time keys; STOP on rc != 0 / FATAL / NaN /
  DIVERGED / NON-CONVERGED > 0; DONE at rst t >= tlim; prune_rst.sh keep 3 EXIT trap as in your BSG run).
- Links 12 h, afterany; your FIN3 link ran 0.93 s/cycle at dt ~18-20 s -> ~1.4e5 cycles ~ 36 h + FIN4 cost (~+15 %) -> queue **4 links**.
- Submit the chain with the FIRST link HELD (`--hold`). Release it only when `docs/handover/GO-2026-10-14-caltech-hegiant-fin4.md` appears on
  bsg-files-1009 (the coordinator writes it after Raven validates your section 3 results). Your 15-min poll is enough.

## 5. Report: NOTE-2026-10-14-caltech-hegiant-fin4.md (push to bsg-files-1009)
- build md5; A/B/C results as listed (rc, BITWISE verdict, budget, bands vs E0, expansion, s/cycle, memory); held chain job ids.
