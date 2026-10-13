# TASK raven -> Caltech: switch the queued BSG production to FIN3 (rotavg-rst2-1013 @ 9972a56d), keep the queue places

From the Raven coordinator (worker PS, prodswitch-1013). USER decision: "we should already switch all production binaries to this
one including other machines". It applies to your queued copy of TASK-2026-10-12-caltech-bsg-requeue.md
(NOTE-2026-10-12-caltech-bsg-requeue.md: build 4404879 -> smoke 4404880 -> links 4404881-83). Nothing has started yet, so this is a
switch BEFORE the start, not a mid-run change. Raven (bsgP, bsgP2n) and Viper (bsgPv) switch the same way.
DeltaAI's running d3bf42d7 copy is NOT affected (it continues unchanged as the old-method comparison run).

## 1. Build (report md5)
- Fork branch **`rotavg-rst2-1013` @ `9972a56d`**. This is the FINAL production hash: Raven and Viper production use the same commit.
  It is fin3-1012 a01c81ad (FIN3) plus five commits:
  - RS's restart commits:
    - a7ec1210: the vet_rotate_avg history and conservative2's cell-F slope in the rst;
    - 1b97d312: the he_star Picard window in the rst;
    - 1fdbe56e: link-end outputs keep the output schedule;
  - 281dd686: the implicit_predictor_frozen fresh default is false (the deck names it anyway);
  - 9972a56d: the averaging history is stored only over its valid radial range (about -55 % rst size).
- Build exactly 9972a56d: the user requires bitwise restarts, and no older fin3 commit has them.
- Same toolchain and flags as your d3bf42d7 build: `-D PROBLEM=he_star_m1 -D Athena_ENABLE_MPI=ON -D Kokkos_ENABLE_CUDA=On
  -D Kokkos_ARCH_HOPPER90=On -D CMAKE_CXX_FLAGS=-ffp-contract=off`, default `Athena_M1_VGD_NMAX`.
- Name `athena_gpu_he_star_m1_9972a56d_nofma`.

## 2. Deck
- `docs/handover/prodswitch-1013/bsgP_seed_fin3.athinput.in` -> `sed "s#@FILES@#<your 10-09 BSG files dir>#g"`.
  The IC files are unchanged (ic_ma2026_f0 246a666e, rosseland 450bc0c1, planck e18ca805).
- It is your current deck `docs/handover/bsgprod-1012/bsgP_seed.athinput.in` (53995cf9) with ONLY these rad_m1 method keys added,
  every one NAMED, so a restart never falls back to an old default:
  `vet_gd_twin_state = trans`, `vet_rotate_avg = 4`, `vet_gd_twin_weight = none`, `vet_interp = auto`,
  `vet_col_lat_cut_thin = 1.0`, `implicit_predictor_frozen = false` (updated 2026-10-13: the frozen predictor was dropped from FIN3), `vet_bc_inner = auto`, `vet_bc_reflect = auto`,
  `implicit_hr_face_gd = none`, `vet_source = relaxed` (`vet_source_noesrc = true` was already named).
- Physics, IC, seed, layout keys (halo_ipc true, lowmem true, 8 blocks/GPU), outputs and tlim 4.96e6 are unchanged. Change no key.

## 3. Switch the queued jobs, keeping their queue places
1. If your link script reads the binary and deck from a config file at job start (as Raven's run.cfg does), edit that file
   (keep a `*.pre_fin3` backup): binary -> the new `_nofma` binary + md5, deck -> the FIN3 deck. The links keep their places.
2. The build job 4404879 and the smoke 4404880 have the old commit / binary inside their scripts:
   - cancel those two only;
   - submit the new build (section 1) and a new smoke (section 4) with the same afterok dependency;
   - re-point the first link: `scontrol update jobid=4404881 dependency=afterok:<new smoke id>`.
   4404881-83 keep their queue priority.
3. If a link script itself names the old binary or deck (not a config file read at start), only then:
   - cancel and resubmit the links with the FIN3 binary and deck;
   - state in the NOTE that the places were lost.

## 4. Smoke (25 cycles), then production
Smoke, as before (`time/nlim=25 time/ndiag=1`, 1 node x 2 H200, 2 ranks x 8 blocks). Check:
- rc 0, no FATAL / NaN / NON-CONVERGED, and dt 87.5 s on every cycle after the first.
- Picard mean (expected ~2.7-3), GPU peak memory, and s/cycle vs your d3bf42d7 numbers if you have any.
  vet_rotate_avg 4 adds the ring of the last 3 closure sets; Raven measures how much on A100.
- The run.log echo lines of `vet_rotate_avg` ("the vet_gd closure is the mean over the last 4 rotated sets"), `vet_gd_twin_state`,
  `vet_interp` and `vet_col_lat_cut_thin`.
- `python3 ic0check.py check <smoke dir> ic0_ref_1009.json`: dens, vel1-3 and m1_f1-3 must be BITWISE. eint and m1_e may differ
  only in 37.6-45.1 Rsun (the seed).
- Restart check (new; required before production).
  - Run the smoke deck continuous to cycle 25 with an extra periodic rst: `output5/dt=1050` on the command line, so rst 00001 is
    written at t = 1050 s = cycle 12.
  - Then, in a NEW srun (a second job is fine), `-r rst/bsg3d.00001.rst time/nlim=25` into a fresh dir.
  - Compare the cycle-25 bin dumps (every variable, byte-for-byte) and the user.hst rows of cycles 13-25 with the continuous run:
    they must be BITWISE.
  - Raven A100 got BITWISE (K=4 and K=0, BSG and AG Car) only with RS's restart fix. d3bf42d7 and fin3 without the fix differ at
    1e-8, because conservative2's cell-F slope was not stored.
  - Raven A100 and Viper MI300A on 1fdbe56e / 281dd686 / 9972a56d: BITWISE.
  - Also check that a chained link keeps the output cadence. Run the smoke to cycle 12 with `output1/dt=200 output2/dt=500`,
    restart from its FINAL rst to cycle 25, and compare with a continuous run with the same keys. The hst rows and the scheduled
    dumps must match one-to-one, with none missing and all BITWISE. Raven and Viper on 1fdbe56e: BITWISE, 11/11 rows.
  - If yours is not bitwise, STOP and report before any production link.

Production: unchanged otherwise. Fresh t = 0 to tlim 4.96e6, chained links (restart from the newest rst with the rst_info.py
last_time keys), STOP on rc != 0 / FATAL / NaN / NON-CONVERGED.

## 4b. Restart-file pruning (user, 2026-10-13): keep only the newest 3 rst files
A full-BSG rst file is about 13 GB with vet_rotate_avg 4, and the run writes about 57 of them (one per day). After each link, your link script
deletes the rst files older than the newest 3 in the run's rst/ dir, under these rules:
- Prune only after a link that ENDED SUCCESSFULLY:
  - srun rc 0;
  - no STOP;
  - no FATAL in the link's run log, and no NaN in its last 60 lines.
- The newest rst must pass a sanity check first. Read its header as rst_info.py does, then require:
  - size > 0;
  - the LEADING part within 1 % of the next-older rst's, unless that one is the cycle-0 file. The leading part is the file size
    minus the vet_rotate_avg history block. That block's size legitimately changes between files: with K-1 vs K slots it is
    about 17 % of a file, and with the trimmed radial range it changes too.
  - How to find the history block: the 8-byte marker `M1ROTAV1` in the header region, then a uint64 nb, then nb/8 int64 h[].
    The block holds h[1] slabs of nmb_total x (nx1+2ng)(nx2+2ng)(nx3+2ng) doubles.
  - The history header must be self-consistent: h[0] = 1, h[1] > 0, 2 <= K = h[2], h[7] in {K-1, K}, and the history smaller
    than the file;
  - time finite, > 0, and > the next-older rst's time; ncycle > 0;
  - the most recently written file must be that newest one;
  - for the rst output block: 0 <= time - last_time < dt.
- Never delete after a failed or STOPped link. Never delete an unreadable file.
- Log every deletion (file, size, t, ncycle) in the run dir, e.g. rst_prune.log.
- Reference implementation: `docs/handover/prodswitch-1013/prune_rst.sh RUN_DIR 3 RC JOBID` (bash + python3, flock), called from an EXIT trap of the
  link script after the srun. Port it, or write the same rule into your link script.
- Test it on a dummy run dir first (hard-linked copies of a smoke run's rst files):
  - success: deletes the oldest file;
  - rc != 0, STOP, FATAL, a truncated newest file, or an unreadable newest file: deletes nothing.
- Report the test in the NOTE.

## 5. Report: NOTE-2026-10-13-caltech-bsg-fin3.md (push to bsg-files-1009)
- SWITCHED: what was edited or resubmitted, the job ids, whether the queue places were kept.
- STARTED: binary md5, smoke numbers, the ic0check output, the restart check verdict.
- Then per link: t, cycle, dt, s/cycle, Picard, counts, newest rst.
