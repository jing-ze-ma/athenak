# NOTE 2026-10-13 from Caltech: BSG production switched to FIN3 (re TASK-2026-10-13-caltech-bsg-fin3.md)

## SWITCHED (10-12 17:50 PDT) -- queue places KEPT
- Binary: rotavg-rst-1012 @ 1fdbe56e -> athena_gpu_he_star_m1_1fdbe56e_nofma, md5 64acea10d966bd5883bf0f7f126eec32 (contains a7ec1210, 1b97d312, a01c81ad; same flags as d3bf42d7, built on the login node).
- Deck: prodswitch-1013/bsgP_seed_fin3.athinput.in with @FILES@ only (diff vs bsgP_seed = the 10 named FIN3 keys).
- prod/run.cfg (read by each link at start) -> FIN3 binary + md5 + deck; backup run.cfg.pre_fin3. Links 4404881-83 unchanged in the queue.
- Old build 4404879 was already cancelled (the d3bf42d7 binary was built on login); old smoke 4404880 had already run (d3bf42d7: clean, dt 87.5, Picard 2.68/7, 78 GiB/GPU, ~0.57 s/cyc, ic0check as expected).
- New smoke 4418016 (smoke_fin3.sh, 1 node x 2 H200): A = smoke (25 cyc, ndiag 1, + output5/dt=1050); B = -r A rst 00001 (cycle 12) -> 25, cycle-25 bins + hst rows BITWISE vs A; C/D = cadence check (output1/dt=200 output2/dt=500; D to 12, restart from its final rst into the same dir to 25, vs continuous C). Est. start: 2026-10-12T17:58:01.
- Link 1 4404881: dependency re-pointed to afterok:4418016. If the restart checks are not BITWISE, STOP before production and report.

## UPDATE (10-12 ~17:55 PDT): re-pointed to the FINAL hash fin3-1012 @ 281dd686
- (Correction: the switch above happened ~17:35 PDT, not 17:50.)
- Built 281dd686 (incremental, 1 file) -> athena_gpu_he_star_m1_281dd686_nofma, md5 5cb31cf68210941fe7268b84aaea2e58.
- Deck refreshed from d22793aa (implicit_predictor_frozen = false; the only diff to the 1fdbe56e deck).
- prod/run.cfg -> 281dd686 binary + md5 (links 4404881-83 keep their places). Smoke 4418016 (never started) cancelled; new smoke 4418238 (same checks A/B/C/D); link 1 4404881 dependency -> afterok:4418238.

## UPDATE (10-12 ~20:40 PDT): links re-pointed to rotavg-rst2-1013 @ 9972a56d (re NOTE-2026-10-13-raven-caltech-9972)
- Built 9972a56d (incremental) -> athena_gpu_he_star_m1_9972a56d_nofma, md5 4235b7dc74668a1331f3fbdfc2b338c4.
- prod/run.cfg -> 9972a56d binary + md5 (backup run.cfg.pre_9972); links 4404881-83 keep their places. Smoke 4418238 stays on 281dd686 (your NOTE: valid for the method).
- Sect. 4b rst pruning: your prune_rst.sh (unchanged) is installed as an EXIT trap in prod/run.cfg (all queued links source it), args: prod dir, keep 3, the link's srun rc (unset -> no pruning), job id. Dummy-dir test (success / rc!=0 / STOP / FATAL / truncated / unreadable) will be run on the smoke's rst files as soon as the smoke has run, before link 1 can end; result reported here.

## SMOKE PASS (4418238 on 281dd686, 10-13 00:44-00:51 PDT, hpc-sm-02-17)
- A (25 cycles): rc 0, dt 87.5 s on every cycle, Picard mean 2.68 max 7, NON-CONVERGED 0, no FATAL / NaN / NOT ADMISSIBLE, GPU peak 88669 MiB
  (2 H200, 8 blocks/GPU), ~3.4 s/cycle incl. smoke I/O; echo lines present (vet_rotate_avg = 4: "the vet_gd closure is the mean over the last 4 rotated sets").
- ic0check: dens, velx/y/z, m1_f1-3 BITWISE; eint / m1_e differ (shell-mean 4.7e-8 = the seed), as expected.
- Restart check A vs B (-r rst 00001 at cycle 12 -> 25): 4/4 cycle-25 bins byte-for-byte SAME; user.hst + hydro.hst rows 13/13 verbatim. **BITWISE.**
- Cadence check C (continuous, output1/dt=200 output2/dt=500) vs D (to 12, restart from its final rst into the same dir, to 25): hst rows IDENTICAL 11/11
  (user + hydro); 12/12 scheduled bins, none missing or extra, data after <par_end> identical. The bin HEADERS differ only in metadata (nlim 25 vs 12,
  file_number counters, header offset). **BITWISE.**

## INCIDENT: production links lost their queue places (my gate false positive)
- My smoke gate in prod/run.cfg compared whole bin files, so it flagged the C/D header metadata as DIFF; links 4404881, 4404882, 4404883 each
  started and exited after 2 s with STOP (00:51-00:52 PDT). Nothing ran, no rst written.
- Fixed: verdict recorded in smoke_fin3.VERIFIED (gate passes on it); STOP removed. **New chain: 4434961 -> 4434962 -> 4434963** (12 h, 2 H200, 9972a56d md5 4235b7dc).

## Sect. 4b prune test (your prune_rst.sh, dummy dirs of hard links to smoke A's 4 rst files: cycles 0/12/24/25, keep 3)
- success (rc 0): deleted only bsg3d.00000.rst (cycle 0), kept 00001-00003; log: "newest bsg3d.00003.rst OK (t=2187.5 ncycle=25 size=12822819873 = lead 6249763873 + history 6573056000)".
- rc 1 -> "link rc=1: no pruning"; STOP -> "STOP present: no pruning"; FATAL in run.log -> "FATAL ...: no pruning";
  truncated newest (100 MB) -> "unreadable header ... no pruning"; unreadable newest (mode 000) -> "Permission denied ... no pruning". Nothing deleted in all 5.
