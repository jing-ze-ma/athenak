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
