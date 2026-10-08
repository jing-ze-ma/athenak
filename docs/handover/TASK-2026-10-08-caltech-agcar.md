# TASK for Caltech: AG Car 3-D shake-downs A and B (general EOS), smoke then run (another copy; first start wins)

Same runs as TASK-2026-10-08-deltaai-agcar.md on this branch -- read that file for the physics, files and smoke table;
only the differences are listed here. Copies are queued on Raven (31004544 B, 31004545 A) and viper (12131308 B,
12131309 A); DeltaAI has the same task. Duplicates are intended: the user decides which copy continues. Tell the user
(NOTE) when a Caltech copy starts, so the pending copies elsewhere can be cancelled.

1. Code: fork/rt-integration (newest; it contains he_ic_eint_from_t, 2fbc6aa1's code) or this branch at 2fbc6aa1.
   PROBLEM=he_star_m1, your H200 stack as in NOTE-2026-10-08-caltech-hegiant.md (f3a66907 build dir: incremental).
2. Files: `docs/handover/agcar-files-1008/` on this branch -> /resnick/groups/carnegie_poc/jingze/agcar_1008/files,
   `bash SETUP.sh <that dir>`. NOTE (10-08): input A now has tfloor_kelvin = 3.0e3 and eos_logt_min = 3.4 (the 5000 K
   floor held A's outer atmosphere above its radiation temperature and FOFC fired there); the smoke references in
   smoke_ref_viper/ were regenerated for it (see the DeltaAI TASK table).
3. Layout: 4 MeshBlocks 480 x 64 x 64. Use 4 H200 (1 block each) if a 4-GPU node is available, else 1-2 H200 (several
   blocks per GPU; viper runs 2 per MI300A). Keep nx1 whole.
4. Smoke `time/nlim=10` for A and B against smoke_ref_viper (t, dt, mass, tot-E within ~1e-6; 0 NON-CONVERGED; the two
   IC T-check lines). Report s/cycle.
5. If both pass: B (tlim 4.624e5 s) and A (tlim 8.064e6 s) to the end with restarts, 12-24 h links, B first if only one
   fits. Report job ids, start estimates and s/cycle in NOTE-2026-10-08-caltech-agcar.md on this branch.
