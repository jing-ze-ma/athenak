# NOTE Caltech -> Raven coordinator: AG Car E1 production (TASK-2026-10-13-caltech-e1)

## QUEUED (10-13 10:52 PDT)
- Binary: `athena_gpu_he_star_m1_dad00933_nofma`, source dad00933acc5a5825137fb661deb402233936086 (dic-prod-1013), incremental from the
  9972a56d tree on the login node (CUDA 12.9, HOPPER90, -ffp-contract=off). **md5 ff7b7f11170340bee7d53dacec851ae9**.
- Files: e1caltech-1013/MD5SUMS all OK; IC e81deb51, rosseland b28e97c5, planck 493f9d01 (one dir `e1_1013/files`).
- Deck: template md5 62115d17 -> `@FILES@` sed only -> `e1_1013/agcA_E1_final_wbrmax.athinput` (read-only) **md5 d29e4615072074392ac492d351433afd**.
  wb_rmax = 2.9715e13 present; vet_gd_halo_ipc = true (kept; smoke Coff checks it).
- Launch (smoke AND production, set as LAUNCH in run.cfg): `srun --mpi=pmix -n 2 -c 8 --cpu-bind=cores e1wrap.sh` where e1wrap.sh = module env +
  `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID` (the BSG production launch). 1 node x 2 H200, 2 ranks x 2 MeshBlocks.
- **Smoke job 4452278** (smoke_e1.sh, 1 h): RUNNING since 10:51 PDT on hpc-sm-01-05.
- **Production links 4452279 (HELD, --hold) -> 4452280 -> 4452281** (afterany, 12 h each, name agcE1c, run dir `e1_1013/prod`,
  E1C_GATE=1, E1C_VERIFIED=e1_1013/smoke/VERIFIED). Released only when GO-e1-caltech.txt exists AND smoke_gate.sh wrote VERIFIED.
- **Start estimate** (`sbatch --test-only`, 1 node x 2 H200, 12 h) at 10:52 PDT: **14:01 PDT today** (~3 h). The smoke started at once.
- BSG fin3 chain 4434961-63 untouched.

## SMOKE PASS (10-13 11:40 PDT) -- links still HELD, waiting for GO-e1-caltech.txt
- First smoke 4452278 died at launch (rc 127, no physics ran): smoke_e1.sh writes `${LAUNCH:+LAUNCH="$LAUNCH"}` into run.cfg through an
  unquoted heredoc, so a multi-word LAUNCH loses its quotes ("--mpi=pmix: command not found"). Workaround: one-word LAUNCH =
  `e1_1013/e1launch.sh` (= `exec srun --mpi=pmix -n 2 -c 8 --cpu-bind=cores e1wrap.sh "$@"`), used for the smoke AND in prod/run.cfg
  (identical launch). Failed dirs kept as smoke_fail_quote. **Fix suggestion for the bundle: quote the value in the heredoc or document one-word LAUNCH.**
- Also: cmplast.py / bin_convert.py import h5py, which was missing here -> `pip install --user h5py` (3.14.0) on login; the in-job
  Coff cmplast line therefore failed, rerun on login below.
- **Smoke 4454178** (hpc-sm-02-09, 1 node x 2 H200, 2 ranks x 2 MeshBlocks), all through link_e1.sh:
  - **C** (fresh 60 cyc): rc 0; FATAL 0, NaN 0, DIVERGED 0, NON-CONVERGED 0, NOT ADMISSIBLE 0; fofc_report at-floor 0.
    `he_ic_equil: done after 1 rounds (deferred he_seed applied)`; round 1: cells 0..348 re-marched, 200 passes, max |drho/rho| 3.165e-02 at r 2.71976e+13.
    **dt 727.2727 s every cycle**. **Picard mean 2.389 max 6** (131 solves; Raven A100 2.42/6).
    `vgdrot-1011 vet_gd_halo_ipc: ON (node-local CUDA IPC, GPU-driven per-shell halo of the exact lists, no host wait)`.
    Init (elapsed at cycle 0) 39.6 s. **s/cycle (cycles 10-60) = 0.600** (44.23 -> 74.23 s).
    Last user.hst row: `4.0000000000000007e+04 7.2727272727272668e+02 2.8424334033570493e+38 9.3538331170191842e+37 2.8423121260304378e+38 ...`
  - **L** (30 + 30, link 2 restarted by the hook from L/rst/agcar3d.00001.rst, t 21818.18, cycle 30; last_time keys 20000/15000/15000/0/0):
    rc 0 both links, clean, Picard 2.08/3; restart trailer line `he_ic_equil: a_ref ... and Phi_eff from the fresh start's equilibrium`.
  - **Coff** (halo_ipc false): rc 0, clean, same Picard; s/cycle 0.618.
- **GATE (smoke_gate.sh): GATE PASS 2026-10-13_11:38:50** -> smoke/VERIFIED written by the gate.
  - CHAINCMP: hydro.hst C rows 11, chain 11, missing 0, extra 0, differing 0; user.hst same; bin: C dumps 8, chain dumps 8, extra 0 -> **CHAINCMP BITWISE**.
  - CMPLAST BITWISE; HST hydro/user BITWISE.
  - **Coff vs C: CMPLAST BITWISE, IPC-OFF HST hydro + user BITWISE** -> vet_gd_halo_ipc = true kept.
- GPU memory: not captured (the end-of-job nvidia-smi ran after the sruns); will read it from link 1's log.
- **Prune test** (hard links of smoke L/rst = 00000 c0, 00001 c30, 00002 c60, KEEP=2 so that success has one to delete):
  success -> deleted only agcar3d.00000.rst; rc=1 -> no pruning; STOP -> no pruning; FATAL in log -> no pruning;
  truncated newest (100 MB head of 00002 as 00003) -> sanity FAILED, no pruning; unreadable newest (4 KB random) -> sanity FAILED, no pruning. 6/6 correct.
- **Link estimate**: 0.60 s/cycle x 33,300 cycles = 5.6 h + outputs -> 1 link if dt holds. 3 links queued (4452279 HELD -> 80 -> 81).
- Next: on GO-e1-caltech.txt -> copy to prod/GO_e1, scontrol release 4452279, push STARTED.

## STARTED (10-13 11:48 PDT)
- GO-e1-caltech.txt (a7b887f5) copied to prod/GO_e1; smoke/VERIFIED present; `scontrol release 4452279` at ~11:45 PDT.
- **Link 1 4452279 RUNNING since 11:48:26 PDT on hpc-sm-02-17**, FRESH start t = 0 (`e1chain: FRESH start`, -t 11:49:55), halo_ipc ON.
  he_star_m1: he_ic_equil: done after 1 rounds (deferred he_seed applied)
- First check (11:53): cycle 650, t 4.727e5 s, dt 727.27, elapsed 305.8 s -> ~0.47 s/cycle so far. Links 4452280 -> 4452281 queued afterany.
- **Raven may cancel its own E1 copy now.**

## LINK 1 STOPPED: Picard DIVERGED at t = 1.24255e7 s, cycle 17085 (10-13 16:19 PDT) -- chain halted, NOT restarted
- Job 4452279 rc 134 (abort): `Picard solve DIVERGED: resid=1.167511e-01 lin_resid=2.701759e-08 > implicit_resid_fatal=1.000000e-02 at time 1.242545e+07 cycle 17085`
  (FATAL at rad_m1_implicit.cpp:13979, both ranks). Hook wrote STOP ("srun rc=134 FATAL DIVERGED"); links 4452280 / 4452281 exited at once on STOP; no pruning.
- **Onset** (no earlier warnings): first hesdirk2 `NOT ADMISSIBLE` at cycle 13622 (t ~9.9e6); per 1000 cycles: 13k 4, 15k 274, 16k 351, 17k 30 (659 total,
  redone with backward Euler). First NC-local line at cycle=15198; 93 `Picard NON-CONVERGED after 200 passes` lines before the divergence
  (resid 1e-7 .. 1e-1). Worst cells: rank 0 (m,k,j,i) (0,66,41,511), (0,66,49,478); rank 1 (0,18,45,347), (0,11,50,350) -> i ~ 350-512 of 536.
  min(E,T) in the NOT ADMISSIBLE lines ~2e-5, once -1 (cycle 17077).
- fofc_report at-floor 0 throughout. dt stayed 727.27 to the end (no dt shrink before the failure).
- Speed: 0.56-0.65 s/cycle until cycle 14000, then 1.3 s (15k) -> 2.2 s (16k) -> 3.4 s/cycle (16k-17k) as the BE redos and 200-pass Picards piled up.
- rst on disk (prune only runs after a clean link): agcar3d.00000 .. 00024 (every 5e5 s; 00024 = t 1.2e7, inside the troubled phase;
  00019 = 9.5e6 is the last before the first NOT ADMISSIBLE). All bins + hst kept. Full log: e1_1013/prod/run.log.4452279.
- **Waiting for the coordinator / DR / user**: do not restart until told which rst and what change. Same binary dad00933 / deck d29e4615.
