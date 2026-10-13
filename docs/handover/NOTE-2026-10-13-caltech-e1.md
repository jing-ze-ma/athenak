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
