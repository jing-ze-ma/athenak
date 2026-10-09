# TASK for DeltaAI: AG Car 3-D production (B, then A) -- URGENT (user GO 10-09)

AG Car (LBV) 3-D runs: he_star_m1 sp wedge 480 x 128 x 128, **4 MeshBlocks 480x64x64 -> 1 node x 4 GH200, 1 rank
per GPU**. B (hot state) is the production you run; A (cool state) is a third copy that also waits in the viper and
Raven queues (first start wins). Push NOTEs on this branch (`bsg-files-1009`) as described below. Account
`bivj-dtai-gh`. Work dir **`/work/nvme/bivj/jma20/agcar_1009/`** only (do not touch bsg_hrdet_1009/, delta_1008/
or older dirs).

## Code and binary
Branch `hrup-bsg-1009` @ **98835d99a** (same code as the viper/Raven AG Car binaries). **Reuse your existing
`athena_hes_gpu_98835d99a16f`** (md5 `c9c6d16705d07dc6e97de6c808d02182`, built with `build_inc_deltaai.sh hes_gpu`,
PROBLEM=he_star_m1, see NOTE-2026-10-09-deltaai-bsg-hrdet.md). Do not rebuild. Use the module stack it was built
with (module reset stack of TASK-2026-10-07-deltaai-hegiant sect. 2) inside every job.

## Files (bundle `docs/handover/agcar-prod-1009/`, 0.86 MB gz)
IC A/B + ext2 Rosseland + Orion-stitched Planck tables (the radial grid is a poly stretch in the input, the EOS
table is built at run time: no other files). Inputs are the Raven production inputs
`agcar_rcx{A,B}_ge_raven_accel_st_mg_hr.athinput` (md5 7cecd99c / d0bc5625) with only the 3 file paths replaced.
```
G=/work/nvme/bivj/jma20/agcar_1009; mkdir -p $G/bundle
git fetch fork bsg-files-1009
git archive fork/bsg-files-1009 docs/handover/agcar-prod-1009 | tar -x -C $G/bundle
P=$G/bundle/docs/handover/agcar-prod-1009
bash $P/SETUP.sh $G/files          # 4 x OK, SETUP_OK; writes $G/files/agcar_rcx{A,B}_ge_accel_st_mg_hr.athinput
```
Run dir per arm `$D` (e.g. `$G/runB`, `$G/smokeB`): copy `$P/rst_info.py` into `$D` (the link script reads
`$D/rst_info.py`, sbatch runs a spool copy of the script) and write `$D/run.cfg`:
```
BIN=<path>/athena_hes_gpu_98835d99a16f
BIN_MD5=c9c6d16705d07dc6e97de6c808d02182
IN=/work/nvme/bivj/jma20/agcar_1009/files/agcar_rcxB_ge_accel_st_mg_hr.athinput     # A: ..._rcxA_...
TLIM=1.3e6          # B (as its input); A: TLIM=8.064e6
XKEYS=""
```
Link script `$P/agcar_dai_link.sh` (copy of the viper/Raven chain logic): reads run.cfg at start; fresh start if
`$D/rst` is empty, else restarts from the newest rst; `$D/CANCEL`/`DONE` -> exit 0, `STOP` -> exit 1; md5 mismatch,
FATAL/NaN or rc != 0 -> writes STOP (no garbage restart chained); GUARD against an older running job of the same
name; srun `-n 4 -c 72 --cpu-bind=cores` + wrapper `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID`; code wall limit = job
limit - 15 min (- 4 min for jobs <= 1 h); nvidia-smi memory log `gpumem.<job>.log`; last line
`rc= fatal= nonconv_last= gpumem_max_MiB=`.
```
sbatch -A bivj-dtai-gh -p <partition> --time=<T> -J <name> -o $D/link.%j.out $P/agcar_dai_link.sh $D [NCYC]
```

## Steps
1. **IMMEDIATELY: NOTE `NOTE-2026-10-09-deltaai-agcar-queue.md`**: wall-time limits (`scontrol show partition
   ghx4-interactive ghx4`: MaxTime, MaxNodes, DefaultTime) and the current start estimate for a 1-node 4-GPU job on
   each (`sbatch --test-only -A bivj-dtai-gh -p <p> -N 1 --gpus-per-node=4 --exclusive --time=<MaxTime> ...`). Push.
2. **Smokes** (ghx4-interactive, `--time=00:20:00`): `smokeA` and `smokeB` run dirs, `NCYC=10` (2nd argument):
   rc 0, no FATAL / NaN, NON-CONVERGED count (the `<rad_m1> implicit transport ... NON-CONVERGED=` lines), steady
   s/cycle (median of cycles without dumps; cycle 0 dumps everything), max GPU memory per GH200 from gpumem.log.
   Report both in the queue NOTE (or a short smoke NOTE). If a smoke fails: STOP, NOTE, no production.
3. **Production B** (fresh start, t = 0) in `$G/runB`, name `agcB_dai`, on regular **ghx4** at its MaxTime, as a
   chain of 3 links (`--dependency=afterany:<prev>`); add links while the run is below TLIM (the script writes
   DONE at tlim). Then **A** in `$G/runA`, name `agcA_dai`, same chain pattern, submitted right after B (A may
   start in parallel on a 2nd node if the queue allows; A is the third copy, so its first link just waits).
   **The moment each production STARTS, push `NOTE-2026-10-09-deltaai-agcar-started.md`** (append: arm, job id,
   start time, node) so viper can cancel its pending copies.
4. **Before submitting every new link** (and on each check): `git fetch fork bsg-files-1009` and look for
   `docs/handover/NOTE-2026-10-09-viper-deltaai-agcar-CANCEL.md`. If it names an arm (A or B): `touch $D/CANCEL`,
   `scancel` that arm's PENDING links (`scancel --state=PENDING`), do not kill a running link unless the note says
   so, and acknowledge in a NOTE.
5. **After each link**: append to `NOTE-2026-10-09-deltaai-agcar-links.md`: arm, job id, rc line, t reached
   (rst_info.py on the newest rst), cycle, s/cycle (steady), last NON-CONVERGED count, max GPU memory, any STOP
   (with its text). Push.

Rules: no destructive git on the fork; never write into the bundle dir; keep all run dirs.
