# TASK for DeltaAI: GPU cost breakdown of the production radiation scheme (prof, user 10-09)

Moved here from Delta (TASK-2026-10-08-delta-prof-batch1 on rad-beam-1008, withdrawn): DeltaAI = short tests.
Question: where does the GPU time of the production radiation scheme go (hydro vs rad_m1 implicit solve vs vet_gd
sweep vs half-range vs vet_col/lat vs halo/MPI vs outputs)? Profiling only, no production, no bitwise gate. Nobody
else runs it.
**Budget: ONE job, 1 node x 4 GH200, ghx4-interactive, `--time=00:30:00` (hard cap 1 h), account `bivj-dtai-gh`,
<= 4 GPU-h charged.** Expected compute: a few minutes (AG Car B 0.29 s/cycle, BSG reduced 0.24 s/cycle on your
smokes / mem rounds). Work dir **`/work/nvme/bivj/jma20/prof_1009/`** (new); read agcar_1009/ and bsg_hrdet_1009/
in place, write nothing there.

**Order:** after the AG Car cancel (NOTE-2026-10-09-viper-deltaai-agcar-CANCEL.md). Interactive allows 1 job per
user: submit this once job 3349057 (or any running AG Car link) has ended.

## 1. Binary (no build)
`athena_hes_gpu_98835d99a16f` (hrup-bsg-1009 98835d99a), md5 `c9c6d16705d07dc6e97de6c808d02182`, the one in
`/work/nvme/bivj/jma20/agcar_1009/smokeB/run.cfg` (prof.sbatch sources BIN from there and checks the md5). Module
stack = the one the AG Car smokes ran with (edit the comment block in prof.sbatch if your jobs load modules).

## 2. Files + tools
```
P=/work/nvme/bivj/jma20/prof_1009; mkdir -p $P/bundle $P/runs
git fetch fork bsg-files-1009 && git archive fork/bsg-files-1009 docs/handover/deltaai-prof-1009 | tar -x -C $P/bundle
```
- AG Car B input: `/work/nvme/bivj/jma20/agcar_1009/files/agcar_rcxB_ge_accel_st_mg_hr.athinput` (SETUP.sh of
  agcar-prod-1009, already done).
- BSG input: `/work/nvme/bivj/jma20/bsg_hrdet_1009/files/bsg_hr_dc5.athinput` (SETUP.sh of bsg-hrdet-1009, done).
- Kokkos Tools simple-kernel-timer, in your existing kokkos-tools master checkout (the memory-events build):
  ```
  cd /work/nvme/bivj/jma20/bsg_hrdet_1009/kt/profiling/simple-kernel-timer && make CXX=g++   # kp_kernel_timer.so + kp_reader
  ```
  If the Makefile is gone: `cmake -B build -DCMAKE_CXX_COMPILER=g++ && make -C build` at the kt root, then set the
  .so path in prof.sbatch and `export KPR=<kp_reader>`. If the checkout is elsewhere, edit the last path in
  prof.sbatch only.

## 3. Job: `sbatch $P/bundle/docs/handover/deltaai-prof-1009/prof.sbatch`
Binding = your standard: `srun -n 4 -c 72 --cpu-bind=cores` + wrapper `CUDA_VISIBLE_DEVICES=$SLURM_LOCALID`
(written by prof_cuda.sh into the run dir). 8 arms in sequence, 4 ranks, 1 MeshBlock per GH200:

| arm | case | mode | nlim |
|---|---|---|---|
| agc_plain | AG Car B production input, FRESH t = 0, 480x128x128 = 4 x (480x64x64), production outputs kept | no tools | 60 |
| agc_kt10 / agc_kt60 | same | simple-kernel-timer | 10 / 60 |
| agc_tmr | same, input copy + `rad_m1/implicit_timers = 3` (code's own fenced category timers) | no tools | 60 |
| bsg_plain | BSG `bsg_hr_dc5.athinput`, reduced 4x1 mesh of the mem rounds (256x128x128 = 4 x (256x64x64)), bin/rst dumps off | no tools | 30 |
| bsg_kt5 / bsg_kt30 | same | simple-kernel-timer | 5 / 30 |
| bsg_tmr | same + implicit_timers | no tools | 30 |

Keys: the command-line keys in prof_cuda.sh (time/nlim, time/ndiag, the reduced-mesh mesh/meshblock keys, output
dcycle/dt) all exist in the inputs; the only NEW key, `implicit_timers`, goes into an input COPY that the script
writes (sed after `<rad_m1>`), never on the command line (your mem_cuda2 lesson). Smoke rule: the first arm
agc_plain is the smoke; the script stops the job if it has rc != 0 or FATAL. A CUDA OOM is a result: report, do not
retry.

## 4. Analysis (login node)
```
R=$P/runs/<jobid>; D=$P/bundle/docs/handover/deltaai-prof-1009; L=$D/labels_98835d99.tsv
python3 $D/prof_group.py $L $R/agc_kt60 $R/agc_kt10 60 10 > $R/agc_prof.txt
python3 $D/prof_group.py $L $R/bsg_kt30 $R/bsg_kt5 30 5  > $R/bsg_prof.txt
```
(steady cycles = long minus short run, mean over the 4 rank files; labels_98835d99.tsv maps each par_for label to its
source file at 98835d99a; groups are the RULES at the top of prof_group.py.) If the parser fails on the kp_reader
format, fix the parser (not the code) and say so.

## 5. NOTE back: `docs/handover/NOTE-2026-10-09-deltaai-prof.md` on THIS branch (bsg-files-1009)
Per case (AG Car B, BSG reduced):
- s/cycle from the plain arm (`s/cycle: mean(6-N)= median=` line; the median excludes dump cycles), plus the kt and
  tmr arms' s/cycle (tool overhead), peak GPU memory, rc/FATAL/NaN/NON-CONVERGED counts.
- the code's counters at the end of the plain arm: `Picard iterations mean= max= NON-CONVERGED=` line, and the
  `<rad_m1> timers` lines of the tmr arm (closure, opacity, solve_pre, tensor, pred, pass_setup, krylov, pass_post,
  solve_end, hyd_c2p, m1_bvals; stages, solves, passes, kry_it).
- from `*_prof.txt`: the top 25 kernels table (time, % of kernel time, % of wall), the groups table (hydro
  fluxes/update/c2p; rad matrix/operator assembly; Krylov/BiCGStab vector ops; mg/line preconditioner; Picard
  update/residual; vet_gd sweep + twin + halo; half-range pass; vet_col/lat; halo/MPI pack-unpack; outputs;
  outside-kernel time = host + MPI waits + launch gaps), the radiation vs hydro fraction line, the top unmapped labels.
- job id, node, queue wait and elapsed (relative), GPU-h (4 x elapsed, x2 interactive charge), binary md5.
Attach the raw files (small) under `docs/handover/deltaai-prof-1009/results/<jobid>/`: both `*_prof.txt`, the kt_*.txt
of the four kt arms, and each arm's out.log tail (`tail -60`). No bins. Keep the run dirs.

Rules: no destructive git on the fork; never write into the bundle dir or other run dirs.
