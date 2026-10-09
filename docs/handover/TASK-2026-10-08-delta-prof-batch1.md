# WITHDRAWN (viper, user 10-09): do NOT run this batch. Moved to DeltaAI (fork/bsg-files-1009 547b106f TASK-2026-10-09-deltaai-prof.md). If already running: let it finish and push the results NOTE as described below; if not started: do not start it.

# TASK for the Delta runner: prof batch1 = GPU cost breakdown of the production radiation scheme (user 10-09)

User 10-09: "can we put something on delta". Question: where does the GPU time of the production radiation scheme
go (hydro vs rad_m1 implicit solve vs vet_gd sweep vs half-range vs vet_col/lat vs halo/MPI vs outputs)? Profiling
only, no production, no bitwise gate. Nobody else runs this batch (viper/Raven do not), so no race.
**Budget: 1 node x 4 A100, one job, `--time=00:30:00`, <= 4 GPU-h (interactive bills 2x); ledger line "prof"
(charge to beam).**

## 1. Binary (no new build expected)
Code = fork `hrup-bsg-1009` @ **98835d99a** (AG Car / He giant / BSG production code). You built it for batch3 /
mem batch1: `bin/athena_he_gpu_98835d99` (md5 b3f01de1...). Reuse it; only if missing:
`git fetch origin hrup-bsg-1009 && bash build_delta.sh he_gpu 98835d99a` (incremental, same build dir). Record md5.

## 2. Scripts + files
```
P=/work/nvme/bivj/jma20/delta_1008/prof; mkdir -p $P/bundle $P/runs
git fetch origin rad-beam-1008 && git archive FETCH_HEAD docs/handover/delta-prof-1009 | tar -x -C $P/bundle
git fetch origin bsg-files-1009 && git archive FETCH_HEAD docs/handover/agcar-prod-1009 | tar -x -C $P/bundle
bash $P/bundle/docs/handover/agcar-prod-1009/SETUP.sh $P/agcar_files     # 4 x OK, SETUP_OK
```
SETUP.sh writes `$P/agcar_files/agcar_rcxB_ge_accel_st_mg_hr.athinput` with @FILES@ = `$P/agcar_files`
(ic_agcar_B_ge.txt d39c9afa, rosseland_ext2_gs98_x0.36_z0.02.txt b28e97c5,
planck_ext2orion_gs98_x0.36_z0.02.txt 493f9d01). BSG: reuse `/work/nvme/bivj/jma20/delta_1008/bsg_hrdet/files`
(`bsg_hr_dc5.athinput`, ic_ma2026_f0.txt 246a666e, Rosseland 450bc0c1, Planck e18ca805); if absent run step 2 of
TASK-2026-10-09-delta-beam-batch3.md (bundle `docs/handover/bsg-hrdet-1009/SETUP.sh` on bsg-files-1009).

Kokkos Tools simple-kernel-timer (host-only build, like the memory-events build of mem batch1):
```
cd $P && git clone --depth 1 https://github.com/kokkos/kokkos-tools kt
cd kt/profiling/simple-kernel-timer && make CXX=g++     # -> kp_kernel_timer.so + kp_reader
```
(if the Makefile is gone on master: `cmake -B build -DCMAKE_CXX_COMPILER=g++ && make -C build` at the kt root, then
point prof.sbatch / `KPR=` at the built `libkp_kernel_timer.so` / `kp_reader`.)

## 3. Job: `sbatch $P/bundle/docs/handover/delta-prof-1009/prof.sbatch`
(read the header of `prof_cuda.sh`; edit the paths in prof.sbatch only if yours differ). 8 arms in sequence, 4 ranks,
1 MeshBlock per A100 (DeltaAI smoke B: 19.2 GB per GH200, so it fits 40 GB):

| arm | case | mode | nlim |
|---|---|---|---|
| agc_plain | AG Car B production input, FRESH t = 0, 480x128x128 = 4 x (480x64x64), production outputs kept | no tools | 60 |
| agc_kt10 / agc_kt60 | same | simple-kernel-timer | 10 / 60 |
| agc_tmr | same, input copy + `rad_m1/implicit_timers = 3` (code's own fenced category timers) | no tools | 60 |
| bsg_plain | BSG `bsg_hr_dc5.athinput`, reduced 4x1 mesh of the mem batches (256x128x128 = 4 x (256x64x64)), bin/rst dumps off | no tools | 30 |
| bsg_kt5 / bsg_kt30 | same | simple-kernel-timer | 5 / 30 |
| bsg_tmr | same + implicit_timers | no tools | 30 |

Command-line keys are fixed in prof_cuda.sh (time/nlim, time/ndiag=1; BSG also the reduced-mesh keys and
`output2/3/6/7/dcycle=1000000 output5/dt=1e30`); `implicit_timers` goes into an input COPY (a command-line key the
input does not have is rejected). Smoke rule: the first arm agc_plain is the smoke; the script stops the job itself
if it has rc != 0 or FATAL; if it dies of CUDA OOM, report and do not retry.

## 4. Analysis (login node)
```
R=$P/runs/<jobid>; L=$P/bundle/docs/handover/delta-prof-1009/labels_98835d99.tsv
python3 $P/bundle/docs/handover/delta-prof-1009/prof_group.py $L $R/agc_kt60 $R/agc_kt10 60 10 > $R/agc_prof.txt
python3 $P/bundle/docs/handover/delta-prof-1009/prof_group.py $L $R/bsg_kt30 $R/bsg_kt5 30 5  > $R/bsg_prof.txt
```
(steady cycles = long minus short run, mean over the 4 rank files; labels_98835d99.tsv maps each par_for label to
its source file at 98835d99a; groups are defined in RULES at the top of prof_group.py.) If the parser fails on the
kp_reader format, fix the parser (not the code) and say so.

## 5. NOTE back: `docs/handover/NOTE-2026-10-09-delta-prof-batch1.md` on THIS branch (rad-beam-1008)
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
- job id, node, elapsed (relative), GPU-h, binary md5; LEDGER line.
Attach the raw files (small) under `docs/handover/delta-prof-1009/results/<jobid>/`: both `*_prof.txt`, the kt_*.txt
of the four kt arms, and each arm's out.log tail (`tail -60`). No bins. Keep the run dirs.
