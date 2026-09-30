# TASK for Caltech: WASP-121b 1x at C256 x nx1 256 from Caltech's own rot-300 restart (remap + benchmark)

**User 09-29:** the 1x fresh start reached rot 300 (viper and Caltech). The science continuation is a high-resolution
run: cubed sphere **C256 (256 x 256 per panel)** and **nx1 256**, 1.0e8 cells. The viper apu queue is long for
multi-node jobs, so Caltech (H200, 4 per node) is the candidate machine. This task: remap + benchmark only.
Do **not** start the production; report the numbers and the user decides the node count.

## What exists
- `problem/remap_file` (pgen key, off by default) is on rt-integration since e804daf1: build any rt-integration
  commit at or after d26b7364 (current tip) with your incremental `build_inc.sh` (dhj target, PROBLEM =
  deep_hot_jupiter_rt, as for w121prod_0928).
- The remap script with **horizontal prolongation** (C_n -> C_(f n), conservative: mass / momenta / eint to
  round-off, limited piecewise-linear, positivity-scaled slopes) is on branch **fork/dhj-remap-h** (9896d452, only
  `docs/handover/scripts/dhj_remap.py` changed). Take the script from there:
  `git show fork/dhj-remap-h:docs/handover/scripts/dhj_remap.py > dhj_remap.py`.
  - viper gates: identity C32 -> C32 bitwise; conservation table below; C256 restart smoke/benchmark on viper
    still queued (8 and 32 nodes), so **your step 3 is also the first C256 smoke**.
- The radial production grid `G_PROD` is in `docs/handover/caltech-2026-09-26/inputs/wasp121_1x/grid_w121_1x.env`.

## Steps
1. Build (rt-integration tip), record the md5.
   **C256 production (user 09-29):** must use rt-integration >= ae767d20 (merge of dhj-ck-conserve) with
   `problem/ck_impl_conserve = 1` (add it to `remap.athinput`, `<problem>` block; solution bitwise unchanged,
   reported ck fluxes consistent with the applied energy in bound cells). NOTE: the same merge sets
   `problem/flux_hst_rkavg = true` by default: the dhj hst flux columns are now RK-weighted fluid fluxes
   (semantics change vs older hst files; set it false to reproduce the old columns).
2. Remap (login node, ~20 s, python3 + numpy):
   ```bash
   source docs/handover/caltech-2026-09-26/inputs/wasp121_1x/grid_w121_1x.env
   RST=/resnick/groups/carnegie_poc/jingze/w121prod_0928/1x/rst/dhj.00600.rst   # your rot-300 rst; read in place
   python3 dhj_remap.py $RST <newdir>/remap_c256 \
     --grid "$G_PROD mesh/nx2=256 mesh/nx3=256 meshblock/nx2=32 meshblock/nx3=32" | tee <newdir>/remap_c256.txt
   ```
   Check the printout as viper's (`/viper/ptmp2/jinma/w121_c256_0929/remap_c256.txt`): mass / eint change
   <= 1e-15 horizontally, radial column changes <= 1e-13; E_kin -0.3..-0.8 % and rho*Phi a few % (expected);
   etot -9e-4 (radial step, as in the 1-D remap). NOTE viper's printout reported **1.38 M cells raised to
   dfloor = 1e-16** (1.4 % of the cells, presumably the top of the new radial grid): report your count and where. **Checked on viper (09-29):** it is INHERITED, not a remap defect: the rot-300 source already has 4214 cells (0.9 %) at dfloor, in the top 5 radial cells (23.6 % / 14.7 % / 12.2 % / 4.4 % / 1.6 % of columns, i = 75..71); radial-only remap 76 -> 256 gives 21441 (the finer top grid splits each floored cell), horizontal-only 292, full C256 = 21441 x 64. All at p <~ 1e-9 bar (outside the accuracy region).
3. Benchmark: fresh start `-i <newdir>/remap_c256/remap.athinput` (NOT `-r`), meshblock 256x32x32 (384 blocks),
   ~70 cycles, `time/ndiag=5 problem/ck_impl_verbose=true`, no bin/rst outputs, on **4 nodes (16 H200)** and
   **8 nodes (32 H200)** (H200 only, exclude hpc-sm-01-09 and hpc-sm-02-16 as before). Template: viper's
   `/viper/ptmp2/jinma/w121_c256_0929/bench.sub` (copied below in essence):
   `srun $BIN -i remap.athinput -d $D time/nlim=70 time/ndiag=5 time/tlim=3.4e7 problem/ck_impl_verbose=true -t 00:24:00`
4. Report: ms/cycle (cycles 20-70), zone-cycles/s per GPU, dt and its limiter (radial vs horizontal; where),
   ck Newton passes and NOT-CONVERGED, GPU memory per GPU, any CUDA local-memory / scratch error (viper's 09-24
   test hit HSA scratch limits at n1 256 with the ck kernels; CUDA limits differ), queue wait of the 8-node job.
   From these: GPU-h per rotation and wall per rotation at 8 / 12 / 16 nodes.

## Viper estimate to compare (not measured at C256)
~40-50 MI300A-GPU-h per rotation (dt ~4 s, radial CFL); H200 expected ~1.6x faster -> ~25-30 H200-h per rotation;
8 nodes ~0.8-1 h/rot, 12 nodes ~0.55-0.65 h/rot.

## Rules
Timing and smoke only; no production, no code changes. Results: `docs/handover/NOTE-2026-09-29-caltech-c256.md`
on a new branch `bench-results-caltech-c256` (push to fork). Pending on viper and not needed for this task:
the ck conservation fix (branch dhj-ck-conserve, energy source of 0.26 % L in the top 2 cells, in gates) will
go into rt-integration before the C256 production.

## UPDATE 09-29 ~21:00: GPU MEMORY (viper C64 x 256 smoke, job 12027341)
The C64 x nx1 256 smoke (24 blocks of 256x32x32 on 2 MI300A, rt-integration ae767d20, ck_impl_conserve = 1) ran clean (281 cycles, 0 ck NOT-CONVERGED, dt 3.7-4.05 s = radial CFL as at C32 x 256; offline estimate at C256: radial 2.6 s vs horizontal 9 s, so dt stays radial ~3.9 s) but used **~103 GB per GPU at 3.1M cells per GPU (~33 KB/cell)**. If that scales with cells, C256 (1.0e8 cells) needs ~3.3 TB of GPU memory: **4 H200 nodes (16 x 141 GB, 6.3M cells/GPU) will NOT fit; 8 nodes (32 GPUs, 3.1M cells/GPU, ~103 GB) should just fit; 12 nodes is comfortable.** So benchmark on **8 and 12 nodes** instead of 4 and 8, and record the peak memory per GPU (nvidia-smi). If 8 nodes runs out of memory, try meshblock 256x16x16 (more, smaller blocks) before going up in nodes. Speed at C64: 3.1e7 zone-cycles/s per MI300A.

## UPDATE 09-30 ~01:15: DEEP-MIXING QUESTION (user) - part of the C256 production analysis
The 10x ledger (viper /viper/ptmp2/jinma/w121prod_0929/budget_rot300_10x, closed to 1e-8 of L) found the deep interior
(10-430 bar) slowly COOLING at C32 (~-0.9 K per 100 rot): the fluid enthalpy flux grows upward from 1.15e28 erg/s at the
face above the wall to 4.4e28 at 10 bar while the bottom supplies only Lrad_bot = 9.2e27 (effective T_int ~1.5x the
input), draining the initial deep adiabat. The 1x deep warms slightly instead. Open: is this upward deep enthalpy flux
resolved circulation or numerical mixing? At C256 (and, for the resolution trend, C64), record from the first rotations on:
per-shell radial energy flux split into mean-meridional <rho v_r>_(theta,phi) h and eddy parts, at 1/10/100 bar and
the wall face; the deep isobar T drift per 10 rot at 10/100 bar and the bottom; compare with C32 (same diagnostics from
w1x/w10x dumps, ana_rot300*/). If the deep flux falls with resolution it is numerical mixing.

## UPDATE 09-30 ~02:30: build >= 793e03c3 for the C256 production
DeltaAI found a cubed-sphere seam flux MPI request leak on uniform meshes (NOTE-2026-09-29-deltaai-seam-mpi-leak.md; one leaked Isend request per off-rank seam neighbour per stage; Cray MPICH aborts after ~8100 cycles) and a CUDA-only dhj MHD startup segfault (NOTE-2026-09-29-deltaai-c2p-track-cuda.md). Both fixed on rt-integration (793e03c3, ad992f7b). The C256 production (hundreds of ranks, long) must use rt-integration >= 793e03c3. The benchmark may use the older build.
