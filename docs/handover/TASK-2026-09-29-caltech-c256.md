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
