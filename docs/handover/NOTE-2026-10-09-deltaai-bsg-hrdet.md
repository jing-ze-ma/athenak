# NOTE DeltaAI: BSG half-range determinism gate -- RESULTS

**bit: PASS  det: PASS**  (DeltaAI GH200, layout B 2 nodes x 8 ranks, job 3347471)

Answer to TASK-2026-10-09-deltaai-bsg-hrdet (bsg-files-1009 ba268b1d).

- **Status: DONE.** Job 3347471 (nodes gh[049,082], ghx4-interactive) started 2026-10-09 09:18 CDT (16:18 CEST),
  elapsed 3 min 38 s (much shorter than the 25 min expected) = 8 GPUs x 218 s = **0.48 GPU-h** (budget 8).
- **Layout (B): 2 nodes x 8 ranks = 16 ranks, 2 MeshBlocks per GH200.** (A) not possible soon: ghx4-interactive
  MaxNodes=2, and a 4-node ghx4 job tests to start on 10-14. Binding: `srun -n 16 --ntasks-per-node=8 -c 36
  --cpu-bind=cores` + wrapper `CUDA_VISIBLE_DEVICES=SLURM_LOCALID/2` (local ranks 2g,2g+1 on Grace socket g and
  GH200 g; Kokkos then sees 1 device per rank). `--exclusive --mem=0 --time=01:00:00` (1 h instead of 30 min,
  because 2 blocks/GPU may run more slowly; only the time used is charged).
- Builds (build_inc_deltaai.sh hes_gpu, PROBLEM=he_star_m1, HOPPER90+ARMV9_GRACE, MPI on, default FMA, Cray CC host):
  - BO `athena_hes_gpu_d385de405fa8` (truerepro2-1009 d385de40) md5 `a66383005fed09b3b4f9409f74f13ab1`
  - BN `athena_hes_gpu_98835d99a16f` (hrup-bsg-1009 98835d99a) md5 `c9c6d16705d07dc6e97de6c808d02182`
- Files: SETUP.sh 3 x OK, SETUP_OK (`/work/nvme/bivj/jma20/bsg_hrdet_1009/files`).
- Gate script = bundle gate_hrdet.sh unchanged except one line after arm 1: abort if bit_old rc != 0 or FATAL (smoke
  rule). GPU memory is sampled every 20 s on both nodes during the whole job (gpumem.log).

## Results (all 4 arms)
| arm | binary | rc | fatal | diverged | nonconv | wall | hst NaN/inf | s/cycle raw (1->last) | s/cycle steady |
|---|---|---|---|---|---|---|---|---|---|
| bit_old (3 cyc) | BO d385de40 | 0 | 0 | 0 | 0 | 21 s | 0 | 3.26 | -- (every cycle writes) |
| bit_new (3 cyc) | BN 98835d99 | 0 | 0 | 0 | 0 | 28 s | 0 | 6.09 | -- (every cycle writes) |
| detH_a (30 cyc) | BN | 0 | 0 | 0 | 0 | 99 s | 0 | 2.84 | 0.86 |
| detH_b (30 cyc) | BN | 0 | 0 | 0 | 0 | 59 s | 0 | 1.36 | 0.86 |

"steady" = median of the cycles without bin output. The raw numbers are dominated by the bin dumps (every 5th cycle):
those cycles take 1.8-4 s, and two of them in detH_a (cycles 25, 30) took 26 s and 24 s, i.e. /work/nvme write
jitter, not compute (the same cycles in detH_b took 1.9 s and 4.0 s). Compute is identical in a and b.

**bit (bit_old vs bit_new):** 14 files compared (12 `bin/*.bin` + 2 `.hst`).
- bitcmp.py: every bin IDENTICAL (cycle 0 and cycle 3, time 262.106), both hst IDENTICAL, **BITWISE_PASS**.
- `cmp`: the 2 hst are byte-identical. The 12 bins differ ONLY in the header: BN dumps 5 more default keys
  (`implicit_flux_faces = x1`, `implicit_flux_beam = closure`, `implicit_blend_alpha = 1`, `implicit_blend_r0 = 1.5`,
  `implicit_blend_xthin = 0`, all "Default value added at run time"), +425 bytes. The bytes after `<par_end>` are
  identical in all 12 files.

**det (detH_a vs detH_b):** 34 files compared (32 `bin/*.bin` at cycles 0,5,...,30 + 2 `.hst`).
- `cmp`: **all 34 byte-identical**, headers included.
- bitcmp.py: **BITWISE_PASS**.

**GPU memory** (nvidia-smi, every 20 s on both nodes): max **38.1 GB per GH200 with 2 MeshBlocks** (about 19 GB per
block, versus the ~33 GB per block estimated from viper's unified memory). 4x4 would therefore fit easily too.

Note: `bitcmp.py` needs a `vis/python` that reads bin format 1.2 (here the 98835d99 tree) and h5py (imported by
bin_convert).

Run dirs kept: `/work/nvme/bivj/jma20/bsg_hrdet_1009/runs/3347471/{bit_old,bit_new,detH_a,detH_b}` (+ gpumem.log,
job log `logs/bsg_hrdet.3347471.out`). Nothing else queued on DeltaAI.
