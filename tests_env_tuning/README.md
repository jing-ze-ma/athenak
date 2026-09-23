# Runtime environment tuning on the APU GPUs (2026-09-23)

Question: does any runtime or environment setting make the production GPU runs faster without
changing the output, as `HSA_NO_SCRATCH_RECLAIM=1` did (22.5 s to 17.0 s per 300 cycles,
`tests_m1/runs_3k_gpu3d/README_HALO.md`)?

**Answer: no.** Nothing beyond the current production setting (`HSA_XNACK=1` plus
`HSA_NO_SCRATCH_RECLAIM=1`) gives a gain above the run-to-run noise with bitwise-identical
output. Keep the submit scripts as they are.

## Setup

- Work dir: `/viper/u2/jinma/ATHENAK/bench/env_tuning_0923`. It has the scripts `common.sh`
  and `job1.sh` to `job6.sh`, the logs `log.out.<jobid>`, and one run dir per arm under
  `runs/<arm>/`.
- **prod4 test case.**
  - Binary: a copy of `bench/cs_mhd_prod4/athena` (md5 `f95130b2...`, commit 395db5bc).
  - Input: a copy of `bench/cs_mhd_prod4/deep_hot_jupiter.athinput`.
  - Run: restart from `bench/cs_mhd_prod3/rst/dhj.00567.rst` for 300 cycles
    (`time/nlim=4430131`), 2 ranks on 2 GPUs on 1 node.
  - Launch: `srun -n 2`, with `ROCR_VISIBLE_DEVICES=$SLURM_LOCALID`. This is the recipe of
    `bench/m1_halo_0923/gate/gpu_v2p4.sh`.
- **M1 3-D box.**
  - Binary: `bench/m1_halo_0923/bin/athena_new2_gpu_box_convection`. This is dccdf506 plus
    halo.patch, which is the src content of d14f96de for this problem; the two_stream change
    in c02ddd34 is default off.
  - Input: `tests_m1/runs_3k_gpu3d/he_slab_m1_3d.athinput` with `implicit_line_solver = pcr`
    and `implicit_bcg_sync = 1` (= `gate/inp/m1_3d.athinput`).
  - Run: 200 cycles on 2 GPUs.
- **Timing:** the `cpu time used` line, which is wall time in the main loop. All runs were on
  apudev, on nodes vipa1327 and vipa1001. Arms were interleaved, with the second repeat in
  reversed order.
- **Identical-output check:**
  - prod4: payload md5 (`bench/m1_halo_0923/tools/pmd5.py`, bytes after `<par_end>`) of the
    final `rst/dhj.00568.rst` and `bin/dhj.mhd_w_bcc.00143.bin`, plus the md5 of
    `dhj.mhd.hst`.
  - Baseline prod4 hashes: rst `8a2b4699`, bin `dfc3d616`, hst `be7ad322`.
  - M1: md5 of both `.hst` files. Baseline: `dba0d078`, `0368269a`.
- **Baseline:** 9 runs over 5 jobs. Mean 17.23 s, standard deviation 0.15 s (0.9 %).
  - Runs: 17.27, 17.20 (job 11944891); 17.15, 17.07 (11944893); 17.32, 17.17 (11944912);
    17.02 (11944926); 17.40, 17.49 (11944939).
  - The same-job baseline moves by up to 0.3 s between jobs, so differences under about 2 %
    are noise.

## Results: prod4, 300 cycles, 2 GPUs

Each row is the baseline plus the one setting shown.

| setting | what it does (source) | runs, s (job) | mean vs baseline | output identical |
|---|---|---|---|---|
| baseline `HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1` | production setting | 9 runs above | 17.23 | reference |
| `HIP_FORCE_DEV_KERNARG=1` | Stores kernel arguments in device memory to cut launch latency. The ROCm doc gives default 1, so this is effectively a no-op. [1][2] | 17.25, 17.07 (11944891) | -0.4 % (noise) | yes |
| `HIP_FORCE_DEV_KERNARG=0` | Kernel arguments stay in host memory. [1][2] | 17.28, 16.88 (11944891) | -0.9 % (noise) | yes |
| `HSA_XNACK=0` | Turns off XNACK page-fault retry, so the GPU can no longer demand-fault host (malloc) pages. [1] | 16.66, 16.74 (11944891); 16.63, 16.84 (11944912) | -3.0 % | prod4 yes. **M1 box CRASHES**: "Memory access fault by GPU", rc 134 (11944926, `runs/m1_xnack0_{a,b}`). Kokkos warns at init that ARCH_AMD_GFX942_APU requires `HSA_XNACK=1`. **REJECT** |
| `HSA_ENABLE_SDMA=0` | Copies use blit kernels instead of the SDMA engines. [1][3] | 17.82, 17.94 (11944891) | **+3.7 % slower** | yes |
| `HSA_ENABLE_INTERRUPT=0` | Completion signals are polled instead of interrupt-driven. [1][3] | 17.36, 17.36 (11944891) | +0.7 % (noise) | yes |
| `GPU_MAX_HW_QUEUES=1` | Maximum number of hardware queues per device; default 4. [1][2] | 16.86, 17.10 (11944893) | -0.8 % against same-job baseline 17.11 (noise) | yes |
| `GPU_MAX_HW_QUEUES=8` | as above | 17.03, 17.02 (11944893) | -0.5 % against same job (noise) | yes |
| `srun --cpu-bind=cores` | Binds each rank to its cores. The default is already bound (see Binding). | 17.12, 17.10 (11944893) | 0 % against same job | yes |
| `srun --cpu-bind=ldoms` | Binds each rank to its NUMA domain, which is also its GPU's domain. | 16.79, 17.05 (11944893) | -1.1 % against same job (noise) | yes |
| `OMP_PROC_BIND` / `OMP_PLACES` | Not applicable: the binary has no OpenMP (no GOMP/omp symbols; the Kokkos host space is Serial). | not run | none | not run |
| `UCX_PROTO_ENABLE=y` | Turns on the UCX v2 protocol selection. The module sets `n`. [4] | 17.28, 17.49 (11944912) | +0.9 % (noise) | yes |
| `UCX_RNDV_THRESH=inf` | Never use the rendezvous protocol, always eager. [4] | 18.08, 17.70 (11944926) | **+3.8 % slower** | yes |
| `UCX_TLS=self,sm,rocm_copy,rocm_ipc` | Restricts UCX to shared-memory and ROCm transports. | fails: `ucp_ep_create ... Destination is unreachable`, rc 14 (11944912) | none | did not run |
| `UCX_TLS=^ud,tcp,rc,dc` | Drops the InfiniBand RC/DC loopback transports. | fails: same error (11944926) | none | did not run |
| `OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm` | Uses Open MPI's own shared-memory path instead of UCX. Device buffers are staged through the rocm accelerator component. | 29.56, 29.89 (11944939; same-job baseline 17.40, 17.49) | **+71 % slower** | yes |
| combo `GPU_MAX_HW_QUEUES=1 HIP_FORCE_DEV_KERNARG=0 --cpu-bind=ldoms` | the three best non-rejected rows together | 17.26, 17.36 (11944912; same-job baseline 17.32, 17.17) | +0.4 % (noise) | yes |

`OMPI_MCA_pml=ob1` without `OMPI_MCA_btl=self,sm` exits with rc 14 at startup (11944912).

## M1 3-D box, 200 cycles, 2 GPUs (job 11944926)

| arm | runs, s | mean | hst identical |
|---|---|---|---|
| baseline | 31.36, 31.09 | 31.23 | reference |
| combo (as above) | 31.99, 31.85 | 31.92 (+2.2 %) | yes |
| `HSA_XNACK=0` | crash, GPU memory access fault | none | rejected |

## MPI / UCX configuration in use (job 11944891, from `ompi_info`, `ucx_info`, `env`)

- **Libraries and default PML:** Open MPI 5.0.8 (openmpi_gpu/5.0) and UCX 1.18.0 (UCX-GPU).
  The PML is ucx (priority 51; the verbose run in `runs/ucxdbg` confirms that ucx is
  selected). The accelerator component is rocm.
- **Environment set by the modules:** `UCX_TLS=^ud,tcp`, `UCX_PROTO_ENABLE=n`,
  `UCX_UNIFIED_MODE=y`, `OMPI_MCA_coll=^hcoll`.
- **UCX devices on the compute node:** rocm_cpy, rocm_ipc, memory, mlx5, and more.
- **Why UCX was tested:** the profile shows the 2-rank exchange (the prod4 idle is MPI waits
  plus fences, per README_HALO.md), so the UCX knobs above were tried. None of them helps.
  The remaining wait is the fences and load imbalance, not the transport.

## Binding (job 11944891, `taskset` in each rank)

- **NUMA layout:** `lscpu` shows 2 NUMA nodes (CPUs 0-23,48-71 and 24-47,72-95).
  `rocm-smi --showtopo` puts GPU0 on NUMA 0 and GPU1 on NUMA 1, linked by XGMI.
- **Default binding is already right:** with `--cpus-per-task=24`, Slurm binds rank 0 to
  cores 0-23 and rank 1 to cores 24-47. Each rank sits on its own GPU's NUMA domain, which
  is why explicit binding changes nothing.

## Recommendation

- Keep the production setting: `HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1` and the default Slurm
  binding. Add nothing.
- Do not use:
  - `HSA_XNACK=0`: it crashes the M1 box, and Kokkos requires XNACK on MI300A.
  - `HSA_ENABLE_SDMA=0` (+3.7 %) and `UCX_RNDV_THRESH=inf` (+3.8 %).
  - pml ob1 (+71 %).
- Further speed has to come from code, for example the 321 `hipDeviceSynchronize` per cycle
  in the MPI path noted in README_HALO.md.

## Sources

1. ROCm documentation, "ROCm environment variables",
   https://rocm.docs.amd.com/en/latest/reference/env-variables.html (fetched 2026-09-23). It
   gives these defaults: `HIP_FORCE_DEV_KERNARG` 1, `GPU_MAX_HW_QUEUES` 4, `HSA_ENABLE_SDMA`
   1, `HSA_NO_SCRATCH_RECLAIM` 0, `HSA_ENABLE_INTERRUPT` 1.
2. The variable names appear in `strings` of ROCm 6.3.4 `lib/libamdhip64.so.6`. That library
   carries no descriptions.
3. The variable names appear in `strings` of ROCm 6.3.4 `lib/libhsa-runtime64.so.1`.
4. UCX 1.18 on this system. The `UCX_*` meanings and defaults come from `ucx_info -cf`
   (for example `UCX_RNDV_THRESH=intra:auto,inter:auto`). No local HTML docs were found under the ROCm or UCX prefixes.
