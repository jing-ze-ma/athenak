# NOTE 2026-11-21gu orion -> caltech (USER via coordinator): are some GPUs faulty? Cheap check alongside the rebuilds (NOTE gt) -- test nodes only

The USER asks whether some Caltech GPUs are simply faulty. Each check takes ~1 min. Run them **alongside** the NOTE gt rebuild; do not delay it.
1. Rerun the **same tiny 1-GPU B_v203 smoke** (NOTE gn / gq: bsg_wngl_smoke, nx2=nx3=32, NR=1, nlim 1-2) on:
   - **hpc-sm-02-03** (GPU 0);
   - **hpc-sm-02-16** (GPU 0);
   - **hpc-sm-02-09 with `CUDA_VISIBLE_DEVICES=3`** (the last GPU of the node used so far).
   Report the `clean_init steady pass 0/1/...` lines for each. **Trapped = pass 1 exactly 3.499e-08**, converging in 4 passes. Healthy = pass 0 9.534e-02 and NOT converged after 50.
2. On each of those nodes (02-03, 02-16, 02-09), report the output of `nvidia-smi -q | grep -iE 'ecc|retired|remapp|xid|errors'` for all 4 GPUs. Also include `nvidia-smi --query-gpu=index,name,driver_version,vbios_version --format=csv`.
Commit everything to caltech_results/H200CHK/gpuhealth.txt.

-- orion (coordinator, for the USER)
