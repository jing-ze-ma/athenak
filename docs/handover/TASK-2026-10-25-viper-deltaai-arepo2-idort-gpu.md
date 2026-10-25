# TASK viper -> deltaai (USER): Arepo2-IDORT GPU on DeltaAI (aarch64): build + gates + g25_red smoke now, g25_full 2 nodes later

Context: Caltech TASK 8513c5b56 + addendum 81b32ab29. On Caltech the GPU IDORT build gives g25_full 2.8x per top step on 1 H200 node
(32 ranks). The comparison with Orion (35.6 s per top step on 560 cores) needs ~576 ranks: Raven does 8 nodes, DeltaAI does 2 nodes.

**Rules**
- No Arepo source, diffs or patches in git.
- Never push to Bitbucket. Arepo2-IDORT is private.
- Code comes as a git bundle that the user copies.
- DeltaAI rules (10-09): ONE job at a time, <= 1 h, <= 2 nodes, all arms in sequence inside that one job.
- Plain `srun -n N --cpu-bind=cores`. `--gpu-bind=closest` hangs multi-rank runs here.

## Stage 1 (now): the 10-24 bundle the user copies
Source on Raven: `/raven/ptmp/jinma/arepo2_gpu/caltech_bundle_1024/` (192 MB + MD5SUMS). It is the same bundle Caltech used:
- `arepo_idort_gpu_1024.bundle`: 5 branches, full history; md5 905ec3cb6ad069aa1e491658956b3da8.
- `data/g25_inputs.tgz`, `data/recipes.tgz`, `data/submodules.tgz`: see MD5SUMS.

Unpack:
1. `git clone -b idort-gpu-s1 arepo_idort_gpu_1024.bundle arepo2`, then `tar xzf submodules.tgz -C arepo2`.
2. Remove the stale `lib/*/.git` gitlink files.

Build (aarch64, GH200):
- HYPRE 2.32.0, CPU only, from recipes.tgz.
- Arepo configs from `gpu1/cfg`, with the CUDA arch set to **sm_90**.
- `-march=native` is fine if you build on the compute node; otherwise use the Grace target.
- Keep `-ffp-contract=off` and nvcc `-fmad=false`: bitwise CPU = GPU depends on them.
- Drop `-lbfd` if it is missing (only BACKTRACE uses it).
- Caltech's other adaptations are in its NOTE 23b200296.
- **Layout:** NumIslandsPerNode 4 (one per GH200); ranks per node divisible by 4 (uneven islands hang); CPUsUsedPerGPU = ranks / 4.

In ONE job (1 node, <= 1 h):
1. The IDORT gates from `gpu1/gate_g1.sh` (cases in `gpu1/tcases`). IdortGpu 0 vs 1, same binary; every gate must say ALL SAME.
   - Caltech list: A1, A4, A4_m4, A6, A2, I2a, B3, A7_rv, I4_fhyb_levels_m4, I1_restart_bitwise.
   - If they do not all fit, run as many as fit and list the rest.
2. A g25_red smoke at 288 ranks, 2 top steps: the CPU arm and the GPU arm of idort-gpu-s1 (`gpu1/time_g1.sh`).
   - Report per full sync and per top step.
   - Reference per full sync: Raven 72 ranks CPU 7.19 s / GPU 3.90 s; Caltech 32 ranks CPU 8.69 s / GPU 4.17 s.

## Stage 2 (later, after viper gates the combined code)
- Code: an incremental bundle with `idort-gpu-nltree` (idort-gpu-s1 + nltree-g25 + Caltech's returned branches + a g25 direct-summation
  fix if confirmed).
- **g25_full on 2 nodes (576 ranks), three arms:**
  1. GPU recipe (GPU IDORT + GPU mesh + GPU node-local tree);
  2. CPU classic tree;
  3. CPU node-local tree.
  - At least 3 top steps after the first.
  - Per arm: s per top step (mean, min, max), split into IDORT, mesh, gravity, CPU partial steps, rest; peak memory per node.
  - Also quote per-core rates: Grace cores differ from Orion's x86.
- If memory allows, the GPU arm on 1 node (288 ranks) as well.
- A separate TASK will follow when the bundle is ready.

## Reply
A NOTE deltaai->viper on bsg-files-1009 with:
- the node facts (cores, GPUs, MPI, CUDA, whether CUDA IPC works between the 4 GPUs of a node);
- the build adaptations;
- the gate list with ALL SAME / differences;
- the smoke timings.

## WITHDRAWN (USER 10-25): DeltaAI is too crowded; do not run this TASK
The g25_full scaling test runs on Raven only (8 nodes, 576 ranks). If anything was started, stop it. No reply is needed.

### Exception (USER 10-25): let job 3454468 run
Keep the already-queued job 3454468: the B3 case3 MPICH-env A/B for the island parallel sort, then its follow-ups.
- It is a diagnosis of the aarch64 sort hang.
- Push its results as a NOTE (which env, if any, fixes the sort).
- Then stop. No stage 2, and no new jobs.
