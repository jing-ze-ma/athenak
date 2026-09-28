# Caltech session 2026-09-25 -> 09-27: state and next steps (START HERE on Caltech)

Read this, then the local memory index (`~/.claude/projects/-resnick-home-jingze-ATHENAK/memory/MEMORY.md`).
A snapshot of those notes is in `docs/handover/claude-memory-2026-09-27-caltech/`. The project-wide handover
is `HANDOVER-2026-09-26.md`: its top blocks list every merge.

The Caltech cluster (Resnick HPC) is for **testing, not production** (user). Production runs on viper.
`rt-integration` is at **50325e8b**, pushed, clean.

## 1. Merged from Caltech (all pushed; viper HIP-gated unless noted)

| merge | what | gates / numbers |
|---|---|---|
| 83876228 | CUDA port: nvcc lambdas, host mirrors (EOSTable), DeepCopyAcross, mb_panel sync_device, kokkos submodule restored | CPU bitwise; GPU vs CPU round-off; restart and 2v1 bitwise; viper HIP bitwise |
| 504d8e30 | defaults: ck_impl_kkt_row = true, ck_impl_xstep = 8 (every > 1) | bitwise gates |
| 48321155 | box_convection + rad_m1 on CUDA (Time2RstSet strided copies) | CPU bitwise; GPU round-off |
| 2e028a38 | Kokkos 4.4.00 -> 4.6.02 (the old pin was 4.4, not 4.6) | all bitwise; ~4 % faster median |
| f8f7232f | WASP-121b inputs: ck_beam_par = true; mg_gf default DROPPED | beam_par -6.2 % on 2 H200; mg_gf +33 % on 2 H200 |
| 587825ef | ck-jlin: cache jlin coefficients, rt_apply j-fastest, fused pass 1 | -13.6 % (2 H200), -16.1 % (1) |
| dff69bc7 | rad_m1 implicit_one_pass_auto (default true): one_pass switches off per solve kind when unpaid | bitwise where it never triggers; restart bitwise |
| 55762ada | ck-store: problem/ck_store_split (storing pass without the chain kernel), lP 9 -> 6 | -7 % more; noise 1.21x |
| 6a9ffb0d | ck-next: beam_tau chord sharing, ck_coef layout, lin_sum, prefetch, re-formed cin/cout | -9.4 % (2 H200), -11 % (1); noise 0.78x. **HIP gate pending on viper** |

Cumulative ck speed-up on H200, nx1 256, production keys:
- 2 GPUs: 32.9 -> 24.0 ms/cycle (-27 %).
- 1 GPU: 59.3 -> 40.5 ms/cycle (-32 %).
- A hydro + cad-only arm is ~14.8 (2 GPUs) / 23.4 (1 GPU), so roughly half of the ck cost is gone.

## 2. Decisions taken (user)

- **Kokkos 5.2.2 NOT adopted.** It is only 1-3 % on dhj and 0 % on M1, and needs a 5096-line DualView rename
  plus C++20. Branch `kokkos5` is kept.
- **one_pass kept.** Its 5-22x round-off deviation (tolerance level) is accepted; the auto-disable is added.
- **Default-on switches that measure as not helping are dropped** (mg_gf so far). A switch that helps on
  MI300A but hurts on H200 gets a backend-dependent default instead.
- **The accuracy bar** is the round-off spread: a noise gate with >= 2-3 1e-14 kick members, deep eint rms
  at 1-100 bar, pass <= ~1.2x.

## 3. Next steps (not started)

**ck, remaining per storing call (1 H200), from the ck-next agent:**
- lin1p: ~8 ms per pass, DRAM-bound. Next idea: do not store ck_ci/ck_co at all.
- ck_beam_tau: ~18.6 ms, compute-bound. Try a shared-memory column tile.
- ck_lin_build: ~15.7 ms, latency-bound.
- ck_coef: ~9.7 ms.

Optional: find the exact GPU contraction that breaks bitwise at ck-next cn10 (diagnostic outputs in
`/resnick/scratch/jingze/cknext/g7`, `g8`; scratch is purged after 14 days).

**M1 on H200:**
- The profile is launch-latency bound: 1344 launches/cycle, 9 % idle.
- ~12 cudaFree + ~520 cudaMalloc per cycle, i.e. per-cycle temporary Views, to be found.
- M1Evt falls back to a full fence on CUDA.
- mg_gc reads pinned memory over PCIe.

## 4. How to work here (details in memory `caltech-cluster-facts`)

- **Slurm**, account `carnegie_poc`, partition `gpu`, 4 GPUs per node (18 H200 nodes, 4 H100).
- **NEVER use node hpc-sm-01-09.** It is a bad H200 (~6.6x slow, out of GPU memory).
  - `run_caltech.sub` already excludes it; add `--exclude=hpc-sm-01-09` to every other GPU job.
  - Pending jobs: `scontrol update JobId=<id> ExcNodeList=hpc-sm-01-09`.
- **QOS**: `debug` (<= 30 min, 1 running job per user, top priority) for tests; `normal` otherwise, with a
  tight `-t` so backfill can start it. H200s are usually full.
  - Correctness runs may use H100 (`--gres=gpu:h100:N`, same sm_90 binary).
  - **Timing only on H200, same binary, interleaved arms in ONE job**: node-to-node differences are ~7 %.
  - Report the median of 8-cycle windows; output stalls ruin means.
- **MPI**: `srun --mpi=pmix` (the site default pmi2 starts singletons). HPC-X is the CUDA-aware MPI.
- **Build**:
  - `docs/handover/caltech-2026-09-26/scripts/build_inc.sh <tag> <cpu|gpu> <commit> [problem]` does
    incremental builds from a git archive snapshot. Problems: `deep_hot_jupiter_rt` (default),
    `box_convection`, `built_in_pgens`.
  - GPU builds go on a compute node: `sbatch -A carnegie_poc -p expansion -c 32 --mem=64G -t 01:00:00
    --wrap "NJ=32 bash build_inc.sh ..."`. Login nodes are overloaded and capped at 8 GB per process.
- **Gates**:
  - `cpugate_caltech.sub` (A/B CPU bitwise, dhj);
  - 2 vs 1 GPU with `/resnick/scratch/jingze/val_p4/run2rank.sub`;
  - noise: `/resnick/scratch/jingze/cknext/noise.py`, acc arms in `cknext/arms.sh`;
  - timing: `/resnick/scratch/jingze/ckjlin/ana2.py`.
  - The hst 3-mom column is a cancelling sum and differs at ~1e-3 between rank counts, so compare rst/bin,
    not hst bytes.
- **Storage**: scratch `/resnick/scratch/jingze` is purged after 14 days idle; keep things worth keeping in
  `/resnick/groups/carnegie_poc/jingze/`.
- **GitHub push**: SSH key `~/.ssh/id_ed25519`; the repo push URL is `git@github.com:jing-ze-ma/athenak.git`.
- **Agents**:
  - Opus 5.5 at medium effort; no haiku/sonnet (user).
  - One agent per branch and worktree. An agent that hands back with jobs in flight can RESUME later, so stop
    it before starting a replacement (memory `agent-collision-lesson`).
  - Agents must not merge or push; the coordinator merges after the user's go-ahead.

## 5. Update 09-27 late: ck-lin2-fma MERGED (8be0ad67)

- -9.9 % on 2 H200 and -10.3 % on 1 H200. Bitwise on CPU and on CUDA.
- **HIP gate pending on viper**: see `TASK-2026-09-27-cklin2-hip-gate.md`.
- Dropped variants (records only): l3 (lin_build per chain pair), CK_BTG 16, team 256.
- Next ck cost per storing call on 1 H200, in ms: lin_build 10.6, jlin 9.7, beam_tau 6.9, lin1p 6.3, coef 3.1.
