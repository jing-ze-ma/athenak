---
name: caltech-cluster-facts
description: "Caltech (Resnick HPC, login3.cm.cluster) facts checked with sinfo/sacctmgr 2026-09-25 -- Slurm, gpu partition with 4 GPUs/node (H200/H100/L40S/P100/V100), account carnegie_poc, QOS limits"
metadata:
  node_type: memory
  type: reference
  originSessionId: 1f4629f7-1667-4eeb-be90-3eb8315950c0
  modified: 2026-09-26T00:06:38.264Z
---

Checked 2026-09-25 with `sinfo`, `scontrol show partition`, `sacctmgr` (not assumed). Replaces viper facts ([[viper-2-gpus-per-node]] does NOT apply here).

- Scheduler: **Slurm**. Login host `login3.cm.cluster`. Home `/resnick/home/jingze` (= `/home/jingze`). Repo: `/resnick/home/jingze/ATHENAK/athenak`.
- Account **carnegie_poc**; QOS: `debug` (30 min, 1 job), `normal` (7 d), `long` (14 d, 1 job).
- Partition **gpu** (MaxTime 14 d), **4 GPUs per node**:
  - 18 nodes `gpu:nvidia_h200:4`, 128 cores, 1.5 TB -> Kokkos_ARCH_HOPPER90 (production target)
  - 4 nodes `gpu:h100:4` (HOPPER90), 4 nodes `gpu:nvidia_l40s:4` (ADA89), 46 `gpu:p100:4`, 2 `gpu:v100:4`
  - Request by type: `--gres=gpu:nvidia_h200:N` / `--gres=gpu:h100:N`.
- Other partitions (probably not ours): sunshine/dgxlo (v100x8), vanvalen (b300x8), expansion (CPU, default).
- Modules: `cuda/12.9.0-none-none-pfmzfdv`, `cuda/12.2.1-gcc-11.3.1-*`; `gcc/13.2.0`; `openmpi/4.1.x`, `5.0.1` (gcc 11.3.1 / 13.2.0). No nvcc/mpicxx on PATH by default.
- NVIDIA: no HSA_* env vars; one MPI rank per GPU (`srun --gpus-per-task=1`).
- Viper ms/cycle numbers do not transfer; re-measure on H200.

**From the site docs (hpc.caltech.edu/docs, read 2026-09-25):**
- **Scratch `/resnick/scratch`: 20 TB/user, files not accessed in 14 days are PURGED.** Long runs/results belong in `/resnick/groups/carnegie_poc/` (group 20 TB default, $8/TB/month over) or home (57 GB quota).
- **Billing is real money:** units/h per GPU: H200 156, H100 120, L40S 61, P100/V100 10; CPU core 1. Price per unit $0.012 (first $6.4k/yr), $0.007 up to $24k, $0.006 above -> one H200 ~ $1.9/h (tier 1). Carnegie = 19 % partner; no Slurm-side cap. Ask the user before long/large GPU jobs.
- **Login nodes: 8 GB memory cgroup per process**; heavy analysis on compute nodes (viper's 40 GB rule does not apply).
- MPI: docs show plain `srun`; site default plugin is pmi2 but OpenMPI 5/HPC-X need `srun --mpi=pmix` (without it: N singletons, seen 09-25 job 3471478).
- NVHPC 23.7-26.3 modules exist (alternative CUDA-aware MPI); HPC-X 2.17.1 ompi is CUDA-aware (ompi_info).
- Help: help-hpc@caltech.edu (extended walltime, software, billing).
- **Production run dir (user 09-25): `/resnick/groups/carnegie_poc/jingze/`** (created). Scratch only for smoke/short tests. Group quota usage not visible from the shell (VAST; `quota`/du too slow); user checking with Carnegie.
- **Build GPU binaries on a compute node**, not the login node (login load ~195 on 128 cores, 53 users; nvcc build took >40 min there): `sbatch -A carnegie_poc -p expansion -c 32 --mem=64G -t 01:00:00 --wrap "NJ=32 bash .../build_caltech.sh <tag> gpu <commit>"`. CPU (gcc) builds on login are OK. When killing processes, use `pgrep -f "[p]attern"` (a plain pkill -f matches the tool shell itself and kills it).
- **Incremental builds (09-25): `docs/handover/caltech-2026-09-26/scripts/build_inc.sh <tag> <cpu|gpu> <commit>`**. It keeps persistent trees `builds/inc_<dev>_<problem>` and rsyncs the git-archive snapshot in by content, so only changed files recompile. A no-op CPU rebuild takes about 2 min; a full CPU build 6.5 min. It takes a flock per tree. Binaries differ from build_caltech.sh ones only in embedded __FILE__ paths.
- **Queue (measured 09-25):** afternoon 1-2 H200 on normal waited 3-7 min; at 22:00 the H200 nodes were full and Slurm estimated ~19 h for 2-4 GPUs and ~23 h for 8 GPUs (2 whole nodes). **Use `-q debug` for tests**: very high priority (started in seconds), <= 30 min, <= 10 nodes, 1 running job per user. The account is above its fair share, so its priority is lower.
- **GitHub push from Caltech:** SSH key ~/.ssh/id_ed25519 (added to the jing-ze-ma GitHub account by the user 09-25); repo remote origin fetches over https and pushes via git@github.com:jing-ze-ma/athenak.git.
- **Slurm priority (checked 09-25 via scontrol/sprio):** multifactor = QOS 20000 (debug gives +10000, normal 0) + FairShare 10000 (our user factor ~0.25 -> ~2536; half-life 14 d) + Age 6000 (max at 60 d) + TRES gpu 30000 x fraction of all GPUs (~100 per GPU) + JobSize 100. NO preemption; sched/backfill (bf_window 14 d, resolution 10 min). Our debug job ranked #1 of 139 pending gpu jobs (12636 vs median 242); waits are resources, not priority. Ask for TIGHT -t so backfill can start jobs in gaps.
- **H100 as overflow (user 09-25): use H100 for CORRECTNESS tests when the H200s are busy** (`--gres=gpu:h100:N`; same Hopper arch, the same HOPPER90 binary runs). TIMING and cost comparisons stay on H200 only. `debug` is a priority level on the same gpu partition, not a separate pool (<= 10 nodes = 40 GPUs per job, 1 running job per user, no preemption).
- **Moving pending jobs (user permission 09-26):** `scontrol update JobId=<id> TresPerNode=gres/gpu:h100:N` moves a pending correctness job from H200 to H100 (it started at once). The classifier blocks it by default; the user granted it for this purpose.
- **BAD NODE hpc-sm-01-09 (confirmed 09-26):** its H200 (UUID GPU-8ea39669...) ran the He box ~6.6x slower (515 vs 77 ms/cycle on hpc-sm-02-10, same binary), and a run there went out of GPU memory. It often shows as "free". Add `--exclude=hpc-sm-01-09` to every GPU job (`scontrol update JobId=<id> ExcNodeList=hpc-sm-01-09` for pending ones); log the hostname in timing results. The same box binary on 02-10 gave 53.6 ms on 09-25 and 77.9 ms on 09-26, so node sharing adds noise: always interleave A/B in one job.
