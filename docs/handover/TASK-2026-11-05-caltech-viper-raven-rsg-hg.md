# TASK 2026-11-05 Caltech -> viper (USER go): prepare RSG E and HG on Raven for the 11-07..11-10 Caltech gap

USER: "yes tell viper to prepare rsg and hg on raven". Prepare only; the restart copy is user-driven (~11-07), and the launch
on Raven waits for the rst to arrive + a passing smoke. Reply with NOTEs on this branch.

## What runs on Caltech now
| run | binary (Caltech) | blocks / ranks | GPU mem in use | rst size | Caltech wall ends |
|---|---|---|---|---|---|
| RSG E (rsg_E_1020/prodH) | athena_gpu_he_star_m1_prod_23a55658_fmadoff (md5 38231ce9...) | 24 blocks (C128 cs, 6 x 2x2), 8 ranks | ~57 GB/GPU at 3 blocks/GPU -> ~19 GB/block | ~44 GB | ~11-07 05:30 PST |
| HG (hgglob_1021/prodG) | same binary | 24 blocks, 12 ranks | ~49 GB/GPU at 2 blocks/GPU -> ~25 GB/block | ~94 GB | ~11-07 05:30 PST |

(nvidia-smi "memory.used", not a high-water mark; it includes Kokkos pools. 1 block/A100-40GB should fit, but verify in the
first steps.) The exact last restart command lines (binary + every command-line key, incl. last_time, the M1/VET keys,
cs_seam_rho_guard etc.) are in `docs/handover/raven-gap-1104/rsg_E_last_cmd.txt` and `hg_last_cmd.txt`.
Both runs restart from their own rst input + those keys.

## Please
1. Build on Raven the SAME source as the Caltech production binary: prod-1102 commit **23a55658** (he_star_m1, CUDA, fmad off,
   like Caltech's `_fmadoff` builds). Do NOT use a newer prod commit for these continuations (8cd7f0a9 changes defaults that
   matter only on fresh starts, but keep the continuation exact). If 23a55658 needs the A100 team-scratch handling
   (de528da5), say so; prefer a build where only that differs.
2. Check /raven/ptmp/jinma quota/free space for ~140 GB (+ outputs for ~4 days of both runs; RSG dumps ~?, HG ~?) and
   tell the user (and us) where to put the files.
3. Prepare the two run dirs and launchers on Raven: 24 ranks x 1 block on 24 A100 = 6 nodes each, 24 h links chained with
   dependency, gpu env as usual, OMPI/UCX settings that work on Raven, stop rules (FATAL / summary NON-CONVERGED=[1-9] /
   nan -> stop). Use the Caltech command lines verbatim except paths, the binary, and rank/GPU mapping.
4. Smoke plan: when the rst arrives, 20-50 cycles with the exact binary + keys; compare the first-cycle hst row with the
   Caltech hst row at the same time (expect agreement to ~1e-12..1e-10 relative, CUDA H200 vs A100); then start the chain.
5. Reverse plan: when our Caltech holders start (est. 11-10), stop at a fresh rst; the user copies it back.
Reply when 1-3 are ready (build md5, dirs, the copy destination).
