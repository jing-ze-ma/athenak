# TASK 2026-10-08 (Delta): sp-blend2 cost arms, batch 1 (AG Car A and B, blend flux, preconditioner screen)

From: the sp-blend2-1008 cost worker (viper). Machine: NCSA Delta, partition gpuA100x4, account bivj-delta-gpu.
Budget: **at most 5 jobs x 1 node x <= 15 min, about 4.5 GPU-h**. Report the GPU-h used.

## Build
- Branch `sp-blend2-1008` on the fork.
  - Commit: **the commit that adds this file** (`git log -1 --format=%h -- docs/handover/TASK-2026-10-08-delta-spb2-batch1.md`).
  - The code is that of sp-blend-1008 @ 4622f424; only docs were added.
- Build with `bash docs/handover/delta-2026-10-08/build_delta.sh he_gpu <sha>`.
  - That script is on fork/rt-integration. Copy it into place if this branch lacks it.
  - CUDA 12.9 and `-ffp-contract=off` for host code, as the Delta bring-up NOTE says.
- Binary: `$B/bin/athena_he_gpu_<sha>`, where `B=/work/nvme/bivj/jma20/delta_1008`.

## Inputs
- The AG Car files already set up by the bring-up smoke: `$B/agcar_files/agcar_shake{A,B}_ge.athinput`
  (agcar-files-1008, MD5SUMS checked).
- These are the production inputs apart from `implicit_flux`. The job script adds
  `rad_m1/implicit_flux=blend rad_m1/implicit_blend=tau_f` to every arm. The `cen` arms set
  `implicit_flux=central` after that, and the later key wins.

## Run (from the checked-out branch)
```bash
cd <worktree of sp-blend2-1008>/docs/handover/delta-spb2-1008
bash run_batch1.sh /work/nvme/bivj/jma20/delta_1008/bin/athena_he_gpu_<sha>
```
This submits 5 jobs (spb2_delta.sh). Output goes to `$B/spb2/`.
- **prepA**: A blend from t = 0 to cycle 236. Wall limit 13 min. Restart at cycle 236 in `$B/spb2/prepA/rst/`.
  This is the same state as Raven gate2/A_blf rst 00001, where the blend costs +41 %.
- **A1, A2, A3** (afterok prepA): each arm restarts from the prepA restart and runs to cycle 276 (40 cycles,
  wall limit 2:40 per arm). A2 is A1 in reverse order, which gives interleaved repeats.
- **B1**: AG Car B arms from t = 0 to cycle 40 (wall limit 3:00 per arm). Do NOT restart a central B state
  with blend: on Raven that diverged in the first step (Picard resid 1.5e6, FATAL).

## Arm table
| job | arm | keys (after the base) | purpose |
| --- | --- | --- | --- |
| A1, A2 | blf | none (rbgs_fwd, current) | reference blend |
| A1, A2 | mg | implicit_precond=mg (levels 3) | transverse semi-coarsening |
| A1, A2 | rbgs | implicit_precond=rbgs | symmetric sweep |
| A1, A2 | mggc | implicit_precond=mg_gc | global band coarse space |
| A1, A2 | cen | implicit_flux=central | cost of central from the same state |
| A3 | mg5 | precond=mg, implicit_mg_levels=5 | deeper mg |
| A3 | pfloat | implicit_precond_float=true | float line solves |
| A3 | ew01 | implicit_lin_ew_max=0.1 | looser adaptive inner tolerance |
| A3 | mgpf | mg + precond_float | combination |
| A3 | blfface | output7 = bin m1_face at cycle 276 | face fluxes for the L(r) check |
| B1 | blf, mg, mggc, cen | as above | B robustness / cost |

## Report: NOTE-2026-10-08-delta-spb2-batch1.md, pushed to sp-blend2-1008
1. Per arm, `bash metrics.sh <rundir>` (in this directory):
   - s/cycle (elapsed= lines, after skipping 4 cycles);
   - the Picard mean / max / NON-CONVERGED line;
   - inner iterations mean / max;
   - breakdowns and line_jacobi fallbacks;
   - vet_col_lat fallbacks;
   - the |F| > cE line;
   - any FATAL;
   - the last hydro.hst row.
2. A physics deviation table against `blf` of the same job:
   - A: `python3 cmp.py $B/spb2/A.blf.<job> $B/spb2/A.<arm>.<job> ...`
   - B: `CASE=B python3 cmp.py ...`
   - Also cmp between the A1 and A2 `blf` runs, which are the same arm twice, i.e. the run-to-run spread.
3. For blfface: shell-mean L from the face flux, 4 pi r_f^2 <m1_fx1>/L0 at r/R_ph = 1.0, 1.02, 1.05, 1.5,
   2.0, 2.9 (R_ph = 2.70140777e13 cm, L0 = 3.41170859e39), next to the cell-flux L that cmp.py prints.
4. The prepA metrics (cycles 0-236) and the GPU-h used.

Keep run trees in `$B/spb2/`, and push only the NOTE. If prepA fails, stop and report; the arms depend on it.
