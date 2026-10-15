# TASK for the Delta session: standing GPU test runner for the radiation work (user 10-08)

User decision: Delta (gpuA100x4, account bivj-delta-gpu, ~99 GPU-h left) runs the short GPU tests of viper's two
top-priority radiation workers, because viper apudev / Raven gpudev are jammed and Delta backfills short jobs fast.
No production on Delta.

## Loop
Every ~15 min (`/loop 15m ...` or a background poll), `git fetch origin` and look for NEW task files named
`docs/handover/TASK-2026-10-08-delta-*-batch*.md` on these branches:
- `sp-blend2-1008` (cost of the sp blend implicit flux: preconditioner / inner solve; worker "spb2")
- `rad-beam-1008` (accuracy for beams in the thin regime: plm_dc, weight window, vet_gd formal solution; worker "beam")
- `cs-floor-diag-1008` (cubed-sphere floor / FOFC counters on WASP-121b; worker "csfloor"; lower priority: run its
  batches only when no spb2/beam batch is waiting; budget <= 10 GPU-h)
A batch file names the sha to build, the inputs (committed next to it or in `agcar-files-1008` on branch
`he-ic-eint-from-t`), the arm table, the metrics and the NOTE file to push back (on the SAME branch).
For each new batch: build that sha incrementally with `docs/handover/delta-2026-10-08/build_delta.sh` (one build dir per
branch; CUDA 12.9; host `-ffp-contract=off`), run the arms (1 node, 4 A100, <= 15 min each unless the batch says
otherwise), extract the metrics, push the NOTE on the branch, then continue polling. Process batches in order of
arrival; one batch at a time per branch.

## Budget and rules
- GPU-hour budget: spb2 <= 40, beam <= 40, csfloor <= 10, keep >= 5 in reserve. Log every job (id, branch, arm, GPU-h) in a ledger
  `docs/handover/delta-2026-10-08/LEDGER.md` pushed with each NOTE. Stop and ask the user when a branch reaches its budget.
- Work dir `/work/nvme/bivj/jma20/delta_1008/` (shared filesystem with DeltaAI: stay inside delta_1008/).
- Never push rt-integration except docs/handover notes addressed to viper; never force-push; never touch DeltaAI's
  directories or jobs.
- Report anomalies (build failure, FATAL, NaN, anything that looks like a code bug) in the NOTE; do not fix code on Delta.

## Update (viper, user 10-08 later)
Batch names now also match `TASK-2026-10-08-delta-*-batch*.md` with suffixes like `batch1b` (re-issues after a
failed batch). spb2 batch1b/2b (b6240ae4) are ready; beam and csfloor batches will appear on their branches.
