# NOTE 2026-11-19k viper -> Orion (USER): MPI deadlock in the clean cell-layout solve (ranks leave the Picard loop at different passes)

The USER asked viper to tell you. It may affect any clean-cc line that shares CleanSolve / CleanBiCGStab.

## What happened
Run: viper pHG 13255125, step .8, RSG C0 smoke from rsg.01658, 3 nodes x 2 MI300A = 6 ranks.
- Binary he_star_m1_c456f3e3d = r9 245ea4418 + clean_xt_cons + clean_rst_k1=1 (import of the restart's 8 basic slope channels,
  following your Q1 hint). Local branch clean-cc-r9-xtcons-viper, not pushed.
- First cycle 750329, now staged thanks to the slope import:
  - dt 1000 and dt 500 fail;
  - dt 250 converges after 2 cuts (5 cuts without the import);
  - its stage 2 ran 20 Picard passes with de swinging between 6.9e-4 and 0.51 while dv stayed ~1e-8 (a limit cycle in E
    only), and was accepted at the pass cap.
- Cycle 750330 (dt 250) then HUNG: run.log silent for 40 min, memory flat.

Stacks, taken with eu-stack on all ranks while hung. All 6 ranks spin at 95-100 % CPU:
- 1 rank: CleanSolve -> CleanCellMarkFaces -> CleanReduce -> MPI_Allreduce. It has LEFT the Picard loop and is in the
  positivity pass.
- 1 rank: CleanBiCGStab (inner lambda) -> MPI_Allreduce.
- 4 ranks: CleanBiCGStab -> CleanCellApply -> CleanGhostFill -> RecvAndUnpackCC, waiting for halo messages.

So one rank decided the Picard loop was done while the others are still in the linear solve: a collective-order mismatch that
never resolves. The reductions we read (CleanReduce, CleanSums) return the same value on every rank, so the split must come
from a RANK-LOCAL decision. Candidates:
- a non-finite (NaN) residual or de/dv compared differently on different ranks (NaN < tol is false; NaN >= tol is also
  false);
- a rank-local fallback in the halo or the lists;
- a rank-local early exit in the convergence estimate (clean_conv_est) or in the loose/accept logic.

## Next on our side
Rerun with a per-rank print of the Picard exit decision (nit, dpic, de, dv, the converged/accept flags) to find the branch.
Then make every loop-exit decision come from a global reduction (or an explicit MPI_Allreduce of the exit flag).

## For you
- Please check whether your lines can take the Picard exit from a rank-local quantity.
- If you have seen unexplained hangs (a run.log that stops, ranks at 100 % CPU, no FATAL), this is a likely cause.
- The rst_k1 / lin_continue keys are local viper keys. The slope import is your Q1 form (no face channels with
  clean_rst_k1=2).

-- viper
