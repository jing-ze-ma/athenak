# NOTE 2026-10-08 (Delta -> spb2 worker): batch 1 STOPPED -- prepA FATAL at startup (job-script key), batch 2 blocked

Answer to TASK-2026-10-08-delta-spb2-batch1. Nothing ran past startup. 0.02 GPU-h used (ledger:
`docs/handover/delta-2026-10-08/LEDGER.md`).

## What happened
- Binaries, both built and ready:
  - batch 1: `athena_he_gpu_nofma_f5b26846`, md5 988bb3d6d3234536eab59729d5e031b2;
  - batch 2: `athena_he_gpu_nofma_6d913aed`, md5 a95d58ee13efaa4c24833de676846fb5.
  - Both: CUDA 12.9, host `-ffp-contract=off`. The target name in build_delta.sh is `he_gpu_nofma`.
- Scripts were run unchanged from `docs/handover/delta-spb2-1008` @ f5b26846:
  `run_batch1.sh` submitted prepA 22760753; A1-A3 22760754-56 (afterok); B1 22760757.
- **prepA 22760753 died 5 s after srun, on all 4 ranks:**
  ```
  ### FATAL ERROR in src/parameter_input.cpp at line 430
  Parameter 'dcycle' in block 'output1' on command line not found
  ```
  - `spb2_delta.sh` BASE sets `output1/dcycle=1 output2/dcycle=$NLIM output3/dcycle=$NLIM`.
  - In `agcar_shake{A,B}_ge.athinput` (agcar-files-1008), output1 (hst), output2 and output3 (bin) are given by
    `dt` only, and AthenaK refuses to override a key that is not in the input.
  - The same BASE is used by every arm, so B1 and all batch-2 arms would fail the same way.
- **The job script exits 0 when the run fails.** The last command is `echo`/`ls`, so Slurm marked prepA
  COMPLETED and released the afterok arms. I cancelled A1-A3 and B1 (22760754-57) before they started, as the
  task says ("If prepA fails, stop and report").
- Batch 2 (`run_batch2.sh`) cannot be submitted: it needs `$B/spb2/prepA/rst`.

## Needed from the worker (Delta does not edit the batch scripts)
1. Make the per-cycle output work with these inputs. Two ways:
   - put `dcycle` keys in a Delta copy of the inputs, e.g. replace `dt` by `dcycle` in output1-3 (dt and dcycle
     both present is legal);
   - or have the script `sed` them into a copy of `$IN` per job.
   - Also check that `output6/file_type=rst` (a block that is not in the input) is accepted from the command
     line; the run never got that far.
2. End `spb2_delta.sh` with the srun exit code, e.g. keep `rc=$?` and `exit $rc` after the loop, or exit at the
   first failing arm. Otherwise afterok does not protect the dependent arms.
3. Push the fix as batch1 (re-issued) or a new batch file. The poller picks up new
   `TASK-2026-10-08-delta-*-batch*.md` files within 15 min; for a re-issued batch1, name it e.g.
   `TASK-2026-10-08-delta-spb2-batch1b.md`.
   - Batch 2 then follows automatically, with binary 6d913aed already built.
