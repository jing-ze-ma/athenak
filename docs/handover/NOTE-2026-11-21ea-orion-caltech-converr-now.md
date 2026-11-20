# NOTE 2026-11-21ea orion -> caltech (USER): run converr.sh NOW on the existing CONVERR arms

Please run `converr.sh` **now**, without waiting for all 14 arms. The main numbers need only CEP / CET3 / CET2 / CET1, and the script skips missing runs.
- It is CPU-only (1 task, -c 8, ~40 GB RAM, no GPU), so it fits next to the FRESHGO arms in any holder (srun --overlap).
- Commit `converr_{rsg,hg}.txt` to `docs/handover/caltech_results/CONVERR/` as soon as they exist. **Keep the bins** for now.
- The remaining CONVERR arms (the CET1_rsg rerun, CEP0 x2, the 3c x4) still wait behind FRESHGO. When they have run, rerun converr.sh and commit again; delete the bins only after that final commit.

-- orion (for the USER)
