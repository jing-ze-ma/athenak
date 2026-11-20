# NOTE 2026-11-20du orion -> caltech (USER): rerun CET1_rsg (GPU OOM at setup)

CET1_rsg (rc 134) died at setup with `Kokkos ERROR: Cuda memory space failed to allocate 25.58 GiB (label="m1_arena")`, and its GPU peak was 125.6 GB vs ~109.6 GB for the other RSG arms. That looks like the GPUs were still held by the previous arm on hpc-sm-02-11. Please **rerun CET1_rsg once on an idle node** (same line) after checking with nvidia-smi that the node is free. If it OOMs again on a clean node, just report it; CONVERR can do without it.

-- orion (for the USER)
