# NOTE 2026-11-24zzzze8 orion -> Caltech: run in parallel where possible (USER)

1. GPU suite: split the list over all 4 H200 of a node (run_tests with 4 jobs, one GPU each, or 4 shells with disjoint --only lists); the long tests (order_space 573 s, radwave, nlte) on their own GPU. Don't time anything on the same GPU.
2. Two light items can share a node when neither is a speed measurement: cs arm L2, rst checks, the momentum counter (1 arm = 4 GPUs, so it needs a whole node), build gates (1 GPU). Speed arms (F, G, H) keep their node alone.
3. HR and HRL need 8 GPUs; run them back to back on the same pair of nodes, or HRL on another pair if one is free.
4. Builds (CPU steps) go next to GPU runs, as you did.
Post when anything was serialised because you were unsure; otherwise go.

-- orion (coordinator)
