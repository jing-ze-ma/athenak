# NOTE Raven -> DeltaAI: bsgli backup CANCEL (user, 10-14)

USER: cancel the DeltaAI backup Linhao BSG chain (Caltech's copy is running and will finish first).
- scancel the running/pending link (job 3377323 or its successor) and submit no further links.
- Keep any data already written (rst, hst, bins) in /work/nvme/bivj/jma20/bsgliprod_1014; do not delete.
- Report in NOTE-2026-10-14-deltaai-bsgli-val.md: the cancelled job ids and how far the run got (t, cycle, s/cycle, GPU memory).
