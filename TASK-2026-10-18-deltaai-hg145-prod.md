# TASK raven -> DeltaAI (2026-10-18): He giant 14p5M PRODUCTION chain on ghx4 (USER GO: "yes try ghx4")

USER explicitly approved a DeltaAI production copy of the He giant 14p5M (overrides the 10-14 "DeltaAI tests only" rule for THIS run).
A Raven copy (4 A100) runs in parallel; keep both unless the USER says otherwise. Reference: Raven /raven/ptmp/jinma/hg145_prod_1018/NOTE.md.

## Setup = your Test A (93246a60), which passed (clean, R/C bitwise)
- Binary: your athena_hes_gpu_4c494982ffb3 (prodmerge-1018-hg @ 4c494982, md5 385860369a7d3d3801df2b4649c0ced7). Do not rebuild.
- Deck: docs/handover/hg145-prod-1018/hg145_prod.athinput.in (d70305a7) with @FILES@ substituted; change no key. Files as Test A.
- Layout: 1 node x 4 GH200, 4 ranks, 1 MeshBlock 704x64x64 per GPU, MPICH_GPU_SUPPORT_ENABLED=1 (as Test A).

## Chain
- Partition **ghx4** (batch), fresh t = 0 in a new run dir (e.g. /work/nvme/bivj/jma20/hg145_prod_1018/run), NOT the Test A dir.
- **tlim 1e9** (long-tlim rule: the run is limited by the links), links = the MAX walltime ghx4 allows; athena `-t` = walltime minus a
  ~15 min margin so the final rst is written. Link 1 + 2 afterany links; each link: FRESH if no rst, else newest valid rst
  (`-r <rst> time/tlim=1e9`).
- Outputs as in the deck (do not reduce). rst prune: keep the newest 3.
- **STOP rule:** stop ONLY on the end-of-run summary `NON-CONVERGED=[1-9]`, FATAL, "Picard solve DIVERGED", NaN in the tail, rc != 0.
  NEVER on per-pass "NON-CONVERGED after N passes" lines or DT_RETRY reason lines.
- Storage: report the dump size per day of sim time after link 1; if /work/nvme quota is a problem, say so (do not delete dumps).

## Report: NOTE-2026-10-18-deltaai-hg145-prod.md
- Job ids + queue ETA immediately; then after link 1 starts: s/cycle, dt, Picard, memory, t reached per link. Post a one-line NOTE at
  each link end (t, cycle, rc, STOP-check result).
