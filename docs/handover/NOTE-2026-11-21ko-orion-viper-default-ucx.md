# NOTE 2026-11-21ko orion -> viper: USER decision — production MPI env = DEFAULT UCX (drop XENV)

USER (2026-11-21 ~15:00 PST): switch viper production to default UCX at the next (re)start of each run.
Drop both XENV variables: UCX_RNDV_THRESH=inf and UCX_TLS=^rocm_ipc. Nothing else changes.
Evidence (Orion DETUCX, v20.7 HIP md5 526fc9a2; full write-up stag2/detucx/RESULTS.md on Orion):
- run-to-run bitwise with default UCX on MI300A: t9h 1 node x 4 ranks (10 runs: default / ^ipc / RNDV / XENV all identical) and
  RSG rsg.684575_d1 2 nodes x 2 MI300A (4 runs identical, 0 CW) - every hst row and bin payload.
- XENV originated from a 10-21 Caltech H200 (CUDA) non-determinism and was never reproduced on MI300A; your own NOTE e1cff5197
  (default bitwise over 3x100 cycles + restart) agrees.
- cost vs default: t9h XENV +160 %, ^rocm_ipc alone +76 %, RNDV=inf alone +151 %; RSG XENV +53 %, ^ipc alone +4 %.
Fallback if multi-node determinism ever breaks: UCX_RNDV_SCHEME=put_ppln (Caltech's fix there).
Also found (owner RSTHIP on Orion): a v20.7 HIP t9h restart that is not bitwise (time 1.5e-9 after the rst) and an hgs 6-ranks-on-1-node
startup hang (with XENV too) — not transport-related.

-- orion (for the USER)
