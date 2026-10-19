# TASK 2026-10-19 raven -> caltech: cs EDGE INSTABILITY arms (CSEDGE), RSG g25 at N=48 on H200 — TOP PRIORITY (RSG prod + global HG stopped/held)

The RSG g25 cs production grew a numerical mode in the single row of cells next to every cubed-sphere panel seam (along-edge
tangential velocity, a one-cell-wide jet toward the cube vertices), at r < 0.045 R_ph, onset ~3.7 d, e-fold 0.5 d at N=128,
KE_int 1e41 -> 1e45. Raven GPU queues start tomorrow; please run these discriminating arms now (diagnostic only, not production).

Binary: your existing H200 build of a2713be0 (e31ea39e) is fine (the old-seed run on a2713be0 showed the same mode).
Deck: docs/handover/rsg-g25-prod/rsg_g25_prod.athinput (b49f142d, old seed) with your local IC/table paths, as for the earlier
H200 determinism task. Mesh keys for every arm: mesh/nx2=48 mesh/nx3=48 meshblock/nx2=24 meshblock/nx3=24 (24 blocks, 1 H200).
Common overrides: time/tlim=1.04e6 time/ndiag=100 output1/dt=2000 output2/dt=43200 output3/dt=43200 output4/dt=1e30
output5/dt=1e30 output6/dt=1e30   (12 d, hydro_w + m1 bins every 0.5 d, no rst).
Arms (1 H200 each, run concurrently):
  B  base (no extra key)
  F  problem/mlt_closure=frozen
  H  hydro/rsolver=hlle
  N  rad_m1/gas_feedback=false          (gas sees gravity - a_ref via the WB potential only; radiation decoupled)
  V  rad_m1/vet_gd=false                (if it refuses to start, report the message and skip)
  B24  base at N=24: mesh/nx2=24 mesh/nx3=24 meshblock/nx2=24 meshblock/nx3=24
Send back (commit to bsg-files-1009, docs/handover/csedge-1019/caltech/<arm>/): rsg.user.hst, rsg.hydro.hst, last 30 lines of run.log,
the exact command line, s/cycle, and the output of
  python3 docs/handover/csedge-1019/csedge_diag.py <athenak>/vis/python <run>/bin rsg  > diag.txt
(a few lines per dump; row0 >> row1 ~ interior = the mode). If an arm is still running at your next check-in, push partial results
(the diag on the dumps so far) — early is more useful than complete. Raven runs the same arms when its queue frees.
