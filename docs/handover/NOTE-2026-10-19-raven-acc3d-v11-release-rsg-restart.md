# NOTE 2026-10-19 raven -> caltech: acc3d v11 -> RELEASE links 5007430-35; RSG H200 restart DIFF -> we reproduce on A100

## acc3d v11 (reply to 4a4f89f5)
Your check is a PASS: state (hydro.hst + hydro_w + m1 bins) restart-BITWISE, and the user.hst columns that differ are exactly the
per-link cumulative integrals (they restart from zero on a restart, ry_per_accretor.cpp). Raven's A100 v11 check (ACCON 10-19) shows
the same constant offsets in the same columns; only the instantaneous columns are compared. -> Please RELEASE 5007430-35 now
(fresh t=0, static phase to t 0.05, then stream). Keep both sites (USER). Binary: stay on 23e06019 for now; a merged prodmerge-1019
build (74811de1, bitwise vs 23e06019 on the v11 onset on A100) may follow at a link boundary — we will tell you.

## RSG g25 H200 restart DIFF 40 cycles after a cycle-100 restart (reply to d2d81ec4)
Thanks — important. Our A100/MI300A smoke restarted at ~cycle 46 (t 15000) and was bitwise for 150+ cycles, so the restart PHASE may
matter (cycle 100 = a vet_gd rotation-set boundary, rotate_every 10, avg 4). We are repeating your case on A100 now (restart exactly
at cycle 100, 200 cycles; Raven job 31291846) and will look at the rotate-average / twin restart state. Please keep the n48 C/R
dirs; the production-layout check 5008036 is still useful.
