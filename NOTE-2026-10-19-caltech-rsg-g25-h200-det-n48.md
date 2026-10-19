# NOTE 2026-10-19 caltech -> raven: RSG g25 H200 determinism (N=48 smoke layout): run-to-run BITWISE, restart DIFF after 40 cycles

Job 5008037: 1 node x 4 H200, binary a2713be0 nofma md5 e31ea39e, deck b49f142d + mesh/nx2,nx3=48, meshblock 24x24 (24 blocks, 6/GPU);
check deck adds only output dcycle lines (hst 1, bin/rst 100, log 50); output6 (m1_vet) on as in the deck.
- C, C2: 200 cycles each, 212 s (~1.06 s/cycle incl. start), peak 30.2 GB/GPU, Picard mean 2.96 max 5, NC 0, no FATAL/NaN/DIVERGED.
- C vs C2: BITWISE (hydro.hst, user.hst, all hydro_w / m1 / m1_vet bins).
- R (from C rsg.00001.rst = cycle 100, t 3.3186e4) vs C: hst rows 100-139 BITWISE; FIRST DIFF at cycle 140 (40 cycles after the restart):
  user.hst [25] P_esrc 9.4500287680638820e34 vs ...046e34 (1e-14 rel), then hydro.hst 2-mom at cycle 145; final hydro_w/m1/m1_vet bins DIFF.
  40 cycles = vet_rotate_avg 4 x vet_gd_rotate_every 10: suspect a rotate-average/twin buffer (or its phase) that is not restored bitwise
  from the rst. Run dir /resnick/groups/carnegie_poc/jingze/rsg_g25_prod/det/n48 (C, C2, R kept).
- Your A100/MI300A restart check was bitwise: was its restart point / run length past 40 cycles after the restart? If not, the A100 run may hit the same.
- 5008036 (production layout, 2 nodes) still queued (est start in ~12 h); will report.
