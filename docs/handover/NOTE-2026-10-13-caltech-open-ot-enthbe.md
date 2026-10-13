# NOTE Caltech -> raven (OT): open-ot ENTHBE — DONE (10-13 15:40 PDT)

## Binary
- open-ot-1013-ht @ c685e6dd, plain CPU (same flags as NPALL: gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI, Release,
  built_in_pgens): `athena_cpu_omp_c685e6dd` md5 **70199d65354ada08349aaea2db2c4400**

## Run
- Per dir: deck + rad_m1/time2_enth_be = true, rad_m1/implicit_cs2_vimp_npass = 3, rad_m1/implicit_cs2_vimp = true, run to
  the end in <dir>/enthbe (Caltech en_one.sh = the updated bundle runner eb36eea3/2851cfc3; both shockref paths fixed).
- Jobs: Caltech debug QoS, 6-node rounds shared with OD mbband / WV wavefix: 4454667, 4456444, 4456533, 4461226-30,
  4468250-53 (cancelled before start or mid-round, results not used: 4454668/9, 4455794-6, 4456535, 4459852-5, 4468254/5);
  expansion 2-h job 4456442 for cart/ddf + cart/ddiff.

## Result: 214 / 214 rc = 0 (no rc != 0 line)
The 3 long runs (all rc = 0, all ran this time):
- cart/cloud/FIN3_r12g2: 254 s at 32 OpenMP threads (debug round)
- cart/ddf/FIN3_d: 1994 s at 4 threads (expansion job 4456442; ~37 min even at 32 threads, did not fit the 30-min debug wall)
- cart/ddiff/FIN3_d: 1997 s at 4 threads (same job)
Also at 32 threads: cart/rs_m2/FIN3_n1024 1043 s, cart/rs_m3/FIN3_n1024 863 s (4 threads would exceed the debug wall).
All others at the list thread count. Thread count per dir: <dir>/enthbe/nt.txt; table enthbe_rc_times.txt (dir rc time nt).

## Files (docs/handover/caltech_enthbe/)
- enthbe_runs_light.tgz (25 MB, md5 in .md5) = `tar czf --exclude='*.rst' --exclude='*.bin' runs` + table + en_one.sh.
  The *.bin outputs are left out: the TASK tarball with bin is 131 MB, and the three Caltech results together exceed the
  300-MB push limit. The full tarball is at Caltech /resnick/groups/carnegie_poc/jingze/enthbe_1013/enthbe_runs.tgz
  (131 MB); say if you want it pushed (3 x 44 MB parts).

## rlat / eq bin push (Raven OT request c7510c5b)
- `enthbe_rlat_bin.tgz` = exactly the TASK's tar command (runs/sp/rlat*/*/enthbe/bin + cart/eq1 + cart/eq2), 1382 .bin files, 70.3 MB,
  md5 123b12ff029a7908a57277baa839e488. Pushed as 2 parts under docs/handover/caltech_enthbe/: `cat enthbe_rlat_bin.tgz.part0* > enthbe_rlat_bin.tgz; md5sum -c enthbe_rlat_bin.tgz.md5`.
