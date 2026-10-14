# TASK raven -> Caltech: FIN4 (merged method candidate) on the IDORT suite — fin4-1014 @ d505d21f

From Raven worker MG, via the coordinator. USER adopted OT's time2_enth_be + implicit_cs2_vimp_npass 3 for the next method
version. FIN4 = FIN3 9972a56d + OT (c685e6dd logic) + OD mbband multi-block fixes + DI he_ic_equil (off). METHOD candidate
test only: production is NOT affected. Small TASK (ENTHBE size, plus 145 ten-cycle BIT pairs).

## 1. Binary
- Branch **`fin4-1014` @ `d505d21f`** on github.com/jing-ze-ma/athenak (pushed). No kokkos change.
- Plain CPU build with EXACTLY the flags of the ENTHBE / mbband builds (gcc 13.2.0 + openmpi 5.0.1, Kokkos OpenMP, MPI on,
  Release, built_in_pgens, no Athena_M1_VGD_NMAX). Report the md5.
- OLD = your `athena_cpu_omp_281dd686` (md5 9ae9f09a...), only for part B.

## 2. KEYS — please read: set NO keys
- Run every deck **AS IS**. Do **NOT** add time2_enth_be, implicit_cs2_vimp_npass, implicit_cs2_vimp or vet_path, and do not
  copy the ENTHBE keys. FIN4's fresh-run DEFAULTS are time2_enth_be = true, implicit_cs2_vimp_npass = 3,
  implicit_cs2_vimp = true, vet_path = none; the point of the test is that the defaults alone reproduce ENTHBE.
- Each hesdirk2 run log prints `time2_enth_be=true implicit_cs2_vimp_npass=3` (cs2 runs) after the time_scheme line; the
  runner records it.

## 3. Run: `bash run_fin4.sh <FIN4> <OLD> 8 AB` (bundle docs/handover/open-1013/caltech_fin4/)
- Roots (edit EB / MB if your trees moved): EB = /resnick/groups/carnegie_poc/jingze/enthbe_1013/caltech_enthbe (the ENTHBE
  tree, it still has <dir>/enthbe with bin + tab), MB = /resnick/groups/carnegie_poc/jingze/mbband_1013/caltech_mbband.
- **Part A** (list_a.txt = the 214 ENTHBE dirs): FIN4 to the end in <EB>/runs/<dir>/fin4, the SAME thread count as
  enthbe/nt.txt (the runner reads it; so cloud and rs n1024 at 32, ddf / ddiff at 4 in an expansion job as before).
  Then cmpout.py enthbe fin4 -> fin4/cmp.txt: BITWISE = every bin (after <par_end>), tab and hst data byte-equal.
  Expected: BITWISE in all 214 (same compiler, same threads; the merged code is the same as c685e6dd on one MeshBlock).
- **Part B** (list_b.txt = the other 145 CVALL dirs): BIT only, OLD vs FIN4, 10 cycles, same threads, cmprst_ra.py (the
  mbband rule) -> <dir>/bit4.txt. Expected BITWISE / RAHIST except about 13 sp Marshak dirs (marsh_r1_T, marsh_r1_u,
  marsh_r1_s, marsh_r1000_u: enth_be acts there at round-off; Raven runs those FULL) and the 2 known p2 FAIL dirs
  (vet_gd too many directions; OLD fails too). No FULL runs in part B.
- Part B first (short), then part A. Finished dirs are skipped on a resubmit (fin4/rc.txt = 0, bit4.txt present).

## 4. Report: docs/handover/NOTE-2026-10-14-caltech-fin4-suite.md
- md5, job ids, every rc != 0 line, counts of fin4_cmpA.txt and fin4_bitB.txt (the runner writes both), and the full
  lines of every part-A DIFF / MISSING and part-B DIFF / FAIL.
- Push fin4_cmpA.txt + fin4_bitB.txt + a light tarball: `tar czf fin4_light.tgz --exclude='*.rst' --exclude='*.bin'
  $(cd $EB; ls -d runs/*/*/*/fin4 runs/*/*/fin4 2>/dev/null)` (run from $EB) under docs/handover/caltech_fin4/.
  The *.bin stay at Caltech unless Raven asks (only needed if part A has DIFF dirs).
- An interim after part B is welcome.

## Coordinator notes
- Priority: normal (after any running production work). The debug QoS rounds as in ENTHBE are fine.
