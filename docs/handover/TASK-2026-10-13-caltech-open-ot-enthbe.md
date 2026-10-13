# TASK raven -> Caltech: open-ot ENTHBE, candidate `time2_enth_be = true` + `implicit_cs2_vimp_npass = 3` on the 214 active IDORT runs

From Raven worker OT (USER GO). METHOD candidate test only; production is NOT affected. Small TASK (same size as NPALL).
OD's mbband TASK has priority if both compete.

## Binary
- Branch `open-ot-1013-ht` @ c685e6dd on github.com/jing-ze-ma/athenak (pushed).
- Equivalently: 281dd686 + `git am patch/0001..0006-*.patch` (no kokkos change).
- Plain CPU build, the same flags as NPALL. Report the md5.

## Run
`bash run_enthbe.sh <binary> 8` (bundle: docs/handover/open-1013/caltech_enthbe/).
- list.txt is the same 214 dirs as NPALL, ordered sp/*_T, lad/*, other sp, cart.
- The 3 long runs come LAST: cart/cloud/FIN3_r12g2, cart/ddf/FIN3_d, cart/ddiff/FIN3_d (about 35 min each at 4 threads).
  - Please run them this time. A longer QoS or the L40S node's CPU cores are fine.
  - If they cannot run within about 3 h, say so; Raven will try them.
- Per dir: deck + the three keys (time2_enth_be, implicit_cs2_vimp_npass = 3, implicit_cs2_vimp = true), run to the
  end into <dir>/enthbe. Finished dirs are skipped on a resubmit.

## Report: docs/handover/NOTE-2026-10-13-caltech-open-ot-enthbe.md
- md5, job ids, rc counts and every rc != 0 line.
- `tar czf enthbe_runs.tgz --exclude='*.rst' runs`, split -b 44M if > 90 MB, plus the md5, under
  docs/handover/caltech_enthbe/.

## Raven OT reply to the 14:50 PDT interim
- Confirmed: all three keys and <dir>/enthbe are what I intended. The bundle script is fixed (2851cfc3).
- Please push the 212 finished runs now as enthbe_runs.tgz (parts + md5), without waiting for ddf / ddiff. Add the two
  dd runs later as a small enthbe_dd.tgz.

## Raven OT: one more small push, please
- Thanks: 214/214 received (light).
- The rlat / rlatt rows are evaluated from the *.bin dumps. Please push ONLY those as enthbe_rlat_bin.tgz (+ md5, split if
  > 44 MB):
  `tar czf enthbe_rlat_bin.tgz runs/sp/rlat*/*/enthbe/bin`
  This covers rlat, rlat_om5, rlatt, rlatt_om5, rlat_T, rlat_om5_T, rlatt_om5_T, rlat_h40 and rlatt_h40.
- Nothing else is needed.
