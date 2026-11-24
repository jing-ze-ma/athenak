# NOTE 2026-11-23zm orion -> Caltech: USER split of N2 (7333073): node A m35mms arm D0 now; node B reserved for a t9l62 test

The USER split N2 4 + 4.

## Node A (first N2 node, 4 H200): m35mms arm D0, start now
- Goal (USER): m35mms to 10 d in ≤ 24 h on one node without loss of accuracy. D0 is the reference: IA12.1 from your real late state.
- `env ATHENA_DT_REPORT=25 <production env> $BIN -r prod_ia11_1123/m35mms/rst/m35mms.00010.rst -d <new dir> $(cat D0_args.txt)`
  - BIN = your athena_ia121_fa59bffd6.
  - D0_args.txt is in `docs/handover/caltech_m35fast/` (md5 in MD5SUMS). It names every hydro / rad_m1 / time key and the 3 mesh layout keys.
  - It runs 2000 cycles (nlim 23945), writes no rst and keeps the production bins.
- If the restart FATALs on `mesh/layout_invariant_mode=tile` (the IA11 rst was written without tile), drop the three `mesh/` keys and say so.
- Post: user + hydro hst; `grep -E "DTREP|DTREG|CLEAN WARNING|CLEAN FALLBACK|FATAL|^CLEAN cycle|elapsed="` of run.log; s/cycle.
- The production m35mms on IA11 stays untouched.
- Later arms D1 / D1h / D2 need a new branch; we'll push it with a NOTE (~21:00-22:00 PST).

## Node B (second N2 node): keep free for a t9l62 restart-transient test arm (NOTE to follow).

-- orion (coordinator)
