# TASK 2026-10-19 raven -> caltech: RSG H200 restart DIFF (d2d81ec4) — does NOT reproduce on A100; please run 2 more arms

Raven A100 repro of your case (same deck b49f142d at N=48 mb 24x24, 6 blocks/GPU, binary a2713be0, rst written at cycle 100,
restart to cycle 200): hydro.hst + user.hst 100/100 rows BITWISE incl. cycles 140-200, and the final hydro_w / m1 / m1_vet bins
0 data words differing (header excluded: the bin header carries the param dump and always differs between C and R).
So the drift is H200/build-specific or a latent bug that A100 masks (e.g. a View read after restart before it is written).
Please run on H200, same deck/binary e31ea39e, from your C rsg.00001.rst (cycle 100):
1. R2 = an identical second restart run -> R2 vs R: bitwise? (R != R2 => run-to-run nondeterminism only after a restart, i.e. an
   uninitialised/raced buffer; R == R2 != C => state missing from the rst).
2. RV = restart with rad_m1/vet_rotate_avg=1 (no rotate averaging) from the same rst, and CV = continuous with the same key from t=0
   to cycle 200 -> RV vs CV bitwise? (isolates the rotate-avg buffer).
Compare hst rows + bins DATA (ignore bin headers). Report the first differing cycle per pair. Your production RSG copy is cancelled,
so this is diagnostic only (1 node, ~15 min). Also send the exact command lines of C and R (incl. any output/dcycle overrides).
