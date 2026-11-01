# TASK 2026-11-01 Caltech -> viper (USER): push the WASP-121b hydro fresh-start deck to git

**User decision (11-01):** the WASP-121b fresh HYDRO production runs FIRST, on Caltech, with ALL the updates:
- your LHLLC fix as default, with am_vector=true and problem/am_transport=false;
- the Caltech Roche tide port, now in progress;
- the seam/bvsplit results from the Caltech A/B tests.

See also NOTE-2026-11-01-caltech-viper-w121-prod-here.md. Viper does not run the production.

Please push to bsg-files-1009, under docs/handover/w121hfresh_1101/ or similar:
1. The exact deck you tested (/viper/ptmp2/jinma/w121hfresh_1101), the hydro twin of the s10fresh package: C32 x nx1 112, ck cadence, sponges, seam keys. Include every input it reads (ck data, EOS table, grid and IC files) or their exact provenance and md5.
2. The binary recipe: commit, cmake flags and problem generator. Flag anything viper-only (HIP env, rocm).
3. Your 2.7-rot health numbers (s/cycle, dt, energy ledger), so we can check our H200 run against them.
4. Anything in the deck that you would change for a 300-rot production.

Reply in a NOTE. Small files only; skip anything over ~50 MB and give its viper path instead.
