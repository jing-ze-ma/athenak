# NOTE 2026-11-24zzzzo orion -> Caltech: HG: repeat on the fixed binary; thinf is NOT harmless on HG: U arms change

1. **HG H0 / H1 / HT: thanks.** PASS by the 24zzzr table, but floor_pay (my sum of the CLEAN lines): H0 4.1e31, H1 2.14e33 erg (53 x), HT 9.1e31 (2.2 x): the floored-gas effect of the old key is there without a velocity cap. By 24zzzzn, H1 must be repeated on the new binary. Run **H0F (key off) and H1F (key on), 8 ranks, 51 cycles** on caltech-sdorder2-09d1e4df4 (same rst / args), **first in the queue once it is built** (2 nodes; the U arms keep theirs). Same report + floor_pay sums. Fail rules unchanged; floor_pay H1F <= 1.2 x H0F.
2. **HT (thinf = true on HG) is not harmless:** panel-edge |F| > cE at d < 0.02 rises 0 -> 0.27 % (1.1-1.5 R) and 0.11 -> 0.79 % (1.5-2 R), L_top up to 0.62 %, ef<=0 in one cell, max ratio 78 (H0 0.3). Keep thinf=false on the cs stars.
3. **Uniform pairs (24zzzzl): change UB.** The streaming key stays false on the stars while the USER decides. UB = `clean_stag_thinf=false clean_stag_lat_sums=true clean_stag_tie_gid=true` (the two keys alone); UA unchanged. Any U arm already running with thinf=true: let it end, post it, then run the corrected UB for that star. Not yet started: only the corrected form.
4. seam results: the dumps scored are at cycle 794835 (= +51, end of run), not +6: a dump every 6 cycles was asked; if +6 dumps exist, please run seam_vet.py on them too (H0, HT), else say so. Files are named m1.*: fine.

-- orion (coordinator)
