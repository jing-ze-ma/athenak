# TASK Caltech: data for the AG Car E1 divergence (DR, e1div-1014) -- analysis only, NO restart, NO new production

From the Raven coordinator (worker DR). E1 link 1 (job 4452279) stopped with Picard DIVERGED at t 1.24255e7 s, cycle 17085
(NOTE-2026-10-13-caltech-e1.md, "LINK 1 STOPPED"). Do NOT restart E1, and do not delete rst/bins. We need profiles and log excerpts.
Worst cells (NOT ADMISSIBLE lines, ghost-inclusive i 347-350 and 478-511) map to r = 1.025-1.029 R_ph (the photosphere/floor junction) and
1.35-1.9 R_ph (the floor layer above wb_rmax = 1.1 R_ph).

## 1. Radial profiles + slices from the bins (python3 + numpy; bin_convert.py from ../e1caltech-1013, which needs h5py -- installed)
```
cd <repo>/docs/handover/e1div-1014
python3 e1div_prof.py <abs path>/e1_1013/prod <abs path>/e1div_out 8      # small CPU job or nice on login; reads every bin pair
```
- Output: prof_<cycle>.npz (~0.12 MB per dump, all dumps 0 .. 1.24e7 s) and slice_<cycle>.npz (~10 MB each, about 10 dumps between 9e6 and 1.24e7 s).
- If memory is tight use NPROC 2 (each worker holds one 536x128x128 dump in float64, ~0.6 GB).

## 2. Files and log excerpts (L = e1_1013/prod/run.log.4452279)
```
O=<abs path>/e1div_out
cp e1_1013/prod/agcar3d.hydro.hst e1_1013/prod/agcar3d.user.hst $O/
grep -n "NOT ADMISSIBLE" $L | head -40 > $O/notadm_first40.txt
grep -n "NOT ADMISSIBLE" $L | awk 'NR%10==0' > $O/notadm_every10.txt
grep -n -i "NON-CONVERGED\|NC-local\|DIVERGED\|FATAL" $L > $O/nc_all.txt
grep -n "cycle=" $L | awk -F'cycle=' '{split($2,a," "); if (a[1]+0 >= 12000) print}' > $O/cycle_lines_12000plus.txt
grep -n -i "fofc\|floor\|clamp\|vceil\|tset" $L | tail -200 > $O/floor_tail200.txt
tail -400 $L > $O/log_tail400.txt
head -250 $L > $O/log_head250.txt
grep -c "" $L > $O/log_nlines.txt; ls -l e1_1013/prod/rst e1_1013/prod/bin | tail -80 > $O/ls_rst_bin.txt
```

## 3. Push back
- `tar czf e1div_out.tgz -C <parent> e1div_out`; if > 45 MB: `split -b 40m e1div_out.tgz e1div_out.tgz.part_` (keep the prof_* + text in the first part
  if you prefer two tars: one with prof_*/hst/txt, one with slice_*).
- Commit to `docs/handover/e1div-1014/caltech_out/` on bsg-files-1009, push, and add a short NOTE (`NOTE-2026-10-14-caltech-e1div.md`): md5 of the parts,
  wall time of the profile run, anything odd.
- Keep everything on disk (rst 00000-00024, bins, logs) until told otherwise.
