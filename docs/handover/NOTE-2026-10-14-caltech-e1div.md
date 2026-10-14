# NOTE Caltech -> Raven DR: E1 divergence data (TASK e1div-1014)

- `e1div_prof.py /resnick/groups/carnegie_poc/jingze/e1_1013/prod e1div_out 6` on the login node (nice): **wall 259 s**, rc 0, no warnings.
  125 prof_<cycle>.npz (every dump 0 .. 1.24e7 s) + 10 slice_<cycle>.npz (9e6 .. 1.24e7 s). Output 71 MB.
- All section-2 files (hst x2, notadm_first40, notadm_every10, nc_all, cycle_lines_12000plus, floor_tail200, log_tail400, log_head250,
  log_nlines, ls_rst_bin) are in the prof tarball.
- Pushed under docs/handover/e1div-1014/caltech_out/ (two tars):
  - `e1div_out_prof.tgz` (11.1 MB, md5 21cebe618b8ec3f62253c8e93395376e): prof_* + hst + all text files.
  - `e1div_out_slice.tgz` (57.5 MB, md5 068d2802abefd4a9cd48676007f8a96a) as parts: `cat e1div_out_slice.tgz.part_* > e1div_out_slice.tgz`.
  - md5s in `e1div_out.md5`.
- Odd: nothing. Note the script needs the user site-packages (h5py via pip --user): `python3 -I` cannot import it.
- All rst 00000-00024, bins and logs stay on disk at e1_1013/prod.
