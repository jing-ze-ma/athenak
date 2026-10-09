# NOTE 2026-10-09 viper: He giant page script make_figs.py pushed (answers TASK-2026-10-08-viper-hegiant-make-figs)

From: viper. To: Caltech. Docs only; no run touched.

## Files: `docs/handover/hegiant-files-1007/page/` (md5 in `page/MD5SUMS`)
- `make_figs.py` (e6bf9e86): the viper page script (he_giant_1006/page/make_figs.py, md5 4cd82316), with three changes.
  - The hard-coded viper paths are replaced by arguments and environment variables.
  - It can now join several run directories.
  - `read_ross` is copied in verbatim from `he_giant_1006/ic/make_ic_mlt_star_ge.py`, so no IC module is imported.
  - Its only non-standard import is `bin_convert` from the AthenaK repo's `vis/python`. The default path is this checkout's
    `vis/python`; override it with `--vis` or `ATHENAK_VIS`. Other imports: numpy, scipy, matplotlib.
- Reference files:
  - Opacity: the script defaults to `../rosseland_tops_hegiant_blend.txt` (already on the branch, c1ef063c). This is the
    run's Rosseland table, used for R_ph = shell-mean tau_R 2/3.
  - MESA profile used for the figures: `M2pt754_Porb100_profile29_cols.data` (31c96bf5). It keeps the MESA format with
    the columns zone, mass, logR, logT, logRho and velocity. This is all make_figs.py reads.
  - Full MESA profile, the input of PROFILE_CHECK: `M2pt754_Porb100_profile29.data.gz` (20458e76; uncompressed
    f0cf56f1). The viper `PROFILE_CHECK.md` is copied as `PROFILE_CHECK_viper_1006.md`.
  - The MESA expansion line (61.8 Rsun + 6.6 km/s t) is hard-coded from that profile header, as before.
- Test output: `test_joined_scout128_smoke897b_numbers.json` and `test_joined_fig_profiles.png`.

## How to run
```
python3 make_figs.py RUNDIR [RUNDIR ...] --out OUTDIR [--table ROSS] [--profile MESA] [--label "He giant N897"] [--vis .../vis/python]
# Caltech, N445 then N897 (time order):
python3 make_figs.py /resnick/groups/carnegie_poc/jingze/hegiant_1007/n445_1g /resnick/groups/carnegie_poc/jingze/hegiant_1007/n897 \
    --out <page dir> --label "He giant N445 -> N897"
```
- Each RUNDIR needs:
  - `bin/*.hydro_w.*.bin` with the matching `bin/*.m1.*.bin`;
  - `hegiant.hydro.hst` and `hegiant.user.hst`;
  - optionally `run*.log` (link stats) and `jobs.txt` (squeue status; skipped if absent).
- Dumps: every dump of every run. The shell means are cached in `OUTDIR/shellmeans_v2.npz`, keyed by absolute path,
  so a rerun only reads new dumps. The slice and shell figures use the latest dump.
- Outputs: fig_rph, fig_hst, fig_profiles, fig_slices and fig_shell (png), plus numbers.json. The script never writes
  into a run dir.
- Environment variables can replace the defaults: HEGIANT_PAGE_OUT, HEGIANT_ROSS_TABLE, HEGIANT_MESA_PROFILE, ATHENAK_VIS.

## Radial-grid change (N445 -> N897 at 10.42 d)
- The old script could not handle the grid change: it took one run dir and cached a single radial grid.
- It now joins runs given in time order. Run i is used only for t < t0 of run i+1, where t0 is the earlier of that
  run's first hst row and first dump.
  - Dump times come from the .bin header, which has 6 significant digits. A dump of the earlier run counts as "before"
    if it is lower by more than half a unit of the 6th digit.
  - If the later run wrote a dump at the restart time, that dump replaces the earlier run's dump.
- Each dump keeps its own r grid. Profile curves are labelled N445/N897, and a profile at each join is added.
- History and log stats are concatenated. A dotted line marks the join in fig_rph and fig_hst.
- numbers.json gains `segments`, `joins_t_d` and a per-dump `nr`.

## Tested on viper (login node, nice)
- Single run, scout128 (N445, 44 dumps, 0-10.422 d):
  - The `Rph_table` and `max_abs_mean_vr_kms` in numbers.json are identical to the old script's page numbers.json
    (all 44 rows).
  - Latest R_ph 73.466 Rsun, as before.
- Joined, scout128 + smoke897b:
  - smoke897b is a viper N897 smoke from the remapped 10.42 d restart, hst from t 900430 s, with 1 dump at 10.428 d.
  - Result: segment 0 has 44 dumps up to 10.4216 d. Its last dump (t 900420.8 s, the remap state) is kept, because
    smoke897b wrote no dump at the restart time. Segment 1 has 1 dump.
  - The latest dump is N897 at 10.428 d, R_ph 73.258 Rsun. All five figures were written.
  - In `test_joined_fig_profiles.png` the N445 profile at 10.42 d and the N897 profile at 10.43 d overlap.
- There is no long N897 run on viper, so a multi-dump N897 leg was not tested.
- Not tested: a run dir without `run*.log` (the code path is guarded) and Caltech's directory layout.
