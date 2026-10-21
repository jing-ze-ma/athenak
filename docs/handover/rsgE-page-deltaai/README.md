# RSG E results-page pipeline (DeltaAI copy) -> handed over to Caltech (USER, 2026-10-21)

Builds the "RSG E DeltaAI" page (https://claude.ai/artifact/JegJZ1mkx7nqHqz3pk2aU9, v6, 0-24 d): time series, radial
profiles, native cubed-sphere slices, T_gas, Mollweide maps + 4 movies (eq, inner, phot, map). DeltaAI run is CANCELLED
(NOTE 8c925750); this is so Caltech can build the same page from its own RSG E run.

Origin: src/ = Raven rsg-page (docs/handover/rsg-page, c542d6b0) + DeltaAI changes: native-cell slices (csplot.py), T_gas from
the code EOS table (eos_gen.py; reproduces Raven tgas.py check values to <1e-4), movies (mk_movies.py), ':'-lists for BINDIR /
EXTRA_RUN / LOGDIR, notes.html. src/upd_page.sh is Raven's original driver (no movies), kept for reference.

## Run it at another site
    ./update.sh          # incremental: only new dumps reduced (~45-60 s, ~3.5 GB RSS each), 4 movies (~3 s per dump)
Output: $OUT/publish/index.html (figures inlined) + $OUT/publish/movies/rsgE_{eq,inner,phot,map}.mp4 -> publish the whole dir.
update.sh has DeltaAI defaults; override by env (all of them are read):
  RUN       dir with rsg.user.hst + run.*.log          BINDIR   ':'-list of dirs with rsg.{hydro_w,m1}.NNNNN.bin
  EXTRA_RUN ':'-list searched for run.*.log            LOGDIR   ':'-list with link.*.out (chain logs)
  OUT       work dir (data/ img/ publish/)             ICF      IC column: docs/handover/rsg-E-1020/ic_rsg_g25_deep.txt
  VISPY     AthenaK vis/python                         PY       python with numpy + matplotlib (+ scipy)
  FFMPEG    ffmpeg binary (default via imageio_ffmpeg) RPH0     5.80019e13 (R_ph of the IC)
  JOBNAME   squeue name used by autoupd.sh             SITE, LAYOUT  text shown on the page
Also: SRC is set to $HERE/src in update.sh - change HERE (line 9) to where you put this dir.
Hard-coded DeltaAI text to edit for your site: src/build_page.py lines ~73-74 (deck md5, start time, output cadence), ~108-117
(title "RSG E DeltaAI", intro), ~151 (footer path); src/mk_movies.py suptitles "(DeltaAI)".
pub.sh  = DeltaAI staging helper for Claude publishing (copies to a dir under Claude's cwd, prints only the text diff); edit H/S.
autoupd.sh = 45-min rebuild loop while the job runs (touch STOP_autoupd to stop); edit H, the squeue user/name and the bin path.
Data reference: DeltaAI run (0-24.5 d, all 1 d dumps) stays at /projects/bivj/jma20/rsg_E_1021/run.
