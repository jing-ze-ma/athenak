# NOTE 2026-11-08 Caltech -> viper/Raven: RSG E + HG gap-dump reductions sent (reply to TASK-2026-11-08-caltech-reduce-gap-dumps)

Done. All requested dumps existed on Caltech; **no dumps missing**.

Reduced in one CPU batch job (Caltech job 6463496, expansion, 32 cores, 2 min 54 s) with the unchanged page pipelines:
- 2-D: `rsgE_page/src/reduce.py` / `hgglob_page/src/reduce.py` (reduction VER 6), called exactly as `catchup_job.sh` does
  (`ONLY=NNNNN RED_THREADS=.. $PY $SRC/reduce.py $RUN data`, env from each page's `update.sh`). RSG via the pageH link view
  (`pageH/sync.sh` run first; 06310-06325 were not linked yet; the links point at prodH/bin). HG from hgglob_1021/prodG/bin.
- 3-D: `vis3d/cs3d.py` `load_reduced` (REDUCE_VERSION c3, same cache key as `dump_list`), reductions only (no clim / frames / encode).
  r_06304 was already cached (and is already on Raven), so not re-sent.

| set | files | Raven destination | md5sum -c on Raven |
|---|---|---|---|
| RSG E red | red_06304..06325 (22, 387 MB) | /raven/ptmp/jinma/pages_1107/cache/rsgE_page/out/data | 22/22 OK |
| RSG E vis3d | r_06305..06325 (21, 153 MB) | /raven/ptmp/jinma/pages_1107/cache/vis3d/rsgE_glob/red | 21/21 OK |
| HG red | red_01730..01740 (11, 202 MB) | /raven/ptmp/jinma/pages_1107/cache/hgglob_page/out/data | 11/11 OK |
| HG vis3d | r_01730..01740 (11, 91 MB) | /raven/ptmp/jinma/pages_1107/cache/vis3d/hggl_glob/red | 11/11 OK |

`MD5SUMS_gap_1108.txt` (computed on Caltech) sits in each destination directory. Times: RSG 6304.0-6325.0 d (06325 at
t = 6325.000 d = the restart you continued from); HG 576.667-580.000 d (01740 at 580.000 d). Raw dumps not sent.
Caltech copies: the same files in $J/rsgE_page/out/data, $J/hgglob_page/out/data, $J/vis3d/{rsgE_glob,hggl_glob}/red; job
script + logs in /resnick/groups/carnegie_poc/jingze/gap_1108.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
