# NOTE 2026-11-07 Caltech -> viper: page caches complete on Raven

Finished 15:34 PST: 185,215,647,382 B in 6772 s (27 MB/s), all 16 streams ok. Location: `/raven/ptmp/jinma/caltech_1107/pages/` (paths relative to Caltech $J; mtimes kept).

Verified against Caltech (sum of file sizes + file counts per dir, IDENTICAL):
| dir | bytes | files |
|---|---|---|
| rsgE_page/out/data | 116,190,691,684 | 6305 |
| hgglob_page/out/data | 33,313,536,430 | 1731 |
| vis3d/rsgE_glob/red | 22,378,745,576 | 2948 |
| vis3d/hggl_glob/red | 8,406,412,836 | 974 |

- `raven_copy_1107/pages_cache.tar` md5 e9aa7fe0b095b18971544cf671c4edcf matches; OK to extract.
- Repoint the builders to these dirs per `caltech-to-raven-1107/pages/*/README.md`. Remember the cs3d.py frame-list patch in `common_global/README.md`, so the Caltech-only `red/r_*.npz` numbers stay in the 3-D movies.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
