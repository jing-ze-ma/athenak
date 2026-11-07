# NOTE 2026-11-07 Caltech -> viper: HG rst complete; page caches being sent

- HG `.part` complete: 97,543,425,853 B (= expected; md5 should be 6109201a71e5cad1d123116733218db5). Please md5 + rename.
- USER approved sending the page reductions. Since 13:41 PST, `raven_copy_1107/bulkpush.sh` is streaming (16 tar streams over the 4 open ssh masters):
  - rsgE_page/out/data, hgglob_page/out/data, vis3d/rsgE_glob/red, vis3d/hggl_glob/red
  - raven_copy_1107/pages_cache.tar + .md5
  - total 185.2 GB
  - Destination `/raven/ptmp/jinma/caltech_1107/pages/<same path relative to $J>`; tar keeps mtimes (figcache keys on name+size+mtime).
- ETA ~2.5 h at ~21 MB/s. Caltech posts a NOTE when every stream reports ok.
- Until then, do not use those dirs. Extract pages_cache.tar only after the md5 matches (e9aa7fe0b095b18971544cf671c4edcf).
- The optional PNG frames are not sent.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
