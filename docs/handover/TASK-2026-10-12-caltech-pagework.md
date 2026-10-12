# TASK DeltaAI -> Caltech: push the BSG true-repro page pipeline (bsg_1009/pagework), scripts only

From DeltaAI (user request 10-12): the "BSG True Reproduction" page (https://claude.ai/artifact/BQ8LZRU8uV1Zw83GtZQQxK)
is to show the new bsgprod run (DeltaAI job chain from 3367085, final method). The page footer says it was regenerated
with `bsg_1009/pagework/make_page.sh`, which is not on GitHub; only viper's base copy is (docs/handover/bsg-hrdet-1009/page).

Please push the pagework SCRIPTS (make_page.sh, the pipe/ copies of bsg_figs.py, gen_status.py, cost_repro.py, ke_compare.py,
compare57.py, cmp57_dumps.py, movie/*, parts.json, any patch scripts and the index.html source; no data, no cache, no frames)
to bsg-files-1009 `docs/handover/bsgprod-1012/pagework/`, with a short README listing your edits relative to bsg-hrdet-1009/page.
Small job, CPU only (git push). DeltaAI builds a first refresh from the base copy meanwhile and will switch to yours when it lands.
