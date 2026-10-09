# TASK 2026-10-09 Caltech -> viper: push the BSG reproduction page scripts (scripts_repro) for the Caltech true repro

The user wants the Caltech BSG true reproduction (chain 4286643-45, /resnick/groups/carnegie_poc/jingze/bsg_1009/prod,
now ~t 3.9 d) on a NEW page compared with Ma+2026 exactly like the "BSG Reproduction Run" page
(https://claude.ai/artifact/KkxYXphcf7mWaajATCd67V: settings table, 57.4-d comparison table, Figs. 2/3/4/6 pairs, KE,
intensity map, movies). That page is built by viper's `bsg_1001/figs/scripts_repro/update_page_repro.sh <pagedir>`
(+ gen_status.py and the figure scripts), which are not in the repo.

Please push to `docs/handover/bsg-hrdet-1009/page/` on bsg-files-1009 (fetch, add, push; never force):
- update_page_repro.sh and every script/module it calls (gen_status.py, the Fig 2/3/4/5/6/7 + KE + intensity + movie
  makers), index.html template if it is generated, and any small reference data they read (paper figure crops
  paper_fig*.png can stay on the artifact; Caltech copies those from it), with md5s;
- one line on arguments/env (run dir(s), hst/bin names for the truerepro2 run: bsg3d.*), and anything viper-specific
  (paths, python env: Caltech has numpy/scipy/matplotlib via spack, no h5py).
Short NOTE when done. Caltech will build the page from its own run and publish it as a new URL.
