# TASK 2026-10-08 for viper: push the He giant page figure script (make_figs.py) for the Caltech N897 run

From: Caltech. The He Giant Expansion page (https://claude.ai/artifact/Wir9zXDuQik8Y2WjCfiKg1) shows viper's scout128 N445
figures at 10.4 d (fig_rph, fig_profiles, fig_slices, fig_shell, fig_hst, made by he_giant_1006/page/make_figs.py on viper).
The user wants all page plots updated from the Caltech N897 continuation (now ~14 d, running to 30 d), but make_figs.py and
its inputs are not in the repo.

Please push to `docs/handover/hegiant-files-1007/page/` on hegiant-opn-1007 (fetch, add, push; never force):
- make_figs.py (+ any helper modules it imports),
- the small reference files it reads (MESA profile / expansion line, PROFILE_CHECK inputs), with md5s,
- one line on how it is run (arguments: run dir(s), which dumps, output dir) and whether it can join a run that changed
  radial grid (N445 -> N897 remap at 10.42 d; Caltech dirs n445_1g and n897 under /resnick/groups/carnegie_poc/jingze/hegiant_1007).
Short NOTE when done. Caltech will regenerate every figure and republish the page.
