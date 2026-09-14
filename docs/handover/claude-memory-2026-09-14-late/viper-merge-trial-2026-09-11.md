---
name: viper-merge-trial-2026-09-11
description: "TRIAL MERGE of origin/polar-average-perf (viper, c26b01ac, 12 commits: rad_tmax_kappa, rt_semi_implicit, time/dt_min, polar_x3_shift, dhj diag/seed) onto the local uncommitted tree: GREEN in worktree regress/wt_merge (merge commit 62109b5f); 3 files/12 hunks resolved (conduction.cpp/.hpp, two_stream_rt.hpp); dhj ck bit-identical to pure viper, red giant bit-identical to local; NOT landed in the main tree yet"
metadata:
  type: project
---
Resolutions to reuse when landing: RadFaceKappa(tmax) carries viper's KappaTemp ceiling into
every kappa_rad site incl. ImplicitRadialUpdate (conduction.cpp:929 `tmax = rad_tmax`); the
semi-implicit branch guard is `if (!explicit_on && semi_imp)` (viper's rt_semi_implicit default
true, local rt_explicit default false; red_giant reads only rt_explicit, dhj only
rt_semi_implicit); viper's dt_min abort runs after the local dt-collapse MinLoc report; hydro_fluxes
auto-merged; polar_x3_shift (sp only, default rotate = behaviour change for sp runs) independent of
cs_corner_poison. Landing sequence: finish the vceil agent (pin10), commit the local work in
pieces, then `git merge origin/polar-average-perf` in the main tree taking the 3 conflicted files
from 62109b5f where the local file is unchanged since the snapshot (two_stream_rt.hpp will need a
re-resolve if the RT NaN guard touched it). Dirs: regress/b_mg_{dhj,rg}, run_ck_{mg,og}, run_rg_mg.
Style: 2219 pre-existing tree-wide violations on BOTH sides (not new).
