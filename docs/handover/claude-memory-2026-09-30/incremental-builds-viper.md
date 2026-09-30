---
name: incremental-builds-viper
description: user 09-28: builds on viper must be incremental (persistent worktree + build dir per target, git checkout of the commit), not fresh git-archive full rebuilds; script /viper/ptmp2/jinma/builds/build_inc_viper.sh (being set up 09-28)
metadata:
  type: feedback
---
Use /viper/ptmp2/jinma/builds/build_inc_viper.sh <target> <commit> (targets dhj_gpu, box_gpu, none_gpu, rg_gpu,
dhj_cpu, box_cpu, none_cpu; README in that dir) instead of git archive + new build dir for every build.
**Why:** every worker was doing 15-25 min full rebuilds (Kokkos included); user 09-28 "I thought we can rebuild by
only changing the files changed". Caltech already does this (build_inc.sh).
**How to apply:** put the script in every worker brief that builds; also reuse an existing binary of the same
commit when one exists (check builds/BUILDS.log and the md5). git archive onto an old tree breaks make (old mtimes).
