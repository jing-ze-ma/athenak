#!/bin/bash
# snapshot of the wt_sph2 working tree (tracked files, uncommitted edits included):
# git archive of `git stash create` (or HEAD when clean); only files that changed are
# rewritten, so the build stays incremental
W=/viper/ptmp2/jinma/sph2_0924; S=$W/src_new; T=/viper/ptmp2/jinma/wt_sph2
c=$(git -C $T stash create); c=${c:-HEAD}
rm -rf $W/snap_tmp; mkdir -p $W/snap_tmp $S
git -C $T archive $c | tar -x -C $W/snap_tmp
rm -rf $W/snap_tmp/kokkos
rsync -rc $W/snap_tmp/ $S/
rm -rf $W/snap_tmp
[ -e $S/kokkos ] || ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos
echo "snap done ($c)"
