#!/bin/bash
# usage: snap.sh <ref|new>: ref = git archive of 7aabf6fd, new = the wt_sporder2 working
# tree (tracked files incl. uncommitted edits); rsync -c keeps the build incremental
W=/viper/ptmp2/jinma/sporder2_0925; T=/viper/ptmp2/jinma/wt_sporder2; a=$1
S=$W/src_$a
if [ $a = ref ]; then c=7aabf6fd; else c=$(git -C $T stash create); c=${c:-HEAD}; fi
rm -rf $W/snap_tmp_$a; mkdir -p $W/snap_tmp_$a $S
git -C $T archive $c | tar -x -C $W/snap_tmp_$a
rm -rf $W/snap_tmp_$a/kokkos
rsync -rc $W/snap_tmp_$a/ $S/
rm -rf $W/snap_tmp_$a
[ -e $S/kokkos ] || ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $S/kokkos
echo "snap $a done ($c)"
