#!/bin/bash
# snapshot of the m1-time2b worktree (tracked + new files) into new/
W=/viper/ptmp2/jinma/wt_time2b; D=/viper/ptmp2/jinma/h2_3x/new
mkdir -p $D
cd $W && git ls-files -co --exclude-standard -- src CMakeLists.txt cmake config.hpp.in inputs vis tst 2>/dev/null | rsync -a --files-from=- $W/ $D/
[ -e $D/kokkos ] || ln -s /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $D/kokkos
