#!/bin/bash
# snapshot of the m1-sp2 worktree (tracked + new files, only src/ and CMake) -> src_new
W=/viper/ptmp2/jinma/s2_0924; T=/viper/ptmp2/jinma/wt_m1sp2
rm -rf $W/src_new; mkdir -p $W/src_new
cd $T && git ls-files -co --exclude-standard -- src CMakeLists.txt config.hpp.in cmake inputs tst \
  | rsync -a --files-from=- $T $W/src_new
ln -sfn /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $W/src_new/kokkos
