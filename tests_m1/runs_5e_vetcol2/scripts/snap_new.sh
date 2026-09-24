#!/bin/bash
W=/viper/ptmp2/jinma/vetcol2_0924; T=/viper/ptmp2/jinma/wt_vetcol2
rm -rf $W/src_new; mkdir -p $W/src_new
cd $T && git ls-files -co --exclude-standard -- src CMakeLists.txt config.hpp.in cmake inputs tst \
  | rsync -a --files-from=- $T $W/src_new
ln -sfn /viper/u2/jinma/ATHENAK/bench/wt_rgbox/kokkos $W/src_new/kokkos
