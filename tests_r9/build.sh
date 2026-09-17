#!/bin/bash -l
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
cd /viper/u2/jinma/ATHENAK/bench/wt_he4/build_gpu_rg
make -j 16 athena 2>&1 | tail -30
cp src/athena /viper/u2/jinma/ATHENAK/bench/wt_he4/tests_r9/athena_v6
echo BUILD_DONE
