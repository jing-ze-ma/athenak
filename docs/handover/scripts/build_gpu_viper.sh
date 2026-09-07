#!/bin/bash -l
module purge
module load gcc/14 rocm/6.3 openmpi_gpu/5.0 cmake/4.0
cd /viper/u2/jinma/ATHENAK/athenak/build_dhj_gpu && make -j 8 2>&1 | grep -v "Failed to get device count" | tail -5 && ls -la src/athena && echo BUILD_DONE
