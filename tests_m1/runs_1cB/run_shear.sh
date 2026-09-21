#!/bin/bash -l
# T4c (1c-B): the sheared advected pulse.  split_vel = recon vs cell at 512, against a
# 4x finer reference.
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_cpu_m1/src/athena
IN=$R/inputs/tests/rad_m1_advect_shear.athinput
O=$R/tests_m1/runs_1cB
run () { local n=$1; shift; local d=$O/shear_$n; rm -rf $d; mkdir -p $d
  (cd $d && $X -i $IN "$@" > run.log 2>&1; echo "exit=$?" >> run.log); }

run ref2048 mesh/nx1=2048 meshblock/nx1=2048 rad_m1/split_vel=recon &
for N in 64 128 256 512; do
  run recon_$N mesh/nx1=$N meshblock/nx1=$N rad_m1/split_vel=recon &
  run cell_$N  mesh/nx1=$N meshblock/nx1=$N rad_m1/split_vel=cell &
done
wait
echo SHEAR_DONE
