#!/bin/bash
# the spatial convergence ladder: dt ~ dx^2 and the same final time, so the operator
# splitting error is the same at every resolution and what converges is the SPACE
# discretisation of the cubed-sphere angular operator.
W=/viper/u2/jinma/ATHENAK/bench/wt_he4_adi
B=$W/build_cpu_chk/src/athena
I=$W/inputs/tests/cubed_sphere_raddiff.athinput
cd "$(dirname "$0")"
run(){ d=$1; shift; rm -rf $d; mkdir -p $d; (cd $d && $B -i $I "$@" > log.txt 2>&1) || { echo "FAIL $d"; tail -6 $d/log.txt; }; }
for p in "16 0.03 20" "32 0.0075 80" "64 0.001875 320"; do set -- $p; n=$1; c=$2; l=$3
  for s in exp sts adi; do
    extra=""
    [ $s = sts ] && extra="hydro/rad_implicit_ang=true hydro/rad_ang_solver=sts"
    [ $s = adi ] && extra="hydro/rad_implicit_ang=true hydro/rad_ang_solver=adi"
    run L_${s}_$n mesh/nx2=$n mesh/nx3=$n meshblock/nx2=$n meshblock/nx3=$n \
        time/cfl_number=$c time/nlim=$l hydro/rad_ang_verbose=true $extra
    echo -n "n=$n "; python3 measure_cs.py L_${s}_$n
    grep -o "|sum V de|/sum V|de| = [0-9.e+-]*" L_${s}_$n/log.txt | tail -1
  done
done
