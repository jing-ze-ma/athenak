#!/bin/bash -l
# WASP-189b CPU gates (login node, nice): parse (-n) of a deck, and a fresh t = 0 start
# (nlim 0, ghost zones written) on a reduced angular mesh (C16, MB 8x8) for the IC hydrostatic
# gate.  CPU binary athena_dhj_cpu_8bcd50c8 (fin51-r-fix descendant of 182477553).
# usage: gate.sh MET (1x|3x)   (run dirs must not exist yet)
M=$1; P=/viper/ptmp2/jinma/dhj_planets_1101/wasp189b; G=$P/gate; N=/viper/ptmp2/jinma/builds/bin/athena_dhj_cpu_8bcd50c8
D=$P/deck/w189_$M.athinput
mkdir $G/parse_$M || exit 1; cd $G/parse_$M || exit 1; ln -s $P/ck ck; ln -s $P/ic ic
nice $N -n -i $D > parse.log 2>&1; echo "parse $M rc=$?"
sed 's/^variable .*hydro_w.*/&\nghost_zones = true/' $D > $G/gz_$M.athinput
K2="mesh/nx2=16 mesh/nx3=16 meshblock/nx2=8 meshblock/nx3=8 time/nlim=0 output1/dt=1.0e9 output2/dt=1.0e9 output3/dt=1.0e9 output4/dt=1.0e9"
mkdir $G/run_gz_$M || exit 1; cd $G/run_gz_$M || exit 1; ln -s $P/ck ck; ln -s $P/ic ic
OMP_NUM_THREADS=1 nice timeout 1200 $N -i $G/gz_$M.athinput $K2 > run.log 2>&1; echo "fresh gz $M rc=$?"
