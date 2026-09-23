#!/bin/bash
# CPU gates of time_scheme = hesdirk2 (runs_3x_hesdirk2), 2-D He slab (slab2d_plm_vimp =
# plm + vimp): G0 off-bitwise vs 2c3c4178, NC count, restart (G8), fallback paths, vet_sc.
W=/viper/ptmp2/jinma/h2_3x; R=$W/cpu_run.sh
V="rad_m1/closure=vet_sc rad_m1/vet_tensor=full"
S="time/tlim=200 output3/dt=100"
$R base g0_base 1 slab2d_plm_vimp $S &
$R new  g0_new  1 slab2d_plm_vimp $S &
$R new  be_200  1 slab2d_plm_vimp_be $S &
$R new  h_200   1 slab2d_plm_vimp_h2 $S &
$R new  h_fail  1 slab2d_plm_vimp_h2 $S rad_m1/time2_dbg_fail=50 &
$R new  hv_200  1 slab2d_plm_vimp_h2 $S $V &
$R new  bv_200  1 slab2d_plm_vimp $S $V &
$R new  h_m2    2 slab2d_plm_vimp_h2 $S meshblock/nx2=16 &
$R new  h_1000  1 slab2d_plm_vimp_h2 time/tlim=1000 output3/dt=500 &
$R new  b_1000  1 slab2d_plm_vimp time/tlim=1000 output3/dt=500 &
wait
# G8: restart the hesdirk2 run from its t=100 file (with the slope) and from the be file
# (no slope: first step backward Euler)
d=$W/cpu/h_rst; rm -rf $d; mkdir -p $d; cd $d
module purge; module load gcc/14 openmpi/5.0
OMP_NUM_THREADS=1 nice -n 10 mpirun -np 1 $W/new/build_boxcpu/src/athena \
  -r $W/cpu/h_200/rst/m1slab.00001.rst -d $d > log.txt 2>&1; echo "rc=$?" >> log.txt
d=$W/cpu/h_rst_noslope; rm -rf $d; mkdir -p $d; cd $d
OMP_NUM_THREADS=1 nice -n 10 mpirun -np 1 $W/new/build_boxcpu/src/athena \
  -r $W/cpu/be_200/rst/m1slab.00001.rst -d $d rad_m1/time_scheme=hesdirk2 > log.txt 2>&1
echo "rc=$?" >> log.txt
echo CPU GATE DONE
