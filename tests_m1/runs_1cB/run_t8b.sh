#!/bin/bash -l
# T8b: 1 rank vs 4 ranks, bitwise.  Needs the MPI build (build_mpi_m1).  NOTE the known
# ROMIO limitation: the output directories must already exist and cannot be on /tmp.
module purge
module load gcc/14 openmpi/4.1
R=/viper/u2/jinma/ATHENAK/athenak
X=$R/build_mpi_m1/src/athena
O=$R/tests_m1/runs_1cB
IP=$R/inputs/tests/rad_m1_advect_pulse.athinput
IB=$O/t8_beam.athinput
TL="time/tlim=4.8e-6 output1/dt=4.8e-6 output2/dt=4.8e-6"
DYN="rad_m1/kappa_p=500 problem/pulse_v=3.0e7 mesh/nx1=512 meshblock/nx1=128"
for NR in 1 4; do
  d=$O/t8b_t4_r$NR; rm -rf $d; mkdir -p $d/tab
  (cd $d && mpirun -np $NR $X -i $IP $TL $DYN > run.log 2>&1; echo "exit=$?" >> run.log)
  d=$O/t8b_beam_r$NR; rm -rf $d; mkdir -p $d/tab $d/rst
  (cd $d && mpirun -np $NR $X -i $IB > run.log 2>&1; echo "exit=$?" >> run.log)
done
echo T8B_DONE
