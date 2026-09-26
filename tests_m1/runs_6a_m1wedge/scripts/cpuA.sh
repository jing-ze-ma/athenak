#!/bin/bash
# gate (A) CPU arms, 4 ranks each, run in parallel; usage: cpuA.sh <bin tag> <outdir> [tlim]
source /etc/profile.d/modules.sh; module purge; module load gcc/14 openmpi/5.0
export OMP_NUM_THREADS=1
W=/viper/ptmp2/jinma/sprhd_0926; B=$W/bin/athena_$1_none_cpu; O=$W/cpu/$2; T=${3:-6.0}
IN=$W/inp/gateA.athinput
r() { local d=$1; shift; mkdir -p $O/$d; (cd $O/$d && timeout 3000 mpirun --oversubscribe -np 4 $B -i $IN -d . time/nlim=-1 time/tlim=$T "$@" > log.txt 2>&1; echo "rc=$?" >> log.txt) & }
r r3
r h9 time/cfl_number=0.9
r v4 rad_m1/vet_col_every=4 rad_m1/time2_vet_col=lag
r h9v4 time/cfl_number=0.9 rad_m1/vet_col_every=4 rad_m1/time2_vet_col=lag
r nowb hydro/wellbalance_dynamic=false rad_m1/force_reference=none problem/wg_phi_eff=false
wait
echo CPUA DONE
