#!/bin/bash
# BSG half-range determinism gate (port of viper /viper/ptmp2/jinma/hrdet_1009/gate8n.sh and Raven gate16.sh).
# usage (inside an allocation of 16 ranks, 1 MeshBlock of 256x64x64 per rank):
#   bash gate_hrdet.sh <BO = d385de40 binary> <BN = 98835d99 binary> <files dir from SETUP.sh> <out dir>
# env: SRUN (default "srun -n $SLURM_NTASKS --cpu-bind=cores")
# arms: bit_old / bit_new = keys-off truerepro2 input, 3 cycles, BO vs BN (must be bitwise identical);
#       detH_a / detH_b   = hr input, 30 cycles, BN twice, bins every 5 cycles (must be bitwise identical).
BO=$1; BN=$2; F=$3; S=$4
SRUN=${SRUN:-srun -n $SLURM_NTASKS --cpu-bind=cores}
export OMP_NUM_THREADS=1 KOKKOS_MAP_DEVICE_ID_BY=mpi_rank MPICH_GPU_SUPPORT_ENABLED=1
export HSA_XNACK=1 HSA_NO_SCRATCH_RECLAIM=1   # ROCm-only (no-op on CUDA), job-env rule
I0=$F/bsg3d_truerepro2.athinput; IH=$F/bsg_hr_dc5.athinput
for x in "$BO" "$BN"; do [ -x "$x" ] || { echo "no binary $x"; exit 1; }; done
for x in "$I0" "$IH"; do [ -f "$x" ] || { echo "no input $x"; exit 1; }; done
go(){ D=$S/$1; mkdir -p $D && cd $D || exit 1
  echo "== $1 $(md5sum $2 | cut -c1-32) $3 start $(date +%T)"; t0=$(date +%s)
  $SRUN $2 -i $3 -d $D $4 < /dev/null > out.log 2>&1; rc=$?
  echo "rc=$rc fatal=$(grep -c FATAL out.log) diverged=$(grep -ci DIVERGED out.log)" \
       "nonconv=$(grep -ci NONCONV out.log) wall=$(( $(date +%s)-t0 ))s $(date +%T)"; }
K3="time/nlim=3 time/ndiag=1 output2/dt=262 output3/dt=262 output6/dt=262 output7/dt=262 output5/dt=1e30"
KH="time/nlim=30 time/ndiag=1 output5/dt=1e30"
go bit_old $BO $I0 "$K3"
go bit_new $BN $I0 "$K3"
go detH_a $BN $IH "$KH"
go detH_b $BN $IH "$KH"
