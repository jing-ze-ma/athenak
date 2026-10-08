#!/bin/bash
#SBATCH -J spb2d
#SBATCH -A bivj-delta-gpu
#SBATCH -p gpuA100x4
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=16
#SBATCH --gpus-per-node=4
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:15:00
#SBATCH -o /work/nvme/bivj/jma20/delta_1008/spb2/j.%j.out
# sp-blend2-1008 cost arms (TASK-2026-10-08-delta-spb2-batchN). 1 node x 4 A100, 4 ranks.
# usage: sbatch spb2_delta.sh <binary> <case A|B> <mode prep|arms> <armsfile> [nlim] [wall]
#  prep: ONE fresh run from t = 0 (the first arm line), restart written at cycle nlim
#        (output6 = rst, dcycle = nlim) into $W/prep<case>/rst
#  arms: every arm line in order; case A restarts from the newest rst in $W/prepA/rst,
#        case B starts fresh from t = 0 (the B solve is hard from the first cycles)
#  arm line: <name> <keys...>   (appended after the base keys; a later key wins)
W=/work/nvme/bivj/jma20/delta_1008/spb2
BIN=$1; C=$2; MODE=$3; ARMS=$4; NLIM=${5:-276}; WALL=${6:-00:03:00}
IN=/work/nvme/bivj/jma20/delta_1008/agcar_files/agcar_shake${C}_ge.athinput
[ -x "$BIN" ] || { echo "no binary $BIN"; exit 1; }
module unload cudatoolkit
module load cuda/12.9
export MPICH_GPU_SUPPORT_ENABLED=1
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank
export OMP_NUM_THREADS=1
export HSA_XNACK=1; export HSA_NO_SCRATCH_RECLAIM=1   # ROCm-only (no-op on CUDA), job-env rule
echo "job $SLURM_JOB_ID node $SLURMD_NODENAME case $C mode $MODE arms $ARMS nlim $NLIM"
md5sum $BIN $IN
BASE="time/nlim=$NLIM time/ndiag=1 output1/dcycle=1 output2/dcycle=$NLIM output3/dcycle=$NLIM"
BASE="$BASE output4/dt=1e30 output5/dt=1e30 rad_m1/implicit_flux=blend rad_m1/implicit_blend=tau_f"
while read -r name keys; do
  [ -z "$name" ] && continue
  case $name in \#*) continue;; esac
  if [ "$MODE" = prep ]; then
    D=$W/prep$C; mkdir -p $D; cd $D || exit 1
    echo "prep $name keys: $keys start $(date +%T)"
    srun -n 4 --cpu-bind=cores $BIN -i $IN -d $D -t $WALL $BASE \
        output6/file_type=rst output6/dcycle=$NLIM $keys < /dev/null > run.log 2>&1
    echo "rc=$? fatal=$(grep -c FATAL run.log) $(date +%T)"; ls -l $D/rst
    break
  fi
  D=$W/$C.$name.$SLURM_JOB_ID; mkdir -p $D; cd $D || exit 1
  echo "arm $name keys: $keys start $(date +%T)"
  if [ "$C" = A ]; then
    RST=$(ls -t $W/prepA/rst/*.rst | head -1); echo "restart $RST"
    srun -n 4 --cpu-bind=cores $BIN -r $RST -d $D -t $WALL $BASE output6/dcycle=100000000 \
        $keys < /dev/null > run.log 2>&1
  else
    srun -n 4 --cpu-bind=cores $BIN -i $IN -d $D -t $WALL $BASE $keys < /dev/null > run.log 2>&1
  fi
  echo "rc=$? fatal=$(grep -c FATAL run.log) $(date +%T)"
done < $ARMS
