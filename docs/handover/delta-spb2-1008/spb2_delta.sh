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
#        (output6 = rst, dcycle = nlim) into $W/prep<case>/rst; exits non-zero on failure
#  arms: every arm line in order; case A restarts from the newest rst in $W/prepA/rst,
#        case B starts fresh from t = 0; exits 2 if any arm failed (rc or FATAL)
#  arm line: <name> <keys...>   (appended after the base keys; a later key wins)
W=/work/nvme/bivj/jma20/delta_1008/spb2
BIN=$1; C=$2; MODE=$3; ARMS=$4; NLIM=${5:-276}; WALL=${6:-00:03:00}
IN0=/work/nvme/bivj/jma20/delta_1008/agcar_files/agcar_shake${C}_ge.athinput
[ -x "$BIN" ] || { echo "no binary $BIN"; exit 1; }
module unload cudatoolkit
module load cuda/12.9
export MPICH_GPU_SUPPORT_ENABLED=1
export KOKKOS_MAP_DEVICE_ID_BY=mpi_rank
export OMP_NUM_THREADS=1
export HSA_XNACK=1; export HSA_NO_SCRATCH_RECLAIM=1   # ROCm-only (no-op on CUDA), job-env rule
echo "job $SLURM_JOB_ID node $SLURMD_NODENAME case $C mode $MODE arms $ARMS nlim $NLIM"
# A per-job COPY of the input that holds every key the command line overrides. A fresh run
# refuses absent keys and absent blocks (and only one rst block is allowed). The copy adds:
# dcycle in output1-3; dcycle in output5 (the input's rst block: prep writes at cycle nlim,
# otherwise never); one extra block output6 = m1_face bin (never, unless an arm names
# output6/dcycle); and the rad_m1 keys the arms set. The A arms restart from the prep copy.
R5=100000000
[ "$MODE" = prep ] && R5=$NLIM
mkdir -p $W/in
IN=$W/in/agcar_shake${C}_ge.$MODE.$SLURM_JOB_ID.athinput
awk -v n=$NLIM -v r5=$R5 '{print} /^<output1>/{print "dcycle = 1"}
  /^<output2>/||/^<output3>/{print "dcycle = " n} /^<output5>/{print "dcycle = " r5}
  /^<rad_m1>/{print "implicit_blend = tau_f"; print "implicit_lin_scaled = false";
              print "implicit_mg_levels = 3"; print "implicit_precond_float = false"}' $IN0 > $IN
{
  echo "<output6>"; echo "file_type = bin"; echo "variable = m1_face"; echo "dcycle = 100000000"
} >> $IN
grep -q "^file_type *= *rst" <(awk '/^<output5>/{f=1;next} /^</{f=0} f' $IN) || { echo "output5 is not rst"; exit 1; }
md5sum $BIN $IN0 $IN
BASE="time/nlim=$NLIM time/ndiag=1 output1/dcycle=1 output2/dcycle=$NLIM output3/dcycle=$NLIM"
BASE="$BASE output4/dt=1e30 rad_m1/implicit_flux=blend rad_m1/implicit_blend=tau_f"
RC=0
while read -r name keys; do
  [ -z "$name" ] && continue
  case $name in \#*) continue;; esac
  if [ "$MODE" = prep ]; then
    D=$W/prep$C; mkdir -p $D; cd $D || exit 1
    echo "prep $name keys: $keys start $(date +%T)"
    srun -n 4 --cpu-bind=cores $BIN -i $IN -d $D -t $WALL $BASE $keys < /dev/null > run.log 2>&1
    rc=$?; nf=$(grep -c FATAL run.log)
    echo "rc=$rc fatal=$nf $(date +%T)"; ls -l $D/rst
    if [ $rc -ne 0 ] || [ $nf -ne 0 ] || ! ls $D/rst/*.rst > /dev/null 2>&1; then exit 2; fi
    exit 0
  fi
  D=$W/$C.$name.$SLURM_JOB_ID; mkdir -p $D; cd $D || exit 1
  echo "arm $name keys: $keys start $(date +%T)"
  if [ "$C" = A ]; then
    RST=$(ls -t $W/prepA/rst/*.rst | head -1); echo "restart $RST"
    srun -n 4 --cpu-bind=cores $BIN -r $RST -d $D -t $WALL $BASE output5/dcycle=100000000 \
        $keys < /dev/null > run.log 2>&1
  else
    srun -n 4 --cpu-bind=cores $BIN -i $IN -d $D -t $WALL $BASE $keys < /dev/null > run.log 2>&1
  fi
  rc=$?; nf=$(grep -c FATAL run.log)
  echo "rc=$rc fatal=$nf $(date +%T)"
  if [ $rc -ne 0 ] || [ $nf -ne 0 ]; then RC=2; fi
done < $ARMS
exit $RC
